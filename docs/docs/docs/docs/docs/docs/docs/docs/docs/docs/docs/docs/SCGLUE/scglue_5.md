# 单细胞样本对齐： 对齐细胞

在上一小节，我们已经获得了训练好的对齐模型，在本小节中，我们将根据模型参数配对细胞。

## 1. 数据准备

### 1.1 导入包

```python
import anndata
import networkx as nx
import scanpy as sc
import scglue
import numpy as np
import pandas as pd
from matplotlib import rcParams
```

### 1.2 导入数据

```python
rna = anndata.read_h5ad("rna_preprocessed.h5ad")
atac = anndata.read_h5ad("atac_preprocessed.h5ad")
glue = scglue.models.load_model("final.dill")
```

## 2. 整合模型

### 2.1 数据预处理

```python
rna.obs['domain']='scRNA-seq'
atac.obs['domain']='scATAC-seq'
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
```

### 2.2 导出细胞参数(特征向量)

```python
rna_loc=pd.DataFrame(rna.obsm['X_glue'], index=rna.obs.index)
atac_loc=pd.DataFrame(atac.obsm['X_glue'], index=atac.obs.index)
```

### 2.3 配对细胞

```python
len1=(len(rna_loc)//5000)+1
xl=[]
for j in range(len1):
  
  c=pd.DataFrame()
  for i in range(len1):
    t1=rna_loc.iloc[5000*(i):5000*(i+1)]
    t2=atac_loc.iloc[5000*(j):5000*(j+1)]
    a=np.corrcoef(t1,t2)[len(t1):,0:len(t1)]
    b=pd.DataFrame(a,index=t2.index,columns=t1.index)  
    
    c=pd.concat([c,b],axis=1)
    del t1
    del t2
    del a
    del b
  for i in range(len(c)):
    xl.append(c.columns[np.where(c.iloc[i]==c.iloc[i].max())[0]].values[0])
  del c
  print('Now epoch is {}'.format(j)) 
res=pd.DataFrame(index=atac_loc.index)
res['pair']=xl
res.to_csv('pair_res.csv')
```

## 3. 过滤重复细胞

### 3.1 导入配对数据

```python
pair=pd.read_csv('pair_res.csv')
pair.columns=['scATAC','scRNA']
```

### 3.2 过滤重复细胞

```python
rna_only_pair=list(set(pair['scRNA']))
atac_only_pair=[]
for i in rna_only_pair:
  atac_only_pair.append(pair[pair['scRNA']==i]['scATAC'].iloc[0])
new_pair=pd.DataFrame()
new_pair['scRNA']=rna_only_pair
new_pair['scATAC']=atac_only_pair
new_name=[]
for i in range(len(new_pair)):
  k='cell_{0}'.format(i)
  new_name.append(k)
new_pair['sample']=new_name
new_pair.to_csv('new_pair.csv')
```

### 3.3 提取scRNA-seq的配对细胞

```python
delli=[]
for i in rna.obs.index:
  if i in rna_only_pair:
    delli.append('tr')
  else:
    delli.append('fa')
rna.obs['delli']=delli
rna_pair=rna[rna.obs['delli']=='tr']
rna_pair.write_h5ad('rna_pair.h5ad',compression="gzip")
```

> View of AnnData object with n_obs × n_vars = 23770 × 35363   

### 3.4 提取scATAC-seq的配对细胞

```python
delli=[]
for i in atac.obs.index:
  if i in atac_only_pair:
    delli.append('tr')
  else:
    delli.append('fa')
atac.obs['delli']=delli
atac_pair=atac[atac.obs['delli']=='tr']
atac_pair.write_h5ad('atac_pair.h5ad',compression="gzip")
```

### 3.5 scRNA-seq与scATAC-seq的细胞类型配对

由于我们配对的细胞仅仅是特征向量相似，而不是真正意义上的相似，可能存在scRNA-seq的男性细胞与scATAC-seq的女性细胞比较相似的情况，为了避免性别与诊断出错，我们对配对后的细胞进行进一步配对

#### 3.5.1 男女性别平衡

由于女性样本比男性样本要多两个，我们随机去除两个女性样本使得男：女=1：1

```python
rna_sex=rna_pair[(rna_pair.obs['Sample.ID']!='Sample-22') & (rna_pair.obs['Sample.ID']!='Sample-27')]
rna_sex.obs.index=rna_pair.obs[(rna_pair.obs['Sample.ID']!='Sample-22') & (rna_pair.obs['Sample.ID']!='Sample-27')].index.values
for i in list(set(rna_sex.obs['cell_type'])):
  print(i,len(rna_sex.obs.loc[rna_sex.obs['cell_type']==i])/len(rna_sex.obs))
```

> ODC 0.622879391860573 
>
> MG 0.0530731435987763 
>
> INH 0.0965977565588208 
>
> EX 0.11754890145545564 
>
> OPC 0.03680355984054881 
>
> ASC 0.06711782701399834 
>
> PER.END 0.0059794196718271995

```python
atac_sex=atac_pair[(atac_pair.obs['Sample.ID']!='Sample-22') & (atac_pair.obs['Sample.ID']!='Sample-27')]
atac_sex.obs.index=atac_pair.obs[(atac_pair.obs['Sample.ID']!='Sample-22') & (atac_pair.obs['Sample.ID']!='Sample-27')].index
for i in list(set(atac_sex.obs['cell_type'])):
  print(i,len(atac_sex.obs.loc[atac_sex.obs['cell_type']==i])/len(atac_sex.obs))
```

> ODC 0.6223660050712599 
>
> MG 0.05993704642825916 
>
> INH 0.10037597272011892 
>
> EX 0.10457287750284165 
>
> OPC 0.03733496546297106 
>
> ASC 0.06933636443123196 
>
> PER.END 0.006076768383317304

```
ret_sex= list(set(rna_sex.obs.index).intersection(atac_sex.obs.index))
len(ret_sex)
```

#### 3.5.2 细胞类型配对-性别

```
rna_sex_F=rna_sex[rna_sex.obs['Sex']=='F']
rna_sex_F.obs.index=rna_sex.obs[rna_sex.obs['Sex']=='F'].index
atac_sex_F=atac_sex[atac_sex.obs['Sex']=='F']
atac_sex_F.obs.index=atac_sex.obs[atac_sex.obs['Sex']=='F'].index
ret_F=list(set(rna_sex_F.obs.index).intersection(atac_sex_F.obs.index))

rna_sex_M=rna_sex[rna_sex.obs['Sex']=='M']
rna_sex_M.obs.index=rna_sex.obs[rna_sex.obs['Sex']=='M'].index
atac_sex_M=atac_sex[atac_sex.obs['Sex']=='M']
atac_sex_M.obs.index=atac_sex.obs[atac_sex.obs['Sex']=='M'].index
ret_M=list(set(rna_sex_M.obs.index).intersection(atac_sex_M.obs.index))

```

#### 3.5.3 细胞类型配对-诊断

```
rna_sex_ad=rna_sex[ret_M+ret_F][rna_sex[ret_M+ret_F].obs['Diagnosis']=='AD']
atac_sex_ad=atac_sex[ret_M+ret_F][atac_sex[ret_M+ret_F].obs['Diagnosis']=='AD']
ret_ad=list(set(rna_sex_ad.obs.index).intersection(atac_sex_ad.obs.index))

rna_sex_ctrl=rna_sex[ret_M+ret_F][rna_sex[ret_M+ret_F].obs['Diagnosis']=='Control']
atac_sex_ctrl=atac_sex[ret_M+ret_F][atac_sex[ret_M+ret_F].obs['Diagnosis']=='Control']
ret_ctrl=list(set(rna_sex_ctrl.obs.index).intersection(atac_sex_ctrl.obs.index))
```

### 3.6 保存新配对列表

```python
pair=pd.read_csv('new_pair.csv')
pair.set_index(pair.columns[0],inplace=True)
pair.set_index(pair.columns[2],inplace=True)
#pair.columns=['scATAC','scRNA']
pair.loc[ret_ad+ret_ctrl].to_csv('mofa_pre_pair.csv')
!cp /content/mofa_pre_pair.csv /content/drive/MyDrive/mofa_gse174367
```

