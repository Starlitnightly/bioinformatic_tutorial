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

```
rna = anndata.read_h5ad("rna_preprocessed.h5ad")
atac = anndata.read_h5ad("atac_preprocessed.h5ad")
glue = scglue.models.load_model("final.dill")
```

## 2. 整合模型

### 2.1 数据预处理

```
rna.obs['domain']='scRNA-seq'
atac.obs['domain']='scATAC-seq'
rna.obsm["X_glue"] = glue.encode_data("rna", rna)
atac.obsm["X_glue"] = glue.encode_data("atac", atac)
```

### 2.2 导出细胞参数

```
rna_loc=pd.DataFrame(rna.obsm['X_glue'], index=rna.obs.index)
atac_loc=pd.DataFrame(atac.obsm['X_glue'], index=atac.obs.index)
```

### 2.3 配对细胞

```
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

```
pair=pd.read_csv('pair_res.csv')
pair.columns=['scATAC','scRNA']
```

### 3.2 过滤重复细胞

```
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

```
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

```
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

