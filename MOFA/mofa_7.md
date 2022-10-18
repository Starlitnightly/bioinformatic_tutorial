# MOFA下游分析: Python环节

## 1. 数据准备

### 1.1 导入依赖

```
import h5py
import pandas as pd
import numpy as np
```

### 1.2 导入数据

```
f_ad = h5py.File('mofa_factor_ad.hdf5','r')   #打开h5文件
# 可以查看所有的主键
for key in f_ad.keys():
    print(f_ad[key].name)
```

> /data /expectations /features /groups /intercepts /model_options /samples /training_opts /training_stats /variance_explained /views

### 1.3 导出权重数据

```
h5ls_out_ad=list(f_ad.keys())
if('views' in h5ls_out_ad):
  view_names=f_ad['views']['views'][:]
  group_names=f_ad['groups']['groups'][:]
  feature_names=np.array([f_ad['features'][i][:] for i in view_names])
  sample_names=np.array([f_ad['samples'][i][:] for i in group_names])
```

### 1.4 标准化函数

```
def normalization(data):
    _range = np.max(abs(data))
    return data / _range
```

### 1.5 获取权重函数

```
def get_weights(f,view,factor,scale=True):
  view_names=f['views']['views'][:]
  f_name=feature_names[np.where(view_names==str.encode(view))[0][0]]
  f_w=f['expectations']['W'][view][factor-1]
  if scale==True:
    f_w=normalization(f_w)
  res=pd.DataFrame()
  res['feature']=f_name
  res['weights']=f_w
  res['abs_weights']=abs(f_w)
  res['sig']='+'
  res.loc[(res.weights<0),'sig'] = '-'

  return res
```

### 1.6 获取特定因子的基因权重

```
res1=get_weights(f_ad,'rna',2)
xl=[]
for i in res1['feature'].values:
  xl.append(bytes.decode(i).replace('Mutation',''))
  
res1['feature']=xl
res1.head()
```

|      |   feature |   weights | abs_weights |  sig |
| :--- | --------: | --------: | ----------: | ---: |
| 0    |      HES4 |  0.158418 |    0.158418 |    + |
| 1    |    ATAD3C | -0.025514 |    0.025514 |    - |
| 2    |      HES5 |  0.244909 |    0.244909 |    + |
| 3    | LINC00982 |  0.620071 |    0.620071 |    + |
| 4    |    PRDM16 |  0.815560 |    0.815560 |    + |

### 1.7 可视化特定因子

```
import matplotlib.pyplot as plt


def plot_high_weight(res,n_feature=10):
  #定义字体
  font1={
    'size':15,
  }
  #排序
  res=res.sort_values('abs_weights',ascending=False)
  res_p=res.iloc[:n_feature].sort_values('weights',ascending=True)
  #获取xy
  x=res_p['weights'].values
  y=range(len(res_p['weights']))

  #定义图像大小
  pp=plt.figure(figsize=(4,6))
  ax=pp.add_subplot(1,1,1)
  #定义+的起点
  max_start=len(res_p)-len(res_p.loc[res_p['sig']=='+'])
  #绘制散点图
  ax.scatter(x[:max_start],y[:max_start],c="#00B0E0",alpha=0.5,linewidths=(np.zeros(len(res_p)-max_start)+2))
  ax.scatter(x[max_start:],y[max_start:],c="#E82089",alpha=0.5,linewidths=(np.zeros(len(res_p)-max_start)+2))
  
  #定义文字位置
  ti=len(res_p)-len(res_p.loc[res_p['sig']=='+']) #+
  fi=0 #-
  ti_max=7 #换行标记
  fi_max=7
  for i in range(10):
    #如果为-
    if res_p['sig'].iloc[i]=='-':
      #计算文字长度
      zl=len(res_p['feature'].iloc[i])
      #计算该行剩余长度
      fi_max-=zl

      ax.annotate(res_p['feature'].iloc[i], xy=(x[i],y[i]), xytext=(-0.7+fi_max/7,fi),weight='heavy',
               #bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5),
               arrowprops=dict(arrowstyle='->',color='grey'),size=12)
      #如果剩余长度不足0，那么换行
      if fi_max<0:
        fi_max=10
        fi+=1
    else:
      zl=len(res_p['feature'].iloc[i])
      ax.annotate(res_p['feature'].iloc[i], xy=(x[i],y[i]), xytext=(1-ti_max/7,ti),weight='heavy',
               #bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5),
               arrowprops=dict(arrowstyle='->',color='grey'),size=12)
      ti_max-=zl
      if ti_max<0:
        ti_max=10
        ti+=1

  plt.xlim(-1.5,1.5)
  plt.yticks(fontsize=15)
  plt.xticks(fontsize=15)
  #设置图注
  plt.legend(["-","+"],loc='best',fontsize=12)
  #设置横纵标题
  ax.set_ylabel('Rank',font1)                   
  ax.set_xlabel('Weights',font1)

plot_high_weight(res1,10)
```

![下载 (3)](mofa_7.assets\下载 (3).png)

### 1.8 Var跟Cor数据导入

```
variance=pd.read_csv('/content/drive/MyDrive/GSE174367/ad/variance_explained.csv')
corr=pd.read_csv('/content/drive/MyDrive/GSE174367/ad/correlate_factors.csv')
variance.index=corr.index
```

### 1.9 Var可视化

```


fig, ax = plt.subplots(figsize=(2,6))         # Sample figsize in inches
g=sns.heatmap(variance,
      cmap=palettable.colorbrewer.sequential.Blues_9.mpl_colors, 
      linewidths=.5,annot_kws={"size": 15},ax=ax
      )
g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15,rotation=45, horizontalalignment='right',)
g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 15)
cbar = g.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=15)

```

![下载 (4)](mofa_7.assets\下载 (4).png)

### 1.10 Cor可视化

```python
fig, ax = plt.subplots(figsize=(3,6))         # Sample figsize in inches
g=sns.heatmap(corr,
      cmap=palettable.colorbrewer.sequential.Reds_9.mpl_colors, 
      linewidths=.5,annot_kws={"size": 15},ax=ax,
      )
g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15,rotation=45, horizontalalignment='right',)
g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 15)
cbar = g.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=15)
```

![下载 (5)](mofa_7.assets\下载 (5).png)

## 2. 单细胞环节

### 2.1 导入数据

```
rna=anndata.read_h5ad("rna_pair_mofa.h5ad")
atac=anndata.read_h5ad("atac_pair_mofa.h5ad")
ret3= list(set(rna.obs.index).intersection(atac.obs.index))
```

#### 2.1.1 rna数据

```
rna_sex=rna[ret3]
for i in list(set(rna_sex.obs['cell_type'])):
  print(i,len(rna_sex.obs.loc[rna_sex.obs['cell_type']==i])/len(rna_sex.obs))
```

> PER.END 0.0039946737683089215 
>
> EX 0.07789613848202397 
>
> ODC 0.7200399467376831 
>
> MG 0.046937416777629824 
>
> INH 0.05858854860186418 
>
> OPC 0.027796271637816245 
>
> ASC 0.06474700399467377

#### 2.1.2 atac数据

```
atac_sex=atac[ret3]
for i in list(set(atac_sex.obs['cell_type'])):
  print(i,len(atac_sex.obs.loc[atac_sex.obs['cell_type']==i])/len(atac_sex.obs))
```

> PER.END 0.003828229027962716 
>
> EX 0.07440079893475367 
>
> ODC 0.7187083888149135 
>
> MG 0.051098535286284955 
>
> INH 0.05858854860186418 
>
> OPC 0.02796271637816245 
>
> ASC 0.06541278295605858

```
ret_sex= list(set(rna_sex.obs.index).intersection(atac_sex.obs.index))
len(ret_sex)
```

### 2.2 scRNA-seq与scATAC-seq细胞类型比较

```
cell_type=list(set(atac_sex[ret_sex].obs['cell_type']))
test_df=pd.DataFrame(columns=cell_type,index=cell_type)
for i in cell_type:
  for j in cell_type:
    test_index=atac_sex[ret_sex].obs[atac_sex[ret_sex].obs['cell_type']==i].index
    test_df.loc[i,j]=len(rna_sex[test_index].obs[rna_sex[test_index].obs['cell_type']==j])/len(test_index)
test_df
```

|         |       ODC |    PER.END |         EX |        INH |      ASC |        OPC |       MG |
| :------ | --------: | ---------: | ---------: | ---------: | -------: | ---------: | -------: |
| ODC     |  0.994908 |          0 | 0.00203666 | 0.00152749 |        0 | 0.00152749 |        0 |
| PER.END | 0.0833333 |   0.916667 |          0 |          0 |        0 |          0 |        0 |
| EX      |         0 | 0.00222222 |   0.991111 | 0.00666667 |        0 |          0 |        0 |
| INH     |         0 |          0 | 0.00510204 |   0.994898 |        0 |          0 |        0 |
| ASC     | 0.0135922 |          0 |  0.0291262 |  0.0038835 | 0.953398 |          0 |        0 |
| OPC     |         0 |          0 | 0.00813008 |          0 |        0 |    0.99187 |        0 |
| MG      |  0.113514 |          0 |  0.0648649 | 0.00540541 | 0.027027 |  0.0108108 | 0.778378 |

```
import seaborn as sns
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4,4))         # Sample figsize in inches
new_blues=sns.color_palette("Reds", 1000)[0:700]
g=sns.heatmap(test_df.astype(float),square=True,cmap=new_blues,
       linewidths=.5,annot_kws={"size": 15},ax=ax)
g.set_xticklabels(g.get_xmajorticklabels(), fontsize = 15,rotation=45, horizontalalignment='right',)
g.set_yticklabels(g.get_ymajorticklabels(), fontsize = 15,rotation=360)
cbar = g.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=15)
plt.savefig("fig1_b_corr.png",dpi=300,bbox_inches = 'tight')
```

<img src="mofa_7.assets\下载 (6).png" alt="下载 (6)" style="zoom:50%;" />

### 2.3 scRNA-seq的AD与Ctrl组聚类

#### 2.3.1 导入数据

```
ad_meta=pd.read_csv('/content/drive/MyDrive/GSE174367/mofa_sex_ad_meta.csv')
ctrl_meta=pd.read_csv('/content/drive/MyDrive/GSE174367/mofa_sex_ctrl_meta.csv')

res_ad=rna_sex[ad_meta['sample']]
res_ctrl=rna_sex[ctrl_meta['sample']]
```

#### 2.3.2 scRNA-seq的factor设置

```
f_ad = h5py.File('mofa_factor_ad.hdf5','r')   #打开h5文件
f_ctrl = h5py.File('mofa_factor_ctrl.hdf5','r')   #打开h5文件

for i in range(f_ad['expectations']['Z']['group0'].shape[0]):
  res_ad.obs['factor{0}'.format(i+1)]=f_ad['expectations']['Z']['group0'][i]
for i in range(f_ctrl['expectations']['Z']['group0'].shape[0]):
  res_ctrl.obs['factor{0}'.format(i+1)]=f_ctrl['expectations']['Z']['group0'][i]
```

#### 2.3.3 绘图设置

```
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white',figsize=14)
scnmt=['#f2535d','#a2e8c4','#e5e24f','#6aafda','#c76fc7','#c4425d','#f48a3d']
```

#### 2.3.4 AD组绘图

```
sc.pp.neighbors(res_ad, metric="cosine")
sc.tl.umap(res_ad,random_state=41822099)
sc.settings.set_figure_params(dpi=80, facecolor='white',figsize=[4,4])
sc.pl.umap(res_ad,color=["cell_type","factor1","factor2","factor4"], wspace=0.2, palette=scnmt,
      cmap='coolwarm',vmin=-5,vmax=5,frameon=False,add_outline=True,outline_color=('white','white'))
```

![下载 (7)](mofa_7.assets\下载 (7).png)

#### 2.3.5 Ctrl组绘图

```
sc.pp.neighbors(res_ctrl, metric="cosine")
sc.tl.umap(res_ctrl,random_state=41822099)
sc.settings.set_figure_params(dpi=80, facecolor='white',figsize=[4,4])
sc.pl.umap(res_ctrl,color=["cell_type","factor1","factor2","factor4"], wspace=0.2, palette=scnmt,
      cmap='coolwarm',vmin=-5,vmax=5,frameon=False,add_outline=True,outline_color=('white','white'))
```

![下载 (8)](mofa_7.assets\下载 (8).png)