# 单细胞样本对齐： 数据预处理

## 1. 数据获取

在本研究中，我们使用了2021年7月在Nature genetic上发表的一篇文章，该文章同时使用了snATAC-seq与snRNA-seq两种不同的技术进行测量阿尔茨海默症晚期病人的皮层。

```python
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/GSE174367_snATAC-seq_cell_meta.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/GSE174367_snRNA-seq_cell_meta.csv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174367/suppl/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5
```

## 2. RNA-seq数据预处理

### 2.1 导入包

```python
import anndata
import networkx as nx
import scanpy as sc
import pandas as pd
import numpy as np
```

### 2.2 导入数据

```python
adata=sc.read_10x_h5('/content/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()
adata.obs.index.name, adata.var.index.name = "cells", "genes"
```

### 2.3 设置meta

```python
meta=pd.read_csv('/content/GSE174367_snRNA-seq_cell_meta.csv.gz')
meta.set_index(meta.columns[0],inplace=True)
```

### 2.4 去除cell_type未知的细胞

```python
xl=[]
for i in adata.obs.index:
  if(i in meta.index):
    xl.append(meta.loc[i]['Cell.Type'])
  else:
    xl.append('Nan')
adata.obs["cell_type"] = xl

adata=adata[adata.obs['cell_type']!='Nan']

Cell=[]
Sample=[]
Batch=[]
Age=[]
Sex=[]
PMI=[]
Tangle=[]
Plaque=[]
Diagnosis=[]
RIN=[]
for i in rna.obs.index:
  if(i in meta.index):
    test=meta.loc[i]
    Cell.append(test['Cell.Type'])
    Sample.append(test['SampleID'])
    Batch.append(test['Batch'])
    Age.append(test['Age'])
    Sex.append(test['Sex'])
    PMI.append(test['PMI'])
    Tangle.append(test['Tangle.Stage'])
    Plaque.append(test['Plaque.Stage'])
    Diagnosis.append(test['Diagnosis'])
    RIN.append(test['RIN'])
  else:
    Cell.append('Nan')
    Sample.append('Nan')
    Batch.append('Nan')
    Age.append('Nan')
    Sex.append('Nan')
    PMI.append('Nan')
    Tangle.append('Nan')
    Plaque.append('Nan')
    Diagnosis.append('Nan')
    RIN.append('Nan')
rna.obs["cell_type"] = Cell
rna.obs["Sample.ID"] = Sample
rna.obs["Batch"] = Batch
rna.obs["Age"] = Age
rna.obs["Sex"] = Sex
rna.obs["PMI"] = PMI
rna.obs["Tangle.Stage"] = Tangle
rna.obs["Plaque.Stage"] = Plaque
rna.obs["Diagnosis"] = Diagnosis
rna.obs["RIN"] = RIN
```

### 2.4 保存数据

```python
rna.write_h5ad('GSE174367rna_61472.h5ad',compression="gzip")
```



## 3. ATAC-seq数据预处理

我们得到的ATAC-seq矩阵较为复杂，所以我们需要对数据进行一定的预处理

### 3.1 导入包

```python
import anndata
import networkx as nx
import scanpy as sc
import pandas as pd
import numpy as np
```

### 3.2 导入数据

```python
atacdata=sc.read_10x_h5('/content/GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5',gex_only=False)
atacdata.obs.index.name, atacdata.var.index.name = "cells", "peaks"
```

### 3.3 染色体位置拆分

```python
atacdata.var["chrom"] = np.vectorize(lambda x: x.split(":")[0])(atacdata.var["gene_ids"])
atacdata.var["chromStart"] = np.vectorize(lambda x: int(x.split(":")[1].split("-")[0]))(atacdata.var["gene_ids"])
atacdata.var["chromEnd"] = np.vectorize(lambda x: int(x.split("-")[1]))(atacdata.var["gene_ids"])
del atacdata.var["gene_ids"]
atacdata.var.head()
```

### 3.4 过滤表达为0的peak

```
sc.pp.filter_genes(atacdata, min_counts=1)
atacdata
```

> AnnData object with n_obs × n_vars = 143401 × 217707    var: 'feature_types', 'genome', 'chrom', 'chromStart', 'chromEnd', 'n_counts'

### 3.5 设置细胞meta

```python
#导入meta
meta1=pd.read_csv('/content/GSE174367_snATAC-seq_cell_meta.csv.gz')
#设置barcode为index
meta1.set_index(meta1.columns[11],inplace=True)
#设置细胞的meta
Cell=[]
Sample=[]
Batch=[]
Age=[]
Sex=[]
PMI=[]
Tangle=[]
Plaque=[]
Diagnosis=[]
RIN=[]
for i in atacdata.obs.index:
  if(i in meta1.index):
    test=meta1.loc[i]
    Cell.append(test['Cell.Type'])
    Sample.append(test['Sample.ID'])
    Batch.append(test['Batch'])
    Age.append(test['Age'])
    Sex.append(test['Sex'])
    PMI.append(test['PMI'])
    Tangle.append(test['Tangle.Stage'])
    Plaque.append(test['Plaque.Stage'])
    Diagnosis.append(test['Diagnosis'])
    RIN.append(test['RIN'])
  else:
    Cell.append('Nan')
    Sample.append('Nan')
    Batch.append('Nan')
    Age.append('Nan')
    Sex.append('Nan')
    PMI.append('Nan')
    Tangle.append('Nan')
    Plaque.append('Nan')
    Diagnosis.append('Nan')
    RIN.append('Nan')
atacdata.obs["cell_type"] = Cell
atacdata.obs["Sample.ID"] = Sample
atacdata.obs["Batch"] = Batch
atacdata.obs["Age"] = Age
atacdata.obs["Sex"] = Sex
atacdata.obs["PMI"] = PMI
atacdata.obs["Tangle.Stage"] = Tangle
atacdata.obs["Plaque.Stage"] = Plaque
atacdata.obs["Diagnosis"] = Diagnosis
atacdata.obs["RIN"] = RIN
atacdata.obs.head()
```

### 3.6 过滤细胞

由于ATAC-seq有130418个细胞，而RNA-seq只有61472个细胞，于是我们就需要随机去除ATAC-seq里面多的细胞，同时要保证细胞类型的比例不变，所以我们按类随机去除细胞

```python
atacdata.obs['ran']=np.zeros(len(atacdata.obs))

cell_type=list(set(rna.obs['cell_type']))
for i in cell_type:
  cell_len=len(rna[rna.obs['cell_type']==i])
  #random select
  a1=atacdata.obs[atacdata.obs['cell_type']==i].sample(n=cell_len).index
  atacdata.obs.loc[a1,'ran']=1

atacdata1=atacdata[atacdata.obs['ran']==1]
atacdata1
```

> View of AnnData object with n_obs × n_vars = 61472 × 217707    obs: 'cell_type', 'ran', 'Sample.ID', 'Batch', 'Age', 'Sex', 'PMI', 'Tangle.Stage', 'Plaque.Stage', 'Diagnosis', 'RIN'    var: 'feature_types', 'genome', 'chrom', 'chromStart', 'chromEnd', 'n_counts'

### 3.7 保存数据

```python
atacdata1.write_h5ad('GSE174367atac_61472.h5ad',compression="gzip")
```

