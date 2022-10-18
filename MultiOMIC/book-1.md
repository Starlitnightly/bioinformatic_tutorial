## 1. 数据预处理

在本章中，我们将介绍如何对数据进行初步的处理，数据来源分别为`Cellranger`,`Cellranger-atac`,`velocyto`三个上游分析结果的数据，格式分别为：

- cellranger：`filtered_feature_bc_matrix.h5`或`filtered_feature_bc_matrix`文件夹
- cellranger-atac: `filtered_peak_bc_matrix`文件夹下的`matrix.mtx`,`barcodes.tsv`,`peaks.bed`三个文件以及`filtered_peak_bc_matrix.h5`文件（如果有的话最好）
- velocyto: `possorted_genome_bam.loom`文件

### 1.1 RNA数据预处理

我们将`filtered_feature_bc_matrix.h5`与`possorted_genome_bam.loom`文件放在同一个目录下，然后使用scvelo的教程读取这两个文件，然后保存成h5ad格式

```python
import scanpy as sc
import scvelo as scv

#读取矩阵文件
adata=sc.read_10x_h5('filtered_feature_bc_matrix.h5')
#使得obs跟var名唯一
adata.var_names_make_unique()
adata.obs_names_make_unique()
#读取velocyto文件
ldata = scv.read('possorted_genome_bam.loom', cache=True)
#合并两文件
adata = scv.utils.merge(adata, ldata)
#保存文件
adata.write_h5ad('rna_raw.h5ad',compression='gzip')
```

### 1.2 ATAC数据预处理

我们找到`filtered_peak_bc_matrix`文件夹下的`matrix.mtx`,`barcodes.tsv`,`peaks.bed`三个文件，使用episcanpy读取后保存

```python
import episcanpy
adata=episcanpy.pp.read_ATAC_10x('matrix.mtx', \
                                cell_names='barcodes.tsv', \
                                var_names='peaks.bed')
adata.write_h5ad('atac_raw.h5ad',compression="gzip")
```

### 1.3 ATAC数据生成Gene-Activity矩阵

这里提供两种不同的方法达到此目的，一个较简单，一个较繁琐。区别在于，简单的方法你不一定有文件，但繁琐的方法你一定能实现
#### 方法1: 使用`filtered_peak_bc_matrix.h5`文件

然后再运行下面的代码生成一个`gene_activity_gene_score.h5`文件
```shell
MAESTRO scatac-genescore \
--format h5 \
--peakcount filtered_peak_bc_matrix.h5 \
--genedistance 10000 \
--species GRCh38 \
--model Enhanced \
-d /data/result \
--outprefix gene_activity
```

#### 方法2：使用前面生成的`atac_raw.h5ad`文件

```python
import scanpy as sc
atac=sc.read('../data/raw_data/atac.h5ad')
atac.to_df().T.to_csv('../data/raw_data/brca_atac.tsv', sep='\t')
```

然后再运行下面的代码生成一个`gene_activity_gene_score.h5`文件
```shell
MAESTRO scatac-genescore \
--format plain \
--peakcount brca_atac.tsv \
--genedistance 10000 \
--species GRCh38 \
--model Enhanced \
-d /data/result \
--outprefix gene_activity
```

### 1.4 Gene-Activity矩阵保存为h5ad

这一步骤稍微有一些繁琐，不过代码照着运行就好了

```python
import os
import collections
import tables
import h5py
import scipy.io
import csv
import gzip
import scipy.sparse as sp_sparse
import argparse as ap
import pandas as pd
FeatureBCMatrix = collections.namedtuple('FeatureBCMatrix', ['ids', 'names', 'barcodes', 'matrix'])
def read_10X_h5(filename):
    """Read 10X HDF5 files, support both gene expression and peaks."""
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, 'matrix')
        except tables.NoSuchNodeError:
            print("Matrix group does not exist in this file.")
            return None
        feature_group = getattr(group, 'features')
        ids = getattr(feature_group, 'id').read()
        names = getattr(feature_group, 'name').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return FeatureBCMatrix(ids, names, barcodes, matrix)
    
```

```python
import re
scatac_count=read_10X_h5('gene_activity_gene_score.h5')
peakmatrix = scatac_count.matrix
features = scatac_count.names.tolist()
features = [re.sub("\W", "_", feature.decode()) for feature in features]
features = [feature.encode() for feature in features]
barcodes = scatac_count.barcodes.tolist()
adata=anndata.AnnData(peakmatrix.T,obs=barcodes,var=features)
adata.obs.index=[i.decode('utf-8') for i in adata.obs[0]]
adata.var.index=[i.decode('utf-8') for i in adata.var[0]]
del adata.obs[0]
del adata.var[0]
adata.write_h5ad('gene_activity_gene_score.h5ad',compression='gzip')
```

到此，我们数据分析所需要的所有基本文件就准备完成了
