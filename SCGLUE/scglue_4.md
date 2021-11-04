# 单细胞样本对齐： 模型训练

在本小节，我们就要开始对齐单细胞了

## 1. 数据准备

### 1.1 导入包

```python
import anndata
import itertools
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scglue
from matplotlib import rcParams
```

### 1.2 导入数据

```python
rna = anndata.read_h5ad("rna_preprocessed.h5ad")
atac = anndata.read_h5ad("atac_preprocessed.h5ad")
graph = nx.read_graphml("prior.graphml.gz")
```

## 2 模型训练

### 2.1 设置模型参数

```python
scglue.plot.set_publication_params()
rcParams["figure.figsize"] = (4, 4)

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="raw", use_rep="X_pca"
)

scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)

graph = graph.subgraph(itertools.chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
))
```

### 2.2 模型训练

```python

glue = scglue.models.SCGLUEModel(
    {"rna": rna, "atac": atac}, sorted(graph.nodes),
    random_seed=0
)

glue.compile()

glue.fit(
    {"rna": rna, "atac": atac},
    graph, edge_weight="weight", edge_sign="sign",
    directory="glue"
)
```

> ```
> [INFO] SCGLUEModel: Setting `graph_batch_size` = 28027
> [INFO] SCGLUEModel: Setting `align_burnin` = 10
> [INFO] SCGLUEModel: Setting `max_epochs` = 38
> [INFO] SCGLUEModel: Setting `patience` = 5
> [INFO] SCGLUEModel: Setting `reduce_lr_patience` = 3
> [INFO] SCGLUETrainer: Using training directory: "glue"
> [INFO] SCGLUETrainer: [Epoch 10] train={'g_nll': 0.383, 'g_kl': 0.004, 'g_elbo': 0.387, 'x_rna_nll': 0.246, 'x_rna_kl': 0.009, 'x_rna_elbo': 0.254, 'x_atac_nll': 0.029, 'x_atac_kl': 0.0, 'x_atac_elbo': 0.03, 'dsc_loss': 0.684, 'gen_loss': 0.286}, val={'g_nll': 0.425, 'g_kl': 0.004, 'g_elbo': 0.429, 'x_rna_nll': 0.246, 'x_rna_kl': 0.009, 'x_rna_elbo': 0.254, 'x_atac_nll': 0.029, 'x_atac_kl': 0.0, 'x_atac_elbo': 0.029, 'dsc_loss': 0.692, 'gen_loss': 0.287}, 24.2s elapsed
> [INFO] SCGLUETrainer: [Epoch 20] train={'g_nll': 0.338, 'g_kl': 0.004, 'g_elbo': 0.342, 'x_rna_nll': 0.244, 'x_rna_kl': 0.009, 'x_rna_elbo': 0.253, 'x_atac_nll': 0.029, 'x_atac_kl': 0.0, 'x_atac_elbo': 0.029, 'dsc_loss': 0.688, 'gen_loss': 0.282}, val={'g_nll': 0.442, 'g_kl': 0.004, 'g_elbo': 0.447, 'x_rna_nll': 0.245, 'x_rna_kl': 0.008, 'x_rna_elbo': 0.253, 'x_atac_nll': 0.029, 'x_atac_kl': 0.0, 'x_atac_elbo': 0.029, 'dsc_loss': 0.691, 'gen_loss': 0.286}, 24.3s elapsed
> Epoch    21: reducing learning rate of group 0 to 2.0000e-04.
> Epoch    21: reducing learning rate of group 0 to 2.0000e-04.
> 2021-10-22 16:24:42,283 ignite.handlers.early_stopping.EarlyStopping INFO: EarlyStopping: Stop training
> [INFO] EarlyStopping: Retoring checkpoint "22"...
> ```

### 2.3 保存模型

```python
glue.save("glue/final.dill")
```