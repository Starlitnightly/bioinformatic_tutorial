# 单细胞样本对齐： 环境配置

单细胞样本对齐，只需要在Python环境即可完成

## 1. Python部分

在Python部分，主要有以下几个包需要被安装：scglue，scanpy，mofapy2以及episcanpy

### 1.1 conda环境

在Python部分需要安装的包，可能与过往需要的包会起到冲突，所以我们新建一个conda环境

```python
conda create -n rna python=3.6
conda activate rna
```

通过上述两行代码，我们现在进入了一个叫rna的python虚拟环境，这个环境是非常干净的，没有什么多余的包，所以我们在下一步中将依次安装需要的依赖

### 1.2 Jupyterlab安装

由于这是一个新的python环境，所以我们需要重新装一下jupyter

```
conda install -c conda-forge jupyterlab
```

### 1.3 scanpy安装

接下来进入正题，我们安装单细胞处理所必需的包-scanpy

```python
conda install seaborn scikit-learn statsmodels numba pytables
conda install -c conda-forge python-igraph leidenalg
pip install scanpy
```

### 1.4 scglue安装

安装完scanpy后，我们安装一下用于单细胞配对的包scglue

```shell
conda install -c defaults -c pytorch -c bioconda -c conda-forge -c scglue scglue --yes
```

### 