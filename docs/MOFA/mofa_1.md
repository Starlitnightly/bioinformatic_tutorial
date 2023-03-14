# MOFA分析： 环境配置

mofa的整个分析，也是分为Python与R两部分，可以说，复杂的生物信息学分析，都与Python跟R密切相关，而对于Mofa分析而言，或许更为复杂。

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

### 1.5 mofapy2安装

我们接下来安装用于多组学因子训练的模型mofa

```
pip install mofapy2
```

### 1.6 episcanpy安装

最后，是表观基因组的单细胞处理的包episcanpy的安装

```shell
pip install git+https://github.com/colomemaria/epiScanpy
```

以上，就是Python环境所需要的全部依赖了

## 2. R环境

### 2.1 安装MOFA2

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MOFA2")
```

### 2.2 安装data.table

```R
install.packages('data.table',type='source')#源码的形式安装
```

### 2.3 安装scater

```R
BiocManager::install("scater",type='source')
```

### 2.4 安装其他

上述三个包的安装比较特别，就单独列了出来，其他的就用下面的命令安装就好

```R
install.packages(c('purrr', 'ggplot2', 'reticulate', 'argparse','RColorBrewer'))
```

到这里，我们复现Nature所需要的R环境就基本配置好了