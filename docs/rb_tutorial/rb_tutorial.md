# Base：环境配置

我们本次教程在Kaggle上运行，需要的数据请向管理员咨询权限

## 1. 视网膜母细胞瘤教程环境配置
!!! note 
    环境一共分三步安装，其中自行选择是否需要装相应的包，包之间有严格的安装顺序，切勿自行调换
    - episcanpy
    - scVelo
    - scvi-tools
    - scglue
    - scrublet
    - pyscenic
    - Pyomic
    - scbasset(可选)
    - pycisTopic(可选)
    - scenicplus（可选）

### 必须安装的包-conda
```Python
!conda install -c conda-forge mamba -y
!mamba install -c annadanese -c conda-forge -c bioconda episcanpy scvelo multivelo scvi-tools scglue scrublet -y --quiet
```

### 必须安装的包-pip
```Python
#--quiet代表静默安装，如果需要看到安装过程去掉即可
#GLUE下游分析涉及，需要安装
!pip install pyscenic --quiet
#Pyomic 我开发的包，需要安装
!pip install Pyomic --quiet
```

### 可选安装的包
```Python
#可选，scATAC-seq的motif估计时需要用
!git clone https://github.com/calico/scBasset.git
!pip install -e scBasset

#可选,pyscenic的多组学套件，整个分析不涉及，所以可以不用安装，节省时间
!git clone https://github.com/aertslab/pycisTopic.git
!git clone https://github.com/aertslab/pycistarget.git
!git clone https://github.com/aertslab/scenicplus.git
!pip install -e pycisTopic
!pip install -e pycistarget
!pip install -e scenicplus
```