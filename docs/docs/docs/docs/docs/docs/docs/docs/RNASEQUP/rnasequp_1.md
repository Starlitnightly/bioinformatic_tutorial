# RNA-seq. 环境配置

我们在前面的教程已经完成了Linux的安装，目前我们的linux可以理解成一个空壳，里面大概什么都没有，于是我们需要安装一些RNA-seq上游分析所必要的包。

## 1. 环境要求

### 1.1 Linux环境

对于Linux环境，我们安装的系统为Ubuntu18.04，满足此系统即可。

### 1.2 Window环境

对于Windows操作系统，我们首先安装WSL即可

## 2. 软件安装

### 2.1 miniconda安装

conda是一个开源的软件包管理系统和环境管理系统，可用于安装多个版本的软件包及其依赖关系，并且可以任意切换（安装conda的目的是为了防止软件版本与包之间互相干扰。）（e.g. 有些软件只能在Python3.6上运行，再新的版本会出现bug，这时就需要用conda来解决这个问题）

我们在这里选择安装miniconda，这是一个轻量级的conda框架，相较于Anaconda的臃肿而言，更加轻便

#### 2.1.1 下载miniconda

```python
# 在linux在使用以下命令下载miniconda
wget-c https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

#### 2.1.2 安装miniconda

```python
# 安装刚刚下载的Miniconda，bash就是运行.sh文件的意思
bash Miniconda3-latest-Linux-x86_64.sh
```

#### 2.1.3 激活conda

```python
#将conda命令添加到环境变量中
source .bashrc
```

#### *2.1.4 国内用户选择清华镜像

```python
# 添加镜像
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config-- add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
conda config-- add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/biocondacondaconfig--setshow_channel_urls yes
```

### 2.2 RNA-seq上游依赖包安装

在这一步中，我们将安装转录组上游分析所需要用到的各个包，我将在最后的一小节中简单介绍一下每个包的用途

#### 2.2.1 创建虚拟环境

为了不干扰Linux系统下其他包的运行，我们将创建一个全新的虚拟环境用来管理RNA-seq分析将用到的包，我们将此虚拟环境命名为rna

```python
#创建名为rna的软件安装环境
conda create -n rna python=3
#查看当前conda环境
conda info --envs
#激活conda的rna环境
source activate rna
```

在每一次退出Linux重进后，都不要忘记了输入source activate rna激活环境

#### 2.2.2 安装RNA-seq上游依赖包

```python
#以下一行命令即可安装完成
conda install -y fastp fastqc multiqc subread bedtools cutadapt trim-galore sra-tools
```

#### *2.2.3 RNA-seq上游依赖包简介

| Package     | Description                                                  |
| ----------- | ------------------------------------------------------------ |
| fastp       | fastq文件质控软件，极其智能                                  |
| fastqc      | 高通量测序数据的高级质控工具                                 |
| multiqc     | 对测序数据进行质量评估（将fastqc生成的多个报告整合成一个文件） |
| subread     | 将reads比对到参考基因组上（速度极快）                        |
| bedtools    | 涵盖各种基因组计算所需要的工具                               |
| cutadapt    | 从高通量测序数据中发现并去除衔接子序列，引物，poly-A尾巴和其他类型的不需要的序列 |
| trim-galore | 是对FastQC和Cutadapt的包装。适用于所有高通量测序，包括RRBS(Reduced Representation Bisulfite-Seq ), Illumina、Nextera 和smallRNA测序平台的双端和单端数据 |
| sra-tools   | 来自NCBI的SRA工具包和SDK是工具和库的集合，这些工具和库用于使用INSDC序列读取档案中的数据 |

