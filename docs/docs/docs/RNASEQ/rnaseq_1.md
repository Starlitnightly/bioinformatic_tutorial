# RNA-seq分析:  环境配置

对于RNA-seq相关的分析，一共分为python跟R两部分。

我们这里处理的是counts矩阵，count矩阵简单来说，可以理解成横坐标为基因，纵坐标为样本的矩阵，如下表所示

| geneid              | sample1 | sample2 |
| ------------------- | ------- | ------- |
| **ENSG00000223972** | 0       | 0       |
| **ENSG00000227232** | 82      | 63      |
| **ENSG00000278267** | 11      | 2       |

我们的目的，就是挖掘这个矩阵中可能包含的信息，各种生物学意义的东西。

## 1. Python部分

对于python部分，在前面分析时所安装的Anaconda，理论上就足够完成全部分析了，但是我们仍有一些包需要安装，比如基因id转换，看到ENSG可能你也不知道这是一个什么东西吧。

```shell
#安装mygene
pip install mygene
```



## 2. R语言部分

对于R语言部分，实际上仅DESeq2这个包是需要的，所以我们就安装这个包就好了

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

