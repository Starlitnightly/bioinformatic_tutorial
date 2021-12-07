# MOFA下游分析：R语言环节

## 1. MOFA数据处理

### 1.1 数据导入

```R
#load library
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
#install.packages('psych')
#install.packages('ggpubr')

#load data
sample_metadata <- fread('mofa_meta.csv')#注意第一列细胞的标头要写为sample
model <- load_model("mofa_factor.hdf5")
samples_metadata(model) <- sample_metadata
```

### 1.2 Var数据导出

```R
#Variance decomposition by Factor
plot_variance_explained(model, max_r2=15)
variance_explained=model@cache[["variance_explained"]][["r2_per_factor"]][["group0"]]
write.table(variance_explained,file='variance_explained.csv',sep=',',row.names =FALSE)

```

### 1.3 Cor数据导出

```R
#Association analysis
p=correlate_factors_with_covariates(model, covariates = c("cell_type","Age","Sex"), plot="log_pval",return_data=TRUE)
write.table(p,file='correlate_factors.csv',sep=',')
```