# 分析3：单细胞样本对齐

前面提到，我们在单细胞多组学中，会得到大量的细胞，每个细胞都有不同的状态。对于mofa而言，2018年的nature使用了scnmt技术，该技术可以同时测定一个细胞的ATAC与RNA情况，但在大部分情况下，我们只能从同一个样本中分别进行ATAC与RNA的测定，而不是对同一个细胞进行测定，但是从理论上来说，应该是有一些细胞是相像的，于是，如何把这些细胞对齐，成为了一个需要解决的问题。