# RNA-seq. 上游分析全教程

## 1. 测序数据下载

在本节中，我们将介绍如何从ncbi上下载原始测序数据SRA文件进行分析

### 1.1 SRA文件下载

#### 1.1.1 单SRA文件下载

 如果你只想分析一个SRA文件，并且已经知道了SRA号，那么我们可以在linux交互环境输入以下命令进行下载

```python
#单SRA下载prefetch
SRR10502962
```

#### 1.1.2 多SRA文件下载

如果你需要批量下载SRA文件，那么你需要得到一个带有多个SRA号的txt文件

e.g. SRR_Acc_List.txt文件内容

```python
SRR5113012
SRR5113013
SRR5113014
```

```python
#多SRA下载
prefetch --option-file SRR_Acc_List.txt
```

### 1.2 sra格式转fastq格式

sra格式的文件一般是经过压缩的测序文件，我们需要转换成原始测序数据fastq格式。

#### 1.2.1 单端测序文件转换

```python
#转换当前目录下全部以.sra结尾的文件
fastq-dump *.sra
```

#### 1.2.2 双端测序文件转换

```python
#--split-files参数可以将其分解为两个fastq文件。
fastq-dump --split-files *.sra
```

## 2. 测序数据质控

为了检测我们测序数据的质量，我们常常需要生成一份质控报告进行直观的观察

#### 2.1 检测数据质量

#### 2.1.1 单份报告生成

在本小节，我们使用fastqc生成质控报告，在fastq文件目录下输入

```python
#批量生成所有fastq文件的质控报告
fastqc *.fastq
```

等待运行结束后，在同目录下有着*_fastqc.html和*_fastqc.zip两个文件，我们可以打开对应的html文件查看该fastq数据的质量。报告的解读见文章

#### 2.1.2 多份报告整合

通过上面的步骤，我们得到了每个fastq文件的质控报告，为了整体进行评估，我们使用multiqc整合报告结果

```python
#整合质控报告结果
multiqc *.zip
```

### 2.2 测序数据过滤

由于我们得到的数据可能包括接头序列，引物，poly-A尾巴和其他类型的不需要的序列，为避免影响下面的分析，我们需要去除这些无关的测序数据。

在这里，相较于其他过滤工具而言，我们选择trim-galore

#### 2.2.1 单端测序

```shell
#新建clean文件夹存放测序结果
mkdir clean
#单端测序数据质量过滤
trim_galore -q 20 \ #设定Phred quality score阈值
	--phred33 \ #选择-phred33或者-phred64，表示测序平台使用的Phred qualityscore
	--stringency 3 \ #设定可以忍受的前后adapter重叠的碱基数，默认为1（非常苛刻）。可以适度放宽，因为后一个adapter几乎不可能被测序仪读到
	--length 20 \ #设定输出reads长度阈值，小于设定值会被抛弃。
	-e 0.1 \ #容错率
	-o /home/seq/clean #输出结果
	/home/seq/SRR10502962_1.fastq #待处理文件
```

#### 2.2.2 双端测序

```shell
#新建clean文件夹存放测序结果
mkdir clean
#双端测序数据质量过滤
trim_galore -q 25 \
	--phred33 \
	--stringency 3 \
	--length 36 \
	-e 0.1 \
	--paired #双端测序
	-o /home/seq/clean
	/home/seq/SRR10502962_1.fastq /home/seq/SRR10502962_2.fastq
```

#### 2.2.3 多组测序数据同时过滤

为了同时完成多组测序数据的同时过滤，我们在这里编写.sh脚本来在Linux系统下批量运行

**config**

```shell
mkdir clean
cd clean
ls /home/seq/*_1.fastq >1
ls /home/seq/*_2.fastq >2
paste 1 2  > config
```

**qc.sh**

```shell
bin_trim_galore=trim_galore
dir='/home/seq/clean'
cat $1 |while read id
do
	arr(${id})
	fq1=${arr[0]}
	fq2=${arr[1]}
$bin_trim_galore -q 25 --phred33 --phred33 --length 36 --stringency 3 --paired -o $dir $fq1 $fq2 
done
```

运行qc.sh

```shell
#config是传递进去的参数
bash qc.sh config
```

## 3. 比对到参考基因组

由于测序仪机器读长的限制，在构建文库的过程中首先需要将DNA片段化，测序得到的序列只是基因组上的部分序列。为了确定测序reads在基因组上的位置，需要将reads比对回参考基因组上，这个步骤叫做mapping

### 3.1 参考基因组文件下载

```shell
#参考基因组下载(hg38)(subread)
wget -O  hg38.fa.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#解压
gunzip hg38.fa.gz
```

### 3.2 subread比对

#### 3.2.1 构建索引

subread有着极其快速的比对效率，这与其构建索引的预处理是不可分开的，所以我们用以下代码来构建索引

```shell
cd /home/seq/index/hg38
subread-buildindex -o hg38 hg38
```

#### 3.2.2 比对

```shell
subread-align -t 0 \ #0代表RNA-seq，1代表DNA-seq
	-T5 \ #线程数
	-i hg38  \ #指定参考基因组的basename
	-r /home/seq/SRR10502962_1.fastq \
	-R /home/seq/SRR10502962_2.fastq \
	-o SRR10502962.bam #输出文件
```

## 4. 统计基因counts数

在这里，我们仅介绍一个工具featureCounts。

featuresCounts软件用于统计基因/转录本上mapping的reads数，也就是用于raw count定量。该软件不仅支持基因/转录本的定量，也支持exon，gene bodies，genomic bins，chromsomal locations等区间的定量

### 4.1 下载gtf基因组注释文件

撰写本教程时的最新版本为103，如有更新可以去官网看看再来下载

```shell
#我们将基因组注释文件存放到refer文件夹中
mkdir /home/seq/refer
cd /home/seq/refer
#下载
wget http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.gtf.gz
```

### 4.2 使用featureCounts进行定量分析

```shell
#选择gtf路径
gtf="/home/seq/refer/Homo_sapiens.GRCh38.103.gtf.gz"
#featureCounts 定量分析
featureCounts -T 5 \ #线程数
	-p \ #针对paired-end数据
	-t exon \ #跟-g一样的意思，其是默认将exon作为一个feature
	-g gene_id \ #从注释文件中提取Meta-features信息用于read count，默认是gene_id
	-a $gtf \ #输入GTF/GFF基因组注释文件
	-o all.id.txt \ #输出文件
	*.bam #待处理数据
	
#去除多余信息，矩阵保存为counts.txt
cat all.id.txt | cut-f1,7- > counts.txt

```

到这里，我们就已经得到一个Counts矩阵了，后续的分析被称为转录组下游分析。将在后面的教程中继续介绍

| Geneid          | SRR10502962.bam |
| --------------- | --------------- |
| ENSG00000223972 | 1               |
| ENSG00000227232 | 132             |
| ENSG00000278267 | 3               |
| ENSG00000243485 | 0               |

