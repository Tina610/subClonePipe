# subClonePipe
通过突变频率预测亚克隆变化

## 使用方法

python getSubClonePipe.py --config configureFile

[configureFile](https://github.com/ChenYuelong/subClonePipe/blob/master/example_sample_config.ini)


---
title: 亚克隆检测算法（程序）研发
notebook: Acorndx_BloodGroup
tags:
---

## 方法

- pyclone: 根据vaf聚类cluster
- ClonEvol:评估clusters之间的从属关系+画图
- fishplot:画图


## 安装

### pyclone

**注意：由于在集群上安装，特别是提交到计算节点进行计算的时候，可能由于matplotlib的设置导致图像画不出来（`RuntimeError: Invalid DISPLAY variable`），需要修改下列文件(如果不用画图功能可以忽略)**

> pyclone/post_process/plot/loci.py
>
> pyclone/post_process/plot/clusters.py
>

添加下列语句到`import`最开始位置，两个文件都需要修改

```
import matplotlib
matplotlib.use('Agg')
```


PyClone is a Python package so it will require the Python language interpreter to be installed. PyClone is designed to work with **Python 2.7.x** series, not Python3.

If you are having a difficult time installing Python or the dependencies I suggest looking at the MiniConda distribution. Thanks to Christof Winter for pointing this out.

The following packages are required

- PyDP >= 0.2.3

- PyYAML >= 3.10

- matplotlib >= 1.2.0 - Required for plotting.

- numpy >= 1.6.2 - Required for plotting and clustering.

- pandas >= 0.11 - Required for multi sample plotting.

- scipy >= 0.11 - Required for plotting and clustering.

- seaborn >= 0.6.0

### CloneEvol

```
install.packages('devtools')
library(devtools)
install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')
```

### fishplot

```
install.packages("devtools")
library(devtools)
install_github("chrisamiller/fishplot")
```

## 使用

### pyclone

1. pyclone第一步


> --in_files 输入文件
>
> --working_dir 分析路径（中间文件及结果生成于何路径）
>
> --samples 输入文件标记


```
PyClone setup_analysis \
	--in_files [1.tsv] [2.tsv] \
	--working_dir analysis/sample1/ \
	--samples [第一个样本标记] [第二个样本标记] \
	--prior total_copy_number

```

> 输出的结果会在working_dir下的trace中，均为bz2文件

2. pyclone第二步

```
PyClone run_analysis --config_file config.yaml
```
3. pyclone第三步

```
PyClone build_table \
	--config_file analysis/sample1/config.yaml \
	--out_file analysis/sample1/table.old_style \
	--table_type old_style
```





### CloneEvol

```
y = infer.clonal.models(variants = x,
	cluster.col.name = 'cluster',
	vaf.col.names = vaf.col.names,
	sample.groups = sample.groups,
	cancer.initiation.model='monoclonal',
	subclonal.test = 'bootstrap',
	subclonal.test.model = 'non-parametric',
	num.boots = 1000,
	founding.cluster = 1,
	cluster.center = 'mean',
	ignore.clusters = NULL,
	clone.colors = clone.colors,
	min.cluster.vaf = 0.01,
	# min probability that CCF(clone) is non-negative
	sum.p = 0.05,
	# alpha level in confidence interval estimate for CCF(clone)
	alpha = 0.05)
```


### fishplot

## 步骤 [github](https://github.com/ChenYuelong/subClonePipe)

- **第一步** ~~根据filter之后的SNP&INDEL结果获得tsv结果。**sd2tsv.pl**：~~ 根据原始multianno.txt，获取tsv结果。**anno2tsv.py**：

> ~~只选取了exonic和splicing的位点~~
>
> 选取所有的点，但是对深度有一个要求，>500X

```
usage: anno2tsv.py [-h] --input_files INFILES [INFILES ...] --samples SAMPLES
                   [SAMPLES ...] --outdir OUTDIR [--driver DRIVER]
                   [--filter FILTER]

转化multianno.txt结果到tsv

optional arguments:
  -h, --help            show this help message and exit
  --input_files INFILES [INFILES ...], -infiles INFILES [INFILES ...]
                        输入文件，一个样本snp和indel结果以"，"分割，不同样本的以" "（空格）分割
  --samples SAMPLES [SAMPLES ...], -s SAMPLES [SAMPLES ...]
                        输入文件标记，空格分割
  --outdir OUTDIR, -o OUTDIR
                        输出目录
  --driver DRIVER, -driver DRIVER
                        driver gene file
  --filter FILTER, -f FILTER
                        需要过滤位点，原始vcf格式
```



~~原始的~~

```

perl sd2tsv.pl [snp.filter] [indel.filter] [driver] [output]

```

- **第二步** 利用pyclone进行位点聚类 **(附加了画图，用来检查无法进行正确推断结果的原因，协助调整)**

```
PyClone setup_analysis --in_files  [s1.tsv] [s2.tsv] --working_dir [example] --samples [s1] [s2] --prior total_copy_number

PyClone run_analysis --config_file [example/config.yaml] --seed 1

PyClone build_table --config_file [example/config.yaml] --out_file [example/table.old_style] --table_type old_style  --max_clusters 6
```





- **第三步** 整理pyclone结果（标记driver）

```
perl getDriver.pl [table.old_style] [driver.txt] [cloneEvaInput.txt]
```

**driver.txt就是我们报出位点的结果，在Leu_interpretation目录下，`HB15xxxxxxx-5-y.txt`就可以**


- **第四步&&第五步** 推断类别从属关系->画图

```
Rscript clonevol.R [cloneEvaInput.txt] [outdir] [sum.p] [alpha]
```

## Pipeline

```
python3 getSubClonePipe.py -config [config.ini] > example.sh
```

### config.ini 说明

文件:<en-media type="application/octet-stream" hash="70b1ef9490a9809073ea988fb7ea6455"/>

利用的python的parserconfig进行的config读取，所以具体规则可参考.

```
[sample] #必须
name=example #必须，样本名称
state=s1,s2 #必须，阶段，可以任意命名（符合python变量规则），例如可以 gangshengbing,shengbinghenjiule。需要逗号分隔

[s1] #必须，上面有多少个阶段，就需要有多少个这个标记，标记名与阶段相同，如果阶段为"gangshengbing,shengbinghenjiule"，那标记应为[gangshengbing]
snp=/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00157-1-61-SNP.txt #必须
indel=/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00157-1-61-INDEL.txt #必须
driver=/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-5-Y.txt #必须
[s2]
snp=/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-1-176-SNP.txt
indel=/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-1-176-INDEL.txt
driver=/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-5-Y.txt
[s3] #如果只有两个阶段，后续的不会读取，所以可保留，也可删除
snp=/home/chenyl/PROJECTS/subClone/rawData/mutation/
indel=/home/chenyl/PROJECTS/subClone/rawData/mutation/
driver=/home/chenyl/PROJECTS/subClone/rawData/mutation/
[s4]
snp=/home/chenyl/PROJECTS/subClone/rawData/mutation/
indel=/home/chenyl/PROJECTS/subClone/rawData/mutation/
driver=/home/chenyl/PROJECTS/subClone/rawData/mutation/
[workdir] #必须
outdir=/home/chenyl/PROJECTS/subClone/analysis/zongxianhua #必须，输出目录
```

### .sh 说明

文件：<en-media type="application/x-sh" hash="bc0d438e465f9171977bc28c908ad800"/>

以上述config.ini生成运行脚本如下，直接运行即可。


```
/usr/bin/python3 /home/chenyl/PROJECTS/subClone/src/anno2tsv.py -infiles   /home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00157-1-61-SNP.txt,/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00157-1-61-INDEL.txt /home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-1-176-SNP.txt,/home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-1-176-INDEL.txt -s  s1 s2 --filter /home/chenyl/PROJECTS/subClone/src/PON.filter -o /home/chenyl/PROJECTS/subClone/analysis/zongxianhua
/home/chenyl/bin/PyClone setup_analysis --in_files  /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/s1.tsv /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/s2.tsv --working_dir /home/chenyl/PROJECTS/subClone/analysis/zongxianhua --samples s1 s2 --prior total_copy_number &&/home/chenyl/bin/PyClone run_analysis --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --seed 1&& /home/chenyl/bin/PyClone build_table --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --out_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/table.old_style --table_type old_style  --max_clusters 6
/usr/bin/perl /home/chenyl/PROJECTS/subClone/src/getDriver.pl /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/table.old_style /home/chenyl/PROJECTS/subClone/rawData/mutation/HB15CK00151-5-Y.txt /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/cloneEvaInput.txt
echo "start plot" && /home/chenyl/bin/PyClone plot_loci --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_loci_density.pdf --plot_type density && /home/chenyl/bin/PyClone plot_loci --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_loci_parallel_coordinates.pdf --plot_type parallel_coordinates && /home/chenyl/bin/PyClone plot_loci --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_loci_scatter.pdf --plot_type scatter && /home/chenyl/bin/PyClone plot_loci --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_loci_similarity_matrix.pdf --plot_type similarity_matrix && /home/chenyl/bin/PyClone plot_loci --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_loci_vaf_parallel_coordinates.pdf --plot_type vaf_parallel_coordinates && /home/chenyl/bin/PyClone plot_loci --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_loci_vaf_scatter.pdf --plot_type vaf_scatter && /home/chenyl/bin/PyClone plot_clusters --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_clusters_density.pdf --plot_type density && /home/chenyl/bin/PyClone plot_clusters --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_clusters_parallel_coordinates.pdf --plot_type parallel_coordinates && /home/chenyl/bin/PyClone plot_clusters --config_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/config.yaml --plot_file /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/plot_clusters_scatter.pdf --plot_type scatter
/usr/bin/Rscript /home/chenyl/PROJECTS/subClone/src/clonevol.R /home/chenyl/PROJECTS/subClone/analysis/zongxianhua/cloneEvaInput.txt /home/chenyl/PROJECTS/subClone/analysis/zongxianhua 0.01 0.01
```


## 结果

**example**

位点结果：<en-media type="text/plain" hash="f249019cc8d36ec7ef68a3b07638ea1f"/>

> 第一列Driver表明这个位点是否为报出位点
>
> id列表示报出位点具体是什么
>
> cluster列表明亚克隆聚类类别

图结果：

模式图：
<en-media type="application/pdf" hash="dca939c34c2dd8488c52b9dfc119a4a6"/>

树图：
<en-media type="application/pdf" hash="1e30bc5c9727cb9a6529d7df4d7ef608"/>

fish图：
<en-media type="application/pdf" hash="bb2682efba228a32a07d6740ca49b22f"/>




## 待改进的问题（手动调整）

1. 并不是每一个样本都可以进行亚克隆分析，pyclone这步不会有什么问题，但是对于亚克隆关系推断的时候，可能会出现无一致关系的可能，通过查看以plot_*开头的pdf文件，查看具体产生不一致的位点，人为去除即可。
2. 对于低频不敏感，而且对于并不是一直存在的位点，可能计算上会有遗漏，具体原因待查。
	1. 一个样本四个阶段，第一个阶段不存在某个位点的突变，可能这个位点就会被忽略，所以这一部分的预处理需要有修改。
	2. 对于这种问题，如果一个样本很重要的一个突变并没有一致性报出，可能会被过滤，这部分需要人工加入。
3. 如果一个样本的多个时期，可能产生不同的亚克隆组分，将导致结果无法得出，fish图无法获得。只有不同的样本能够检测出相同的亚克隆，流程才能顺利进行
