# 1、RRBS简介

RRBS即reduced representation bisulfite sequencing，其测序原理图如下：

![](https://ws1.sinaimg.cn/large/8a8b2c15gy1fui9ewkokoj20mm0d374w.jpg)

第一步：使用MSPI内切酶处理，其识别序列是CCGG，切完之后，双链序列如下：

CGG--------------------------------C

      C---------------------------------GGC

第二步：自动补齐

CGG--------------------------------CCG

GCC--------------------------------GGC

第三步：bisulfite conversion

![](https://ws1.sinaimg.cn/large/8a8b2c15gy1fui9w38py3j20kk0b3408.jpg)



将双链中没有甲基化的C变为U，甲基化的C保持不变，这样转化之后，原来互补的两条双链就不互补了，这样的两条序列相当于一对duplex

# 2、数据预处理

## 2.1 fastqc质量检测

```shell
fastqc -t 2 -o ./00_fastqc R1.fastq R2.fastq
```

主要看Q30的百分比，以及reads中有无接头

![](https://ws1.sinaimg.cn/large/8a8b2c15gy1fuia0cy0krj20q60eumxy.jpg)

从图中可看出，reads中只有illumina universal adapter，因此在后续的处理中只需处理这一种接头

## 2.2 trim_galore去除接头

```shell
trim_galore --phred33 --illumina --stringency 3 -e 0.1  \
            --gzip --length 35 --rrbs --fastqc -o ./cut_adapter \
            --paired R1.fastq R2.fastq
```

详细参数使用`trim_galore -h`查看，这里特别解释一下--rrbs这个参数，是专门用于用MspI处理的RRBS数据，加上这个参数，在切除接头之后会自动切除2bp的碱基，即上面说的自动补齐的那两个碱基。另外--illumina这个参数是用来切除illumina universal adapter（AGATCGGAAGAGC）的，另外参数--nextera切除nextera adapter（CTGTCTCTTATA），--small_rna切除Illumina Small RNA 3' Adapter（GATCGTCGGACT）

# 3、使用bismark转换基因组序列

参考基因组下载：

https://support.illumina.com/sequencing/sequencing_software/igenome.html

```shell
bismark_genome_preparation --path_to_bowtie dir_bowtie --genomic_composition \
           --verbose dir_ref dir_hg38.fa
```

# 4、序列比对

```shell
bismark dir_ref -1 $R1 -2 $R2 --path_to_bowtie $bowtie \
        -o ./01_align --rg_tag --rg_id $fname --rg_sample $fname \
        --unmapped --prefix $fname --basename $fname --nucleotide_coverage
```

-------------------------------------------------------------------------------------

上述代码耗时较长，改为

```shell
bismark dir_ref -1 $R1 -2 $R2 --path_to_bowtie $bowtie \
		--prefix $fname --basename $fname \
        -o ./01_align --non_directional --temp_dir $tmp
```



比对好的结果为bam格式的，其内容与BWA比对后的bam文件稍微有些不同

![](https://ws1.sinaimg.cn/large/8a8b2c15gy1fujmvp4y19j210t0270sq.jpg)

上面是一条bam文件的记录，共17列，下面分别解释一下各列的含义：

**第1列**：

```shell
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG
#reads编号，_1代表来自原来的R1，后面的ATCACG为index序列
```

**第2列**:

此列为一数值，具体解释如下图

![](https://ws1.sinaimg.cn/large/8a8b2c15gy1fujn0hqjnzj207u05vwea.jpg)

	相同数值对于R1和R2来说代表不同的链，此外，这些数值也是由多个值相加所得，下面这几个标签都是2的n次方，随机挑选其中的几个，他们的和是唯一的：

|      |                                                              |
| ---- | ------------------------------------------------------------ |
| 1    | 代表这个序列采用的是PE双端测序                               |
| 2    | 代表这个序列和参考序列完全匹配，没有插入缺失                 |
| 4    | 代表这个序列没有mapping到参考序列上                          |
| 8    | 代表这条序列的另一端序列没有比对到参考序列上，如果这条是R1，它对应的R2端没有比对上 |
| 16   | 代表这个序列比对到参考序列的负链上                           |
| 32   | 代表这个序列对应的另一端序列比对到参考序列的负链上           |
| 64   | 代表这个序列是R1端序列                                       |
| 128  | 代表这个序列是R2端序列                                       |
| 256  | 代表这个序列不是主要的比对，一条序列可能比对到了多个位置，只有一个是首要的比对位置 |
| 512  | 代表这个序列在QC时失败了，被过滤不掉了                       |
| 1024 | 代表这个序列是PCR重复序列                                    |
| 2048 | 代表这个序列是补充的比对，不常用                             |

**第3-4列**：

```shell
chr4    184651014
# 代表这条reads比对到4号染色体的184651014这个位置
```

**第5列**：

为一数值，Mapping quality，比对的质量分数，越高说明该read比对到参考基因组上的位置越唯一

**第6列**：

CIGAR值，碱基匹配上的碱基数。

|      |                                           |      |      |                                             |
| ---- | ----------------------------------------- | ---- | ---- | ------------------------------------------- |
| M    | match/mismatch                            |      | I    | insert                                      |
| D    | deletion                                  |      | N    | skipped（跳过这段区域）                     |
| S    | soft clipping（被剪切的序列存在于序列中） |      | H    | hard clipping（被剪切的序列不存在于序列中） |
| P    | padding(填充)                             |      |      |                                             |



```shell
135M
#代表135个碱基在比对时完全比对
3S6M1P1I4M
#代表前三个碱基被剪切去除了，然后6个比对上了，然后打开了一个缺口，有一个碱基插入，最后是4个比对上了，是按照顺序的
```

**第7列**：

	MRNM(chr)，mate的reference sequence name，实际上就是mate比对到的染色体号，若是没有mate，则是*，若为=则代表与此条比对到相同的染色体上。

**第8列**：

mate position，mate比对到参考序列上的第一个碱基位置，若无mate,则为0

**第9列**：

insert size

**第10列**：

Sequence，就是read的碱基序列，如果是比对到负链上则对read进行了reverse completed

```shell
GTGGGGGGTTGTGTGGTTAGGATGGAGGAGTAGGTTATTAAGGTTGAGGTTGAGTTGGGTATTTTTTTTAGGGTGTTTAAGAAGTTTTTGTATATGATGAATGTG
```

**第11列**：

ASCII，read质量的ASCII编码

```shell
JJJJJJJJJJJFJJJJFJJJJJJJJJJJAJFJFJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJAJJJJJJJFJJJJJFJJJJJJJJJAJJJFJJJJJJJJJ
```

**第12列**：

NM-tag，与bowtie中NM:i相同。read string转换成reference string需要的最少核苷酸的edits:插入/缺失/替换

```
NM:i:48
```

**第13列**：

XX-tag，对错配的描述，不包括indel(插入和缺失)

```shell
XX:Z:2C4C2C6C1C1AC18C3C15C4C2CC14
#表示2个碱基完全匹配，一个C替换，接着4个碱基完全匹配，一个C发生替换
```

**第14列**：

XM-tag，methylation call string

```shell
XM:Z:........xz.z....hx............x....h.hh....xz....xz...x....h.hh..hhhx....z...h......h
```

![](https://ws1.sinaimg.cn/large/8a8b2c15gy1fujnnga752j209q06daa0.jpg)

CHG、CHH、CG代表甲基化的三种模式

CG：C后面接的是G

CHG：C后面跟了任意一个碱基，再接上一个G

CHH：C后面跟的是除了G以外的任何碱基

**第15列**：

XR-tag，read conversion state for the alignment 。

共两种转换：CT和GA，GA就是指将read里的所有G转换成A

```shell
XR:Z:GA
```

**第16列**：

XG-tag，genome conversion state for the alignment。

共两种：GA和CT。CT是指将全基因组上所有的C转换成T

```shell
XG:Z:CT
```

**第17列**：

表示样本名称，是在做比对时通过--rg_tag自己设置的

```
RG:Z:FF01N_ATCACG_S9_L007
```



# 5、去除重复序列

这个步骤对于RRBS来说不需要，WGBS需要

```shell
deduplicate_bismark --paired --bam bam_file --output_dir dir_bamfile
```

# 6、抽提出甲基化的统计信息

```shell
bismark_methylation_extractor --paired-end --report --output dir --gzip --bedGraph bam_file
```

产生的结果文件如下：

```shell
CHG_OB_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz
CHG_OT_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz
CHH_OB_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz
CHH_OT_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz
CpG_OB_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz
CpG_OT_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz
FF01N_ATCACG_S9_L007_pe.deduplicated.bedGraph.gz
FF01N_ATCACG_S9_L007_pe.deduplicated.bismark.cov.gz
FF01N_ATCACG_S9_L007_pe.deduplicated.M-bias.txt
FF01N_ATCACG_S9_L007_pe.deduplicated_splitting_report.txt
```

现在来看看每个结果文件中的内容：

**CpG_OT_FF01N_ATCACG_S9_L007_pe.deduplicated.txt.gz**

****

前六个文件中的格式是一样的：

```shell
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG      +       chr4    184651148       Z
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG      +       chr4    184651086       Z
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG      -       chr4    184651072       z
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG      +       chr4    184651056       Z
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG      +       chr4    184651018       Z
E00500:282:HN3M3CCXY:7:1101:8613:2628_1:N:0:ATCACG      +       chr4    184651016       Z
E00500:282:HN3M3CCXY:7:1101:11576:2628_1:N:0:ATCACG     -       chr22   46508093        z
E00500:282:HN3M3CCXY:7:1101:11576:2628_1:N:0:ATCACG     -       chr22   46508039        z
E00500:282:HN3M3CCXY:7:1101:13078:2628_1:N:0:ATCACG     +       chr11   8730303 Z
```

一共5列，第一列是read的ID，第二列表示是否甲基化，第三、四列是比对到的染色体上的位置，第五列表达的信息同第二列是一致的，Z表示CpG中的C甲基化，z表示CpG中的C未甲基化

**FF01N_ATCACG_S9_L007_pe.deduplicated.bedGraph.gz**

```shell
chr18   10795   10796   100
chr18   10815   10816   100
chr18   10816   10817   83.3333333333333
chr18   10823   10824   100
chr18   10824   10825   83.3333333333333
chr18   10829   10830   40
chr18   10830   10831   50
chr18   10831   10832   80
chr18   10832   10833   83.3333333333333
```

一共4列，前三列代表位置信息，间隔为1bp，第三列是该位置甲基化的百分比

**FF01N_ATCACG_S9_L007_pe.deduplicated.bismark.cov.gz**

```shell
chr18   10796   10796   100     5       0
chr18   10816   10816   100     5       0
chr18   10817   10817   83.3333333333333        5       1
chr18   10824   10824   100     5       0
chr18   10825   10825   83.3333333333333        5       1
chr18   10830   10830   40      2       3
chr18   10831   10831   50      3       3
chr18   10832   10832   80      4       1
chr18   10833   10833   83.3333333333333        5       1
chr18   10840   10840   50      5       5
```

	一共6列，前三列代表位置，第四列代表该位置的甲基化的百分比，第五列表示该位置上甲基化的coverage，第六列表示该位置上未甲基化的coverage。第4列的值=第5列的值/（第5列+第6列的值）

**FF01N_ATCACG_S9_L007_pe.deduplicated.M-bias.txt**

```shell
CpG context (R1)
================
position        count methylated        count unmethylated      % methylation   coverage
1       1152112 795693  59.15   1947805
2       12334   7034    63.68   19368
3       14795   8508    63.49   23303
4       33622   47811   41.29   81433
5       87472   52611   62.44   140083
6       37551   43574   46.29   81125
7       75471   57671   56.68   133142
```

这个文件展示了，该样本每条reads（原长151bp）每个位置上的甲基化的百分比及覆盖度

**FF01N_ATCACG_S9_L007_pe.deduplicated_splitting_report.txt**

```shell
Processed 3280866 lines in total
Total number of methylation call strings processed: 6561732

Final Cytosine Methylation Report
=================================
Total number of C's analysed:   122221758

Total methylated C's in CpG context:    8297224
Total methylated C's in CHG context:    228899
Total methylated C's in CHH context:    493695

Total C to T conversions in CpG context:        6289349
Total C to T conversions in CHG context:        30742632
Total C to T conversions in CHH context:        76169959

C methylated in CpG context:    56.9%
C methylated in CHG context:    0.7%
C methylated in CHH context:    0.6%
```

利用上述文件可以生成一个html格式的综合性报告

```shell
bismark2report -o S11_graphical.html \
               --alignment_report FF01N_TTAGGC_S11_L007_PE_report.txt \
               --dedup_report FF01N_TTAGGC_S11_L007_pe.deduplication_report.txt \
               --splitting_report FF01N_TTAGGC_S11_L007_pe.deduplicated_splitting_report.txt \
               --mbias_report FF01N_TTAGGC_S11_L007_pe.deduplicated.M-bias.txt \
               --nucleotide_report FF01N_TTAGGC_S11_L007_pe.nucleotide_stats.txt
```



# 7、利用methylkit（R包）进行可视化

## 7.1 安装

```shell
source("https://bioconductor.org/biocLite.R")
biocLite("methylKit")
```

To view documentation for the version of this package installed in your system, start R and enter:

```shell
browseVignettes("methylKit")
```

## 7.2 数据读入

sam/bam格式的都可以，但是要注意读入的文件必须是经过samtools排序后的，如果读入的是bam文件，则还需要是index后的。

# 8、自己写脚本进行后续的分析

主要利用*.bismark.cov.gz这个结果文件。

## 8.1 将每个文件加上title

```shell
#!/bin/bash
for i in *.gz
do
        fname=`echo $i|sed 's/.cov.gz//'`
        less $i|awk '{printf("%s\t%d\t%d\t%.2f\t%d\t%d\n",$1,$2,$3,$4,$5,$6) }'|awk '{print $1,$2,$3,$4,$5,$6,$5+$6}'>$fname.cov.txt
        awk 'BEGIN{OFS="\t"}NR==1{gsub(/.cov.txt$/,"",FILENAME);print "Chr","Start","End",FILENAME"_Methy%",FILENAME"_Met",FILENAME"_Unmet",FILENAME"_Depth"}{NF=NF;print}' $fname.cov.txt>../05_sample_merge/$fname.met.txt
done
```

除了原本的5列之外，还计算第六列，即该位点的深度。

或者也可以使用这个python脚本

[E:\学习笔记\数据分析\RRBS甲基化测序\bismark_results_handle.py]()

## 8.2 合并所有文件

[E:\学习笔记\数据分析\RRBS甲基化测序\rrbs_met_merge.py]()

合并之后文件格式如下：

```
Chr,Start,End,108ZB_Methy,108ZB_Met,108ZB_Unmet,108ZB_Total,109ZB_Methy,109ZB_Met,109ZB_Unmet,109ZB_Total,118ZB_Methy,118ZB_Met,118ZB_Unmet,118ZB_Total,17ZB_Methy,17ZB_Met,17ZB_Unmet,17ZB_Total,21ZB_Methy,21ZB_Met,21ZB_Unmet,21ZB_Total,29ZB_Methy,29ZB_Met,29ZB_Unmet,29ZB_Total,35ZB_Methy,35ZB_Met,35ZB_Unmet,35ZB_Total,37ZB_Methy,37ZB_Met,37ZB_Unmet,37ZB_Total,48ZB_Methy,48ZB_Met,48ZB_Unmet,48ZB_Total,67ZB_Methy,67ZB_Met,67ZB_Unmet,67ZB_Total,82ZB_Methy,82ZB_Met,82ZB_Unmet,82ZB_Total,83ZB_Methy,83ZB_Met,83ZB_Unmet,83ZB_Total,89ZB_Methy,89ZB_Met,89ZB_Unmet,89ZB_Total,92ZB_Methy,92ZB_Met,92ZB_Unmet,92ZB_Total,96ZB_Methy,96ZB_Met,96ZB_Unmet,96ZB_Total
chr1 ,10469,10469,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,100,1,0,1,NA,NA,NA,NA,NA,NA,NA,NA,0,0,1,1,NA,NA,NA,NA,NA,NA,NA,NA
chr1 ,10470,10470,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0,1,1,NA,NA,NA,NA,100,1,0,1,100,1,0,1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA
chr1 ,10471,10471,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,0,0,1,1,NA,NA,NA,NA,NA,NA,NA,NA,0,0,1,1,NA,NA,NA,NA,NA,NA,NA,NA
```

如果样本过多，则合并的时候可能会出现内存错误，此时可以试着先将每个样本进行初步的过滤，例如，以下条件可适当利用：

①每个位点的深度至少为5

[E:\学习笔记\数据分析\RRBS甲基化测序\prefilter_by_depth.py]()

②先将每一对N-T样本合并，如果某一个位点在N-T中都测到，则保留

之后再将所有样本合并

## 8.3 过滤(找差异甲基化位点)

接下来则对合并的样本进行过滤，每一个位点都要满足以下条件：

①平均测序深度至少为10

②该位点无论是在正常还是肿瘤中都至少被4个样本测到

按上述两个条件过滤之后，对每个位点做wilcoxon test，并筛选p≤0.05的位点，及差异甲基化位点（DMP）

[E:\学习笔记\数据分析\RRBS甲基化测序\wilcoxon_test.py]()





