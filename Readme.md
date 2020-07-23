## Getting started
#### Requirements
    gcc: 4.7 or higher
    zlib: 1.2.3.5 or higher
    htslib: 1.9 or higher
    pthread library

#### Install
    git clone https://github.com/BGI-flexlab/SOAPnuke.git
    cd SOAPnuke
    make

#### QuickStart

    filter:
    
    SOAPnuke filter -1 test.r1.fq.gz -2 test.r2.fq.gz -C clean_1.fq.gz -D clean_2.fq.gz -o result -T 8
    
    filterHts:
    
    SOAPnuke filterHts --ref chr21.fa -1 input.bam -2 output.cram  -o result SOAPnuke filterHts -1 input.bam -2 output.bam  -o result

    filterStLFR:

    SOAPnuke filterStLFR -1 fq1.list -2 fq2.list -C clean1.gz -D clean2.gz -o result -T 8 -c config

**Note:
We don't recommend using 1.X anymore, since 2.X outweighs it much in performance and preprocessing effect. This README is basically applied to 2.X as well.**
	
## Introduction

As a novel analysis tool developed for quality control and preprocessing of FASTQ and SAM/BAM data, SOAPnuke includes three modules for different usage scenarios: **filter**, **filterHts** and **filterStLFR**. 

All usages start with executable file **SOAPnuke**, and different modules are invoked with different sub-commands like this:

	Program: SOAPnuke
	Version: 2.1.3
	Contact: GongChun<gongchun@genomics.cn>  ChenYuxin<chenyuxin@genomics.cn>
	Command:
         filter        preprocessing sequences
         filterHts     preprocessing HTS sequences
         filterStLFR   preprocessing stLFR sequences

## Parameter
Module **filter** has following paramters for running:
	
	-1, --fq1		FILE		fq1 file(required)
	-2, --fq2		FILE		fq2 file, used when pe
	-C, --cleanFq1		STR		clean fq1 file name(required,gz format)
	-D, --cleanFq2		STR		clean fq2 file name
	-o, --outDir		STR		output directory, directory must exists
	-8, --outFileType	STR		output file format: fastq or fasta[fastq]
	-0, --log		STR		log file

	-f, --adapter1		STR		adapter sequence of fq1 file (5' adapter when filtersRNA mode)
	-r, --adapter2		STR		adapter sequence of fq2 file (for PE reads or 3' adapter when filtersRNA)
	-Z, --contam1		STR		contaminant sequence(s) for fq1 file, split by comma
	-z, --contam2		STR		contaminant sequence(s) for fq2 file, split by comma
	-Y, --ctMatchR		FLOAT/STR	contam's shortest consistent matching ratio [default:0.2]
	-c, --global_contams	STR		global contaminant sequences which need to be detected, split by comma if more than 1
	-d, --glob_cotm_mR	STR		minimum match ratio in global contaminant sequences detection
	-k, --glob_cotm_mM	STR		maximum mismatch number in global contaminant sequences detection
	-5, --seqType		INT		Sequence fq name type, 0->old fastq name, 1->new fastq name [0]
							old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2
							new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC
	-R, --trimFq1		STR		trim fq1 file name(gz format)
	-W, --trimFq2		STR		trim fq2 file name
	-K, --tile		STR		tile number to ignore reads, such as [1101-1104,1205]
	-F, --fov		STR		fov number to ignore reads (only for zebra-platform data), such as [C001R003,C003R004]
	
	Adapter related:
	-J, --ada_trim				trim read when find adapter[default:discard]
	-a, --contam_trim			trim read when find contam[default:discard]
	
	find 5' adapter
	-S, --adaRCtg		INT		mini 5' adapter continuous alignment length (default: 6)
	-s, --adaRAr		FLOAT		mini alignment rate when find 5' adapter: alignment/tag (default: 0.8)
	
	find 3' adapter
	-U, --adaRMa		INT		mini alignment length when find 3' adapter (default: 5)
	-u, --adaREr		FLOAT		Max error rate when find 3' adapter (mismatch/match) (default: 0.4)
	-b, --adaRMm		INT		Max mismatch number when find 3' adapter (default: 4)

	-l, --lowQual		INT		low quality threshold  [default:5]
	-q, --qualRate		FLOAT		low quality rate  [default:0.5]
	-n, --nRate		FLOAT		N rate threshold  [default:0.05]
	-m, --mean		FLOAT		filter reads with low average quality
	-p, --highA		FLOAT		filter reads if ratio of A in a read exceed [FLOAT]
	-g, --polyG_tail	FLOAT		filter reads if found polyG in tail [INT]
	-X, --polyX		INT		filter reads if a read contains polyX [INT]
	-i, --index				remove index
	-L, --totalReadsNum	INT/FLOAT	number/fraction of reads you want to keep in the output clean fq file(cannot be assigned when -w is given).
						It will extract reads randomly through the total clean fq file by default, you also can get the head reads
						for save time by add head suffix to the integer(e.g. -L 10000000head)
	-t, --trim		INT,INT,INT,INT	trim some bp of the read's head and tail, they means: (PE type:read1's head and tail and read2's head and tail  [0,0,0,0]; SE type:read head and tail [0,0])
	-x, --trimBadHead	INT,INT		Trim from head ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]
	-y, --trimBadTail	INT,INT		Trim from tail ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]

	-O, --overlap		INT		filter the small insert size.Not filter until the value exceed 1[-1]
	-P, --mis		FLOAT		the maximum mismatch ratio when find overlap between PE reads(depend on -O)[0.1]

	-e, --patch		INT		reads number of a patch processed[400000]
	-T, --thread		INT		threads number used in process[6]

	-Q, --qualSys		INT		quality system 1:64, 2:33[default:2]
	-G, --outQualSys	INT		out quality system 1:64, 2:33[default:2]
	-3, --maxReadLen	INT		read max length,default 49 for filtersRNA
	-4, --minReadLen	INT		read min length,default 18 for filtersRNA,30 for other modules
	-w, --output_clean	INT		max reads number in each output clean fastq file

	-7, --pe_info				Add /1, /2 at the end of fastq name.[default:not add]
	-B, --baseConvert	STR		convert base when write data,example:TtoU

	-h, --help				help
	-v, --version				show version


## Availability

SOAPnuke is released under [GPLv3][1]. The latest source code is [freely
available at github][2]. 

## Frequently asked questions (FAQs)

1. [What types of data format does SOAPnuke work with?](#dataf)
2. [What types of data does SOAPnuke work with?](#data)
3. [Can I choose a part of quality-control statistics instead of all?](#qcn)

#### <a name="dataf"></a>1. What types of data format does SOAPnuke work with?

SOAPnuke only works with sequence data of FASTQ format temporarily, but we are considering supporting SAM/BAM to provide operations on aligned datasets.

#### <a name="data"></a>2. What types of data does SOAPnuke work with?

**Filter** module are designed to deal with sequences from most experiments. Compared to it, the other three modules with edited parameter-sets and tailored functions deal with specific data. These designs have refered to existing pipelines of data preprocessing.

 
#### <a name="qcn"></a>3. Can I choose a part of quality-control statistics instead of all?

We will realize this feature in the following versions, but all those statistics will be computed for now.
 




[1]: http://en.wikipedia.org/wiki/GNU_General_Public_License
[2]: https://github.com/BGI-flexlab/SOAPnuke
