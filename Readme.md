##Getting started

	# Make sure libraries like boost, zlib, log4cplus, openssl have been installed.
	git clone https://github.com/BGI-flexlab/SOAPnuke.git
	cd SOAPnuke/src
	# Make sure those libraries required are accessible in ENV{LD_LIBRARY_PATH}
	cmake .
	make
	./SOAPnuke
	
##Introduction

SOAPnuke is a novel analysis tool developed for ultrafast quality control and preprocessing of high throughput sequencing (HTS) data. It consists of four modules: Filter, Filter-DGE, Filter-sRNA and Filter-meta. The first module is designed for general usage, while the rest three with suitable parameters or tailored functions are for datasets from specific experiments as their names suggest. All of these modules combines quality control and preprocessing to speed up the report on statistics graphs of raw datasets, preprocessed datasets and preprocessing status. Moreover, Hadoop MapReduce framework is introduced to provide great acceleration for parallel large datasets performing with good scalability, while it also supports standalone running.

All usages start with executable file **SOAPnuke**, and different modules are invoked with different sub-commands: **filter** for Filter,
**filterdge** for Filter-DGE, **filtersrna** for Filter-sRNA and **filtermeta** for Filter-meta module.

##Availability

SOAPnuke is released under [GPLv3][1]. The latest source code is [freely
available at github][2]. Before compiling our program, please make sure that you have installed these following libraries: boost, zlib, log4cplus and openssl. After you acquire the source code, modify the **CmakeLists.txt** according to locations of libraries mentioned above. Then, you can use `cmake .` and `make` to compile.

##Frequently asked questions (FAQs)

1. [What types of data format does SOAPnuke work with?](#dataf)
2. [What types of data does SOAPnuke work with?](#data)
3. [Can I choose a part of quality-control statistics instead of all?](#qcn)

####<a name="dataf"></a>1. What types of data format does SOAPnuke work with?

SOAPnuke only works with sequence data of FASTQ format temporarily, but we are considering supporting SAM/BAM to provide operations on aligned datasets.

####<a name="data"></a>2. What types of data does SOAPnuke work with?

**Filter** module are designed to deal with sequences from most experiments. Compared to it, the other three modules with edited parameter-sets and tailored functions deal with specific data. These designs have refered to existing pipelines of data preprocessing.

**Filter-DGE** module deals with DGE datasets. DGE reads are all single-end so only one input file is available. Each read has two adapters at both ends and a ‘CATG’ segment neighboring targeted sequences, which is 17bp in length, in direction of 5’ ends. We also consider those features when designing parameters.

**Filter-sRNA** module deals with small-RNA reads. In sRNA experiment, mRNA sequences must be filtered out, so a parameter of mRNA_FilePath is available. In addition, functions are reduced with reference to existing pipelines.

**Filter-meta** module deals with metagenomics data. It does not have new functions since meta data do not have unique features compared to other data. However, we still modify its functions and parameters according to existing pipelines.
 
####<a name="qcn"></a>3. Can I choose a part of quality-control statistics instead of all?

We will make it happen in the following versions, but all those statistics will be computed for now.
 




[1]: http://en.wikipedia.org/wiki/GNU_General_Public_License
[2]: https://github.com/BGI-flexlab/SOAPnuke