#ifndef _GLOBAL_PARAMETER_H
#define _GLOBAL_PARAMETER_H
#include <string>
#include <functional>
using namespace::std;
class C_global_parameter{
//common parameter
public:
	/*
	PROCESS_THREAD_NUM(2), IS_STREAMING(false), filterTile_(false), tileIsFov_(false), misMatch_(1), matchRatio_(0.5), lowQual_(5),
                                         qualRate_(0.5), nRate_(0.05), highA_(0), polyX_(0), minMean_(0.0), filterIndex_(false),
                                         rmdup_(false), dupRateOnly_(false), cutReadNum_(0), headTrim_(0), tailTrim_(0), headTrim2_(0), tailTrim2_(0),
                                         memLimit_(700 * MEM_UNIT), qualSys_(ILLUMINA_), isFilterSmallInsertSize_(false), overlap_(10),
                                         mis_(0.1), readLen_(0), readLen2_(0), outDir_("."), onlyStat_(false), isPE_(true),minReadLength(50),cutAdaptor(false),cutBasesNumber(0),
                                         isAdptList_(true), isFull_(false), size_(0), cleanQualSys_(ILLUMINA_), filterAdapter_(true),seqType_(0),outType_(0),highAType_(0)
	*/
	C_global_parameter():is_streaming(false),seq_type("0"),index_remove(false),qualityPhred(33),outputQualityPhred(33),adapter_discard_or_trim("discard"),contam_discard_or_trim("discard"),adapter_method("hd"),whether_add_pe_info(false),output_file_type("fastq"),lowQual(5),lowQualityBaseRatio(0.5),meanQuality(-1),n_ratio(0.05),highA_ratio(-1),polyG_tail(-1),polyX_num(-1),overlap_length(-1),peMismatchRatio(0.1),max_read_length(-1),min_read_length(30),output_clean(0),have_output1(0),have_output2(0),adaMis(2),adaMR(0.5),ctMatchR("0.2"),adaEdge(6),adaRCtg(6),adaRAr(0.8),adaRMa(5),adaREr(0.4),adaRMm(4),threads_num(6),patchSize(0),split_line(10000000),mode(""),total_reads_num(0),f_total_reads_ratio(0),l_total_reads_num(0),total_reads_num_random(true),clean_file_reads(0){};
	C_global_parameter(int argc,char* argv[]);

	string mode;
	string module_name;
	//input and output file
	string fq1_path,fq2_path;
	string trim_fq1,trim_fq2,clean_fq1,clean_fq2;
	string output_dir;
	string log;

	bool is_streaming;
	//input and output file type
	string seq_type;	//Sequence fq name type, 0->old fastq name, 1->new fastq name [0]\n";
       // old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
        // new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
	string output_file_type;	//fasta or fastq
	//adapter
	string adapter_discard_or_trim;
	string adapter_method;
	string adapter1_seq,adapter2_seq;
	string contam_discard_or_trim;
	string contam1_seq,contam2_seq;
	string ctMatchR;
	string global_contams,g_mrs,g_mms;
	//tile to removed
	string tile;
	string fov;
	//index remove
	bool index_remove;
	//base quality
	int qualityPhred,outputQualityPhred;
	int lowQual;
	float lowQualityBaseRatio;
	int meanQuality;
	string trimBadHead,trimBadTail;
	//base content
	float n_ratio,highA_ratio,polyG_tail;
	int polyX_num;
	string trim;	//pe trim and se trim
	string base_convert;
	//PE reads
	//string whether_remove_overlap_peReads,overlap_length,peMismatchRatio;
	int overlap_length;
	float peMismatchRatio;
	bool whether_add_pe_info;
	//computer resource
	int threads_num;
	int patchSize;
	int split_line;
	//whether only stat
	//string whether_only_stat;
	//read length limit
	int max_read_length,min_read_length;
	//reads number limit
	float total_reads_num;
	float f_total_reads_ratio;
	unsigned long long l_total_reads_num;
	bool total_reads_num_random;
	//output reads number,the only variable
	//int output_reads_num;
	unsigned long long output_clean;
	unsigned long long have_output1,have_output2;
	unsigned long long clean_file_reads;
//special parameter
	//filtersRNA module
		//adapter find method2 parameter
	int adaRCtg;
	float adaRAr;
	int adaRMa;
	float adaREr;
	int adaRMm;
	//filter module
		//adapter find method1 parameter
	int adaMis;
	float adaMR;
	int adaEdge;
	//filterMeta module
		//adapter find method is same as filter module
};

#endif
