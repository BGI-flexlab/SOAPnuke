#ifndef _GLOBAL_VARIABLE_H
#define _GLOBAL_VARIABLE_H

#include <string>
#include <set>
//#include <map>
using namespace::std;
#define READ_MAX_LEN 500
#define MAX_QUAL 42
#define MIN_QUAL 0
//#define REAL_MAX 50
class C_filter_stat{
public:
	C_filter_stat(){
        in_adapter_list_num=0;
        include_adapter_seq_num=0;
        include_contam_seq_num=0;
        n_ratio_num=0;
        highA_num=0;
        polyX_num=0;
        tile_num=0;
        fov_num=0;
        low_qual_base_ratio_num=0;
        mean_quality_num=0;
        short_len_num=0;
        long_len_num=0;
        over_lapped_num=0;
        no_3_adapter_num=0;
        include_adapter_seq_num1=0;
        include_adapter_seq_num2=0;
        int_insertNull_num=0;
        include_adapter_seq_num_overlap=0;
        include_contam_seq_num1=0;
        include_contam_seq_num2=0;
        include_contam_seq_num_overlap=0;
        n_ratio_num1=0;
        n_ratio_num2=0;
        n_ratio_num_overlap=0;
        highA_num1=0;
        highA_num2=0;
        highA_num_overlap=0;
        polyX_num1=0;
        polyX_num2=0;
        polyX_num_overlap=0;
        low_qual_base_ratio_num1=0;
        low_qual_base_ratio_num2=0;
        low_qual_base_ratio_num_overlap=0;
        mean_quality_num1=0;
        mean_quality_num2=0;
        mean_quality_num_overlap=0;
        short_len_num1=0;
        short_len_num2=0;
        short_len_num_overlap=0;
        long_len_num1=0;
        long_len_num2=0;
        long_len_num_overlap=0;
		include_global_contam_seq_num=include_global_contam_seq_num1=include_global_contam_seq_num2=include_global_contam_seq_num_overlap=0;
		polyG_num=polyG_num1=polyG_num2=polyG_num_overlap=0;
        readsNumWithstLFRbarcode=0;
        dupReadsNum=0;
	};
	//int output_reads_num;
	int in_adapter_list_num;
	int include_adapter_seq_num,include_adapter_seq_num1,include_adapter_seq_num2,include_adapter_seq_num_overlap;
	int include_contam_seq_num,include_contam_seq_num1,include_contam_seq_num2,include_contam_seq_num_overlap;
	int include_global_contam_seq_num,include_global_contam_seq_num1,include_global_contam_seq_num2,include_global_contam_seq_num_overlap;
	int n_ratio_num,n_ratio_num1,n_ratio_num2,n_ratio_num_overlap;
	int highA_num,polyX_num;
	int highA_num1,highA_num2,highA_num_overlap;
	int polyX_num1,polyX_num2,polyX_num_overlap;
	int polyG_num;
	int polyG_num1,polyG_num2,polyG_num_overlap;
	int tile_num,fov_num;
	int low_qual_base_ratio_num,low_qual_base_ratio_num1,low_qual_base_ratio_num2,low_qual_base_ratio_num_overlap;
	int mean_quality_num,mean_quality_num1,mean_quality_num2,mean_quality_num_overlap;
	int short_len_num,long_len_num;
	int short_len_num1,short_len_num2,short_len_num_overlap;
	int long_len_num1,long_len_num2,long_len_num_overlap;
	int over_lapped_num;
	int no_3_adapter_num,int_insertNull_num;
	int readsNumWithstLFRbarcode;
	set<string> stLFRbarcodeNum;
	int dupReadsNum;
};
class C_general_stat{
public:
	C_general_stat();
	int read_max_length;
	int read_length;
	long long reads_number;
	long long base_number;
	long long a_number,c_number,g_number,t_number,n_number;
	//unsigned long long a_ratio,c_ratio,g_ratio,t_ratio,n_ratio;
	long long q20_num,q30_num;
	//unsigned long long q20_ratio,q30_ratio;
};

class C_reads_pos_base_stat{
public:
	C_reads_pos_base_stat();
	long long position_acgt_content[READ_MAX_LEN][5];
	//map<int,map<char,int> > position_acgt_content_ratio;
};

class C_reads_pos_qual_stat{
public:
	C_reads_pos_qual_stat();
	long long position_qual[READ_MAX_LEN][MAX_QUAL];
	//map<int,map<int,int> > position_qual_ratio;
};
class C_reads_trim_stat{
public:
	C_reads_trim_stat();
	long long hlq[READ_MAX_LEN],ht[READ_MAX_LEN];
	long long ta[READ_MAX_LEN],tlq[READ_MAX_LEN],tt[READ_MAX_LEN];
};
class C_fastq_file_stat{
public:
	C_fastq_file_stat();
	C_general_stat gs;
	C_reads_pos_base_stat bs;
	C_reads_pos_qual_stat qs;
	C_reads_trim_stat ts;
	void clear();
};
class C_global_variable{
public:
	C_global_variable();
	C_filter_stat fs;
	C_fastq_file_stat raw1_stat,raw2_stat,trim1_stat,trim2_stat,clean1_stat,clean2_stat;
};
#endif
