#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include "global_parameter.h"
#include "global_variable.h"
#include "zlib.h"
#include "process_argv.h"
#include <string>
#include <vector>

class C_sequence_stat_result{	//basic sequence stat information
public:
	C_sequence_stat_result():seq_len(0),a_num(0),c_num(0),g_num(0),t_num(0),n_num(0),a_ratio(0),c_ratio(0),g_ratio(0),t_ratio(0),n_ratio(0),contig_base(0),gc_ratio(0){};
	int seq_len;
	int a_num,c_num,g_num,t_num,n_num;
	float a_ratio,c_ratio,g_ratio,t_ratio,n_ratio;
	int contig_base;
	float gc_ratio;
};
class C_fastq_stat_result:public C_sequence_stat_result{	//fastq stat information
public:
	C_fastq_stat_result():in_adapter_list(0),include_adapter_seq(-1),include_5_adapter(false),include_3_adapter(-1),include_contam(-1),include_global_contam(-1),qual_len(0),low_qual_base_num(0),low_qual_base_ratio(0),mean_quality(0),q20_num(0),q30_num(0),dup(false){};
	bool in_adapter_list;
	int include_adapter_seq;
	bool include_5_adapter;	//only used in sRNA
	int include_3_adapter;	//only used in sRNA
	int include_contam,include_global_contam;
	string read_tile,read_fov;
	int qual_len;
	int low_qual_base_num;
	float low_qual_base_ratio;
	float mean_quality;
	int q20_num,q30_num;
	bool dup;
};
class C_pe_fastq_stat_result{	//pe fastq stat information
public:
	C_pe_fastq_stat_result():over_lapped(1),containStLFRbarcode(false),stLFRbarcode(""),dup(false){};
	C_fastq_stat_result fastq1_result;
	C_fastq_stat_result fastq2_result;
	bool over_lapped;	//pe fastq whether overlapped
	bool containStLFRbarcode;
	string stLFRbarcode;
	bool dup;
};

class C_sequence{	//sequence
public:
	//void output();
	string sequence;
};

class C_fasta:public C_sequence{	//fasta
public:
	//void output();
	string seq_name;
};

class C_fastq:public C_sequence{	//fastq,include adapter seq and trim information
public:
	string seq_id;
	string qual_seq;
//	string adapter_seq;
//	string adapter_seq2;
//	string contam_seq;
	//string global_contams;
	string adapter_orientation,adapter_orientation2;
	string head_trim_len,tail_trim_len;
	int head_hdcut,head_lqcut,tail_hdcut,tail_lqcut,adacut_pos;
	int contam_pos;
	int global_contam_5pos,global_contam_3pos;
	int raw_length;
	int pairPos;
	string toString();
};
class C_pe_fastq{	//pe fastq
public:
	//void output();
	C_fastq fastq1;
	C_fastq fastq2;
};
class C_single_fastq_filter{	//SE fastq filter
public:
	C_single_fastq_filter(C_fastq& a,C_global_parameter& gp);
	void se_trim(C_global_parameter& gp);
	int se_discard(C_filter_stat* fs,C_global_parameter& gp);
	int sRNA_discard(C_filter_stat* fs,C_global_parameter& gp);
	C_fastq read;
	C_fastq_stat_result read_result;
};
class C_pe_fastq_filter{	//PE fastq filter
public:
	C_pe_fastq_filter(C_fastq& a,C_fastq& b,C_global_parameter& gp);
	void pe_trim(C_global_parameter& gp);
	int pe_discard(C_filter_stat* fs,C_global_parameter& gp);
	int pe_dis(bool a,bool b);
	C_fastq fq1;
	C_fastq fq2;
	C_pe_fastq_stat_result reads_result;
	//added for stLFR module
	void setStLFRbarcode(string bcs){
        reads_result.stLFRbarcode=bcs;
	}
	string getStLFRbarcode(){
	    return reads_result.stLFRbarcode;
	}
};
class C_sRNA_fastq_filter:public C_single_fastq_filter{	//sRNA fastq filter
public:
	C_sRNA_fastq_filter(C_fastq& a,C_global_parameter& gp);
	void sRNA_trim(C_global_parameter& gp);
	int sRNA_discard(C_filter_stat* fs,C_global_parameter& gp);
};



#endif
