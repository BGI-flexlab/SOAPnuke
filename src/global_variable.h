#ifndef _GLOBAL_VARIABLE_H
#define _GLOBAL_VARIABLE_H

#include <string>
#include <set>
#include "global_parameter.h"
#include <cstdint>
// #include <map>
using namespace ::std;
#define READ_MAX_LEN 1000
#define MAX_QUAL 42
#define MIN_QUAL 0
// #define REAL_MAX 50
class C_filter_stat
{
public:
	C_filter_stat()
	{
		in_adapter_list_num = 0;
		include_adapter_seq_num = 0;
		include_contam_seq_num = 0;
		n_ratio_num = 0;
		highA_num = 0;
		polyX_num = 0;
		tile_num = 0;
		fov_num = 0;
		low_qual_base_ratio_num = 0;
		mean_quality_num = 0;
		short_len_num = 0;
		long_len_num = 0;
		over_lapped_num = 0;
		no_3_adapter_num = 0;
		include_adapter_seq_num1 = 0;
		include_adapter_seq_num2 = 0;
		int_insertNull_num = 0;
		include_adapter_seq_num_overlap = 0;
		include_contam_seq_num1 = 0;
		include_contam_seq_num2 = 0;
		include_contam_seq_num_overlap = 0;
		n_ratio_num1 = 0;
		n_ratio_num2 = 0;
		n_ratio_num_overlap = 0;
		highA_num1 = 0;
		highA_num2 = 0;
		highA_num_overlap = 0;
		polyX_num1 = 0;
		polyX_num2 = 0;
		polyX_num_overlap = 0;
		low_qual_base_ratio_num1 = 0;
		low_qual_base_ratio_num2 = 0;
		low_qual_base_ratio_num_overlap = 0;
		mean_quality_num1 = 0;
		mean_quality_num2 = 0;
		mean_quality_num_overlap = 0;
		short_len_num1 = 0;
		short_len_num2 = 0;
		short_len_num_overlap = 0;
		long_len_num1 = 0;
		long_len_num2 = 0;
		long_len_num_overlap = 0;
		include_global_contam_seq_num = include_global_contam_seq_num1 = include_global_contam_seq_num2 = include_global_contam_seq_num_overlap = 0;
		polyG_num = polyG_num1 = polyG_num2 = polyG_num_overlap = 0;
		readsNumWithstLFRbarcode = 0;
		dupReadsNum = 0;
	};
	// uint64_t output_reads_num;
	uint64_t in_adapter_list_num;
	uint64_t include_adapter_seq_num, include_adapter_seq_num1, include_adapter_seq_num2, include_adapter_seq_num_overlap;
	uint64_t include_contam_seq_num, include_contam_seq_num1, include_contam_seq_num2, include_contam_seq_num_overlap;
	uint64_t include_global_contam_seq_num, include_global_contam_seq_num1, include_global_contam_seq_num2, include_global_contam_seq_num_overlap;
	uint64_t n_ratio_num, n_ratio_num1, n_ratio_num2, n_ratio_num_overlap;
	uint64_t highA_num, polyX_num;
	uint64_t highA_num1, highA_num2, highA_num_overlap;
	uint64_t polyX_num1, polyX_num2, polyX_num_overlap;
	uint64_t polyG_num;
	uint64_t polyG_num1, polyG_num2, polyG_num_overlap;
	uint64_t tile_num, fov_num;
	uint64_t low_qual_base_ratio_num, low_qual_base_ratio_num1, low_qual_base_ratio_num2, low_qual_base_ratio_num_overlap;
	uint64_t mean_quality_num, mean_quality_num1, mean_quality_num2, mean_quality_num_overlap;
	uint64_t short_len_num, long_len_num;
	uint64_t short_len_num1, short_len_num2, short_len_num_overlap;
	uint64_t long_len_num1, long_len_num2, long_len_num_overlap;
	uint64_t over_lapped_num;
	uint64_t no_3_adapter_num, int_insertNull_num;
	uint64_t readsNumWithstLFRbarcode;
	set<string> stLFRbarcodeNum;
	uint64_t dupReadsNum;
};
class C_general_stat
{
public:
	C_general_stat();
	uint64_t read_max_length;
	uint64_t read_length;
	uint64_t reads_number;
	uint64_t base_number;
	uint64_t a_number, c_number, g_number, t_number, n_number;
	// unsigned uint64_t a_ratio,c_ratio,g_ratio,t_ratio,n_ratio;
	uint64_t q20_num, q30_num;
	// unsigned uint64_t q20_ratio,q30_ratio;
};

class C_reads_pos_base_stat
{
public:
	C_reads_pos_base_stat();
	uint64_t position_acgt_content[READ_MAX_LEN][5];
	// map<int,map<char,int> > position_acgt_content_ratio;
};

class C_reads_pos_qual_stat
{
public:
	C_reads_pos_qual_stat();
	C_reads_pos_qual_stat(C_global_parameter &gp);
	uint64_t *position_qual[READ_MAX_LEN];
	// map<int,map<int,int> > position_qual_ratio;
};
class C_reads_trim_stat
{
public:
	C_reads_trim_stat();
	uint64_t hlq[READ_MAX_LEN], ht[READ_MAX_LEN];
	uint64_t ta[READ_MAX_LEN], tlq[READ_MAX_LEN], tt[READ_MAX_LEN];
};
class C_fastq_file_stat
{
public:
	C_fastq_file_stat();
	C_fastq_file_stat(C_global_parameter &gp);
	C_general_stat gs;
	C_reads_pos_base_stat bs;
	C_reads_pos_qual_stat qs;
	C_reads_trim_stat ts;
	void clear();
};
class C_global_variable
{
public:
	C_global_variable();
	C_global_variable(C_global_parameter &gp);
	C_filter_stat fs;
	C_fastq_file_stat raw1_stat, raw2_stat, trim1_stat, trim2_stat, clean1_stat, clean2_stat;
};
#endif
