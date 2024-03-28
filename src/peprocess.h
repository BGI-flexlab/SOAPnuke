#ifndef _PEPROCESS_H
#define _PEPROCESS_H

#include <iostream>
#include <string>
#include <vector>
#include <mutex>
#include <queue>
#include <thread>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <algorithm>
#include <string.h>
#include "process_argv.h"
#include "global_parameter.h"
#include "global_variable.h"
#include "sequence.h"
#include "BloomFilter.h"
#include "ReverseBloomFilter.h"
#include "rmdup.h"
// #include "Md5.h"
// 0:bf 1:rbf 2:normal
#define RMDUP 2
#define max_thread 72
struct PEstatOption
{
	vector<C_fastq> *fq1s, *fq2s;
	C_fastq_file_stat *stat1, *stat2;
	PEstatOption()
	{
		fq1s = fq2s = NULL;
		stat1 = stat2 = NULL;
	}
};
struct PEcalOption
{
	C_filter_stat *local_fs;
	vector<C_fastq> *fq1s, *fq2s;
	vector<C_fastq> *trim_result1, *trim_result2, *clean_result1, *clean_result2;
	PEcalOption()
	{
		local_fs = NULL;
		fq1s = NULL;
		fq2s = NULL;
		trim_result1 = trim_result2 = clean_result1 = clean_result2 = NULL;
	}
	~PEcalOption()
	{
	}
};
class peProcess
{
public:
	explicit peProcess(C_global_parameter m_gp);
	void process();
	void print_stat();
	void update_stat(C_fastq_file_stat &fq1s_stat, C_fastq_file_stat &fq2s_stat, C_filter_stat &fs_stat, string type);
	void *stat_pe_fqs(PEstatOption opt, string dataType);
	virtual void filter_pe_fqs(PEcalOption *opt);
	virtual void filter_pe_fqs(PEcalOption *opt, int index);
	int read(vector<C_fastq> &pe1, vector<C_fastq> &pe2, ifstream &infile1, ifstream &infile2);
	void peWrite(vector<C_fastq> &pe1, vector<C_fastq> &pe2, gzFile out1, gzFile out2);
	void peWrite(vector<C_fastq> &pe1, vector<C_fastq> &pe2, FILE *out1, FILE *out2);
	void preOutput(int type, C_fastq &a);
	void output_fastqs(string type, vector<C_fastq> &fq1, gzFile outfile);
	void output_fastqs(string type, vector<C_fastq> &fq1, FILE *outfile);
	void merge_stat();
	void C_fastq_init(C_fastq &a, C_fastq &b);

	virtual void *sub_thread(int index);
	virtual void *sub_thread_rmdup_step1(int index);
	void catRmFile(int index, int cycle, string type, bool gzFormat);
	void catRmFile(vector<int> indexes, int cycle, string type, bool gzFormat);
	void *smallFilesProcess();
	void peStreaming_stat(C_global_variable &local_gv);
	void make_tmpDir();
	void remove_tmpDir();
	void create_thread_read(int index);
	void create_thread_smalltrimoutputFile(int index, int cycle);
	void create_thread_smallcleanoutputFile(int index, int cycle);
	void thread_process_reads(int index, int cycle, vector<C_fastq> &fq1s, vector<C_fastq> &fq2s);
	void closeSmallTrimFileHandle(int index);
	void closeSmallCleanFileHandle(int index);
	void run_cmd(string cmd);
	void run_extract_random();
	void *sub_extract(string in, int mo, string out);
	void rmTmpFiles();
	void check_disk_available();
	void extractReadsToFile(int cycle, int thread_index, int reads_number, string position, int &output_index, bool gzFormat);
	void extractReadsToFile(int cycle, int thread_index, int reads_number, string position, bool gzFormat);
	void reArrangeReads(int cycle, bool gzFormat, bool split, int &splitIndex, int &preNum);
	void addCleanList(int tmp_cycle, int index);
	int dupNum;
	mutex checkDup;
	BloomFilter *dupDB;
	ReverseBloomFilter *RdupDB;
	//    set<string> checkDupMap;
	gzFile dupOut1, dupOut2;
	gzFile *dupThreadOut1, *dupThreadOut2;
	mutex logLock;

public:
	C_global_parameter gp;
	C_global_variable gv;
	gzFile *multi_gzfq1, *multi_gzfq2;
	FILE **multi_Nongzfq1, **multi_Nongzfq2;
	gzFile *gz_trim_out1, *gz_trim_out2;
	gzFile *gz_clean_out1, *gz_clean_out2;
	FILE **nongz_trim_out1;
	FILE **nongz_trim_out2;
	FILE **nongz_clean_out1;
	FILE **nongz_clean_out2;
	ofstream of_log;
	off_t *t_start_pos;
	off_t *t_end_pos;
	int fq1fd, fq2fd;
	string tmp_dir;
	int limit_end;

	mutex pe_cal_m;
	mutex read_m, write_m;
	int buffer;
	string new_fq1_path, new_fq2_path;
	C_filter_stat *local_fs;
	C_fastq_file_stat *local_raw_stat1, *local_raw_stat2, *local_trim_stat1, *local_trim_stat2, *local_clean_stat1, *local_clean_stat2;
	string *sticky_reads1, *sticky_reads2;
	gzFile gz_fq1, gz_fq2;
	map<int, int> pe1_out, pe2_out;
	// int queue_buffer;
	int bq_check, pair_check;
	int cur_cat_cycle;																		 // current cycle of cat small files
	vector<string> *readyTrimFiles1, *readyTrimFiles2, *readyCleanFiles1, *readyCleanFiles2; // list of small files in threads
	vector<int> *clean_file_readsNum;														 // reads number of small files in threads
	int *sub_thread_done;
	int end_sub_thread;
	int patch;
	uint64_t *threadReadsNum;
	uint64_t *totalData;
	uint64_t *threadCurReadReadsNumIdx;
	vector<vector<uint64_t *>> threadData;
	vector<vector<size_t>> threadDataNum;
	bool *dupFlag;

private:
	int used_threads_num;
};

#endif
