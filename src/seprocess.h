#ifndef _SEPROCESS_H
#define _SEPROCESS_H

#include <iostream>
#include <string>
#include <vector>
#include <mutex>
#include <queue>
#include <thread>
#include "process_argv.h"
#include "global_parameter.h"
#include "global_variable.h"
#include "sequence.h"

#define max_thread 48
struct SEstatOption
{
	vector<C_fastq>* fq1s;
	C_fastq_file_stat* stat1;
};
struct SEcalOption
{
	C_filter_stat* se_local_fs;
	vector<C_fastq>* fq1s;
	vector<C_fastq> *trim_result1,*clean_result1;
};
struct SEthreadOpt
{
	vector<string> files;
	int index;
};
class seProcess{
public:
	seProcess(C_global_parameter m_gp);
	void q_clear(queue<C_fastq>& a);
	void process();
	void brother();
	void print_stat();
	void update_stat(C_fastq_file_stat& fq1s_stat,C_filter_stat& fs_stat,string type);
	void* stat_se_fqs(SEstatOption opt);
	void* filter_se_fqs(SEcalOption opt);
	void* sub_thread(int index);
	int read(vector<C_fastq>& se1,ifstream& infile1);
	void seWrite(vector<C_fastq>& se1,string type,gzFile out1);
	void seWrite_split(vector<C_fastq>& se1);
	//void peRead();
	void* doCal(int index);
	void  preOutput(int type,C_fastq& a);
	void output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile);
	void output_fastqs2(int type,vector<C_fastq> &fq1,ofstream& outfile);
	void output_split_fastqs(string type,vector<C_fastq> &fq1);
	void seWrite(int num);
	void run_pigz_split(int type);
	void get_line_number(int* line_num);
	void merge_stat();
	void merge_data();
	void C_fastq_init(C_fastq& a);
	void process_nonssd();
	void* sub_thread_nonssd(int index);
	void* sub_thread_nonssd_multiOut(int index);
	int read_gz(vector<C_fastq>& se1);
	void merge_stat_nonssd();
	void seStreaming_stat(C_global_variable& local_gv);
	void remove_tmpDir();
	void make_tmpDir();
	void run_pigz();
	void thread_process_reads(int index,vector<C_fastq> &fq1s);
	void create_thread_outputFile(int index);
	void* sub_thread_nonssd_realMultiThreads(int index);
	void create_thread_read(int index);
	//void peOutput(outputOption opt);
public:
	C_global_parameter gp;
	C_global_variable gv;
	queue<C_fastq> read_queue1,cal_trim_queue1,cal_clean_queue1;
	int read_num;
	int cal_num;
	int write_num;
	int read_end,cal_end;
	string brother_stat,write_stat;
	mutex read_cal_m,cal_write_m;
	long long processed_reads;
	vector<C_fastq> trim_merge1[max_thread],clean_merge1[max_thread];
	gzFile gzfp1;
	ifstream nongzfp1;
	gzFile gz_trim_out1[max_thread];
	gzFile gz_clean_out1[max_thread];
	gzFile gz_trim_out1_nonssd,gz_clean_out1_nonssd;
	ofstream of_log;
	string tmp_dir;
	off_t t_start_pos[max_thread];
	off_t t_end_pos[max_thread];
	char* src1;
	int fq1fd;
	gzFile multi_gzfq1[max_thread];
private:
	vector<C_fastq> fq1s;
	vector<C_fastq> trim_output_fq1;
	vector<C_fastq> clean_output_fq1;
	
	//C_fastq fastq1;
	int used_threads_num;
	string random_num;
};
#endif