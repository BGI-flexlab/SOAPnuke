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
#include "process_argv.h"
#include "global_parameter.h"
#include "global_variable.h"
#include "sequence.h"

#define max_thread 48
struct C_reads_trim_stat_2
{
	C_reads_trim_stat stat1;
	C_reads_trim_stat stat2;
};
struct PEstatOption
{
	vector<C_fastq>* fq1s,*fq2s;
	C_fastq_file_stat* stat1,*stat2;
};
struct PEcalOption
{
	C_filter_stat* local_fs;
	vector<C_fastq>* fq1s,*fq2s;
	vector<C_fastq> *trim_result1,*trim_result2,*clean_result1,*clean_result2;
};
struct PEthreadOpt
{
	vector<string> files;
	int index;
};
class peProcess{
public:
	peProcess(C_global_parameter m_gp);
	void q_clear(queue<C_fastq>& a);
	void process();
	void brother();
	void print_stat();
	void update_stat(C_fastq_file_stat& fq1s_stat,C_fastq_file_stat& fq2s_stat,C_filter_stat& fs_stat,string type);
	void* stat_pe_fqs(PEstatOption opt);
	void filter_pe_fqs(PEcalOption opt);
	void* sub_thread(int index);
	int read(vector<C_fastq>& pe1,vector<C_fastq>& pe2,ifstream& infile1,ifstream& infile2);
	void peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,string type,gzFile out1,gzFile out2);
	void peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,string type,FILE* out1,FILE* out2);
	//void add_raw_trim(C_fastq_file_stat& a,C_fastq_file_stat& a2,C_reads_trim_stat& b,C_reads_trim_stat& b2);
	void peWrite_split(vector<C_fastq>& pe1,vector<C_fastq>& pe2);
	//void peRead();
	void* doCal(int index);
	void  preOutput(int type,C_fastq& a);
	void output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile);
	void output_fastqs(string type,vector<C_fastq> &fq1,FILE* outfile);
	void output_fastqs2(int type,vector<C_fastq> &fq1,ofstream& outfile);
	void output_split_fastqs(string type,vector<C_fastq> &fq1);
	//void peWrite(int num);
	void run_pigz(int type);
	void get_line_number(int* line_num);
	void merge_stat();
	void merge_data();
	void merge_trim_data();
	void merge_clean_data();
	void merge_clean_data(int index);
	void C_fastq_init(C_fastq& a,C_fastq& b);
	void process_nonssd();
	void* sub_thread_nonssd(int index);
	void* sub_thread_nonssd_multiOut(int index);
	void* sub_thread_nonssd_realMultiThreads(int index);
	int read_gz(vector<C_fastq>& pe1,vector<C_fastq>& pe2);
	void merge_stat_nonssd();
	void peStreaming_stat(C_global_variable& local_gv);
	void* gzread_(int index,gzFile a);
	void make_tmpDir();
	void remove_tmpDir();
	void* sub_thread_nonssd_threadBuffer(int index);
	void* monitor_read_thread();
	void* monitor_write_thread();
	void create_thread_read(int index);
	void create_thread_outputFile(int index);
	void create_thread_trimoutputFile(int index);
	void create_thread_cleanoutputFile(int index);
	void thread_process_reads(int index,vector<C_fastq> &fq1s,vector<C_fastq> &fq2s);
	void run_cmd(string cmd);
	void run_extract_random();
	void process_some_reads(int index,int out_number);
	void merge_stat(int index);
	void merge_data(int index);
	void limit_process_reads(int index,vector<C_fastq> &fq1s,vector<C_fastq> &fq2s,gzFile gzfq1,gzFile gzfq2);
	//void peOutput(outputOption opt);
public:
	C_global_parameter gp;
	C_global_variable gv;
	//queue<C_fastq> read_queue1,read_queue2,cal_trim_queue1,cal_trim_queue2,cal_clean_queue1,cal_clean_queue2;
	int read_num;
	int cal_num;
	int write_num;
	int read_end,cal_end;
	string brother_stat,write_stat;
	mutex read_cal_m,cal_write_m;
	long long processed_reads;
	//vector<C_fastq> trim_merge1[max_thread],trim_merge2[max_thread],clean_merge1[max_thread],clean_merge2[max_thread];
	gzFile gzfp1,gzfp2;
	gzFile multi_gzfq1[max_thread],multi_gzfq2[max_thread];
	ifstream nongzfp1,nongzfp2;
	gzFile gz_trim_out1[max_thread],gz_trim_out2[max_thread];
	gzFile gz_clean_out1[max_thread],gz_clean_out2[max_thread];
	FILE* nongz_clean_out1[max_thread];
	FILE* nongz_clean_out2[max_thread];
	gzFile gz_trim_out1_nonssd,gz_trim_out2_nonssd,gz_clean_out1_nonssd,gz_clean_out2_nonssd;
	ofstream of_log;
	off_t t_start_pos[max_thread];
	off_t t_end_pos[max_thread];
	char* src1,*src2;
	int fq1fd,fq2fd;
	string tmp_dir;
	int limit_end;
	mutex thread_read_m[max_thread],thread_write_m[max_thread];
	bool file_end;
private:
	vector<C_fastq> fq1s,fq2s;
	vector<C_fastq> trim_output_fq1,trim_output_fq2;
	vector<C_fastq> clean_output_fq1,clean_output_fq2;
	
	//C_fastq fastq1,fastq2;
	int used_threads_num;
	string random_num;
};
/*
class peProcess{
public:
	peProcess(C_global_parameter gp,C_global_variable global_gv);
	void process();

public:
	C_global_parameter gp;
	C_global_variable gv;
};
*/
/*
vector<string> get_pe_hard_trim(string a);
vector<string> get_se_hard_trim(string a);
void peProcess(C_global_parameter gp,C_global_variable& gv);
void peProcess2(C_global_parameter gp,C_global_variable& gv);
*/
void output_fastqs(int type,vector<C_fastq> fq1,gzFile outfile,C_global_parameter gp);

#endif