#ifndef _PEPROCESS_H
#define _PEPROCESS_H

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
	void* filter_pe_fqs(PEcalOption opt);
	void* sub_thread(PEthreadOpt opt);
	int read(vector<C_fastq>& pe1,vector<C_fastq>& pe2,ifstream& infile1,ifstream& infile2);
	void peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,string type,gzFile out1,gzFile out2);
	//void peRead();
	void* doCal(int index);
	void  preOutput(int type,C_fastq& a);
	void output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile);
	void output_fastqs2(int type,vector<C_fastq> &fq1,ofstream& outfile);
	void peWrite(int num);
	void run_pigz_split(int type);
	void get_line_number(int* line_num);
	void merge_stat();
	void merge_data();
	void C_fastq_init(C_fastq& a,C_fastq& b);
	void process_nonssd();
	void* sub_thread_nonssd(int index);
	int read_gz(vector<C_fastq>& pe1,vector<C_fastq>& pe2);
	void merge_stat_nonssd();
	void peStreaming_stat(C_global_variable& local_gv);
	//void peOutput(outputOption opt);
public:
	C_global_parameter gp;
	C_global_variable gv;
	queue<C_fastq> read_queue1,read_queue2,cal_trim_queue1,cal_trim_queue2,cal_clean_queue1,cal_clean_queue2;
	int read_num;
	int cal_num;
	int write_num;
	int read_end,cal_end;
	string brother_stat,write_stat;
	mutex read_cal_m,cal_write_m;
	long long processed_reads;
	vector<C_fastq> trim_merge1[max_thread],trim_merge2[max_thread],clean_merge1[max_thread],clean_merge2[max_thread];
	gzFile gzfp1,gzfp2;
	ifstream nongzfp1,nongzfp2;
	gzFile gz_trim_out1[max_thread],gz_trim_out2[max_thread];
	gzFile gz_clean_out1[max_thread],gz_clean_out2[max_thread];
	gzFile gz_trim_out1_nonssd,gz_trim_out2_nonssd,gz_clean_out1_nonssd,gz_clean_out2_nonssd;
	ofstream of_log;
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