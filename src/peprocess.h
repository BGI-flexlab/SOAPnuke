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
struct PEstatOption
{
	vector<C_fastq>* fq1s,*fq2s;
	C_fastq_file_stat* stat1,*stat2;
    PEstatOption(){
        fq1s=fq2s=NULL;
        stat1=stat2=NULL;
    }
};
struct PEcalOption
{
	C_filter_stat* local_fs;
	vector<C_fastq>* fq1s,*fq2s;
	vector<C_fastq> *trim_result1,*trim_result2,*clean_result1,*clean_result2;
	PEcalOption(){
	    local_fs=NULL;
	    fq1s=NULL;
	    fq2s=NULL;
        trim_result1=trim_result2=clean_result1=clean_result2=NULL;
	}
};
class peProcess{
public:
    explicit peProcess(C_global_parameter m_gp);
	void process();
	void print_stat();
	void update_stat(C_fastq_file_stat& fq1s_stat,C_fastq_file_stat& fq2s_stat,C_filter_stat& fs_stat,string type);
	void* stat_pe_fqs(PEstatOption opt,string dataType);
	void filter_pe_fqs(PEcalOption* opt);
	void* sub_thread(int index);
	int read(vector<C_fastq>& pe1,vector<C_fastq>& pe2,ifstream& infile1,ifstream& infile2);
	void peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,gzFile out1,gzFile out2);
	void peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,FILE* out1,FILE* out2);
	//void add_raw_trim(C_fastq_file_stat& a,C_fastq_file_stat& a2,C_reads_trim_stat& b,C_reads_trim_stat& b2);
	void peWrite_split(vector<C_fastq>& pe1,vector<C_fastq>& pe2);
	//void peRead();
	void  preOutput(int type,C_fastq& a);
	void output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile);
	void output_fastqs(string type,vector<C_fastq> &fq1,FILE* outfile);
	void output_split_fastqs(string type,vector<C_fastq> &fq1);
	//void peWrite(int num);
	void run_pigz(int type);
	void merge_stat();
	void merge_trim_data();
	void merge_clean_data();
	void merge_clean_data(int index);
	void C_fastq_init(C_fastq& a,C_fastq& b);
	void process_nonssd();
	void* sub_thread_nonssd_realMultiThreads(int index);
	void merge_stat_nonssd();
	void peStreaming_stat(C_global_variable& local_gv);
	void make_tmpDir();
	void remove_tmpDir();
	//void* monitor_read_thread();
	//void* monitor_write_thread();
	void create_thread_read(int index);
	void create_thread_trimoutputFile(int index);
	void create_thread_cleanoutputFile(int index);
	void thread_process_reads(int index,vector<C_fastq> &fq1s,vector<C_fastq> &fq2s);
	void run_cmd(string cmd);
	void run_extract_random();
	void process_some_reads(int index,int out_number);
	void merge_stat(int index);
	void limit_process_reads(int index,vector<C_fastq> &fq1s,vector<C_fastq> &fq2s,gzFile gzfq1,gzFile gzfq2);
	void check_disk_available();
	//void peOutput(outputOption opt);
public:
	C_global_parameter gp;
	C_global_variable gv;
	//queue<C_fastq> read_queue1,read_queue2,cal_trim_queue1,cal_trim_queue2,cal_clean_queue1,cal_clean_queue2;
	//vector<C_fastq> trim_merge1[max_thread],trim_merge2[max_thread],clean_merge1[max_thread],clean_merge2[max_thread];
//	gzFile gzfp1,gzfp2;
	gzFile *multi_gzfq1,*multi_gzfq2;
	//ifstream nongzfp1,nongzfp2;
	gzFile* gz_trim_out1,*gz_trim_out2;
	gzFile* gz_clean_out1,*gz_clean_out2;
	FILE** nongz_clean_out1;
	FILE** nongz_clean_out2;
	ofstream of_log;
	off_t* t_start_pos;
	off_t* t_end_pos;
//	char* src1,*src2;
	int fq1fd,fq2fd;
	string tmp_dir;
	int limit_end;
	//mutex* thread_read_m,*thread_write_m;
	//bool file_end;

	mutex pe_cal_m;
	mutex read_m,write_m;
	int buffer;
	string new_fq1_path,new_fq2_path;
	C_filter_stat *local_fs;
	C_fastq_file_stat* local_raw_stat1,*local_raw_stat2,*local_trim_stat1,*local_trim_stat2,*local_clean_stat1,*local_clean_stat2;
	//queue<vector<C_fastq> > thread_read_buffer1[max_thread],thread_read_buffer2[max_thread];
	string *sticky_reads1,*sticky_reads2;
	gzFile gz_fq1,gz_fq2;
	map<int,int> pe1_out,pe2_out;
	//int queue_buffer;
	int bq_check,pair_check;
private:
	int used_threads_num;
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