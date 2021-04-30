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
#include "BloomFilter.h"
#include <algorithm>
#include <fstream>
#include "rmdup.h"

#define max_thread 72
#define RMDUP 2
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
class seProcess{
public:
	explicit seProcess(C_global_parameter m_gp);
	void process();
	void print_stat();
	void update_stat(C_fastq_file_stat& fq1s_stat,C_filter_stat& fs_stat,string type);
	void* stat_se_fqs(SEstatOption opt,string dataType);
	//void add_raw_trim(C_fastq_file_stat& a,C_reads_trim_stat& b);
	void filter_se_fqs(SEcalOption opt);
	void* sub_thread(int index);
	int read(vector<C_fastq>& se1,ifstream& infile1);
	void seWrite(vector<C_fastq>& se1,gzFile out1);
	void seWrite(vector<C_fastq>& se1,FILE* out1);
	//void peRead();
	void  preOutput(int type,C_fastq& a);
	void output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile);
	void output_fastqs(string type,vector<C_fastq> &fq1,FILE* outfile);
	void merge_stat();
	void C_fastq_init(C_fastq& a);
	void seStreaming_stat(C_global_variable& local_gv);
	void remove_tmpDir();
	void make_tmpDir();
	void thread_process_reads(int index,int& cycle,vector<C_fastq> &fq1s);
	void create_thread_read(int index);
	void run_extract_random();
	void limit_process_reads(int index,vector<C_fastq> &fq1s,gzFile gzfq1);
	void check_disk_available();
    void catRmFile(int index,int cycle,string type,bool gzFormat);
    void catRmFile(vector<int> indexes,int cycle,string type,bool gzFormat);
    void* smallFilesProcess();

    void create_thread_smalltrimoutputFile(
            int index
            , int cycle);

    void create_thread_smallcleanoutputFile(
            int index
            , int cycle);

    void closeSmallTrimFileHandle(int index);

    void closeSmallCleanFileHandle(int index);

    void run_cmd(string cmd);

    void *sub_extract(
            string in
            , int mo
            , string out);

    void rmTmpFiles();

    void extractReadsToFile(
            int cycle
            , int thread_index
            , int reads_number
            , string position
            , int &output_index
            , bool gzFormat);

    void extractReadsToFile(
            int cycle
            , int thread_index
            , int reads_number
            , string position
            , bool gzFormat);

    void addCleanList(
            int tmp_cycle
            , int index);

    void *sub_thread_rmdup_step1(int index);

    void filter_se_fqs(
            SEcalOption opt
            , int index);

    int
            dupNum;
    mutex
            checkDup;
    BloomFilter
            *dupDB;
    gzFile
            dupOut1;
    gzFile
            *dupThreadOut1;
    set<string>
            checkDupMap;
    mutex
            logLock;
    //void peOutput(outputOption opt);
public:
    C_global_parameter
            gp;
    C_global_variable
            gv;
    gzFile
            *gz_trim_out1;
    gzFile
            *gz_clean_out1;
    FILE
            **nongz_clean_out1;
    FILE
            **nongz_trim_out1;
	ofstream of_log;
	string tmp_dir;
	gzFile* multi_gzfq1;
    FILE** multi_Nongzfq1;

    mutex
            se_stat_m,
            se_write_m;
    C_filter_stat
            *se_local_fs;
    C_fastq_file_stat
            *se_local_raw_stat1,
            *se_local_trim_stat1,
            *se_local_clean_stat1;
    int
            se_bq_check;
    int
            cur_cat_cycle;   //current cycle of cat small files
    vector<string>
            *readyTrimFiles1,
            *readyCleanFiles1; //list of small files in threads
    vector<int>
            *clean_file_readsNum; //reads number of small files in threads
    int
            *sub_thread_done;
    int
            end_sub_thread;
    int
            patch;
    uint64_t
            *threadReadsNum;
    uint64_t
            *totalData;
    uint64_t
            *threadCurReadReadsNumIdx;
    vector<vector<uint64_t *> >
            threadData;
    vector<vector<size_t> >
            threadDataNum;
    bool
            *dupFlag;
private:
    vector<C_fastq>
            fq1s;
    vector<C_fastq>
            trim_output_fq1;
    vector<C_fastq>
            clean_output_fq1;

    //C_fastq fastq1;
    int
            used_threads_num;
    string
            random_num;
};
#endif