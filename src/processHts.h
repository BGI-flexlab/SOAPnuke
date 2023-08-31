//
// Created by berry on 2020-04-08.
//
#ifndef SOAPNUKE_PROCESSHTS_H
#define SOAPNUKE_PROCESSHTS_H
#include <htslib/sam.h>
#include <htslib/cram.h>

//#include <cram/cram_samtools.h>
#include <htslib/bgzf.h>
#include <string.h>
#include <thread>
#include <limits.h>
#include "global_parameter.h"
#include "seprocess.h"
#include "peprocess.h"
#include "peFilterTmpOut.h"

using namespace::std;
#define READBUF 1000
//同时支持bam、cram的读写
//todo
//以后支持写入fastq
class processHts{
public:
    bool pe;
    C_global_parameter gp;
    seProcess *seP;
    peProcess *peP;
    string inputCram;
    htsFile *in;
    string outputFormat;
    //多线程相关数据结构
    htsFile *cleanOut,*filteredOut;
    htsFile** threadIn;
    htsFile** threadCleanOut,**threadFilteredOut;
    string tmpDir;
    int* threadProgress;
    int* threadWriteProgress;
    mutex writeLock;
    bool* threadDone;
    BGZF *cleanBam,*filteredBam;

//    cram_fd* cramOut;
    bam_hdr_t *header;
    bam1_t *aln;
    ofstream log;
    int8_t* seq_comp_table;
    int readsNumInPatch;
    //主线程读取数据模式下的每个子线程读取数量
    int lineNumPerThread;
    int keysNumber;
    const htsFormat* inputFormat;
    string inputFormatString;
    string outFormat;
    mutex* blockDataMutex;
    mutex writeMutex;
//    bamBlock** dataPool;
//    bamBlock** dataPool2;
    int pool1Stat;
    int pool2Stat;
    string* seLastID;
public:
    explicit processHts(C_global_parameter m_gp);

//    void process();

    void processSE();
    void processPE();

    void openHts();

    void *se_sub_thread(int i);
    void *peSingThread(int i);


    bam1_t** readSEData(vector<C_fastq> &pVector);
    char *get_read(bam1_t*);
    char *reverse(char *str);

    char *get_quality(bam1_t*);

    void filter_se_fqs(SEcalOption opt, bool *filterFlag);

    void writeBackToCram(bam1_t**,bool *filterFlag,int size);

    void writeBam();
    void writeCram();

    void closeHts(htsFile *fp);

    void prepare();

    void clean();

    bam1_t** readPEData(vector<C_fastq> &fq1s, vector<C_fastq> &fq2s, int *IDstat,bam1_t* lastLine,int &readLineNum);

    void writeBackToCram(bam1_t** data, bool *filterFlag, int size, int *IDstat);

    void check_disk_available();

    void filter_pe_fqs(PEcalOption *opt,bool* filterFlag);

    void writeSam();

    void pe_sortedByPos_sub_thread(int index);
    
    void writeFile();

    int readPEData(bam1_t** data,int patch);

    void parseRead(bam1_t *data, C_fastq& read);

    void multiThreadsPEprocess();
    void multiThreadsSEprocess();

    void* peSubThread(vector<C_fastq> fq1s,vector<C_fastq> fq2s,int start,int end,bool* threadFilter,int index);

    void* seSubThread(bam1_t** data,int start,int end,bool* threadFilter,int index);

    void writeFile(int index,int patch);

    void writeSam(int index, int patch);

    void writeBam(int index, int patch);

    void writeCram(int index, int patch);

    void closeThreadFd(int index);

    void *peDependentIOSubThread(int index);

    bam1_t **
    readPEData(vector<C_fastq> &fq1s, vector<C_fastq> &fq2s, int *IDstat,int &readLineNum, int index, long long &totalNum);

    void openHts(int index);

    void writeBackToCram(bam1_t **data, bool *filterFlag, int size, int index);

    void catSmallFiles();

    void catCram(vector<string> smallFiles, htsFile *out);

    void catBam(vector<string> smallFiles, BGZF* fp);

    void *seDependentIOSubThread(int index);

    bam1_t **readSEData(vector<C_fastq> &fq1s, int *IDstat, int &readLineNum, int index, long long &totalNum);

//    bamBlock *readPEData(bam1_t *lastLine);

    void* peSubThread2(int index);

    bam1_t **readSEData(int &lineNum);
};
#endif //SOAPNUKE_PROCESSHTS_H