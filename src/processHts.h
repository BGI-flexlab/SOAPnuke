//
// Created by berry on 2020-04-08.
//

#ifndef SOAPNUKE_PROCESSHTS_H
#define SOAPNUKE_PROCESSHTS_H
#include <htslib/sam.h>
#include <htslib/cram.h>
#include "global_parameter.h"
#include "seprocess.h"
#include "peprocess.h"
#include <string.h>
using namespace::std;
#define READBUF 500
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
    htsFile *in,*out;
//    cram_fd* cramOut;
    bam_hdr_t *header;
    bam1_t *aln;
    ofstream log;
    int8_t* seq_comp_table;
    int readsNuminPatch;
    const htsFormat* inputFormat;
    string inputFormatString;
    string outFormat;
    explicit processHts(C_global_parameter m_gp);

//    void process();

    void processSE();
    void processPE();

    void openHts();

    void *se_sub_thread(int i);
    void *pe_sub_thread(int i);


    bam1_t** readSEData(vector<C_fastq> &pVector);
    char *get_read();
    char *reverse(char *str);

    char *get_quality();

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
};
#endif //SOAPNUKE_PROCESSHTS_H
