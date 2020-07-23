//
// Created by berry on 2020-07-14.
//

#ifndef SOAPNUKE_PROCESSSTLFR_H
#define SOAPNUKE_PROCESSSTLFR_H
#include<map>
#include "peprocess.h"

class processStLFR: public peProcess{
    //    void process();
public:
    explicit processStLFR(C_global_parameter m_gp);
    void filter_pe_fqs(PEcalOption* opt) override;
    void *sub_thread(int index) override;
private:
    std::map<string,int> barcodeList;
    vector<int> barcodeStartPos;
    vector<int> barcodeLength;

    string stLFRprocessBarcode(C_fastq &fastq1, C_fastq &fastq2);
};
#endif //SOAPNUKE_PROCESSSTLFR_H
