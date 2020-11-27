//
// Created by berry on 2020-07-14.
//

#include <fstream>
#include "processStLFR.h"
#include "read_filter.h"
#define READBUF 500
processStLFR::processStLFR(C_global_parameter m_gp) : peProcess(m_gp) {
/*set<string> barcodeList;
    vector<int> barcodeStartPos;
    vector<int> barcodeLength;*/
    ifstream bl(m_gp.barcodeListPath);
    string line;
    int count=0;
    while(getline(bl,line)){
        count++;
        vector<string> eles;
        line_split(line,eles);

        string barcodeSeq=eles[0];
        for(int i=0;i<barcodeSeq.length();i++){
            barcodeSeq[i]=toupper(barcodeSeq[i]);
        }
        string barcodeIndex=eles[1];
//        string rcBarcodeSeq=reversecomplementary(barcodeSeq);

        //set 0/1 mismatch for barcode sequence
        string tmp="ACGT";
        for(int i=0;i<barcodeSeq.length();i++){
            for(int j=0;j<tmp.length();j++){
                string tmpBarcodeSeq=barcodeSeq;
                tmpBarcodeSeq[i]=tmp[j];
                int v=atoi(eles[1].c_str());
                barcodeList[tmpBarcodeSeq]=v;
            }
        }
    }
    gp.barcodeNumInList=count;
    bl.close();
    vector<string> eles;
    line_split(m_gp.barcodeRegionStr,',',eles);
    if(eles.size()!=3){
        cerr<<"Error:barcode region format error"<<endl;
        exit(1);
    }
    map<int,int> tmpMap;
    for(vector<string>::iterator ix=eles.begin();ix!=eles.end();ix++){
        vector<string> eles2;
        line_split(*ix,'_',eles2);
        if(eles2.size()!=2){
            cerr<<"Error:barcode region format error"<<endl;
            exit(1);
        }
        tmpMap[atoi(eles2[0].c_str())-1]=atoi(eles2[1].c_str());
        barcodeStartPos.push_back(atoi(eles2[0].c_str())-1);
//        barcodeLength.push_back(atoi(eles2[1].c_str()));
    }
    sort(barcodeStartPos.begin(),barcodeStartPos.end());
    for(int i=0;i<3;i++){
        if(barcodeStartPos[i]<0){
            cerr<<"Error:barcode region format error, barcode start pos should be positive integer"<<endl;
            exit(1);
        }
        barcodeLength.push_back(tmpMap[barcodeStartPos[i]]);
    }
    if(barcodeStartPos.size()!=barcodeLength.size()){
        cerr<<"Error:code error"<<endl;
        exit(1);
    }
//    if(gp.rmdup) {
//        //estimate total reads number
//
//        long long guessedReadsNum = 0;
//        if (gp.approximateReadsNum == 0) {
//            string fqPath=gp.fq1_path;
//            if(gp.inputAsList){
//                ifstream inList(gp.fq1_path);
//                while(getline(inList,fqPath)){
//                    guessedReadsNum += guessReadsNum(fqPath);
//                }
//                inList.close();
//            }
//        } else {
//            guessedReadsNum = gp.approximateReadsNum;
//        }
//        int multiple = 20;
//        while (multiple * guessedReadsNum > bfSize) {
//            multiple -= 5;
//            if (multiple < 10) {
//                cerr << "Error:reads number maybe is too large to do remove duplication" << endl;
//                exit(1);
//            }
//        }
//        dupDB = new BloomFilter(guessedReadsNum, multiple);
//        dupNum=0;
//    }
}

void processStLFR::filter_pe_fqs(PEcalOption *opt,int index) {
    vector<C_fastq>::iterator i2=opt->fq2s->begin();
    vector<C_fastq>::iterator i_end=opt->fq1s->end();
    bool* dupFilter=new bool[opt->fq1s->size()];
    if(gp.rmdup){
        if(RMDUP!=2)
            checkDup.lock();
        memset(dupFilter,false,opt->fq1s->size());
        int iter=0;
        for(vector<C_fastq>::iterator i=opt->fq1s->begin();i!=i_end;i++){
            string checkSeq=(*i).sequence+(*i2).sequence;
//            if(checkDupMap.find(checkSeq)!=checkDupMap.end()){
//                dupNum++;
//            }else{
//                checkDupMap.insert(checkSeq);
//            }
            if(RMDUP==0) {
                if (dupDB->query(checkSeq)) {
                    dupNum++;
                    dupFilter[iter] = true;
                    gzwrite(dupOut1, (*i).toString().c_str(), (*i).toString().size());
                    gzwrite(dupOut2, (*i2).toString().c_str(), (*i2).toString().size());
                } else {
                    dupDB->add();
                }
            }else if(RMDUP==1){
                if(RdupDB->query(checkSeq)){
                    dupNum++;
                    dupFilter[iter] = true;
                    gzwrite(dupOut1, (*i).toString().c_str(), (*i).toString().size());
                    gzwrite(dupOut2, (*i2).toString().c_str(), (*i2).toString().size());
                }else{
                    RdupDB->add();
                }
            }else{
                if(dupFlag[threadCurReadReadsNumIdx[index]-opt->fq1s->size()+iter]){
//                    dupNum++;
                    dupFilter[iter] = true;
                    gzwrite(dupThreadOut1[index], (*i).toString().c_str(), (*i).toString().size());
                    gzwrite(dupThreadOut2[index], (*i2).toString().c_str(), (*i2).toString().size());
                }
            }
            iter++;
            i2++;
            if(i2==opt->fq2s->end()){
                break;
            }
        }
        if(RMDUP!=2)
            checkDup.unlock();
    }
    i2=opt->fq2s->begin();
    i_end=opt->fq1s->end();
    int iter=0;
    for(vector<C_fastq>::iterator i=opt->fq1s->begin();i!=i_end;i++){
        string barcodeCombine=stLFRprocessBarcode(*i,*i2);
        C_pe_fastq_filter pe_fastq_filter=C_pe_fastq_filter(*i,*i2,gp);
        if(dupFilter[iter]){
            pe_fastq_filter.reads_result.dup=true;
        }
        iter++;
        pe_fastq_filter.setStLFRbarcode(barcodeCombine);
        /*int head_hdcut,head_lqcut,tail_hdcut,tail_lqcut,adacut_pos;
    int contam_pos;
    int global_contam_pos;
    int raw_length;*/
        pe_fastq_filter.pe_trim(gp);
        if(gp.adapter_discard_or_trim=="trim" || gp.contam_discard_or_trim=="trim" || !gp.trim.empty() || !gp.trimBadHead.empty() || !gp.trimBadTail.empty()){
            (*i).head_hdcut=pe_fastq_filter.fq1.head_hdcut;
            (*i).head_lqcut=pe_fastq_filter.fq1.head_lqcut;
            (*i).tail_hdcut=pe_fastq_filter.fq1.tail_hdcut;
            (*i).tail_lqcut=pe_fastq_filter.fq1.tail_lqcut;
            (*i).adacut_pos=pe_fastq_filter.fq1.adacut_pos;
            //(*i).contam_pos=pe_fastq_filter.fq1.contam_pos;
            //(*i).global_contam_pos=pe_fastq_filter.fq1.global_contam_pos;
            //(*i).raw_length=pe_fastq_filter.fq1.raw_length;
            (*i2).head_hdcut=pe_fastq_filter.fq2.head_hdcut;
            (*i2).head_lqcut=pe_fastq_filter.fq2.head_lqcut;
            (*i2).tail_hdcut=pe_fastq_filter.fq2.tail_hdcut;
            (*i2).tail_lqcut=pe_fastq_filter.fq2.tail_lqcut;
            (*i2).adacut_pos=pe_fastq_filter.fq2.adacut_pos;
            //(*i2).contam_pos=pe_fastq_filter.fq2.contam_pos;
            //(*i2).global_contam_pos=pe_fastq_filter.fq2.global_contam_pos;
            //(*i2).raw_length=pe_fastq_filter.fq2.raw_length;
        }
        if(!gp.trim_fq1.empty()){
            preOutput(1,pe_fastq_filter.fq1);
            preOutput(2,pe_fastq_filter.fq2);
            opt->trim_result1->emplace_back(pe_fastq_filter.fq1);
            opt->trim_result2->emplace_back(pe_fastq_filter.fq2);
        }
        if(pe_fastq_filter.pe_discard(opt->local_fs,gp)!=1){
            if(!gp.clean_fq1.empty()){
                preOutput(1,pe_fastq_filter.fq1);
                preOutput(2,pe_fastq_filter.fq2);
                opt->clean_result1->emplace_back(pe_fastq_filter.fq1);
                opt->clean_result2->emplace_back(pe_fastq_filter.fq2);
            }
        }
        i2++;
        if(i2==opt->fq2s->end()){
            break;
        }
    }
    delete[] dupFilter;
}

void *processStLFR::sub_thread(int index) {
    if(gp.inputAsList) {
        logLock.lock();
        of_log << get_local_time() << "\tthread " << index << " start" << endl;
        logLock.unlock();
//        create_thread_read(index);
        int thread_cycle = -1;

        ifstream inputList(gp.fq1_path);
        ifstream inputList2(gp.fq2_path);
        string fq1Path,fq2Path;
        char buf1[READBUF], buf2[READBUF];
        C_fastq fastq1, fastq2;
        C_fastq_init(fastq1, fastq2);
        int thread_read_block = 4 * gp.patchSize * patch;
        vector<C_fastq> fq1s, fq2s;
//        bool inputGzformat = true;
        long long file1_line_num(0), file2_line_num(0);
        long long block_line_num1(0), block_line_num2(0);
        while(1){
            getline(inputList,fq1Path);
            getline(inputList2,fq2Path);
            if(fq1Path=="" || fq2Path=="") {
                if (!fq1s.empty()) {
                    if (end_sub_thread == 1) {
                        break;
                    }
                    int tmp_cycle = file1_line_num / (thread_read_block * gp.threads_num);
                    if (tmp_cycle != thread_cycle && tmp_cycle > 0) {
                        addCleanList(thread_cycle, index);
                    }
                    thread_cycle = tmp_cycle;
                    threadCurReadReadsNumIdx[index]=file1_line_num/4;
                    thread_process_reads(index, thread_cycle, fq1s, fq2s);
                    if (limit_end > 0) {
                        break;
                    }
                }
                gzclose(multi_gzfq1[index]);
                gzclose(multi_gzfq2[index]);
                break;
            }
            multi_gzfq1[index] = gzopen(fq1Path.c_str(), "rb");
            if (!multi_gzfq1[index]) {
                cerr << "Error:cannot open the file," << fq1Path << endl;
                exit(1);
            }
            gzsetparams(multi_gzfq1[index], 2, Z_DEFAULT_STRATEGY);
            gzbuffer(multi_gzfq1[index], 2048 * 2048);
            multi_gzfq2[index] = gzopen(fq2Path.c_str(), "rb");
            if (!multi_gzfq2[index]) {
                cerr << "Error:cannot open the file," << fq2Path << endl;
                exit(1);
            }
            gzsetparams(multi_gzfq2[index], 2, Z_DEFAULT_STRATEGY);
            gzbuffer(multi_gzfq2[index], 2048 * 2048);

            gzFile tmpRead = gzopen(fq1Path.c_str(), "rb");
            int spaceNum = 0;
            if (gzgets(tmpRead, buf1, READBUF) != NULL) {
                string tmpLine(buf1);
                while (isspace(tmpLine[tmpLine.size() - 1])) {
                    spaceNum++;
                    tmpLine.erase(tmpLine.size() - 1);
                }
            }
            gzclose(tmpRead);
            while (1) {
                if (gzgets(multi_gzfq1[index], buf1, READBUF) != NULL) {
                    if ((file1_line_num / thread_read_block) % gp.threads_num == index) {
                        block_line_num1++;
                        if (block_line_num1 % 4 == 1) {
                            fastq1.seq_id.assign(buf1);
                            fastq1.seq_id.erase(fastq1.seq_id.size() - spaceNum, spaceNum);
                        }
                        if (block_line_num1 % 4 == 2) {
                            fastq1.sequence.assign(buf1);
                            fastq1.sequence.erase(fastq1.sequence.size() - spaceNum, spaceNum);
                        }
                        if (block_line_num1 % 4 == 0) {
                            fastq1.qual_seq.assign(buf1);
                            fastq1.qual_seq.erase(fastq1.qual_seq.size() - spaceNum, spaceNum);
                        }
                    }
                    file1_line_num++;
                }
                if (gzgets(multi_gzfq2[index], buf2, READBUF) != NULL) {
                    if ((file2_line_num / thread_read_block) % gp.threads_num == index) {
                        block_line_num2++;
                        if (block_line_num2 % 4 == 1) {
                            fastq2.seq_id.assign(buf2);
                            fastq2.seq_id.erase(fastq2.seq_id.size() - spaceNum, spaceNum);
                        } else if (block_line_num2 % 4 == 2) {
                            fastq2.sequence.assign(buf2);
                            fastq2.sequence.erase(fastq2.sequence.size() - spaceNum, spaceNum);
                        } else if (block_line_num2 % 4 == 0) {
                            fastq2.qual_seq.assign(buf2);
                            fastq2.qual_seq.erase(fastq2.qual_seq.size() - spaceNum, spaceNum);
                            fq1s.emplace_back(fastq1);
                            fq2s.emplace_back(fastq2);

                            if (fq1s.size() == gp.patchSize) {
                                if (end_sub_thread == 1) {
                                    break;
                                }
                                int tmp_cycle = file2_line_num / (thread_read_block * gp.threads_num);
                                if (tmp_cycle != thread_cycle && tmp_cycle > 0) {
                                    addCleanList(thread_cycle, index);
                                }
                                thread_cycle = tmp_cycle;
                                threadCurReadReadsNumIdx[index]=file1_line_num/4;
                                thread_process_reads(index, thread_cycle, fq1s, fq2s);
                                if (index == 0) {
                                    of_log << get_local_time() << " processed_reads:\t" << file1_line_num / 4 << endl;
                                }
                            }
                        }
                    }
                    file2_line_num++;
                } else {
                    break;
                }
            }
        }
        if (thread_cycle >= 0)
            addCleanList(thread_cycle, index);
        check_disk_available();
        sub_thread_done[index] = 1;
        logLock.lock();
        of_log << get_local_time() << "\tthread " << index << " done\t" << endl;
        logLock.unlock();
        return &bq_check;
    }else{
        peProcess::sub_thread(index);
        return &bq_check;
    }
}

string processStLFR::stLFRprocessBarcode(C_fastq &fastq1, C_fastq &fastq2) {
    if(fastq1.seq_id.find("/1")==string::npos || fastq2.seq_id.find("/2")==string::npos){
        cerr<<"Error:Reads1 and Reads2 ID error in /1 or /2,"<<fastq1.seq_id<<endl;
        exit(1);
    }
    string fastq1_id=fastq1.seq_id.substr(0,fastq1.seq_id.size()-2);
    string fastq2_id=fastq2.seq_id.substr(0,fastq2.seq_id.size()-2);
    if(fastq1_id!=fastq2_id){
        cerr<<"Error:Fastq reads ID unequal at reads,"<<fastq1_id<<"\t"<<fastq2_id<<endl;
        exit(1);
    }
    bool find=true;
    stringstream combineValue;
    for(int i=0;i<3;i++) {
        if(fastq2.sequence.size()<barcodeStartPos[i]+barcodeLength[i]){
            cerr<<"Error:given position and length exceeds the read sequence("<<fastq2.sequence.size()<<"), please check barcodeRegionStr parameter,"<<barcodeStartPos[i]<<"_"<<barcodeLength[i]<<endl;
            exit(1);
        }
        string expectedB = fastq2.sequence.substr(barcodeStartPos[i],barcodeLength[i]);
        if(barcodeList.find(expectedB)!=barcodeList.end()){
            combineValue<<barcodeList[expectedB];
        }else{
            find=false;
            break;
        }
        if(i!=2) {
            combineValue<< "_";
        }
    }
    if(find){
        if(!gp.tenX) {
            fastq1.seq_id = fastq1_id + "#" + combineValue.str() + "/1";
            fastq2.seq_id = fastq2_id + "#" + combineValue.str() + "/2";
        }else{
            fastq1.seq_id=fastq1_id+"_1\tBX:Z:"+combineValue.str();
            fastq2.seq_id=fastq2_id+"_2\tBX:Z:"+combineValue.str();
        }
        fastq2.sequence=fastq2.sequence.substr(0,barcodeStartPos[0]);
        fastq2.qual_seq=fastq2.qual_seq.substr(0,barcodeStartPos[0]);
        return combineValue.str();
    }else{
        if(!gp.tenX) {
            fastq1.seq_id = fastq1_id + "#" + "0_0_0" + "/1";
            fastq2.seq_id = fastq2_id + "#" + "0_0_0" + "/2";
        }else{
            fastq1.seq_id=fastq1_id+"_1\tBX:Z:"+"0_0_0";
            fastq2.seq_id=fastq2_id+"_2\tBX:Z:"+"0_0_0";
        }
        if(!gp.notCutNoLFR){
            fastq2.sequence=fastq2.sequence.substr(0,barcodeStartPos[0]);
            fastq2.qual_seq=fastq2.qual_seq.substr(0,barcodeStartPos[0]);
        }
        return "0_0_0";
    }
}

void *processStLFR::sub_thread_rmdup_step1(int index) {
    if(gp.inputAsList) {
        logLock.lock();
        of_log << get_local_time() << "\tthread " << index << " start" << endl;
        logLock.unlock();
//        create_thread_read(index);
        int thread_cycle = -1;

        ifstream inputList(gp.fq1_path);
        ifstream inputList2(gp.fq2_path);
        string fq1Path,fq2Path;
        char buf1[READBUF], buf2[READBUF];
        C_fastq fastq1, fastq2;
        C_fastq_init(fastq1, fastq2);
        int thread_read_block = 4 * gp.patchSize * patch;
        vector<C_fastq> fq1s, fq2s;
//        bool inputGzformat = true;
        long long file1_line_num(0), file2_line_num(0);
        long long block_line_num1(0), block_line_num2(0);
        vector<string> seqs;
        while(1){
            getline(inputList,fq1Path);
            getline(inputList2,fq2Path);
            string fq1seq,fq2seq;
            if(fq1Path=="" || fq2Path=="") {
                if (!seqs.empty()) {
                    uint64_t* curData=new uint64_t[seqs.size()];
                    memset(curData,0,sizeof(uint64_t)*seqs.size());
                    for(int i=0;i<seqs.size();i++){
                        curData[i]=hash<string>()(seqs[i]);
//                        curData[i]=new unsigned char[16];
//                        MDString(seqs[i].c_str(),curData[i]);
                    }
                    threadData[index].emplace_back(curData);
                    threadReadsNum[index]+=seqs.size();
                    threadDataNum[index].emplace_back(seqs.size());
                    seqs.clear();
                    if (limit_end > 0) {
                        break;
                    }
                }
                gzclose(multi_gzfq1[index]);
                gzclose(multi_gzfq2[index]);
                break;
            }
            multi_gzfq1[index] = gzopen(fq1Path.c_str(), "rb");
            if (!multi_gzfq1[index]) {
                cerr << "Error:cannot open the file," << fq1Path << endl;
                exit(1);
            }
            gzsetparams(multi_gzfq1[index], 2, Z_DEFAULT_STRATEGY);
            gzbuffer(multi_gzfq1[index], 2048 * 2048);
            multi_gzfq2[index] = gzopen(fq2Path.c_str(), "rb");
            if (!multi_gzfq2[index]) {
                cerr << "Error:cannot open the file," << fq2Path << endl;
                exit(1);
            }
            gzsetparams(multi_gzfq2[index], 2, Z_DEFAULT_STRATEGY);
            gzbuffer(multi_gzfq2[index], 2048 * 2048);
            gzFile tmpRead = gzopen(fq1Path.c_str(), "rb");
            int spaceNum = 0;
            if (gzgets(tmpRead, buf1, READBUF) != NULL) {
                string tmpLine(buf1);
                while (isspace(tmpLine[tmpLine.size() - 1])) {
                    spaceNum++;
                    tmpLine.erase(tmpLine.size() - 1);
                }
            }
            gzclose(tmpRead);

            while (1) {
                if (gzgets(multi_gzfq1[index], buf1, READBUF) != NULL) {
                    if ((file1_line_num / thread_read_block) % gp.threads_num == index) {
                        block_line_num1++;
                        if (block_line_num1 % 4 == 2) {
                            fq1seq.assign(buf1);
                            fq1seq.erase(fq1seq.size() - spaceNum,spaceNum);
                        }
                    }
                    file1_line_num++;
                }
                if (gzgets(multi_gzfq2[index], buf2, READBUF) != NULL) {
                    if ((file2_line_num / thread_read_block) % gp.threads_num == index) {
                        block_line_num2++;
                        if (block_line_num2 % 4 == 2) {
                            fq2seq.assign(buf2);
                            fq2seq.erase(fq2seq.size() - spaceNum,spaceNum);
                            string ligatedStr=fq1seq+fq2seq;
                            seqs.emplace_back(ligatedStr);
                        }
                        if (seqs.size() == gp.patchSize) {
                            uint64_t* curData=new uint64_t[seqs.size()];
                            memset(curData,0,sizeof(uint64_t)*seqs.size());
                            for(int i=0;i<seqs.size();i++){
                                curData[i]=hash<string>()(seqs[i]);
//                                curData[i]=new unsigned char[16];
//                                MDString(seqs[i].c_str(),curData[i]);
                            }
                            threadData[index].emplace_back(curData);
                            threadDataNum[index].emplace_back(seqs.size());
                            threadReadsNum[index]+=seqs.size();
                            seqs.clear();
                            if (index == 0) {
                                of_log << get_local_time() << " processed_reads:\t" << file1_line_num / 4 << endl;
                            }
                        }
                    }
                    file2_line_num++;
                } else {
                    break;
                }
            }
        }
        if (thread_cycle >= 0)
            addCleanList(thread_cycle, index);
        check_disk_available();
//        sub_thread_done[index] = 1;
        logLock.lock();
        of_log << get_local_time() << "\tthread " << index << " done\t" << endl;
        logLock.unlock();
        return &bq_check;
    }else{
        peProcess::sub_thread_rmdup_step1(index);
        return &bq_check;
    }
}
void processStLFR::filter_pe_fqs(PEcalOption *opt) {
    vector<C_fastq>::iterator
            i2
            = opt->fq2s
                 ->begin();
    vector<C_fastq>::iterator
            i_end
            = opt->fq1s
                 ->end();
//    bool* dupFilter=new bool[opt->fq1s->size()];
//    if(gp.rmdup){
//        checkDup.lock();
//        memset(dupFilter,false,opt->fq1s->size());
//        int iter=0;
//        for(vector<C_fastq>::iterator i=opt->fq1s->begin();i!=i_end;i++){
//            string checkSeq=(*i).sequence+(*i2).sequence;
////            if(checkDupMap.find(checkSeq)!=checkDupMap.end()){
////                dupNum++;
////            }else{
////                checkDupMap.insert(checkSeq);
////            }
//            if(dupDB->query(checkSeq)){
//                dupFilter[iter]=true;
//                gzwrite(dupOut1,(*i).toString().c_str(),(*i).toString().size());
//                gzwrite(dupOut2,(*i2).toString().c_str(),(*i2).toString().size());
//            }else{
//                dupDB->add();
//            }
//            iter++;
//            i2++;
//            if(i2==opt->fq2s->end()){
//                break;
//            }
//        }
//        checkDup.unlock();
//    }
//    i2=opt->fq2s->begin();
//    i_end=opt->fq1s->end();
    int
            iter
            = 0;
    for (
            vector<C_fastq>::iterator
            i
            = opt->fq1s
                 ->begin();
            i
            != i_end;
            i++
            )
    {
        string
                barcodeCombine
                = stLFRprocessBarcode(
                *i
                , *i2
                                     );
        C_pe_fastq_filter
                pe_fastq_filter
                = C_pe_fastq_filter(
                        *i
                        , *i2
                        , gp
                                   );
//        if(dupFilter[iter]){
//            pe_fastq_filter.reads_result.dup=true;
//        }
        iter++;
        pe_fastq_filter.setStLFRbarcode(barcodeCombine);
        /*int head_hdcut,head_lqcut,tail_hdcut,tail_lqcut,adacut_pos;
    int contam_pos;
    int global_contam_pos;
    int raw_length;*/
        pe_fastq_filter.pe_trim(gp);
        if (gp.adapter_discard_or_trim
            == "trim"
            || gp.contam_discard_or_trim
               == "trim"
            || !gp.trim
                  .empty()
            || !gp.trimBadHead
                  .empty()
            || !gp.trimBadTail
                  .empty())
        {
            (*i).head_hdcut
                            = pe_fastq_filter.fq1
                                             .head_hdcut;
            (*i).head_lqcut
                            = pe_fastq_filter.fq1
                                             .head_lqcut;
            (*i).tail_hdcut=pe_fastq_filter.fq1.tail_hdcut;
            (*i).tail_lqcut=pe_fastq_filter.fq1.tail_lqcut;
            (*i).adacut_pos=pe_fastq_filter.fq1.adacut_pos;
            //(*i).contam_pos=pe_fastq_filter.fq1.contam_pos;
            //(*i).global_contam_pos=pe_fastq_filter.fq1.global_contam_pos;
            //(*i).raw_length=pe_fastq_filter.fq1.raw_length;
            (*i2).head_hdcut=pe_fastq_filter.fq2.head_hdcut;
            (*i2).head_lqcut=pe_fastq_filter.fq2.head_lqcut;
            (*i2).tail_hdcut=pe_fastq_filter.fq2.tail_hdcut;
            (*i2).tail_lqcut=pe_fastq_filter.fq2.tail_lqcut;
            (*i2).adacut_pos=pe_fastq_filter.fq2.adacut_pos;
            //(*i2).contam_pos=pe_fastq_filter.fq2.contam_pos;
            //(*i2).global_contam_pos=pe_fastq_filter.fq2.global_contam_pos;
            //(*i2).raw_length=pe_fastq_filter.fq2.raw_length;
        }
        if(!gp.trim_fq1.empty()){
            preOutput(1,pe_fastq_filter.fq1);
            preOutput(2,pe_fastq_filter.fq2);
            opt->trim_result1->emplace_back(pe_fastq_filter.fq1);
            opt->trim_result2->emplace_back(pe_fastq_filter.fq2);
        }
        if(pe_fastq_filter.pe_discard(opt->local_fs,gp)!=1){
            if(!gp.clean_fq1.empty()){
                preOutput(1,pe_fastq_filter.fq1);
                preOutput(2,pe_fastq_filter.fq2);
                opt->clean_result1->emplace_back(pe_fastq_filter.fq1);
                opt->clean_result2->emplace_back(pe_fastq_filter.fq2);
            }
        }
        i2++;
        if(i2==opt->fq2s->end()){
            break;
        }
    }
//    delete[] dupFilter;
}
