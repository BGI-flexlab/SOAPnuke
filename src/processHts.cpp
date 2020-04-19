//
// Created by berry on 2020-04-08.
//

#include "processHts.h"
#define bam_is_read1(b) (((b)->core.flag&BAM_FREAD1)!=0)
#define bam_is_read2(b) (((b)->core.flag&BAM_FREAD2)!=0)
processHts::processHts(C_global_parameter m_gp) {

    gp=m_gp;
    //屏蔽trim功能
    gp.adapter_discard_or_trim="discard";
    gp.contam_discard_or_trim="discard";
    gp.trim="";
    gp.trimBadHead="";
    gp.trimBadTail="";

    //每次读取的数据量，影响内存占用量
    readsNuminPatch=100000;
    //inputCram=m_gp.fq1_path; //fq1_path沿用了之前处理fastq的命名，代表了输入文件
    openHts();

    seq_comp_table= new int8_t[16]{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };//整数转ATCG编码值
    //判断是SE还是PE
    while(sam_read1(in,header,aln)>=0){
        if(aln->core.flag&BAM_FPAIRED){
            pe=true;
        }else{
            pe=false;
        }
        break;
    }
    bam_destroy1(aln);
    if(pe){
        peP=new peProcess(gp);
    }else{
        seP=new seProcess(gp);
    }


    closeHts(in);


}
void processHts::check_disk_available(){
    if(access(gp.fq1_path.c_str(),0)==-1){
        cerr<<"Error:input raw fastq not exists suddenly, please check the disk"<<endl;
        exit(1);
    }
    if(access(gp.output_dir.c_str(),0)==-1){
        cerr<<"Error:output directory cannot open suddenly, please check the disk"<<endl;
        exit(1);
    }
}

void processHts::processSE() {
    prepare();
    se_sub_thread(0);
    seP->check_disk_available();
    seP->merge_stat();
    seP->print_stat();
    seP->check_disk_available();
    clean();
}
void processHts::processPE() {
    prepare();
    pe_sub_thread(0);
    check_disk_available();
    peP->merge_stat();
    peP->print_stat();
    check_disk_available();
    clean();
}
void* processHts::pe_sub_thread(int index){

    C_fastq fastq1,fastq2;
    peP->C_fastq_init(fastq1,fastq2);
    vector<C_fastq> fq1s,fq2s;
//    string inputFormat="cram";
    if(gp.fq2_path.rfind(".bam")==gp.fq2_path.size()-4){
        writeBam();
    }else {
        if(gp.fq2_path.rfind(".sam")==gp.fq2_path.size()-4){
            writeSam();
        }else {
            writeCram();
        }
    }
    bam1_t* startAln=bam_init1();
    while(1) {
        int* readNameGroupCount=new int[readsNuminPatch];
        memset(readNameGroupCount,0,readsNuminPatch);
        int readLineNum=0;
        bam1_t** bamData=readPEData(fq1s,fq2s,readNameGroupCount,startAln,readLineNum);
        int readsNum=fq1s.size();
        if(readsNum==0){
            delete[] readNameGroupCount;
            break;
        }
        bool* filterFlag=new bool[readsNum];
        memset(filterFlag,0,readsNum);

        PEcalOption* opt2=new PEcalOption();
        vector<C_fastq>  clean_result1;
        vector<C_fastq>  clean_result2;
        opt2->local_fs=&(peP->local_fs[index]);
        opt2->fq1s=&fq1s;
        opt2->fq2s=&fq2s;
        opt2->clean_result1=&clean_result1;
        opt2->clean_result2=&clean_result2;
        filter_pe_fqs(opt2,filterFlag);
        PEstatOption opt_raw,opt_clean;
        opt_raw.fq1s=&fq1s;
        opt_raw.stat1=&(peP->local_raw_stat1[index]);
        opt_raw.fq2s=&fq2s;
        opt_raw.stat2=&(peP->local_raw_stat2[index]);
        peP->stat_pe_fqs(opt_raw,"raw");		//statistic raw fastqs
        //add_raw_trim(local_raw_stat1[index],local_raw_stat2[index],raw_cut.stat1,raw_cut.stat2);
        fq1s.clear();
        fq2s.clear();
        opt_clean.stat1=&(peP->local_clean_stat1[index]);
        opt_clean.stat2=&(peP->local_clean_stat2[index]);
        opt_clean.fq1s=&clean_result1;
        opt_clean.fq2s=&clean_result2;
        peP->stat_pe_fqs(opt_clean,"clean");

        bool* peFilterFlag=new bool[readLineNum];
        memset(peFilterFlag,0,readLineNum);
        //assgin filter flag to each read line
        int iter=0;
        for(int i=0;i<readsNum;i++){
            for(int j=0;j<readNameGroupCount[i];j++){
                peFilterFlag[iter]=filterFlag[i];
                iter++;
            }
        }
        writeBackToCram(bamData,peFilterFlag,iter);
        for(int i=0;i<iter;i++){
            bam_destroy1(bamData[i]);
            bamData[i]=NULL;
        }
        delete[] bamData;
        clean_result1.clear();
        clean_result2.clear();
        delete[] peFilterFlag;
        delete[] filterFlag;
        delete[] readNameGroupCount;
    }

    return &seP->se_bq_check;
}



void* processHts::se_sub_thread(int index){
//    seP.of_log<<get_local_time()<<"\tthread "<<index<<" start"<<endl;
//	}
//    seP->create_thread_read(index);

    C_fastq fastq1;
    seP->C_fastq_init(fastq1);
    vector<C_fastq> fq1s;
//    string inputFormat="cram";
    if(gp.fq2_path.rfind(".bam")==gp.fq2_path.size()-4){
        writeBam();
    }else {
        writeCram();
    }
    while(1) {
        bam1_t** bamData=readSEData(fq1s);
        int readsNum=fq1s.size();
        if(readsNum==0){
            break;
        }
        bool* filterFlag=new bool[readsNum];
        memset(filterFlag,0,readsNum);
        vector<C_fastq> clean_result1;
        SEcalOption opt2;
        opt2.se_local_fs = &(seP->se_local_fs[index]);
        opt2.fq1s = &fq1s;
        opt2.clean_result1=&clean_result1;
        filter_se_fqs(opt2,filterFlag);        //filter raw fastqs by the given parameters
        SEstatOption opt_raw;
        opt_raw.fq1s = &fq1s;
        opt_raw.stat1 = &(seP->se_local_raw_stat1[index]);
        seP->stat_se_fqs(opt_raw, "raw");        //statistic raw fastqs
        fq1s.clear();
        SEstatOption opt_clean;
        opt_clean.stat1=&(seP->se_local_clean_stat1[index]);
        opt_clean.fq1s=&clean_result1;
        seP->stat_se_fqs(opt_clean,"clean");
        writeBackToCram(bamData,filterFlag,readsNum);
        for(int i=0;i<readsNum;i++){
            bam_destroy1(bamData[i]);
            bamData[i]=NULL;
        }
        delete[] filterFlag;
        delete[] bamData;
        clean_result1.clear();
    }

    return &seP->se_bq_check;
}
void processHts::clean(){
    log << get_local_time() << "\tAnalysis accomplished!" << endl;
    log.close();
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    closeHts(out);
    closeHts(in);
}
void processHts::prepare(){
    string mkdir_str="mkdir -p "+gp.output_dir;
    if(system(mkdir_str.c_str())==-1){
        cerr<<"Error:mkdir fail"<<endl;
        exit(1);
    }
    log.open(gp.log);
    if(!log){
        cerr<<"Error, cannot write to such file,"<<gp.log<<endl;
        exit(-1);
    }
    openHts();
    string format;
    switch(inputFormat->format){
        case 4:format="bam";break;
        case 6:format="cram";break;
        default:{
            cerr<<"Error, only support BAM/CRAM in this module"<<endl;
            exit(1);
        }
    }
    inputFormatString=format;
    string whetherPE=pe?"PE":"SE";
    log<<"input file format:"<<format<<endl;
    log<<"reads in file are "<<whetherPE<<endl;
    log<<get_local_time()<<"\tAnalysis start!"<<endl;
}
void processHts::writeBackToCram(bam1_t** data,bool* filterFlag,int size){
    for(int i=0;i<size;i++){
        if(filterFlag[i]){
            if((data[i]->core.flag & 512)!=1) {
                data[i]->core.flag += 512;
            }
        }
        if(sam_write1(out,header,data[i])<0){
            cerr<<"Error, write file error"<<endl;
            exit(1);
        }
    }
}
void processHts::filter_pe_fqs(PEcalOption* opt,bool* filterFlag){
    //C_reads_trim_stat_2 cut_pos;
    vector<C_fastq>::iterator i2=opt->fq2s->begin();
    vector<C_fastq>::iterator i_end=opt->fq1s->end();
    int iter=0;
    for(vector<C_fastq>::iterator i=opt->fq1s->begin();i!=i_end;i++){
        C_pe_fastq_filter pe_fastq_filter=C_pe_fastq_filter(*i,*i2,gp);
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
            peP->preOutput(1,pe_fastq_filter.fq1);
            peP->preOutput(2,pe_fastq_filter.fq2);
            opt->trim_result1->emplace_back(pe_fastq_filter.fq1);
            opt->trim_result2->emplace_back(pe_fastq_filter.fq2);
        }
        if(pe_fastq_filter.pe_discard(opt->local_fs,gp)!=1){
//            peP->preOutput(1,pe_fastq_filter.fq1);
//            peP->preOutput(2,pe_fastq_filter.fq2);
            opt->clean_result1->emplace_back(pe_fastq_filter.fq1);
            opt->clean_result2->emplace_back(pe_fastq_filter.fq2);
        }else{
            filterFlag[iter]=true;
        }
        iter++;
        i2++;
        if(i2==opt->fq2s->end()){
            break;
        }
    }
    //return cut_pos;
}

void processHts::filter_se_fqs(SEcalOption opt,bool* filterFlag){
    //C_reads_trim_stat cut_pos;
    int iter=0;
    for(vector<C_fastq>::iterator i=opt.fq1s->begin();i!=opt.fq1s->end();i++){
        C_single_fastq_filter se_fastq_filter=C_single_fastq_filter(*i,gp);
        se_fastq_filter.se_trim(gp);
        if(gp.adapter_discard_or_trim=="trim" || gp.contam_discard_or_trim=="trim" || !gp.trim.empty() || !gp.trimBadHead.empty() || !gp.trimBadTail.empty()){
            (*i).head_hdcut=se_fastq_filter.read.head_hdcut;
            (*i).head_lqcut=se_fastq_filter.read.head_lqcut;
            (*i).tail_hdcut=se_fastq_filter.read.tail_hdcut;
            (*i).tail_lqcut=se_fastq_filter.read.tail_lqcut;
            (*i).adacut_pos=se_fastq_filter.read.adacut_pos;
            //(*i).contam_pos=se_fastq_filter.read.contam_pos;
            //(*i).global_contam_pos=se_fastq_filter.read.global_contam_pos;
            //(*i).raw_length=se_fastq_filter.read.raw_length;
        }
        //*i=se_fastq_filter.read;
        if(!gp.trim_fq1.empty()){
            seP->preOutput(1,se_fastq_filter.read);
            opt.trim_result1->emplace_back(se_fastq_filter.read);
        }
        int whether_discard(0);
        if(gp.module_name=="filtersRNA"){
            whether_discard=se_fastq_filter.sRNA_discard(opt.se_local_fs,gp);
        }else{
            whether_discard=se_fastq_filter.se_discard(opt.se_local_fs,gp);
        }
        if(whether_discard!=1){
//            seP->preOutput(1,se_fastq_filter.read);
            opt.clean_result1->emplace_back(se_fastq_filter.read);
        }else{
            filterFlag[iter]=true; //set filter flag
        }
        iter++;
    }
    //return cut_pos;
}
bam1_t** processHts::readPEData(vector<C_fastq> &fq1s,vector<C_fastq> &fq2s,int* IDstat,bam1_t* lastLine,int &readLineNum) {


    C_fastq fastq1,fastq2;
    bam1_t** data=new bam1_t*[readsNuminPatch*2];
    for(int i=0;i<readsNuminPatch*2;i++){
        data[i]=bam_init1();
    }
    int iter=0;
    int readIDNum=0;
    peP->C_fastq_init(fastq1,fastq2);
    aln=bam_init1();
    string lastQname="";
    bool scanRead1=false;
    bool scanRead2=false;
    set<string> read1ID,read2ID;
    int lineNum=0;
    if(bam_get_qname(lastLine)!=NULL){
        readIDNum++;
        if(bam_copy1(aln,lastLine)==NULL){
            cerr<<"Error, bam copy error"<<endl;
            exit(1);
        }
        char* seq=get_read();
        char* qual=get_quality();
        char* qname=bam_get_qname(aln);
        if(seq==NULL || qual==NULL || qname==NULL){
            cerr<<"Error, parse "<<inputFormatString<<" file error"<<endl;
            exit(1);
        }
        string qnameS=qname;

        if(bam_is_read1(aln)){
            fastq1.seq_id.assign(qname);
            fastq1.sequence.assign(seq);
            fastq1.qual_seq.assign(qual);
            scanRead1=true;
        }else{
            fastq2.seq_id.assign(qname);
            fastq2.sequence.assign(seq);
            fastq2.qual_seq.assign(qual);
            scanRead2=true;
        }
        delete[] seq;
        delete[] qual;
        if(bam_copy1(data[lineNum],aln)==NULL){
            cerr<<"Error, copy data in memory failed"<<endl;
            exit(1);
        }
        lastQname=qnameS;
        lineNum++;
    }
    while(sam_read1(in,header,aln)>=0){
//        aln->core.seq;
        char* seq=get_read();
        char* qual=get_quality();
        char* qname=bam_get_qname(aln);
        if(seq==NULL || qual==NULL || qname==NULL){
            cerr<<"Error, parse "<<inputFormatString<<" file error"<<endl;
            exit(1);
        }
        string qnameS=qname;

        if(bam_copy1(data[lineNum],aln)==NULL){
            cerr<<"Error, copy data in memory failed"<<endl;
            exit(1);
        }


        if(lastQname!=qnameS && lastQname!="" && lineNum>1){
            if(fastq1.sequence.size()==0 || fastq2.sequence.size()==0){
                cerr<<"code error"<<endl;
                exit(1);
            }
            fq1s.emplace_back(fastq1);
            fq2s.emplace_back(fastq2);
            IDstat[iter]=readIDNum;
            readIDNum=0;
            if(lineNum+10>=readsNuminPatch*2){
                if(bam_copy1(lastLine,aln)==NULL){
                    cerr<<"Error, bam copy error"<<endl;
                    exit(1);
                }
                return data;
            }


            scanRead1=false;
            scanRead2=false;
            iter++;
        }
        lastQname=qnameS;

        if(bam_is_read1(aln) && !scanRead1){
            fastq1.seq_id.assign(qname);
            fastq1.sequence.assign(seq);
            fastq1.qual_seq.assign(qual);
            scanRead1=true;
        }else{
            if(!scanRead2) {
                fastq2.seq_id.assign(qname);
                fastq2.sequence.assign(seq);
                fastq2.qual_seq.assign(qual);
                scanRead2 = true;
            }
        }
        delete[] seq;
        delete[] qual;
        lineNum++;
        readLineNum=lineNum;
        readIDNum++;

    }
    if(fq1s.size()>0) {
        if (fastq1.sequence.size() == 0 || fastq2.sequence.size() == 0) {
            cerr << "code error" << endl;
            exit(1);
        }
        fq1s.emplace_back(fastq1);
        fq2s.emplace_back(fastq2);
        IDstat[iter] = readIDNum;
    }
    return data;
}
bam1_t** processHts::readSEData(vector<C_fastq> &fq1s) {


    C_fastq fastq1;
    bam1_t** data=new bam1_t*[readsNuminPatch];
    for(int i=0;i<readsNuminPatch;i++){
        data[i]=bam_init1();
    }
    int iter=0;
    seP->C_fastq_init(fastq1);
    aln=bam_init1();
    while(sam_read1(in,header,aln)>=0){
//        aln->core.seq;
        char* seq=get_read();
        char* qual=get_quality();
        char* qname=bam_get_qname(aln);
        if(seq==NULL || qual==NULL || qname==NULL){
            cerr<<"Error, parse cram file error"<<endl;
            exit(1);
        }
        fastq1.seq_id.assign(qname);
        fastq1.sequence.assign(seq);
        fastq1.qual_seq.assign(qual);
        delete[] seq;
        delete[] qual;
//        delete qname;
        fq1s.emplace_back(fastq1);
        if(bam_copy1(data[iter],aln)==NULL){
            cerr<<"Error, copy data in memory failed"<<endl;
            exit(1);
        }
        if (fq1s.size() == readsNuminPatch) {
            return data;
            break;
        }
        iter++;
    }
    for(int i=fq1s.size();i<readsNuminPatch;i++){
        data[i]=NULL;
    }
    return data;
}
void processHts::openHts(){
    in=hts_open(gp.fq1_path.c_str(),"r");
    if(in==NULL){
        cerr<<"Error, cannot open such file,"<<inputCram<<endl;
        exit(1);
    }
    if(gp.fq1_path.rfind(".cram")==gp.fq1_path.size()-5 || gp.fq2_path.rfind(".cram")==gp.fq2_path.size()-5) {
        string refPath = gp.reference;
        string refFai = refPath + ".fai";
        if (hts_set_fai_filename(in, refFai.c_str()) < 0) {
            cerr << "Error, reference is needed, cannot open such file," << refFai << endl;
            exit(1);
        }
    }
    inputFormat=hts_get_format(in);
    if(inputFormat->category==1){
        if(inputFormat->format!=4 && inputFormat->format!=6){
            cerr<<"Error, only support BAM/CRAM in this module"<<endl;
            exit(1);
        }
    }else{
        cerr<<"Error, only support BAM/CRAM in this module"<<endl;
        exit(1);
    }
    header=sam_hdr_read(in);
    if(header==NULL){
        cerr<<"Error, get header fail"<<endl;
        exit(1);
    }
    aln=bam_init1();
    if(aln==NULL){
        cerr<<"Error, init error,"<<__FILE__<<":line"<<__LINE__<<endl;
        exit(1);
    }
}
void processHts::writeCram(){
    string outPath=gp.output_dir+"/"+gp.fq2_path;
    out=sam_open(outPath.c_str(),"wc");
    if(out==NULL){
        cerr<<"Error, cannot write to such file,"<<gp.fq2_path<<endl;
        exit(1);
    }
//    if(cram_set_header(out->fp.cram,header)==-1){
//        cerr<<"Error, set header for cram failed"<<endl;
//        exit(1);
//    }
    string refPath=gp.reference;
    string refFai=refPath+".fai";
    if(hts_set_fai_filename(out,refFai.c_str())<0){
        cerr<<"Error, cannot open such file,"<<refFai<<endl;
        exit(1);
    }
    if(sam_hdr_write(out,header)==-1){
        cerr<<"Error, write header to cram failed"<<endl;
        exit(1);
    }
}
void processHts::writeBam(){
    string outPath=gp.output_dir+"/"+gp.fq2_path;
    out=sam_open(outPath.c_str(),"wb");
    if(out==NULL){
        cerr<<"Error, cannot write to such file,"<<gp.fq2_path<<endl;
        exit(1);
    }
    if(sam_hdr_write(out,header)<0){
        cerr<<"Error, write file error"<<endl;
        exit(1);
    }
}
void processHts::writeSam(){
    string outPath=gp.output_dir+"/"+gp.fq2_path;
    out=sam_open(outPath.c_str(),"w");
    if(out==NULL){
        cerr<<"Error, cannot write to such file,"<<gp.fq2_path<<endl;
        exit(1);
    }
    if(sam_hdr_write(out,header)<0){
        cerr<<"Error, write file error"<<endl;
        exit(1);
    }
}
void processHts::closeHts(htsFile* fp){
    if(hts_close(fp)<0){
        cerr<<"Error, cannot close cram file"<<endl;
        exit(1);
    }
}
char *processHts::reverse(char *str)
{
    int i = strlen(str)-1,j=0;
    char ch;
    while (i>j) {
        ch = str[i];
        str[i]= str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
}
char *processHts::get_read()
{
    int len = aln->core.l_qseq + 1;
    char *read = new char[len];
    memset(read,0,len);
    char *seq = (char *)bam_get_seq(aln);
    int n;

    if (!read) return NULL;

    for (n=0; n < aln->core.l_qseq; n++) {
        if (aln->core.flag & BAM_FREVERSE) read[n] = seq_nt16_str[seq_comp_table[bam_seqi(seq,n)]];
        else                               read[n] = seq_nt16_str[bam_seqi(seq,n)];
    }
    if (aln->core.flag & BAM_FREVERSE) reverse(read);
    return read;
}
char* processHts::get_quality()
{
    char* qual_out = new char[ aln->core.l_qseq + 1];
    memset(qual_out,0,aln->core.l_qseq + 1);
    char *q = (char *)bam_get_qual(aln);
    int n;

    if (!qual_out) return NULL;

    if (*q == '\xff') {
        free(qual_out);
        qual_out = NULL;
        return NULL;
    }

    for (n=0; n < aln->core.l_qseq; n++) {
        qual_out[n] = q[n]+33;
    }
    if (aln->core.flag & BAM_FREVERSE) reverse(qual_out);
    return qual_out;
}
