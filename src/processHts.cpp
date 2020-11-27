//
// Created by berry on 2020-04-08.
//
#include "processHts.h"
#define bam_is_read1(b) (((b)->core.flag&BAM_FREAD1)!=0)
#define bam_is_read2(b) (((b)->core.flag&BAM_FREAD2)!=0)
#define BAMCPERR "bam copy error"
//if MultiThreadsMethod is set to 1, IO independent
//if is set to 2, IO is locked
#define MultiThreadsMethod 2
processHts::processHts(C_global_parameter m_gp) {

    gp=m_gp;
    //屏蔽trim功能
    gp.adapter_discard_or_trim="discard";
    gp.contam_discard_or_trim="discard";
    gp.trim="";
    gp.trimBadHead="";
    gp.trimBadTail="";

    if(gp.fq2_path.rfind(".bam")==gp.fq2_path.size()-4){
        outputFormat="bam";
    }else if(gp.fq2_path.rfind(".sam")==gp.fq2_path.size()-4){
        outputFormat="sam";
    }else if(gp.fq2_path.rfind(".cram")==gp.fq2_path.size()-5){
        outputFormat="cram";
    }else{
        cerr<<"Error:only support sam/bam/cram suffix format output"<<endl;
        exit(1);
    }
    //每次读取的数据量，影响内存占用量
    //readsNumInPatch=gp.threads_num>1?100000/gp.threads_num:100000;
    readsNumInPatch=100000;

    lineNumPerThread=100000;
    //key number limit in map
    keysNumber=50000000;
    //inputCram=m_gp.fq1_path; //fq1_path沿用了之前处理fastq的命名，代表了输入文件
    if(gp.threads_num>1) {

//        dataPool=new bamBlock*[gp.threads_num];
//        dataPool2=new bamBlock*[gp.threads_num];
        pool1Stat=0;
        pool2Stat=0;
        //clean的时候删除
        threadIn=new htsFile*[gp.threads_num];
        threadProgress=new int[gp.threads_num];
        threadDone=new bool[gp.threads_num];
        threadWriteProgress=new int[gp.threads_num];
        for(int i=0;i<gp.threads_num;i++){
            openHts(i);
            threadProgress[i]=-2;
            threadDone[i]=0;
            threadWriteProgress[i]=0;
        }
        //threadProgress 代表处理批次进度
        //threadDone 代表是否完成
        threadCleanOut = new htsFile *[gp.threads_num];
//        threadFilteredOut = new htsFile *[gp.threads_num];
        //由于要输出小文件，线程的写句柄是变化的，所以无法在构造函数中提前写好
//        for(int i=0;i<gp.threads_num;i++){
//            writeFile(i);
//        }
    }
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
        seP=new seProcess(gp);
    }else{
        seP=new seProcess(gp);
    }
    bam_hdr_destroy(header);

    closeHts(in);
    blockDataMutex=new mutex[gp.threads_num];
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

    if(gp.threads_num>1){
        seLastID = new string[gp.threads_num];
        for (int i = 0; i < gp.threads_num; i++) {
            seLastID[i] = "";
        }

        if (outputFormat == "bam") {
            if(MultiThreadsMethod==2){
                writeFile();
            }else {
                string cleanOut = gp.output_dir + "/" + gp.fq2_path;
                cleanBam = bgzf_open(cleanOut.c_str(), "w");
                if (bam_hdr_write(cleanBam, header) < 0) {
                    cerr << "Error:couldn't write header" << __FILE__ << " " << __LINE__ << endl;
                    exit(1);
                }
            }
//            string filterOut=gp.output_dir+"/_failQC."+outputFormat;
//            filteredBam=bgzf_open(filterOut.c_str(), "w");
//            if (bam_hdr_write(filteredBam, header) < 0) {
//                cerr<<"Error:couldn't write header"<<__FILE__<<" "<<__LINE__<<endl;
//                exit(1);
//            }
        } else {
            writeFile();
        }
        if(MultiThreadsMethod==1) {
            tmpDir = gp.output_dir + "/TMP";
            if (access(tmpDir.c_str(), 0) != -1) {
                string rmCmd = "rm -rf " + tmpDir;
                if (system(rmCmd.c_str()) == -1) {
                    cerr << "Error:when running " << rmCmd << endl;
                    exit(1);
                }
            }
//        mkdir(tmpDir.c_str(),0755);
            mkDir(tmpDir.c_str(), 0755);
        }
        thread t_array[gp.threads_num];
        for (int i = 0; i < gp.threads_num; i++) {
            t_array[i] = thread(bind(&processHts::seDependentIOSubThread, this, i));
        }
        if(MultiThreadsMethod==1){
            thread catFiles = thread(bind(&processHts::catSmallFiles, this));
            catFiles.join();
        }
        for (int i = 0; i < gp.threads_num; i++) {
            t_array[i].join();
        }
    }else {
        se_sub_thread(0);
    }
    seP->check_disk_available();
    seP->merge_stat();
    seP->print_stat();
    seP->check_disk_available();
    clean();
}
void processHts::processPE() {
    prepare();

//    处理按照基因组位置排序的PE bam/cram
//    pe_sortedByPos_sub_thread(0);
//  读取一批数据，分发给其它线程处理
    if(gp.threads_num>1){
        if (outputFormat == "bam") {
            if(MultiThreadsMethod==2) {
                writeFile();
            }else {
                string cleanOut = gp.output_dir + "/" + gp.fq2_path;
                cleanBam = bgzf_open(cleanOut.c_str(), "w");
                if (bam_hdr_write(cleanBam, header) < 0) {
                    cerr << "Error:couldn't write header" << __FILE__ << " " << __LINE__ << endl;
                    exit(1);
                }

//                string filterOut = gp.output_dir + "/_failQC." + outputFormat;
//                filteredBam = bgzf_open(filterOut.c_str(), "w");
//                if (bam_hdr_write(filteredBam, header) < 0) {
//                    cerr << "Error:couldn't write header" << __FILE__ << " " << __LINE__ << endl;
//                    exit(1);
//                }
            }
        } else {
            writeFile();
        }
//            writeFile();
        if(MultiThreadsMethod==1) {
            tmpDir = gp.output_dir + "/TMP";
            if (access(tmpDir.c_str(), 0) != -1) {
                string rmCmd = "rm -rf " + tmpDir;
                if (system(rmCmd.c_str()) == -1) {
                    cerr << "Error:when running " << rmCmd << endl;
                    exit(1);
                }
            }
//            mkdir(tmpDir.c_str(), 0755);
            mkDir(tmpDir.c_str(), 0755);
        }
        thread t_array[gp.threads_num];
        for (int i = 0; i < gp.threads_num; i++) {
            t_array[i] = thread(bind(&processHts::peDependentIOSubThread, this, i));
        }
        if(MultiThreadsMethod==1) {
            thread catFiles = thread(bind(&processHts::catSmallFiles, this));
            catFiles.join();
        }
        for (int i = 0; i < gp.threads_num; i++) {
            t_array[i].join();
        }
//        }else{
//            multiThreadsPEprocess();
//        }
       //
    }else{
        //由于该函数内涉及peProcess类中的线程函数，需要给1个线程索引作为参数，这里只设置为0即为单线程
        peSingThread(0);
    }
    check_disk_available();
    peP->merge_stat();
    peP->print_stat();
    check_disk_available();
    clean();
}
//合并小文件,需要一个存储各个线程的输出进度
void processHts::catSmallFiles(){
    //threadProgress[gp.threads_num] 存储了各个线程已完成的patch数目
    //首先判断是否完成
    int curCatCycle=0;
    while(1) {
        bool allDone = true;
        for (int i = 0; i < gp.threads_num; i++) {
            if (threadDone[i] == 0) {
                allDone = false;
                break;
            }
        }
        vector<string> cleanList, filteredList;
        if (allDone) {
            int availableCycle = -1;
            for (int i = 0; i < gp.threads_num; i++) {
                if (threadProgress[i] > availableCycle) {
                    availableCycle = threadProgress[i];
                }
            }
            for (int i = curCatCycle; i <=  availableCycle; i++) {
                for (int j = 0; j < gp.threads_num; j++) {
                    if (threadProgress[j] >= i) {
                        string cleanFileName = tmpDir + "/" + gp.fq2_path + ".t" + to_string(j) + "." + to_string(i);
                        cleanList.push_back(cleanFileName);
//                        string filteredFileName =
//                                tmpDir + "/_failQC." + outputFormat + ".t" + to_string(j) + "." + to_string(i);
//                        filteredList.push_back(filteredFileName);
                    }
                }
            }
        } else {
            int availableCycle = INT_MAX;
            for (int i = 0; i < gp.threads_num; i++) {
                if (threadProgress[i] < availableCycle) {
                    availableCycle = threadProgress[i];
                }
            }
            if (availableCycle != INT_MAX && availableCycle>=0) {
                for (int i = curCatCycle; i <= availableCycle; i++) {
                    for (int j = 0; j < gp.threads_num; j++) {
                        //按顺序拿出可合并列表
                        string cleanFileName = tmpDir + "/" + gp.fq2_path + ".t" + to_string(j) + "." + to_string(i);
                        cleanList.push_back(cleanFileName);
//                        string filteredFileName =
//                                tmpDir + "/_failQC." + outputFormat + ".t" + to_string(j) + "." + to_string(i);
//                        filteredList.push_back(filteredFileName);
                    }
                }
                curCatCycle=availableCycle+1;
            }
        }
        //将列表中的文件cat到指定的输出文件中
        if (cleanList.size() > 0) {
            if(outputFormat=="cram") {
                catCram(cleanList, cleanOut);
            }else if(outputFormat=="bam"){
                catBam(cleanList,cleanBam);
            }else{
                cerr<<"Error:not support such format output"<<endl;
                exit(1);
            }
            for (int i = 0; i < cleanList.size(); i++) {
                remove(cleanList[i].c_str());
            }
            cleanList.clear();
        }
//        if(filteredList.size()>0){
//            if(outputFormat=="cram") {
//                catCram(filteredList,filteredOut);
//            }else if(outputFormat=="bam"){
//                catBam(filteredList,filteredBam);
//            }else{
//                cerr<<"Error:not support such format output"<<endl;
//                exit(1);
//            }
//            for(int i=0;i<filteredList.size();i++){
//                remove(filteredList[i].c_str());
//            }
//            filteredList.clear();
//        }
        if(allDone){
            break;
        }
        sleep(2);
    }
};
void processHts::catBam(vector<string> smallFiles,BGZF* fp){
    BGZF *in = NULL;
    uint8_t *buf = NULL;
    int BUF_SIZE=0x10000;
    int BGZF_EMPTY_BLOCK_SIZE=28;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;
    int nfn=smallFiles.size();


    buf = (uint8_t*) malloc(BUF_SIZE);
    if (!buf) {
        fprintf(stderr, "[%s] Couldn't allocate buffer\n", __func__);
        exit(1);
    }
    for(i = 0; i < nfn; ++i){
        int len,j;

        in = bgzf_open(smallFiles[i].c_str(), "r");


        if (in == 0) {
            cerr<<"Error:fail to open file,"<<smallFiles[i]<<endl;
            exit(1);
        }
        if (in->is_write) exit(1);
        sam_hdr_t *old = bam_hdr_read(in);
        if (old == NULL) {
            cerr<<"Error:couldn't read header for "<<smallFiles[i]<<endl;
            exit(1);
        }
        if (in->block_offset < in->block_length) {
            if (bgzf_write(fp, (char *)in->uncompressed_block + in->block_offset, in->block_length - in->block_offset) < 0) {
                cerr<<"Error:bgzf write error "<<__FILE__<<" "<<__LINE__<<endl;
                exit(1);
            };
            if (bgzf_flush(fp) != 0){
                cerr<<"Error:cannot flush bgzf"<<endl;
                exit(1);
            }
        }
        j=0;
        while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0) {
            if(len<es){
                int diff=es-len;
                if(j==0) {
                    fprintf(stderr, "[%s] ERROR: truncated file?: '%s'.\n", __func__, smallFiles[i].c_str());
                    exit(1);
                }
                if (bgzf_raw_write(fp, ebuf, len) < 0){
                    cerr<<"Error:cannot write bgzf "<<__FILE__<<" "<<__LINE__<<endl;
                    exit(1);
                }

                memcpy(ebuf,ebuf+len,diff);
                memcpy(ebuf+diff,buf,len);
            } else {
                if(j!=0) {
                    if (bgzf_raw_write(fp, ebuf, es) < 0){
                        cerr<<"Error:cannot write bgzf "<<__FILE__<<" "<<__LINE__<<endl;
                        exit(1);
                    }
                }
                len-= es;
                memcpy(ebuf,buf+len,es);
                if (bgzf_raw_write(fp, buf, len) < 0){
                    cerr<<"Error:cannot write bgzf "<<__FILE__<<" "<<__LINE__<<endl;
                    exit(1);
                }
            }
            j=1;
        }
        int GZIPID1=31;
        int GZIPID2=139;
        /* check final gzip block */
        {
            const uint8_t gzip1=ebuf[0];
            const uint8_t gzip2=ebuf[1];
            const uint32_t isize=*((uint32_t*)(ebuf+es-4));
            if(((gzip1!=GZIPID1) || (gzip2!=GZIPID2)) || (isize!=0)) {
                fprintf(stderr, "[%s] WARNING: Unexpected block structure in file '%s'.", __func__, smallFiles[i].c_str());
                fprintf(stderr, " Possible output corruption.\n");
                exit(1);
                if (bgzf_raw_write(fp, ebuf, es) < 0){
                    cerr<<"Error:error in "<<__FILE__<<" "<<__LINE__<<endl;
                    exit(1);
                }
            }
        }
        bam_hdr_destroy(old);
        bgzf_close(in);
        in = NULL;
        if(bgzf_flush(fp)==-1){
            cerr<<"Error:bgzf flush error"<<endl;
            exit(1);
        }
    }
    free(buf);
}
void processHts::catCram(vector<string> smallFiles,htsFile* out){
    int nfn=smallFiles.size();
//    samFile *out;
    cram_fd *out_c;
    int i;
    out_c = out->fp.cram;
//    cram_set_option(out_c, CRAM_OPT_VERSION, vers);
    //fprintf(stderr, "Creating cram vers %s\n", vers);
    for (i = 0; i < nfn; ++i) {
        samFile *in;
        cram_fd *in_c;
        cram_container *c;
        sam_hdr_t *old_h;

        in = sam_open(smallFiles[i].c_str(), "rc");
        if (in == 0) {
            cout<<"Error:cannot open such file,"<<smallFiles[i]<<endl;
            exit(1);
        }
        in_c = in->fp.cram;

        old_h = sam_hdr_read(in);
        if (!old_h) {
            cout<<"Error:fail to read the header of file "<<smallFiles[i]<<endl;
            exit(1);
        }
        // Copy contains and blocks within them
        while ((c = cram_read_container(in_c))) {
            cram_block *blk;

            if (cram_container_is_empty(in_c)) {
                if (cram_write_container(out_c, c) != 0)
                    exit(1);

                // Container compression header
                if (!(blk = cram_read_block(in_c)))
                    exit(1);
                if (cram_write_block(out_c, blk) != 0) {
                    cram_free_block(blk);
                    exit(1);
                }
                cram_free_block(blk);
                cram_free_container(c);

                continue;
            }
            int32_t num_slices;

            // Not switching rg so do the usual read/write loop
            if (cram_write_container(out_c, c) != 0)
                exit(1);

            // Container compression header
            if (!(blk = cram_read_block(in_c)))
                exit(1);
            if (cram_write_block(out_c, blk) != 0) {
                cram_free_block(blk);
                exit(1);
            }
            cram_free_block(blk);


            // Container num_blocks can be invalid, due to a bug.
            // Instead we iterate in slice context instead.
            (void)cram_container_get_landmarks(c, &num_slices);
            cram_copy_slice(in_c, out_c, num_slices);

            cram_free_container(c);
        }

        bam_hdr_destroy(old_h);
        sam_close(in);
    }
}
//该方法是主线程读写，其它线程计算的模式
//这是主线程方法
//void processHts::multiThreadsPEprocess(){
//    writeFile();
//    thread t_array[gp.threads_num];
//    for (int i = 0; i < gp.threads_num; i++) {
//        t_array[i] = thread(bind(&processHts::peSubThread2, this, i));
//    }
//    vector<C_fastq> fq1s,fq2s;
//    bam1_t* startAln=bam_init1();
////    thread *subThreads=new thread[gp.threads_num];
//    long long totalReadsNum=0;
//
//    int blockNum=0;
//    while(1) {
//        int *readNameGroupCount = new int[readsNumInPatch];
//        memset(readNameGroupCount, 0, readsNumInPatch);
//        int readLineNum = 0;
//        //主线程读取数据，然后放到内存池供其它线程使用
//        //判断内存池是否可用，可用的话就往里面写新读的数据，不然就等待
//        clock_t start=clock();
//        bamBlock* bamData = readPEData(startAln);
//        clock_t finish=clock();
//        cout<<finish-start<<endl;
//        if(bamData->getRecordsNum()==0){
//
//            break;
//        }
//        int assignThreadIndex=blockNum%gp.threads_num;
//        if(pool1Stat==0){
//            if(pool2Stat!=0) {
//                dataPool[assignThreadIndex] = bamData;
//                if (assignThreadIndex == gp.threads_num - 1) {
//                    pool1Stat = 1;
//                }
//            }else{
//                blockDataMutex[assignThreadIndex].lock();
//                dataPool2[assignThreadIndex]=bamData;
//                if(assignThreadIndex==gp.threads_num-1){
//                    pool2Stat=1;
//                }
//                blockDataMutex[assignThreadIndex].unlock();
//            }
//        }else{
//            blockDataMutex[assignThreadIndex].lock();
//            dataPool2[assignThreadIndex]->~bamBlock();
//            dataPool2[assignThreadIndex]=dataPool[assignThreadIndex];
//            if(assignThreadIndex==gp.threads_num-1){
//                pool2Stat=0;
//            }
//            blockDataMutex[assignThreadIndex].unlock();
//        }
//
//        //将过滤flag合并到filterFlag
//    }
//}
void processHts::multiThreadsSEprocess(){
    writeFile();
//    vector<C_fastq> fq1s;
//    thread *subThreads=new thread[gp.threads_num];
    thread t_array[gp.threads_num];

    bool breakFlag=0;
    while(1) {
        if(breakFlag){
            break;
        }
        int readLineNum=0;
        bam1_t** bamData=readSEData(readLineNum);
        if(readLineNum<lineNumPerThread*gp.threads_num){
            breakFlag=1;
        }
        if(readLineNum==0){
            break;
        }
        bool *filterFlag = new bool[readLineNum];
        memset(filterFlag, 0, readLineNum);
        int meanAssigned=lineNumPerThread;
        int* assignedNum=new int[gp.threads_num];
        memset(assignedNum,0,gp.threads_num);
        bool** threadFilter=new bool*[gp.threads_num];
        memset(threadFilter,0,gp.threads_num);
        int realUseThreadsNum=0;
        int lastEndIndex=-1;
        string lastReadID="";
        bool lastReadFilter=false;
        int startIndex=-1;
        int endIndex=-1;
        C_fastq fastq1;
        seP->C_fastq_init(fastq1);
        int checkStart=0;
        for(int i=0;i<gp.threads_num;i++){
            if(i==0){
                while(checkStart<readLineNum){
                    parseRead(bamData[0],fastq1);
                    if(fastq1.seq_id==lastReadID){
                        checkStart++;
                    }else{
                        break;
                    }
                }
                startIndex=checkStart;
            }else{
                startIndex=lastEndIndex+1;
            }
            endIndex=(i+1)*meanAssigned;

            realUseThreadsNum=i;
            if(endIndex>readLineNum){
                endIndex=readLineNum;
                if(startIndex>=endIndex){
                    break;
                }
            }
            if(endIndex<readLineNum) {
                int checkEnd = endIndex;
                parseRead(bamData[checkEnd - 1], fastq1);
                string checkReadID = fastq1.seq_id;
                while (checkEnd--) {
                    parseRead(bamData[checkEnd - 1], fastq1);
                    if (fastq1.seq_id != checkReadID) {
                        break;
                    } else {
                        checkReadID = fastq1.seq_id;
                    }
                }
                endIndex=checkEnd+1;
                lastEndIndex=checkEnd;
            }
            threadFilter[i]=new bool[endIndex-startIndex];
            memset(threadFilter[i],0,endIndex-startIndex);
            assignedNum[i]=endIndex-startIndex;
            t_array[i]=thread(bind(&processHts::seSubThread,this,bamData,startIndex,endIndex,threadFilter[i],i));
        }
        for(int i=0;i<=realUseThreadsNum;i++){
            t_array[i].join();
        }
//        vector<C_fastq>().swap(fq1s);
        //将过滤flag合并到filterFlag
        int iter=0;
        for(int i=0;i<realUseThreadsNum;i++){
            for(int j=0;j!=assignedNum[i];j++){
                if(iter<checkStart){
                    filterFlag[iter]=lastReadFilter;
                }else {
                    filterFlag[iter] = threadFilter[i][j];
                }
                iter++;
            }
            delete[] threadFilter[i];
        }
        lastReadFilter=filterFlag[iter-1];
        parseRead(bamData[readLineNum-1],fastq1);
        lastReadID=fastq1.seq_id;
        delete[] threadFilter;
        delete[] assignedNum;

        writeBackToCram(bamData,filterFlag,readLineNum);
        delete[] filterFlag;
        for(int i=0;i<readLineNum;i++){
            bam_destroy1(bamData[i]);
            bamData[i]=NULL;
        }
        delete[] bamData;
    }
}

void* processHts::seSubThread(bam1_t** data,int start,int end,bool* threadFilter,int index){
    //parse fastq from data(bam1_t**) to fq1s
    vector<C_fastq> fq1s;
    C_fastq fastq1;
    seP->C_fastq_init(fastq1);
    string lastReadID="";
    int* readsGroupNum=new int[end-start];
    memset(readsGroupNum,0,end-start);
    int iter=0;
    for(int i=start;i<end;i++){
        parseRead(data[i],fastq1);
        if(fastq1.seq_id!=lastReadID) {
            fq1s.emplace_back(fastq1);
            readsGroupNum[iter]=1;
            iter++;
        }else{
            readsGroupNum[iter-1]++;
        }
        lastReadID=fastq1.seq_id;
    }
//    int testTotalNum=0;
//    for(int i=0;i<fq1s.size();i++){
//        if(readsGroupNum[i]>1){
//            cout<<readsGroupNum[i]<<"\there"<<endl;
//        }
//        testTotalNum+=readsGroupNum[i];
//    }
//    cout<<testTotalNum<<endl;
    int readsNum=fq1s.size();
    bool* readsNumFilter=new bool[readsNum];
    memset(readsNumFilter,0,readsNum);
    vector<C_fastq> clean_result1;
    SEcalOption opt2;
    opt2.se_local_fs = &(seP->se_local_fs[index]);
    opt2.fq1s = &fq1s;
    opt2.clean_result1=&clean_result1;
    filter_se_fqs(opt2,readsNumFilter);        //filter raw fastqs by the given parameters
    for(int i=0;i<readsNum;i++){

    }
    SEstatOption opt_raw;
    opt_raw.fq1s = &fq1s;
    opt_raw.stat1 = &(seP->se_local_raw_stat1[index]);
    seP->stat_se_fqs(opt_raw, "raw");        //statistic raw fastqs
    vector<C_fastq>().swap(fq1s);
    SEstatOption opt_clean;
    opt_clean.stat1=&(seP->se_local_clean_stat1[index]);
    opt_clean.fq1s=&clean_result1;
    seP->stat_se_fqs(opt_clean,"clean");
    vector<C_fastq>().swap(clean_result1);
    iter=0;
    for(int i=0;i<readsNum;i++){
        for(int j=0;j<readsGroupNum[i];j++){
            threadFilter[iter]=readsNumFilter[i];
            iter++;
        }
    }
    delete[] readsGroupNum;
    delete[] readsNumFilter;
    return (void*)NULL;
}

//void* processHts::peSubThread2(int index){
//    //获得fq1s和fq2s
//    vector<C_fastq> fq1s,fq2s;
//    while(1){
//        blockDataMutex[index].lock();
//        bamBlock* threadBlock=dataPool2[index];
//        int* readGroupIndex=threadBlock->getReadGroupIndex();
//        if(threadBlock->getRecordsNum()==0){
//            blockDataMutex[index].unlock();
//            sleep(1);
//        }else{
//            int lineNum=0;
//            for(int i=0;i<threadBlock->getReadsNum();i++){
//                C_fastq fastq1,fastq2;
//                char* seq=get_read((threadBlock->getRecords())[lineNum]);
//                char* qual=get_quality((threadBlock->getRecords())[lineNum]);
//                char* qname=bam_get_qname((threadBlock->getRecords())[lineNum]);
//                if(seq==NULL || qual==NULL || qname==NULL){
//                    cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
//                    exit(1);
//                }
//                for(int j=0;j<readGroupIndex[i];j++){
//                    if(bam_is_read1((threadBlock->getRecords())[lineNum])){
//                        fastq1.seq_id.assign(qname);
//                        fastq1.sequence.assign(seq);
//                        fastq1.qual_seq.assign(qual);
//                    }else{
//                        if(bam_is_read2((threadBlock->getRecords())[lineNum])) {
//                            fastq2.seq_id.assign(qname);
//                            fastq2.sequence.assign(seq);
//                            fastq2.qual_seq.assign(qual);
//                        }
//                    }
//                    lineNum++;
//                }
//                fq1s.emplace_back(fastq1);
//                fq2s.emplace_back(fastq2);
//            }
//            bool* threadFilter=new bool[threadBlock->getReadsNum()];
//            memset(threadFilter,0,threadBlock->getReadsNum());
//            PEcalOption* opt2=new PEcalOption();
//            vector<C_fastq>  clean_result1;
//            vector<C_fastq>  clean_result2;
//            opt2->local_fs=&(peP->local_fs[index]);
//            opt2->fq1s=&fq1s;
//            opt2->fq2s=&fq2s;
//            opt2->clean_result1=&clean_result1;
//            opt2->clean_result2=&clean_result2;
//            //过滤，将过滤结果传给threadFilter
//            filter_pe_fqs(opt2,threadFilter);
//            //统计原始fastq和过滤后fastq的信息
//            PEstatOption opt_raw,opt_clean;
//            opt_raw.fq1s=&fq1s;
//            opt_raw.stat1=&(peP->local_raw_stat1[index]);
//            opt_raw.fq2s=&fq2s;
//            opt_raw.stat2=&(peP->local_raw_stat2[index]);
//            peP->stat_pe_fqs(opt_raw,"raw");		//statistic raw fastqs
//            vector<C_fastq>().swap(fq1s);
//            vector<C_fastq>().swap(fq2s);
//            opt_clean.stat1=&(peP->local_clean_stat1[index]);
//            opt_clean.stat2=&(peP->local_clean_stat2[index]);
//            opt_clean.fq1s=&clean_result1;
//            opt_clean.fq2s=&clean_result2;
//            peP->stat_pe_fqs(opt_clean,"clean");
//            vector<C_fastq>().swap(clean_result1);
//            vector<C_fastq>().swap(clean_result2);
//            return (void*)NULL;
//        }
//    }
//
//}
//这是分线程
//处理分配给某个线程的fastq数据，传递的是总的数据和索引，由于共同操作总数据，需要在线程里面拷贝出来，避免线程间冲突
//此方法用于主线程读取数据，其它线程只处理数据（占用cpu，不占IO），如果改成线程读写分离，那么此方法就用不到
void* processHts::peSubThread(vector<C_fastq> fq1s,vector<C_fastq> fq2s,int start,int end,bool* threadFilter,int index){
    if(fq1s.size()<end || fq2s.size()<end){
        cerr<<"Error:code error"<<endl;
    }
    vector<C_fastq> realFq1,realFq2;
    for(int i=start-1;i!=end;i++){
        realFq1.push_back(fq1s[i]);
        realFq2.push_back(fq2s[i]);
    }

    PEcalOption* opt2=new PEcalOption();
    vector<C_fastq>  clean_result1;
    vector<C_fastq>  clean_result2;
    opt2->local_fs=&(peP->local_fs[index]);
    opt2->fq1s=&realFq1;
    opt2->fq2s=&realFq2;
    opt2->clean_result1=&clean_result1;
    opt2->clean_result2=&clean_result2;
    //过滤，将过滤结果传给threadFilter
    filter_pe_fqs(opt2,threadFilter);
    //统计原始fastq和过滤后fastq的信息
    PEstatOption opt_raw,opt_clean;
    opt_raw.fq1s=&realFq1;
    opt_raw.stat1=&(peP->local_raw_stat1[index]);
    opt_raw.fq2s=&realFq2;
    opt_raw.stat2=&(peP->local_raw_stat2[index]);
    peP->stat_pe_fqs(opt_raw,"raw");		//statistic raw fastqs
    vector<C_fastq>().swap(realFq1);
    vector<C_fastq>().swap(realFq2);
    opt_clean.stat1=&(peP->local_clean_stat1[index]);
    opt_clean.stat2=&(peP->local_clean_stat2[index]);
    opt_clean.fq1s=&clean_result1;
    opt_clean.fq2s=&clean_result2;
    peP->stat_pe_fqs(opt_clean,"clean");
    vector<C_fastq>().swap(clean_result1);
    vector<C_fastq>().swap(clean_result2);
    return (void*)NULL;
}
//此方法用于处理按基因组位置排序的数据，写了一部分，待完成
void processHts::pe_sortedByPos_sub_thread(int index){
    map<size_t ,peFilterTmpOut*> dataInfo;
    writeFile();
    //open a tmp file to write, which store information like this:
    // readID   peFilterTmpOut
    // information in peFilterTmpOut:
    //  read1 or 2,multi-aligned or not,filter or not,line number from raw bam/cram file
    string tmpOut=gp.output_dir+"/_tmpFile";
    ofstream ofTmpOut(tmpOut.c_str());
    if(!ofTmpOut){
        cerr<<"Error:cannot write to such file,"<<tmpOut<<endl;
        exit(1);
    }
    //todo
    //open another tmp file to store supplemantary align record

    //first cycle to scan bam/cram file
    map<size_t,vector<long long> > supAlignInfo;
    long long lineNumber=0;
    C_filter_stat* fq1FilterStat=new C_filter_stat();
    C_filter_stat* fq2FilterStat=new C_filter_stat();
    while(1){
        //read a number of reads
        bam1_t** bamData=new bam1_t*[readsNumInPatch];
        for(int i=0;i<readsNumInPatch;i++){
            bamData[i]=bam_init1();
        }
        int readLineNumber=readPEData(bamData,readsNumInPatch);
        if(readLineNumber==0){
            break;
        }
        hash<string> h;
        for(int i=0;i<readLineNumber;i++){
            lineNumber++;
            C_fastq fastq1;
            seP->C_fastq_init(fastq1);
            parseRead(bamData[i],fastq1);
            bool readType=bam_is_read1(bamData[i]);
            size_t intKey=h(fastq1.seq_id);
            if((bamData[i]->core.flag & BAM_FSUPPLEMENTARY)!=0){
                supAlignInfo[intKey].push_back(lineNumber);
            }else{
                C_single_fastq_filter fastq_filter=C_single_fastq_filter(fastq1,gp);
                int whether_discard=0;
                if(readType){
                    whether_discard=fastq_filter.se_discard(fq1FilterStat,gp);
                }else{
                    whether_discard=fastq_filter.se_discard(fq2FilterStat,gp);
                }
                peFilterTmpOut* info=new peFilterTmpOut(readType,whether_discard,lineNumber);
                if(dataInfo.find(intKey)!=dataInfo.end()){
                    peFilterTmpOut* mateReadInfo=dataInfo[intKey];
                    if(mateReadInfo->readType != readType){
                        //pe reads both get
                        int peFilter=whether_discard || mateReadInfo->filter;
                        ofTmpOut<<intKey<<"\t"<<peFilter<<"\t"<<mateReadInfo->lineNumber<<endl;
                        ofTmpOut<<intKey<<"\t"<<peFilter<<"\t"<<lineNumber<<endl;
                        dataInfo.erase(intKey);
                    }
                }else{
                    dataInfo[intKey]=info;
                }
            }
        }
        if(readLineNumber<readsNumInPatch){
            break;
        }
    }
    ofTmpOut.close();
    if(dataInfo.size()>0){
        //将未处理的reads输出到这个文件，看看有哪些例外情况
        ofstream ofNotProcessedReads(gp.output_dir+"/_unProcessReads");
        if(!ofNotProcessedReads){
            cerr<<"Error:cannot write to such file,"<<gp.output_dir+"/_unProcessReads"<<endl;
            exit(1);
        }
        for(map<size_t,peFilterTmpOut*>::iterator ix=dataInfo.begin();ix!=dataInfo.end();ix++){
            ofNotProcessedReads<<(*ix).first<<"\t"<<(*ix).second->toString()<<endl;
        }
        ofNotProcessedReads.close();
    }
    //处理含有2048的reads
    ifstream ifPairedFile(gp.output_dir+"/_tmpFile");
    ofstream ofWholeReadsFile(gp.output_dir+"/_tmpFile2");
    if(!ofWholeReadsFile){
        cerr<<"Error:cannot write to such file,"<<gp.output_dir+"/_tmpFile2";
        exit(1);
    }
    if(!ifPairedFile){
        cerr<<"Error:cannot open such file,"<<gp.output_dir+"/_tmpFile"<<endl;
        exit(1);
    }
    string lineInfo;

    while(getline(ifPairedFile,lineInfo)){
        vector<string> eles;
        line_split(lineInfo,'\t',eles);
        size_t intKey=atoi(eles[0].c_str());
        ofWholeReadsFile<<eles[2]<<"\t"<<eles[1]<<endl;
        if(supAlignInfo.find(intKey)!=supAlignInfo.end()){
//            int filterResult=atoi(eles[1].c_str());
            for(vector<long long>::iterator ix=supAlignInfo[intKey].begin();ix!=supAlignInfo[intKey].end();ix++){
                ofWholeReadsFile<<*ix<<"\t"<<eles[1]<<endl;
            }
            supAlignInfo.erase(intKey);
            if(supAlignInfo.size()==0){
                break;
            }
        }else{
//            cerr<<"Error:code error or raw data error"<<endl;
//            exit(1);
        }
    }
    ofWholeReadsFile.close();
    ifPairedFile.close();
    //todo
    //unaccomplished
}
//解析bam1_t数据结构成C_fastq
void processHts::parseRead(bam1_t* data,C_fastq& read){
    char* seq=get_read(data);
    char* qual=get_quality(data);
    char* qname=bam_get_qname(data);
    if(seq==NULL || qual==NULL || qname==NULL){
        cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
        exit(1);
    }
    read.seq_id.assign(qname);
    read.sequence.assign(seq);
    read.qual_seq.assign(qual);
    delete[] seq;
    delete[] qual;
}
//最初设计目的是用于单线程，因此会有个记录分块读取后的尾巴数据的变量startAln，以此来保留上一个读取数据
//若设计成多线程，需要在readPEData中额外添加线程索引，并且保留startAln（防止线程数设置为1）
//多线程情况下，此方法将拥有独立IO和计算功能，最后由主线程将输出结果合并
//为了保证输出结果顺序，将每次处理结果都输出，然后由主线程同时做合并
void* processHts::peSingThread(int index){

    vector<C_fastq> fq1s,fq2s;
//    string inputFormat="cram";
    writeFile();
    bam1_t* startAln=bam_init1();
    while(1) {
        int* readNameGroupCount=new int[readsNumInPatch];
        memset(readNameGroupCount,0,readsNumInPatch);
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
        vector<C_fastq>().swap(fq1s);
        vector<C_fastq>().swap(fq2s);
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
        vector<C_fastq>().swap(clean_result1);
        vector<C_fastq>().swap(clean_result2);
        delete[] peFilterFlag;
        delete[] filterFlag;
        delete[] readNameGroupCount;
        delete opt2;
    }
    return &peP->bq_check;
}
//这是多线程PE独立IO方法
void* processHts::peDependentIOSubThread(int index){

    vector<C_fastq> fq1s,fq2s;
//    string inputFormat="cram";
    int patch=0;
//    bam1_t* startAln=bam_init1();
    long long curTotalNum=0;
//    writeFile(index,0);
    while(1) {
        //打开写小文件的句柄
        if(MultiThreadsMethod==1) {
            writeFile(index, patch);
        }
        int* readNameGroupCount=new int[readsNumInPatch];
        memset(readNameGroupCount,0,readsNumInPatch);
        int readLineNum=0;
        bam1_t** bamData=readPEData(fq1s,fq2s,readNameGroupCount,readLineNum,index,curTotalNum);
        int readsNum=fq1s.size();
        if(readsNum==0){
            delete[] readNameGroupCount;
            if(MultiThreadsMethod==1) {
                closeHts(threadCleanOut[index]);
//            closeHts(threadFilteredOut[index]);
            }
            threadProgress[index]=patch;
            threadDone[index]=true;
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
        vector<C_fastq>().swap(fq1s);
        vector<C_fastq>().swap(fq2s);
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
        //线程间排队输出到同一文件
//        while(true) {
//            if (index == 0) {
//                if (patch == 0 || threadProgress[gp.threads_num - 1] == patch-1) {
//                    writeMutex.lock();
//                    threadProgress[index]=patch;
//                    writeBackToCram(bamData, peFilterFlag, iter);
//                    patch++;
//                    writeMutex.unlock();
//                    break;
//                }
//            } else {
//                if (threadProgress[index - 1] == patch) {
//                    writeMutex.lock();
//                    threadProgress[index]=patch;
//                    writeBackToCram(bamData, peFilterFlag, iter);
//                    patch++;
//                    writeMutex.unlock();
//                    break;
//                }
//            }
//            usleep(20000);
//        }
        if(MultiThreadsMethod==1){
            writeBackToCram(bamData,peFilterFlag,iter,index);
        }else {
            while (1) {
                int lastIndex = index == 0 ? gp.threads_num - 1 : index - 1;
                if (threadWriteProgress[lastIndex] - threadWriteProgress[index] == 1 ||
                    (index == 0 && threadWriteProgress[index] == threadWriteProgress[lastIndex])) {
                    writeLock.lock();
                    writeBackToCram(bamData, peFilterFlag, iter);
                    writeLock.unlock();
                    threadWriteProgress[index]++;
                    break;
                } else {
                    usleep(100000);
                }
            }
        }
        if(index==0){
            log<<get_local_time()<<"\tprocessed reads: "<<curTotalNum<<endl;
        }
        if(MultiThreadsMethod==1) {
            closeHts(threadCleanOut[index]);
        }
//        closeHts(threadFilteredOut[index]);
        int bamSize=int(readsNumInPatch*2.3);
        for(int i=0;i<bamSize;i++){
            bam_destroy1(bamData[i]);
            bamData[i]=NULL;
        }
        delete[] bamData;
        vector<C_fastq>().swap(clean_result1);
        vector<C_fastq>().swap(clean_result2);
        delete[] peFilterFlag;
        delete[] filterFlag;
        delete[] readNameGroupCount;
        delete opt2;
        if(readsNum<readsNumInPatch){
            threadProgress[index]=patch;
            threadDone[index]=true;
            break;
        }
        threadProgress[index]=patch;
        patch++;
    }
    return &peP->bq_check;
}
void* processHts::seDependentIOSubThread(int index){

    vector<C_fastq> fq1s;
//    string inputFormat="cram";
    int patch=0;
//    bam1_t* startAln=bam_init1();
    long long curTotalNum=0;
    long long processedReads=0;
    while(1) {
        //打开写小文件的句柄
        if(MultiThreadsMethod==1) {
            writeFile(index,patch);
        }
        int readLineNum=0;
        int* readNameGroupCount=new int[readsNumInPatch];
        memset(readNameGroupCount,0,readsNumInPatch);
        bam1_t** bamData=readSEData(fq1s,readNameGroupCount,readLineNum,index,curTotalNum);
        int readsNum=fq1s.size();
        processedReads+=readsNum;
        if(readLineNum<readsNum){
            cerr<<"Error:code error"<<__FILE__<<"\t"<<__LINE__<<endl;
            exit(1);
        }
        if(readsNum==0){
            delete[] readNameGroupCount;
            if(MultiThreadsMethod==1) {
                closeHts(threadCleanOut[index]);
            }
//            closeHts(threadFilteredOut[index]);
            threadProgress[index]=patch;
            threadDone[index]=true;
            break;
        }
        bool* filterFlag=new bool[readsNum];
        memset(filterFlag,0,readsNum);

        SEcalOption opt2;
        vector<C_fastq>  clean_result1;
        opt2.se_local_fs=&(seP->se_local_fs[index]);
        opt2.fq1s=&fq1s;
        opt2.clean_result1=&clean_result1;
        filter_se_fqs(opt2,filterFlag);
        SEstatOption opt_raw,opt_clean;
        opt_raw.fq1s=&fq1s;
        opt_raw.stat1=&(seP->se_local_raw_stat1[index]);
        seP->stat_se_fqs(opt_raw,"raw");		//statistic raw fastqs
        //add_raw_trim(local_raw_stat1[index],local_raw_stat2[index],raw_cut.stat1,raw_cut.stat2);
        vector<C_fastq>().swap(fq1s);
        opt_clean.stat1=&(seP->se_local_clean_stat1[index]);
        opt_clean.fq1s=&clean_result1;
        seP->stat_se_fqs(opt_clean,"clean");
        bool* seFilterFlag=new bool[readLineNum];
        memset(seFilterFlag,0,readLineNum);
        //assgin filter flag to each read line
        int iter=0;
        for(int i=0;i<readsNum;i++){
            for(int j=0;j<readNameGroupCount[i];j++){
                seFilterFlag[iter]=filterFlag[i];
                iter++;
            }
        }
//        cout<<"index:\t"<<index<<"\titer:\t"<<iter<<"\treadLineNum:\t"<<readLineNum<<endl;
        if(iter!=readLineNum){
            cerr<<"Error:code error"<<__FILE__<<"\t"<<__LINE__<<endl;
            exit(1);
        }
        if(MultiThreadsMethod==1){
            writeBackToCram(bamData,seFilterFlag,iter,index);
        }else {
            while (1) {
                int lastIndex = index == 0 ? gp.threads_num - 1 : index - 1;
                if (threadWriteProgress[lastIndex] - threadWriteProgress[index] == 1 ||
                    (index == 0 && threadWriteProgress[index] == threadWriteProgress[lastIndex])) {
                    writeLock.lock();
                    writeBackToCram(bamData, seFilterFlag, iter);
                    writeLock.unlock();
                    threadWriteProgress[index]++;
                    break;
                } else {
                    usleep(100000);
                }
            }
        }
//        writeBackToCram(bamData,seFilterFlag,iter,index);
        if(index==0){
            log<<get_local_time()<<"\tprocessed reads: "<<curTotalNum<<endl;
        }
        if(MultiThreadsMethod==1) {
            closeHts(threadCleanOut[index]);
        }
//        closeHts(threadFilteredOut[index]);
        threadProgress[index]=patch;
        int bamSize=int(readsNumInPatch*1.5);
        for(int i=0;i<bamSize;i++){
            bam_destroy1(bamData[i]);
            bamData[i]=NULL;
        }
        delete[] bamData;
        vector<C_fastq>().swap(clean_result1);
        delete[] seFilterFlag;
        delete[] filterFlag;
        if(readsNum<readsNumInPatch){
            delete[] readNameGroupCount;
            threadProgress[index]=patch;
            threadDone[index]=true;
            break;
        }
        patch++;
    }
    return &seP->se_bq_check;
}
int processHts::readPEData(bam1_t** data,int patch){
    bam1_t* aln=bam_init1();
    int i=0;
    for(;i<patch;i++){
        if(sam_read1(in,header,aln)>=0){
            if(bam_copy1(data[i], aln)==NULL){
                cerr<<BAMCPERR<<endl;
                exit(1);
            }
        }else{
            break;
        }
    }
    bam_destroy1(aln);
    return i;
}
void processHts::writeFile(){
    if(gp.fq2_path.rfind(".bam")==gp.fq2_path.size()-4){
        writeBam();
    }else if(gp.fq2_path.rfind(".sam")==gp.fq2_path.size()-4){
        writeSam();
    }else if(gp.fq2_path.rfind(".cram")==gp.fq2_path.size()-5){
        writeCram();
    }else{
        cerr<<"Error:only support sam/bam/cram suffix format output"<<endl;
        exit(1);
    }
}
void processHts::writeFile(int index,int patch){
    if(gp.fq2_path.rfind(".bam")==gp.fq2_path.size()-4){
        writeBam(index,patch);
    }else if(gp.fq2_path.rfind(".sam")==gp.fq2_path.size()-4){
        writeSam(index,patch);
    }else if(gp.fq2_path.rfind(".cram")==gp.fq2_path.size()-5){
        writeCram(index,patch);
    }else{
        cerr<<"Error:only support sam/bam/cram suffix format output"<<endl;
        exit(1);
    }
}
void processHts::closeThreadFd(int index){
    closeHts(threadCleanOut[index]);
    closeHts(threadFilteredOut[index]);
}
void* processHts::se_sub_thread(int index){
//    seP->of_log<<get_local_time()<<"\tthread "<<index<<" start"<<endl;
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
        vector<C_fastq>().swap(fq1s);
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
        vector<C_fastq>().swap(clean_result1);
    }

    return &seP->se_bq_check;
}
void processHts::clean(){
    log << get_local_time() << "\tAnalysis accomplished!" << endl;
    log.close();
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    delete[] threadWriteProgress;
    if(gp.threads_num>1){
        for (int i = 0; i < gp.threads_num; i++) {
            closeHts(threadIn[i]);
        }
        if(outputFormat=="bam") {
            if(MultiThreadsMethod==1) {
                if (bgzf_close(cleanBam) < 0) {
                    cerr << "Error:fail to close clean bam" << endl;
                    exit(1);
                }
            }else{
                closeHts(cleanOut);
//                closeHts(in);
            }
//            if (bgzf_close(filteredBam) < 0) {
//                cerr << "Error:fail to close filtered bam" << endl;
//                exit(1);
//            }
//        }else{
//            closeHts(cleanOut);
//            closeHts(filteredOut);
        }else{
            if(MultiThreadsMethod==2){
                closeHts(cleanOut);
            }
        }
        if(MultiThreadsMethod==1) {
            string rmCmd = "rm -rf " + tmpDir;
            if (system(rmCmd.c_str()) == -1) {
                cerr << "Error:when running " << rmCmd << endl;
                exit(1);
            }
        }
//        closeHts(cleanOut);
    }else {
        closeHts(cleanOut);
//        closeHts(filteredOut);
        closeHts(in);
    }
}
void processHts::prepare(){
    mkDir(gp.output_dir);
//    string mkdir_str="mkdir -p "+gp.output_dir;
//    if(system(mkdir_str.c_str())==-1){
//        cerr<<"Error:mkdir fail"<<endl;
//        exit(1);
//    }
    log.open(gp.log);
    if(!log){
        cerr<<"Error:cannot write to such file,"<<gp.log<<endl;
        exit(-1);
    }
    openHts();
    string format;
    switch(inputFormat->format){
        case 4:format="bam";break;
        case 6:format="cram";break;
        default:{
            cerr<<"Error:only support BAM/CRAM in this module"<<endl;
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
            if (sam_write1(cleanOut, header, data[i]) < 0) {
                cerr << "Error:write file error" << endl;
                exit(1);
            }
        }else {
            if (sam_write1(cleanOut, header, data[i]) < 0) {
                cerr << "Error:write file error" << endl;
                exit(1);
            }
        }
    }
}
void processHts::writeBackToCram(bam1_t** data,bool* filterFlag,int size,int index){
    for(int i=0;i<size;i++){
        if(filterFlag[i]){
            if((data[i]->core.flag & 512)!=1) {
                data[i]->core.flag += 512;
            }
            if (sam_write1(threadCleanOut[index], header, data[i]) < 0) {
                cerr << "Error:write file error" << endl;
                exit(1);
            }
        }else {
            if (sam_write1(threadCleanOut[index], header, data[i]) < 0) {
                cerr << "Error:write file error" << endl;
                exit(1);
            }
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
        //pe_fastq_filter.pe_trim(gp);
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
//        if(!gp.trim_fq1.empty()){
//            peP->preOutput(1,pe_fastq_filter.fq1);
//            peP->preOutput(2,pe_fastq_filter.fq2);
//            opt->trim_result1->emplace_back(pe_fastq_filter.fq1);
//            opt->trim_result2->emplace_back(pe_fastq_filter.fq2);
//        }
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
    bam1_t** data=new bam1_t*[readsNumInPatch*2];
    for(int i=0;i<readsNumInPatch*2;i++){
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
            cerr<<"Error:bam copy error"<<endl;
            exit(1);
        }
        char* seq=get_read(aln);
        char* qual=get_quality(aln);
        char* qname=bam_get_qname(aln);
        if(seq==NULL || qual==NULL || qname==NULL){
            cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
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
        seq=NULL;
        delete[] qual;
        qual=NULL;
        if(bam_copy1(data[lineNum],aln)==NULL){
            cerr<<"Error:copy data in memory failed"<<endl;
            exit(1);
        }
        lastQname=qnameS;
        lineNum++;
    }
    while(sam_read1(in,header,aln)>=0){
//        aln->core.seq;
        char* seq=get_read(aln);
        char* qual=get_quality(aln);
        char* qname=bam_get_qname(aln);
        if(seq==NULL || qual==NULL || qname==NULL){
            cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
            exit(1);
        }
        string qnameS=qname;

        if(bam_copy1(data[lineNum],aln)==NULL){
            cerr<<"Error:copy data in memory failed"<<endl;
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
            if(lineNum+10>=readsNumInPatch*2){
                if(bam_copy1(lastLine,aln)==NULL){
                    cerr<<"Error:bam copy error"<<endl;
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
            if(bam_is_read2(aln) && !scanRead2) {
                fastq2.seq_id.assign(qname);
                fastq2.sequence.assign(seq);
                fastq2.qual_seq.assign(qual);
                scanRead2 = true;
            }
        }
        delete[] seq;
        seq=NULL;
        delete[] qual;
        qual=NULL;
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
//bamBlock* processHts::readPEData(bam1_t* lastLine) {
//    bamBlock* bB=new bamBlock(readsNumInPatch);
//    bam1_t** data=bB->getRecords();
//    int* readIDGroupCount=bB->getReadGroupIndex();
//    int readReadsNum=0;//最后实际读取的reads数目
//    int linesOfOneReadID=0;
//    aln=bam_init1();
//    string lastQname="";
//    bool scanRead1=false;
//    bool scanRead2=false;
//    set<string> read1ID,read2ID;
//    //读取的行数，即records数目
//    int lineNum=0;
//    //处理上次未处理的最后一行，如果htslib支持回滚就简单了
//    if(bam_get_qname(lastLine)!=NULL){
//        linesOfOneReadID++;
//        if(bam_copy1(aln,lastLine)==NULL){
//            cerr<<"Error:bam copy error"<<endl;
//            exit(1);
//        }
//        char* qname=bam_get_qname(aln);
//        if(qname==NULL){
//            cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
//            exit(1);
//        }
//        string qnameS=qname;
//        if(bam_copy1(data[lineNum],aln)==NULL){
//            cerr<<"Error:copy data in memory failed"<<endl;
//            exit(1);
//        }
//        lastQname=qnameS;
//        lineNum++;
//    }
//    while(sam_read1(in,header,aln)>=0){
////        aln->core.seq;
//        char* qname=bam_get_qname(aln);
//        if(qname==NULL){
//            cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
//            exit(1);
//        }
//        string qnameS=qname;
//        if(lastQname!=qnameS && lastQname!="" && lineNum>1){
//            readIDGroupCount[readReadsNum]=linesOfOneReadID;
//            linesOfOneReadID=0;
//            if(lineNum+10>=readsNumInPatch*2){
//                if(bam_copy1(lastLine,aln)==NULL){
//                    cerr<<"Error:bam copy error"<<endl;
//                    exit(1);
//                }
//                lineNum++;
//                readReadsNum++;
//                bB->setRecordsNum(lineNum);
//                bB->setReadsNum(readReadsNum);
//                return bB;
//            }
//            readReadsNum++;
//        }
//        if(bam_copy1(data[lineNum],aln)==NULL){
//            cerr<<"Error:copy data in memory failed"<<endl;
//            exit(1);
//        }
//        lastQname=qnameS;
//        lineNum++;
//        linesOfOneReadID++;
//    }
//    bB->setRecordsNum(lineNum);
//    bB->setReadsNum(readReadsNum);
//    return bB;
//}
bam1_t** processHts::readPEData(vector<C_fastq> &fq1s,vector<C_fastq> &fq2s,int* IDstat,int &readLineNum,int index,long long& totalNum) {
    C_fastq fastq1,fastq2;
    //考虑到supplementary比对情况，实际数据量要大于2被reads数量，如果碰到极端情况（supplementary比对非常多），则在后面申请个更大的数组
    int bamSize=int(readsNumInPatch*2.3);
    bam1_t** data=new bam1_t*[bamSize];
    for(int i=0;i<bamSize;i++){
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
    bam1_t *threadAln=bam_init1();
    while(sam_read1(threadIn[index],header,threadAln)>=0){
        char* qname=bam_get_qname(threadAln);
        if(qname==NULL){
            cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
            exit(1);
        }
        string qnameS=qname;
        if(lastQname!=qnameS && lastQname!="") {
            totalNum++;
        }
        //判断是否属于该线程需要处理的数据
        if((totalNum/readsNumInPatch)%gp.threads_num==index){
            //判断上一对reads是否读取完毕
            char* seq=get_read(threadAln);
            char* qual=get_quality(threadAln);
            if(seq==NULL || qual==NULL){
                cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
                exit(1);
            }
            if(bam_copy1(data[lineNum],threadAln)==NULL){
                cerr<<"Error:copy data in memory failed"<<endl;
                exit(1);
            }
            if(lastQname!=qnameS && lastQname!="" && lineNum>1){
                fq1s.emplace_back(fastq1);
                fq2s.emplace_back(fastq2);
                IDstat[iter]=readIDNum;
                readIDNum=0;
                if(iter==readsNumInPatch-1){
                    bam_destroy1(threadAln);
                    return data;
                }
                scanRead1=false;
                scanRead2=false;
                iter++;
            }
            //还未读到一条新read
            if(bam_is_read1(threadAln) && !scanRead1){
                fastq1.seq_id.assign(qname);
                fastq1.sequence.assign(seq);
                fastq1.qual_seq.assign(qual);
                scanRead1=true;
            }else{
                if(bam_is_read2(threadAln) && !scanRead2) {
                    fastq2.seq_id.assign(qname);
                    fastq2.sequence.assign(seq);
                    fastq2.qual_seq.assign(qual);
                    scanRead2 = true;
                }
            }
            delete[] seq;
            seq=NULL;
            delete[] qual;
            qual=NULL;
            lineNum++;//bam1_t数组的索引
            readLineNum=lineNum;
            readIDNum++; //每一个readID对应的比对记录数目
        }else{
            //若不是要处理的数据
            if(fq1s.size()>0){
                bam_destroy1(threadAln);
                break;
            }
        }
        lastQname=qnameS;
    }
    if(fastq1.sequence.size()>0 && fastq2.sequence.size()) {
        fq1s.emplace_back(fastq1);
        fq2s.emplace_back(fastq2);
        IDstat[iter] = readIDNum;
    }
//    bam_destroy1(threadAln);
    return data;
}
bam1_t** processHts::readSEData(vector<C_fastq> &fq1s,int* IDstat,int &readLineNum,int index,long long& totalNum) {
    C_fastq fastq1;
    //考虑到supplementary比对情况，实际数据量要大于2被reads数量，如果碰到极端情况（supplementary比对非常多），则在后面申请个更大的数组
    //todo
    int bamSize=int(readsNumInPatch*1.5);
    bam1_t** data=new bam1_t*[bamSize];
    for(int i=0;i<bamSize;i++){
        data[i]=bam_init1();
    }
    int iter=0;
    seP->C_fastq_init(fastq1);
    aln=bam_init1();
    string lastQname=seLastID[index];
//    if(totalNum>0){
//        lastQname.assign(bam_get_qname(lastLine));
//    }
//    bool scanRead1=false;
    int lineNum=0;
    bam1_t *threadAln=bam_init1();
//    bool test=true;
//    cout<<"index:\t"<<index<<"\tcurrent totalNum:\t"<<totalNum<<"\t"<<(totalNum/readsNumInPatch)%gp.threads_num<<endl;
    while(sam_read1(threadIn[index],header,threadAln)>=0){
        char* qname=bam_get_qname(threadAln);
        if(qname==NULL){
            cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
            exit(1);
        }
        string qnameS=qname;
        if(lastQname!=qnameS) {
            totalNum++;
        }
        //判断是否属于该线程需要处理的数据
        if(((totalNum-1)/readsNumInPatch)%gp.threads_num==index){
            //判断上一对reads是否读取完毕
            char* seq=get_read(threadAln);
            char* qual=get_quality(threadAln);
            if(seq==NULL || qual==NULL){
                cerr<<"Error:parse "<<inputFormatString<<" file error"<<endl;
                exit(1);
            }
            if(bam_copy1(data[lineNum],threadAln)==NULL){
                cerr<<"Error:copy data in memory failed"<<endl;
                exit(1);
            }
            lineNum++;//bam1_t数组的索引
            readLineNum=lineNum;
            fastq1.seq_id.assign(qname);
            fastq1.sequence.assign(seq);
            fastq1.qual_seq.assign(qual);
            if(lastQname!=qnameS){
                fq1s.emplace_back(fastq1);
                IDstat[iter] = 1;
                iter++;
            }else{
                IDstat[iter-1]++;
            }
            //还未读到一条新read
            delete[] seq;
            seq=NULL;
            delete[] qual;
            qual=NULL;
             //每一个readID对应的比对记录数目
        }else{
            //若不是要处理的数据
            if(fq1s.size()>0){
                bam_destroy1(threadAln);
                lastQname=qnameS;
                break;
            }
        }
        lastQname=qnameS;
    }
    if(fq1s.size()>0) {
        if (fastq1.sequence.size() == 0) {
            cerr << "code error" << endl;
            exit(1);
        }
    }
//    bam_destroy1(threadAln);
    seLastID[index]=lastQname;
    return data;
}
bam1_t** processHts::readSEData(int& lineNum) {
    bam1_t** data=new bam1_t*[lineNumPerThread*gp.threads_num];
    for(int i=0;i<lineNumPerThread*gp.threads_num;i++){
        data[i]=bam_init1();
    }
    int iter=0;
//    seP->C_fastq_init(fastq1);
    aln=bam_init1();
    while(sam_read1(in,header,aln)>=0){
        lineNum++;
//        aln->core.seq;

//        char* seq=get_read(aln);
//        char* qual=get_quality(aln);
//        char* qname=bam_get_qname(aln);
//        if(seq==NULL || qual==NULL || qname==NULL){
//            cerr<<"Error:parse cram file error"<<endl;
//            exit(1);
//        }
//        fastq1.seq_id.assign(qname);
//        fastq1.sequence.assign(seq);
//        fastq1.qual_seq.assign(qual);
//        delete[] seq;
//        delete[] qual;
//        delete qname;
//        fq1s.emplace_back(fastq1);
        if(bam_copy1(data[iter],aln)==NULL){
            cerr<<"Error:copy data in memory failed"<<endl;
            exit(1);
        }
        if (lineNum == lineNumPerThread*gp.threads_num) {
            return data;
            break;
        }
        iter++;
    }
//    for(int i=fq1s.size();i<readsNumInPatch;i++){
//        data[i]=NULL;
//    }
    return data;
}
bam1_t** processHts::readSEData(vector<C_fastq> &fq1s) {
    bam1_t** data=new bam1_t*[lineNumPerThread*gp.threads_num];
    for(int i=0;i<lineNumPerThread*gp.threads_num;i++){
        data[i]=bam_init1();
    }
    int iter=0;
    C_fastq fastq1;
    seP->C_fastq_init(fastq1);
    aln=bam_init1();
    while(sam_read1(in,header,aln)>=0){
//        aln->core.seq;

        parseRead(aln,fastq1);
        fq1s.emplace_back(fastq1);
        if(bam_copy1(data[iter],aln)==NULL){
            cerr<<"Error:copy data in memory failed"<<endl;
            exit(1);
        }
        if (fq1s.size() == readsNumInPatch) {
            return data;
            break;
        }
        iter++;
    }
    for(int i=fq1s.size();i<readsNumInPatch;i++){
        data[i]=NULL;
    }
    return data;
}
void processHts::openHts(){
    in=hts_open(gp.fq1_path.c_str(),"rc");
    if(in==NULL){
        cerr<<"Error:cannot open such file,"<<inputCram<<endl;
        exit(1);
    }
    if(gp.fq1_path.rfind(".cram")==gp.fq1_path.size()-5 || gp.fq2_path.rfind(".cram")==gp.fq2_path.size()-5) {
        string refPath = gp.reference;
        string refFai = refPath + ".fai";
        if (hts_set_fai_filename(in, refFai.c_str()) < 0) {
            cerr << "Error:reference is needed, cannot open such file," << refFai << endl;
            exit(1);
        }
    }
    inputFormat=hts_get_format(in);
    if(inputFormat->category==1){
        if(inputFormat->format!=4 && inputFormat->format!=6){
            cerr<<"Error:only support BAM/CRAM in this module"<<endl;
            exit(1);
        }
    }else{
        cerr<<"Error:only support BAM/CRAM in this module"<<endl;
        exit(1);
    }
    header=sam_hdr_read(in);
    if(header==NULL){
        cerr<<"Error:get header fail"<<endl;
        exit(1);
    }

    aln=bam_init1();

    if(aln==NULL){
        cerr<<"Error:init error,"<<__FILE__<<":line"<<__LINE__<<endl;
        exit(1);
    }
}
void processHts::openHts(int index){
    threadIn[index]=hts_open(gp.fq1_path.c_str(),"r");
    if(threadIn[index]==NULL){
        cerr<<"Error:cannot open such file,"<<inputCram<<endl;
        exit(1);
    }
    if(gp.fq1_path.rfind(".cram")==gp.fq1_path.size()-5 || gp.fq2_path.rfind(".cram")==gp.fq2_path.size()-5) {
        string refPath = gp.reference;
        string refFai = refPath + ".fai";
        if (hts_set_fai_filename(threadIn[index], refFai.c_str()) < 0) {
            cerr << "Error:reference is needed, cannot open such file," << refFai << endl;
            exit(1);
        }
    }
    inputFormat=hts_get_format(threadIn[index]);
    if(inputFormat->category==1){
        if(inputFormat->format!=4 && inputFormat->format!=6){
            cerr<<"Error:only support BAM/CRAM in this module"<<endl;
            exit(1);
        }
    }else{
        cerr<<"Error:only support BAM/CRAM in this module"<<endl;
        exit(1);
    }
    header=sam_hdr_read(threadIn[index]);
    if(header==NULL){
        cerr<<"Error:get header fail"<<endl;
        exit(1);
    }
    aln=bam_init1();
    if(aln==NULL){
        cerr<<"Error:init error,"<<__FILE__<<":line"<<__LINE__<<endl;
        exit(1);
    }
}
void processHts::writeCram(){
    string outPath=gp.output_dir+"/"+gp.fq2_path;
    cleanOut=sam_open(outPath.c_str(),"wc");
    if(cleanOut==NULL){
        cerr<<"Error:cannot write to such file,"<<outPath<<endl;
        exit(1);
    }
    string refPath=gp.reference;
    string refFai=refPath+".fai";
    if(hts_set_fai_filename(cleanOut,refFai.c_str())<0){
        cerr<<"Error:cannot open such file,"<<refFai<<endl;
        exit(1);
    }
    if(sam_hdr_write(cleanOut,header)==-1){
        cerr<<"Error:write header to cram failed"<<endl;
        exit(1);
    }
//    string outPath2=gp.output_dir+"/_failQC.cram";
//    filteredOut=sam_open(outPath2.c_str(),"wc");
//    if(filteredOut==NULL){
//        cerr<<"Error:cannot write to such file,"<<outPath2<<endl;
//        exit(1);
//    }
//    if(hts_set_fai_filename(filteredOut,refFai.c_str())<0){
//        cerr<<"Error:cannot open such file,"<<refFai<<endl;
//        exit(1);
//    }
//    if(sam_hdr_write(filteredOut,header)==-1){
//        cerr<<"Error:write header to cram failed"<<endl;
//        exit(1);
//    }
}
void processHts::writeBam(){
    string outPath=gp.output_dir+"/"+gp.fq2_path;
    cleanOut=sam_open(outPath.c_str(),"wb");
    if(cleanOut==NULL){
        cerr<<"Error:cannot write to such file,"<<outPath<<endl;
        exit(1);
    }
    if(sam_hdr_write(cleanOut,header)<0){
        cerr<<"Error:write file error"<<endl;
        exit(1);
    }
//    string outPath2=gp.output_dir+"/_failQC.bam";
//    filteredOut=sam_open(outPath2.c_str(),"wb");
//    if(filteredOut==NULL){
//        cerr<<"Error:cannot write to such file,"<<outPath2<<endl;
//        exit(1);
//    }
//    if(sam_hdr_write(filteredOut,header)<0){
//        cerr<<"Error:write file error"<<endl;
//        exit(1);
//    }
}
void processHts::writeSam(){
    string outPath=gp.output_dir+"/"+gp.fq2_path;
    cleanOut=sam_open(outPath.c_str(),"w");
    if(cleanOut==NULL){
        cerr<<"Error:cannot write to such file,"<<gp.fq2_path<<endl;
        exit(1);
    }
    if(sam_hdr_write(cleanOut,header)<0){
        cerr<<"Error:write file error"<<endl;
        exit(1);
    }
//    string outPath2=gp.output_dir+"/_failQC.sam";
//    filteredOut=sam_open(outPath.c_str(),"w");
//    if(filteredOut==NULL){
//        cerr<<"Error:cannot write to such file,"<<outPath2<<endl;
//        exit(1);
//    }
//    if(sam_hdr_write(filteredOut,header)<0){
//        cerr<<"Error:write file error"<<endl;
//        exit(1);
//    }
}

void processHts::writeCram(int index,int patch){
    string outPath=tmpDir+"/"+gp.fq2_path+".t"+to_string(index)+"."+to_string(patch);
    threadCleanOut[index]=sam_open(outPath.c_str(),"wc");
    if(threadCleanOut[index]==NULL){
        cerr<<"Error:cannot write to such file,"<<outPath<<endl;
        exit(1);
    }
    string refPath=gp.reference;
    string refFai=refPath+".fai";
    if(hts_set_fai_filename(threadCleanOut[index],refFai.c_str())<0){
        cerr<<"Error:cannot open such file,"<<refFai<<endl;
        exit(1);
    }
    if(sam_hdr_write(threadCleanOut[index],header)==-1){
        cerr<<"Error:write header to cram failed"<<endl;
        exit(1);
    }
//    string outPath2=tmpDir+"/_failQC.cram"+".t"+to_string(index)+"."+to_string(patch);
//    threadFilteredOut[index]=sam_open(outPath2.c_str(),"wc");
//    if(threadFilteredOut[index]==NULL){
//        cerr<<"Error:cannot write to such file,"<<outPath2<<endl;
//        exit(1);
//    }
//    if(hts_set_fai_filename(threadFilteredOut[index],refFai.c_str())<0){
//        cerr<<"Error:cannot open such file,"<<refFai<<endl;
//        exit(1);
//    }
//    if(sam_hdr_write(threadFilteredOut[index],header)==-1){
//        cerr<<"Error:write header to cram failed"<<endl;
//        exit(1);
//    }
}
void processHts::writeBam(int index,int patch){
    string outPath=tmpDir+"/"+gp.fq2_path+".t"+to_string(index)+"."+to_string(patch);
    threadCleanOut[index]=sam_open(outPath.c_str(),"wb");
    if(threadCleanOut[index]==NULL){
        cerr<<"Error:cannot write to such file,"<<outPath<<endl;
        exit(1);
    }
    if(sam_hdr_write(threadCleanOut[index],header)<0){
        cerr<<"Error:write file error"<<endl;
        exit(1);
    }
//    string outPath2=tmpDir+"/_failQC.bam"+".t"+to_string(index)+"."+to_string(patch);
//    threadFilteredOut[index]=sam_open(outPath2.c_str(),"wb");
//    if(threadFilteredOut[index]==NULL){
//        cerr<<"Error:cannot write to such file,"<<outPath2<<endl;
//        exit(1);
//    }
//    if(sam_hdr_write(threadFilteredOut[index],header)<0){
//        cerr<<"Error:write file error"<<endl;
//        exit(1);
//    }
}
void processHts::writeSam(int index,int patch){
    string outPath=tmpDir+"/"+gp.fq2_path+".t"+to_string(index)+"."+to_string(patch);
    threadCleanOut[index]=sam_open(outPath.c_str(),"w");
    if(threadCleanOut[index]==NULL){
        cerr<<"Error:cannot write to such file,"<<gp.fq2_path<<endl;
        exit(1);
    }
    if(sam_hdr_write(threadCleanOut[index],header)<0){
        cerr<<"Error:write file error"<<endl;
        exit(1);
    }
//    string outPath2=tmpDir+"/_failQC.sam"+".t"+to_string(index)+"."+to_string(patch);
//    threadFilteredOut[index]=sam_open(outPath.c_str(),"w");
//    if(threadFilteredOut[index]==NULL){
//        cerr<<"Error:cannot write to such file,"<<outPath2<<endl;
//        exit(1);
//    }
//    if(sam_hdr_write(threadFilteredOut[index],header)<0){
//        cerr<<"Error:write file error"<<endl;
//        exit(1);
//    }
}

void processHts::closeHts(htsFile* fp){
    if(hts_close(fp)<0){
        cerr<<"Error:cannot close cram file"<<endl;
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
char *processHts::get_read(bam1_t *aln)
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
char* processHts::get_quality(bam1_t *aln)
{
    char* qual_out = new char[ aln->core.l_qseq + 1];
    memset(qual_out,0,aln->core.l_qseq + 1);
    char *q = (char *)bam_get_qual(aln);
    int n;

    if (!qual_out) return NULL;

    if (*q == '\xff') {
        delete[] qual_out;
        qual_out = NULL;
        return NULL;
    }

    for (n=0; n < aln->core.l_qseq; n++) {
        qual_out[n] = q[n]+33;
    }
    if (aln->core.flag & BAM_FREVERSE) reverse(qual_out);
    return qual_out;
}
