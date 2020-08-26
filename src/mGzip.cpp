//
// Created by berry on 2020-08-11.
//

#include "mGzip.h"
using namespace::std;
bool mGzip::check_mGzip(string filePath) {
    bool whether_mGzip=false;
//    gzFile in=gzopen(filePath.c_str(),"rb");
    FILE* in=fopen(filePath.c_str(),"rb");
    if(in==NULL){
        cerr<<"Error, fail to open file:"<<filePath<<endl;
        exit(1);
    }
    int bufSize=100;
    char buf[bufSize];
    if(fread(buf,1,14,in)==14){
        if(!(Byte(buf[0])==0x1f && Byte(buf[1])==0x8b)){
            fclose(in);
            return whether_mGzip;
        }
        if(buf[12]=='I' && buf[13]=='G'){
            whether_mGzip=true;
        }else{
            whether_mGzip=false;
        }
    }
    fclose(in);

//    gzFile input=gzopen(filePath.c_str(),"rb");
//    char* test=new char[1000];
//    int iter=0;
//    while(gzgets(input,test,1000)!=NULL){
//        cout<<test<<endl;
//        iter++;
//    }
//    cout<<iter<<endl;
    return whether_mGzip;
}

vector<threadDataInfo> mGzip::allocate(int threadsNum, vector<string> filesPath) {
    vector<threadDataInfo> threadNeededInfo;
    for(int fileIndex=0;fileIndex<filesPath.size();fileIndex++){
        string filePath=filesPath[fileIndex];
        FILE* in=fopen(filePath.c_str(),"rb");
        if(in==NULL){
            cerr<<"Error, fail to open file:"<<filePath<<endl;
            exit(1);
        }
        //span file
        long curOffSet=0;
        int totalBlocks=0;

        vector<long> blockStarts;
        vector<long> spanLength;
//        vector<long> compressSize;
        vector<long> rawSize;
        vector<u_int> crc32;
        string filename;
        int method,flag;
        while(getOneBlock(in,method,flag,curOffSet,blockStarts,spanLength,rawSize,crc32,filename)>0){
            totalBlocks++;
        }
        //allocate blocks to thread
        int* blocksNumToThreads=new int[threadsNum];
        for(int i=0;i<threadsNum;i++){
            blocksNumToThreads[i]=0;
        }
        for(int i=0;i<totalBlocks;i++){
            int thread_idx=i%threadsNum;
            blocksNumToThreads[thread_idx]++;
        }
        int curBlockIndex=0;
        for(int i=0;i<threadsNum;i++){
            threadDataInfo curThreadInfo;
            curThreadInfo.index=i;
            vector<u_int> curFileCrc32;
            vector<long> curFileRawBlockSize;
            vector<long> curFileSeekOffSetStart;
            vector<long> curFileSpanLength;
            for(int j=0;j<blocksNumToThreads[i];j++){
                if(curThreadInfo.fileList.size()<=fileIndex){
                    curThreadInfo.rawFileName.push_back(filename);
                    curThreadInfo.method.push_back(method);
                    curThreadInfo.flag.push_back(flag);
                    curThreadInfo.fileList.push_back(filePath);
                }
                curFileCrc32.push_back(crc32[curBlockIndex]);
                curFileRawBlockSize.push_back(rawSize[curBlockIndex]);
                curFileSeekOffSetStart.push_back(blockStarts[curBlockIndex]);
                curFileSpanLength.push_back(spanLength[curBlockIndex]);
                curBlockIndex++;
                if(curBlockIndex>=totalBlocks){
                    break;
                }
            }
            curThreadInfo.crc32.push_back(curFileCrc32);
            curThreadInfo.rawBlockSize.push_back(curFileRawBlockSize);
            curThreadInfo.seekOffSetStart.push_back(curFileSeekOffSetStart);
            curThreadInfo.spanLength.push_back(curFileSpanLength);
            threadNeededInfo.push_back(curThreadInfo);
        }
    }
    return threadNeededInfo;
}

int mGzip::getOneBlock(FILE *pFile,  int& method, int& flag, long &curOffSet, vector<long> &blockStarts, vector<long> &spanLength, vector<long> &rawSize,
                        vector<u_int> &crc32, string &fileName) {
//    uLongf testSize=441+30+4;
//    unsigned char* readBuf=new unsigned char[testSize];
//    fread(readBuf,1,testSize,pFile);
//    uLongf decomLen=1000000;
//    Bytef* decom=new Bytef[decomLen];
//    uncompress(decom,&decomLen,readBuf,testSize);
    int firstReadSize=14;
    char* buf=new char[firstReadSize];
    if(fread(buf,1,firstReadSize,pFile)==firstReadSize){
        if(!(buf[12]=='I' && buf[13]=='G')){
            cerr<<"Error, input is not with mGzip format"<<endl;
            exit(1);
        }
    }else{
        return -1;
    }
    method=buf[2];
    flag=buf[3];
    curOffSet+=firstReadSize;
    int extraLen=(Byte(buf[11])<<8)+buf[10];
    delete[] buf;
    int secReadSize=extraLen-2;
    buf=new char[secReadSize];
    if(fread(buf,1,secReadSize,pFile)==secReadSize){
        long blockSize=0;
        for(int i=secReadSize;i>2;i--){
            blockSize+=Byte(buf[i-1])<<(8*(i-3));
        }
        delete[] buf;
        buf=new char[1];
        //get raw file name
        int nextLen=0;
        curOffSet+=secReadSize;
        bool getFileName=fileName==""?true:false;
        while(fread(buf,1,1,pFile)==1){
            nextLen++;
            if(Byte(buf[0])==0x00){
                break;
            }else{
                if(getFileName)
                    fileName += char(buf[0]);
            }
        }
        curOffSet+=nextLen;
        delete[] buf;
        long compressBlockSize=blockSize-firstReadSize-secReadSize-nextLen-8;
        blockStarts.push_back(curOffSet);
        spanLength.push_back(compressBlockSize);
        //
        const char* hello="@CL100012597L1C001R002_0/1\n"
                          "TGAGAACCCAGGGGTGAGCATTGATTCTTTGTTCCCCAACAACCTTCTTC\n"
                          "+\n"
                          "FEFBFGFFFGGFFFEFAFFFFEFAFAFF:FFBDDGFFFEFEEDEF2FF-F\n"
                          "@CL100012597L1C001R002_1/1\n"
                          "TATGAAACCATATGGGAAGGGAGGAAGCAGGTAGAGTATAGATTAACCAA\n"
                          "+\n"
                          "FGFFFFFFFGGFFFEFFFCFDAE:DBFEFFFDFF=GB9'FF1E7D(*F-,\n"
                          "@CL100012597L1C001R002_4/1\n"
                          "GACATTTACTGTGGATCCAGCACTATTCCAAGTGCTCTTTAGGCCCCGAG\n"
                          "+\n"
                          "GBFGFFFFBFEFGFDFFFBFFFF<EE?G?B;E/;F+:F>/59B'104*2&\n"
                          "@CL100012597L1C001R002_6/1\n"
                          "GGCCGCAGCTGGGACTTCACTGCTGAGTGGAGGGAGGGCAGCGGGCCACG\n"
                          "+\n"
                          "?FFFFFFFFFF:FEFFFF;AEABFECE;E=2B86,E2DA:C7+'F1F'F;\n"
                          "@CL100012597L1C001R002_7/1\n"
                          "TGAATCACTAAGACAACTCCAATTCTGTTACTAGAAGGAATTAGTTTGGG\n"
                          "+\n"
                          "FGFFFFFFFGFGFDGFGFGEF0EFF;GFEGEF+FFD4/9EE.39FEF-F&\n"
                          "@CL100012597L1C001R002_8/1\n"
                          "CTTCCTGTTTTGTAGATCTTGGTCTTCTCATTGTATTCTCACATAGCAGG\n"
                          "+\n"
                          "FDFF;FFFFFFFFFFEFFEFF>>FF7FE?@F4FC,E,9>@A=>D1EF>DB\n"
                          "@CL100012597L1C001R002_9/1\n"
                          "GGTGCCCTTATAAGAAGCAGCAGCAACACAAAACATCTCTCTCCTTCATG\n"
                          "+\n"
                          "FFFFFFEFFFFFFFFFFFEFFFFEFFFE@EFBEA@FFF<F57(:';F&0@\n"
                          "@CL100012597L1C001R002_10/1\n"
                          "TTTCAAAATATATTATAAAGCTACTGTAATCAGGACAGCATAGTATGGGA";
        Bytef* com=new Bytef[compressBlockSize];
        Bytef* decom=new Bytef[compressBlockSize];
        uLongf comprLen=compressBlockSize;
        uLong len =  1;
        Bytef* dest3=new Bytef[compressBlockSize];
        compress(com, &comprLen, (const Bytef*)hello, len);
        uncompress(decom, &len, com, comprLen);
        Bytef* source=new Bytef[compressBlockSize];
        if(fread(source,1,compressBlockSize,pFile)==compressBlockSize){
            uLong destSize=1000000;
            unsigned char *dest=new unsigned char[destSize];
            uLong sourceSize=compressBlockSize;
            uncompress(dest, &destSize,com, comprLen);
            uncompress(dest, &destSize,source, sourceSize);

//            cout<<"here"<<endl;

        }
//        fseek(pFile,compressBlockSize,SEEK_CUR);
        curOffSet+=compressBlockSize;
        int blockTailLen=8;
        buf=new char[blockTailLen];
        if(fread(buf,1,blockTailLen,pFile)==blockTailLen){
            u_int crc=0;
            for(int i=3;i>=0;i--){
                crc+=Byte(buf[i])<<(blockTailLen*i);
            }
            crc32.push_back(crc);
            int rawLength=0;
            for(int i=7;i>=4;i--){
                rawLength+=Byte(buf[i])<<((i-4)*blockTailLen);
            }
            rawSize.push_back(rawLength);
        }else{
            delete[] buf;
            return -1;
        }
        delete[] buf;
        curOffSet+=blockTailLen;
    }else{
        return -1;
    }
    return 1;
}
