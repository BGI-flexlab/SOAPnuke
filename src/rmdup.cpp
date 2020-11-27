//
// Created by berry on 2020-10-20.
//


#include "rmdup.h"


rmdup::rmdup(uint64_t* data, uint64_t size) {
    this->data=data;
    this->size=size;
}

void rmdup::markDup(bool *dupFlag) {
    int prime=getPrime();
    int* count=new int[prime];
    memset(count,0,sizeof(int)*prime);
//    unordered_set<string> tmpset;
    int contentSize=0;
    int multisize=0;
//    vector<uint64_t> dupIdx,dupIdx2;
//    unordered_set<uint64_t> tmpset2,dupset;
    for(uint64_t i=0;i<size;i++){
//        string curStr=md52string(data[i]);
//        if(tmpset.find(curStr)!=tmpset.end()){
//            checkDupNum++;
////            dupset.insert(curStr);
//            cout<<"here\t"<<curStr<<endl;
//            dupIdx.push_back(i);
//        }else {
//            tmpset.insert(curStr);
//        }
        int idx=data[i]%prime;
        count[idx]++;
    }
//    cout<<checkDupNum<<endl;
    uint32_t* multiIdx=new uint32_t[prime];
    memset(multiIdx,0,sizeof(uint32_t)*prime);
    int tmpIter=0;
    for(int i=0;i<prime;i++){
        if(count[i]>1){
            multiIdx[i]=tmpIter;
            tmpIter++;
            contentSize++;
            multisize+=count[i];
        }
    }
    uint64_t** dataContent=new uint64_t*[contentSize];
//    memset(dataContent,0,sizeof(uint64_t*)*contentSize);
    int contentIdx=0;
    for(int i=0;i<prime;i++){
        if(count[i]>1){
            dataContent[contentIdx]=new uint64_t[count[i]];
            contentIdx++;
        }
    }


    int* curCount=new int[contentSize];
    memset(curCount,0,sizeof(int)*contentSize);

    for(uint64_t i=0;i<size;i++){
        int idx=data[i]%prime;
        if(count[idx]>1){
            if(multiIdx[idx]>=contentSize){
                cerr<<"Error:code error"<<endl;
                exit(1);
            }
            dataContent[multiIdx[idx]][curCount[multiIdx[idx]]]=data[i];
            curCount[multiIdx[idx]]++;
            if(curCount[multiIdx[idx]]==count[idx]){
                unordered_set<uint64_t> exists;
//                string curStr;
                for(int j=0;j<count[idx];j++){
                    uint64_t curDa=dataContent[multiIdx[idx]][j];
//                    curStr=md52string(dataContent[multiIdx[idx]][j]);
                    if(exists.find(curDa)==exists.end()){
                        exists.insert(curDa);
                    }else{
//                        out<<dataContent[multiIdx[idx]][j]<<endl;
                        dataContent[multiIdx[idx]][j]=-1;
                    }
                }
            }
        }
    }
    memset(curCount,0,sizeof(int)*contentSize);
    int iter=0;
    for(uint64_t i=0;i<size;i++) {
        int idx = data[i] % prime;
        if (count[idx] > 1) {
            if (multiIdx[idx] >= contentSize) {
                cerr << "Error:code error" << endl;
                exit(1);
            }
            if (dataContent[multiIdx[idx]][curCount[multiIdx[idx]]] == -1) {
                dupFlag[i]=true;
//                dupIdx2.push_back(i);
                iter++;
            }
            curCount[multiIdx[idx]]++;
        }
    }
//    for(vector<uint64_t>::iterator ix=dupIdx.begin();ix!=dupIdx.end();ix++){
//        cout<<*ix<<endl;
//    }
    for(uint64_t i=0;i<size;i++){
        int idx=data[i]%prime;
        if(count[idx]>1){
            if(dataContent[multiIdx[idx]]!=NULL){
//                for(int j=0;j<count[idx];j++){
//                    delete[] dataContent[multiIdx[idx]][j];
//                    dataContent[multiIdx[idx]][j]=NULL;
//                }
                delete[] dataContent[multiIdx[idx]];
                dataContent[multiIdx[idx]]=NULL;
            }
        }
    }
    delete[] dataContent;
    delete[] multiIdx;
    delete[] curCount;
    delete[] count;
}

int rmdup::getPrime() {
    uint32_t realSize=size>pow(2,32)-1?pow(2,32)-1:size;
    int prime=0;

    float low=49/50;
    float high=1;
    while(1) {
        for (int i = realSize * low; i < realSize * high; i++) {
            bool isPrime = true;
            for (int j = 2; j <= sqrt(i); j++) {
                if (i % j == 0) {
                    isPrime = false;
                    break;
                }
            }
            if (isPrime) {
                prime = i;
            }
        }
        if(prime!=0){
            break;
        }else{
            high=low;
            low-=1/50;
            if(low<0.6){
                cerr<<"code error, report the bug to gongchun@genomics.cn"<<endl;
                exit(1);
            }
        }
    }
    return prime;
}

rmdup::~rmdup() {
//    for(uint64_t i=0;i<size;i++){
//        delete[] data[i];
//        data[i]=NULL;
//    }
    delete[] data;
    data=NULL;
}

uint64_t rmdup::getlong(unsigned char *str) {
    uint64_t value=0;
    for(int i=0;i<8;i++){
        value+=(uint64_t)str[i]<<(8*(7-i));
    }
    return value;
}

string rmdup::md52string(unsigned char *str) {
    string value;
    unsigned char *tmpStr=str;
    value=to_string(getlong(str));
    tmpStr+=8;
    value+=to_string(getlong(tmpStr));
    return value;
}
