//
// Created by berry on 2020-10-19.
//

#include "ReverseBloomFilter.h"

ReverseBloomFilter::ReverseBloomFilter(long readsNum, float multiple) {
    arrSize=readsNum*multiple;
    arr=new long[arrSize];
}

ReverseBloomFilter::ReverseBloomFilter(long readsNum, float multiple, long memSizeUsedInRmdup) {
    arrSize=readsNum*multiple;
    while (arrSize  > maxRBfSize) {
        multiple -= 0.5;
        if (multiple < 1) {
            cerr << "Error:reads number maybe is too large to do remove duplication" << endl;
            exit(1);
        }
    }
    long actualMem=arrSize*8;
    long memG=actualMem/(1024*1024);
    if(actualMem>memSizeUsedInRmdup){
        arrSize=memSizeUsedInRmdup/8;
        cerr<<"Error:given memSize is small, maybe it should be at least "<<memG<<"G"<<endl;
        exit(1);
    }
    while(actualMem<minRBfSize){
        multiple+=0.5;
        arrSize=readsNum*multiple;
        actualMem=arrSize*8;
    }
    arr=new long[arrSize];
    memset(arr,-1,arrSize);
}

bool ReverseBloomFilter::query(string &a) {
    curHash=hash<string>()(a);
    while(curHash<arrSize){
        curHash*=pow(2,10);
    }
    idx=curHash%arrSize;
    if(arr[idx]==-1){
        return false;
    }else{
        if(arr[idx]==curHash){
            return true;
        }else{
            return false;
        }
    }
}

int ReverseBloomFilter::add() {
    if(idx>=arrSize){
        return -1;
    }else{
        arr[idx]=curHash;
        return 0;
    }
}
