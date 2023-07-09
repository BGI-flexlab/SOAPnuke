//
// Created by berry on 2020-10-09.
//

#include "BloomFilter.h"
#include <cstdint>

BloomFilter::BloomFilter(long long sampleSize, int multiple):multiple(20) {
    if(sampleSize==0){
        cerr<<"Error:no reads found in input file"<<endl;
        exit(1);
    }
    hashNum=round(log(2)*multiple);
    realUseBitSize=sampleSize*multiple;
    if(realUseBitSize>maxBfSize || realUseBitSize>maxBfSize){
        cerr<<"Error:reads number maybe is too large to do remove duplication"<<endl;
        exit(1);
    }
    realUseByteSize=realUseBitSize/8+10;
    arr=new char[realUseByteSize];
    memset(arr,0,realUseByteSize);
    indexs=new long[hashNum];
//    cout<<"memSize used in rmdup:"<<realUseByteSize/(1024*1024)<<"M"<<endl;
}

BloomFilter::~BloomFilter() {
    delete[] indexs;
    delete[] arr;
}

float BloomFilter::expectedFP() {
    return pow(0.6185,multiple);
}

bool BloomFilter::query(string& seq) {
    bool exist=true;
    for(int i=0;i<hashNum;i++) {
        unsigned long value = createHash2(seq,i);
            long index=value%realUseBitSize;
        indexs[i]=index;
        if(exist && !getPosStatus(index)){
            exist=false;
        }
    }
    return exist;
}

void BloomFilter::add() {
    for(int i=0;i<hashNum;i++){
        setPosStatus(indexs[i]);
    }
}

unsigned long BloomFilter::createHash(string &str,int j) {
    switch(j){
        case 0:{
            unsigned long v=5381;
            for(int i = 0; i < str.length(); i++){
                v = ((v << 5) + v) + str[i];
            }
            return v;
        }
        case 1:{
            unsigned long hash = str.length();
            for(int i = 0; i < str.length(); i++)
            {
                hash = ((hash << 5) ^ (hash >> 27)) ^ str[i];
            }
            return hash;
        }
        case 2:{
            unsigned long hash = 0;
            for(int i = 0; i < str.length(); i++)
            {
                hash = str[i] + (hash << 6) + (hash << 16) - hash;
            }
            return hash;
        }
        case 3:{    //  BKDRHash
            unsigned long seed = 131; // 31 131 1313 13131 131313 etc..
            unsigned long hash = 0;
            for(int i = 0; i < str.length(); i++)
            {
                hash = (hash * seed) + str[i];
            }
            return hash;
        }
        case 4:{
            unsigned long hash = 0;
            unsigned long x    = 0;
            for(int i = 0; i < str.length(); i++)
            {
                hash = (hash << 4) + str[i];
                if((x = hash & 0xF0000000L) != 0)
                {
                    hash ^= (x >> 24);
                }
                hash &= ~x;
            }
            return hash;
        }
        case 5:{
            unsigned long BitsInUnsignedInt = (unsigned long)(4 * 8);
            unsigned long ThreeQuarters     = (unsigned long)((BitsInUnsignedInt  * 3) / 4);
            unsigned long OneEighth         = (unsigned long)(BitsInUnsignedInt / 8);
            unsigned long HighBits          = (unsigned long)(0xFFFFFFFF) << (BitsInUnsignedInt - OneEighth);
            unsigned long hash              = 0;
            unsigned long test              = 0;
            for(int i = 0; i < str.length(); i++)
            {
                hash = (hash << OneEighth) + str[i];
                if((test = hash & HighBits)  != 0)
                {
                    hash = (( hash ^ (test >> ThreeQuarters)) & (~HighBits));
                }
            }
            return hash;
        }
        case 6:{
            unsigned long hash = 1315423911;
            for(int i = 0; i < str.length(); i++)
            {
                hash ^= ((hash << 5) + str[i] + (hash >> 2));
            }
            return hash;
        }
        case 7:{    //RSHash
            int b     = 378551;
            int a     = 63689;
            unsigned long hash = 0;
            for(int i = 0; i < str.length(); i++)
            {
                hash = hash * a + str[i];
                a    = a * b;
            }
            return hash;
        }
        default:{
            cerr<<"Error:code error"<<endl;
            exit(1);
        }
    }

}

bool BloomFilter::getPosStatus(long idx) {
    uint64_t bigIdx=idx/8;
    int smallIdx=idx%8;
    char c=arr[bigIdx];
    return (c<<smallIdx) & 0b10000000;
}

void BloomFilter::setPosStatus(long idx) {
    uint64_t bigIdx=idx/8;
    int smallIdx=idx%8;
    arr[bigIdx]=arr[bigIdx] | (0b10000000>>smallIdx);
}

unsigned long BloomFilter::createHash2(string &str, int i) {
    string newStr=to_string(i)+str;
    return hash<string>()(newStr);
}
