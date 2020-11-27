//
// Created by berry on 2020-10-09.
//

#ifndef SOAPNUKE_BLOOMFILTER_H
#define SOAPNUKE_BLOOMFILTER_H

#include<string>
#include<bitset>
#include<iostream>
#include<math.h>
#include<string.h>
using namespace::std;
#define maxBfSize 1024L*1024*1024*200
class BloomFilter {
private:
    int multiple;
    int hashNum;
    long realUseBitSize;
    long* indexs;
    char* arr;
public:
    BloomFilter(long long sampleSize,int multiple);
    bool query(string& seq);
    void add();
    bool getPosStatus(long idx);
    void setPosStatus(long idx);
    static unsigned long createHash(string& seq,int i);
    float expectedFP();
    long realUseByteSize;
    ~BloomFilter();

    unsigned long createHash2(string &str, int i);
};


#endif //SOAPNUKE_BLOOMFILTER_H
