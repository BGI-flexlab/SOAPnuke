//
// Created by berry on 2020-10-20.
//

#ifndef SOAPNUKE_RMDUP_H
#define SOAPNUKE_RMDUP_H
#include<iostream>
#include<fstream>
#include<unordered_set>
#include<vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include "gc.h"
using namespace::std;
#define AUTHOREMAIL "gongchun@genomics.cn"
class rmdup {
private:
    uint64_t* data;
    uint64_t size;
    uint32_t getPrime();
public:
    rmdup(uint64_t* data,uint64_t size);
    ~rmdup();
    void markDup(bool* dupFlag);
    uint64_t getlong(unsigned char* str);
    string md52string(unsigned char* str);
};


#endif //SOAPNUKE_RMDUP_H
