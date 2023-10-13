//
// Created by berry on 2020-10-19.
//

#ifndef SOAPNUKE_REVERSEBLOOMFILTER_H
#define SOAPNUKE_REVERSEBLOOMFILTER_H
#include <string>
#include <iostream>
#include <string.h>
#include <math.h>
#include <cstdint>
using namespace ::std;
#define maxRBfSize 1024L * 1024 * 1024 * 4
#define minRBfSize 1024L * 1024 * 1000
class ReverseBloomFilter
{
	long idx;
	uint64_t curHash;
	long *arr;

public:
	ReverseBloomFilter(long readsNum, float multiple);
	ReverseBloomFilter(long readsNum, float multiple, long memSizeUsedInRmdup);
	bool query(string &a);
	int add();
	long arrSize;
};

#endif // SOAPNUKE_REVERSEBLOOMFILTER_H
