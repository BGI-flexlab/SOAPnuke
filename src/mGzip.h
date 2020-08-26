//
// Created by berry on 2020-08-11.
//

#ifndef SOAPNUKE_MGZIP_H
#define SOAPNUKE_MGZIP_H
#include<string>
#include <vector>
#include <zlib.h>
#include<iostream>
using namespace::std;
//only support SE data.
class threadDataInfo{
public:
    //implicit construct function
    threadDataInfo(){
        index=-1;
    }
    int index;
    //If input is a single file, the size would all be 1
    vector<string> fileList;
    //store the compress and raw size of each block, one thread may process several blocks
    vector<string> rawFileName;
    vector<int> method;
    vector<int> flag;
//    vector<long> mtime;
//    vector<vector<long> > compressedBlockSize;
    vector<vector<u_int> > crc32;
    vector<vector<long> > rawBlockSize;
    vector<vector<long> > seekOffSetStart;
    vector<vector<long> > spanLength;
};

class mGzip {
public:
    static bool check_mGzip(string filePath);
    static vector<threadDataInfo> allocate(int threadsNum, vector<string> filesPath);

    static int getOneBlock(FILE *pFile, int& method, int& flag, long &curOffSet, vector<long> &blockStarts, vector<long> &length, vector<long> &size,
                            vector<u_int> &crc32, string &fileName);
};


#endif //SOAPNUKE_MGZIP_H
