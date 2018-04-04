/*
 * DNAProcessor.h
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#ifndef PREPROCESSTOOL_3_LOCAL_DNAPROCESSOR_H_
#define PREPROCESSTOOL_3_LOCAL_DNAPROCESSOR_H_

#include "Common.h"
#include <openssl/md5.h>
#include <boost/thread/shared_mutex.hpp>

namespace PreProcessTool {


//线程任务函数的参数
    typedef struct TaskParam
    {
        Read *reads1;
        Read *reads2;
        int left, right;
        FqInfo info1, info2;

        TaskParam()
        {
            reads1 = NULL;
            reads2 = NULL;
        }
    } TaskParam;

//统计read截断了两端后的信息
    typedef struct StatisInfo
    {
        int q20, q30;
        int a, c, g, t, n;
        int lowQual;
        int readLen;

        StatisInfo()
        {
            bzero(this, sizeof(StatisInfo));
        }
    }StatisInfo;

//保存统计结果
    typedef struct StatisResult
    {
        bool hasAdpt;
        bool isLowQual;
        bool nExceed;
        bool isPolyA;
        bool isPolyX;
        //由于pe时，需要read1和read2的平均质量值
        int  sumQuality;

        StatisResult()
        {
            hasAdpt = false;
            isLowQual = false;
            nExceed = false;
            isPolyA = false;
            isPolyX = false;
            sumQuality = 0;
        }
    }StatisResult;


//memory limit's unit is MB
    const int MEM_UNIT = 1024 * 1024;

    class FilterProcessor
    {
    public:
        FilterProcessor();

        static void printUsage();
        static void printVersion();

        int filter(int argc, char **argv);

    private:
        /**
         * initial the member variables
         */
        int processParams(int argc, char **argv);

        /**
         * get the adapter list files' read name into readsName_
         */
        int getReadsNameFromFile(string filename, set<string> &readsName);

        /**
         * thread's task
         */
        void task(TaskParam *param);

        /**
         * output clean data
         */
        void outputRawDataTask(gzFile &file,Read *reads, unsigned int size);

        /**
         * ouput clean data
         */
        void outputCleanDataTaskStreamingSE(Read *reads);
        void outputCleanDataTaskStreamingPE(Read *reads,Read *reads2);
        void outputCleanDataTask(gzFile &file, Read *reads);
        void outputCleanDataTask2(gzFile &file, Read *reads, const char* name);

        void outputCleanData(gzFile &file, const StrRead &read);

        void outputTempData(ofstream& file, Read* reads1, Read* reads2, int outTypeMark);
        void outputTempData(ofstream& file, Read* reads);

        void outputDupData(ofstream& file, Read* reads1, Read* reads2, unsigned int size, int outTypeMark);
        void outputDupData(ofstream& file, Read* reads, unsigned int size);

        //merge all sorted files and remove duplication
        void mergeSortedFiles(int num, FqInfo *info1, FqInfo *info2, gzFile &outFile1, gzFile &outFile2);
        void mergeSortedFiles(int num, FqInfo *info, gzFile &outFile);

        int adaptorIndex(Read *read1,string adapter,int adpLen,set<string> &readsName,StatisResult &sr);

        bool statisticsPE(Read *read1, Read *read2, int index, FqInfo *info, FqInfo *info2);

        bool statisticsSE(Read *read, int index, FqInfo *info);


        /**
         * called by statistics
         */
        StatisInfo auxStatistics(Read *read, int headTrim, int tailTrim, string adapter, int adptLen, set<string> &readsName, FqInfo &info, StatisResult &sr);

        //statistics the clean read
        StatisInfo auxStatistics(StrRead *read);

        /**
         * 用于adapter list情况
         */
        bool hasAdapter(set<string> &readsName, const char *seqName);

        /**
         * 用于adapter为序列的情况
         */
        //bool hasAdapter(const char *sequence, int readLen, const char *adapter, int adptLen);
        int hasAdapter(const char *sequence, int readLen, const char *adapter, int adptLen);

        /**
         * 判断insert size是否过小
         */
        bool isSmallSize(const char *sequence1, int readLen1, const char *sequence2, int readLen2);

        /**
         * 统计碱基的分布
         */
        void calculateBaseDistribute(Read *read, FqInfo &info, int readLen);
        void calculateBaseDistribute(StrRead *read, FqInfo &info, int readLen);

        //void getTiles(string tiles, set<int> &tileSet);

    private:
        unsigned int PROCESS_THREAD_NUM;
        bool IS_STREAMING;
        char* streamingInput;
        gzFile fqStreaming;

        int maskLowQual; // mask low quality base with N

        int trimBadHeadMaxLength;
        int trimBadHeadQuility;
        int trimBadTailMaxLength;
        int trimBadTailQuility;

        int TtoU;
        int UtoT;

        set<string> readsName1_;   //store the read's name which in adapter list file
        set<string> readsName2_;

        bool filterTile_;
        set<string> tiles_;

        bool tileIsFov_;

        string adapter1_;
        string adapter2_;

        string fqFile1_;
        string fqFile2_;

        int misMatch_;  //接头允许错配数
        float matchRatio_;  //接头比对比列

        int lowQual_;    //low quality
        float qualRate_;
        float nRate_;

        float polyA_;
        int polyX_;

        float minMean_;

        bool filterIndex_;    //indicate whether filter the read name's index or not

        bool rmdup_;     //indicate whether remove duplication or not
        bool dupRateOnly_; // calculate duplicaitons rate only

        unsigned long cutReadNum_;  //截取指定数据量clean data的碱基量,单位M, 0表示不截取

        int headTrim_;   //trim the 5' end of read some bp
        int tailTrim_;   //trim the 3' end of read some bp
        int headTrim2_;  //
        int tailTrim2_;

        long memLimit_;
        QualitySystem qualSys_;  //the fq file's quality system


        bool isFilterSmallInsertSize_;  //是否过滤掉insert size过短的
        //判断为small insert size的条件
        int overlap_;
        float mis_;

        int readLen_;    //read length in fq1 file
        int readLen2_;   //read length in fq2 file
        string outDir_;  //all output file's output path
        bool onlyStat_;

        bool isPE_;
        int minReadLength;
        bool cutAdaptor;
        int cutAdaptorOri;
        unsigned long cutBasesNumber;
        bool isAdptList_;
        //用于指示clean data是否已经足够了
        bool isFull_;
        //如果adapter为序列，则只是adapter序列的长度
        int adapterLen1_;
        int adapterLen2_;

        boost::mutex dupMutex_;
        map<string, int> duplications_;    //存放fq1和fq2的read碱基序列连接起来的字符串,用于判断是否有重复
        unsigned char md5Seq_[MD5_DIGEST_LENGTH];   //存放碱基序列MD5值

        boost::mutex sizeMutex_;
        unsigned int *cleanDataIndexs_;  //存放处理线程判断为clean read的那些read在buffer中的下标
        boost::uint32_t size_;   //cleanDataIndexs_里面的元素个数

        boost::uint32_t doneNum_;  //处理完的线程数

        int encodSize_;  //read的碱基序列编码所需长度

        string lanID_;
        QualitySystem cleanQualSys_;

        bool filterAdapter_;   //是否需要过滤adapter

        //待输出的文件名
        string cleanFq1_;
        string cleanFq2_;
        string rawFq1_;
        string rawFq2_;

        int seqType_;
        int outType_;
        int polyAType_;
    };

}  // namespace PreProcessTool


#endif /* PREPROCESSTOOL_3_LOCAL_DNAPROCESSOR_H_ */
