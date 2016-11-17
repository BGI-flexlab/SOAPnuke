/*
 * Common.h
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#ifndef PREPROCESSTOOL_3_LOCAL_METACOMMON_H_
#define PREPROCESSTOOL_3_LOCAL_METACOMMON_H_

#include "CommonInclude.h"
#include <string.h>
#include "Logger.h"

namespace MetaPreProcessTool {

//本工具支持的最大读长
const static int MAX_LENGTH = 512;
//本工具支持的最大质量值
const static int MAX_QUALITY = 100;

//本工具支持的质量体系
enum QualitySystem
{
	ILLUMINA_ = 64,
	SANGER_ = 33,
};

//read 信息
typedef struct {
	char *readName;
	char *baseSequence;
	char *optionalName;
	char *baseQuality;
}Read;


/**
 * 存放read1和read2的readname连接后的字符串的编码
 */
typedef struct ReadSeq
{
	ReadSeq(int num)
	{
		size = num;
		readName = new uint64_t[size + 1];
		for (int i=0; i<=size; i++)
		{
			readName[i] = 0;
		}
	}

	ReadSeq(const ReadSeq &rs)
	{
		size = rs.size;
		readName = new uint64_t[size + 1];
		for (int i=0; i<=size; i++)
		{
			readName[i] = rs.readName[i];
		}
	}

	~ReadSeq()
	{
		if (readName != NULL)
		{
			delete []readName;
			readName = NULL;
		}
	}

	uint64_t *readName;
	int size;

}ReadSeq;

struct Compare
{
    bool operator()(const ReadSeq &a, const ReadSeq &b)
    {
    	int len = a.size;
        for (int i=0;i<len; i++ )
        {
            if (a.readName[i] < b.readName[i])
            {
                return true;
            }
            else if (a.readName[i] > b.readName[i])
            {
                return false;
            }
        }
        return false;
    }
};

/////////////////////////////////////
/////default output files name///////
/////////////////////////////////////
//statistics files name
const static string LOG_FILE = "out.log";
const static string SEQUENCING_QUALITY = "Basic_Statistics_of_Sequencing_Quality";
const static string FILTERED_READS = "Statistics_of_Filtered_Reads";
const static string BASE_DISTRIBUTIONS = "Base_distributions_by_read_position";
const static string DISTRIBUTION_OF_Q20_Q30 = "Distribution_of_Q20_Q30_bases_by_read_position";
const static string BASE_QUALITY_VALUE_DISTRIBUTION = "Base_quality_value_distribution_by_read_position";
//clean fq file name's prefix
const static string CLEAN_FQ_PREFIX = "Clean_";
//raw fq file prefix
const static string RAW_FQ_PREFIX = "Raw_";

//代表FqInfo接头体中成员base二维数组的每行的下标
enum Base
{
	A_ = 0, C_, G_, T_, N_,
};

//FqInfo结构用于存储过滤的所有统计信息。
struct FqInfo
{
	unsigned int rawReadLength;     //raw data 读长
	unsigned int cleanReadLength;   //clean data 读长
	unsigned int cleanReadLengthTmp;
	unsigned long rawTotalReadNum;  //raw data read个数
	unsigned long cleanTotalReadNum;//clean data read 个数
	unsigned long rawTotalBaseNum;  //raw data total base number
	unsigned long cleanTotalBaseNum; //clean data total base number
	unsigned long rawBaseA;     //raw data base A number
	unsigned long cleanBaseA;   //clean data base A number
	unsigned long rawBaseC;     //raw data base C number
	unsigned long cleanBaseC; //clean data base C number
	unsigned long rawBaseG; //raw data base G number
	unsigned long cleanBaseG; //clean data base G number
	unsigned long rawBaseT; //raw data base T number
	unsigned long cleanBaseT; //clean data base T number
	unsigned long rawBaseN; //raw data base N number
	unsigned long cleanBaseN; //clean data base N number
	unsigned long rawQ20; //rawfq文件中碱基质量>=20的碱基总数
	unsigned long cleanQ20; //cleanfq文件中碱基质量>=20的碱基总数
	unsigned long rawQ30; //rawfq文件中碱基质量>=30的碱基总数
	unsigned long cleanQ30; //cleanfq文件中碱基质量>=30的碱基总数

	unsigned long duplicationNum; //重复个数
	unsigned long adapterNum;  //the number of read which contain adapter in raw data
	unsigned long nExceedNum;  //the number of read which n rate was exceed in raw data
	unsigned long lowQualNum;  //low qualtiy read number in raw data
	unsigned long lowMeanNum;  //low mean quality read number in raw data
	unsigned long smallInsertNum;  //samll inert number in raw data
	unsigned long polyNNum;    //polyA number in raw data

	unsigned long totalDuplicationNum;
	unsigned long totalAdapterNum;
	unsigned long totalNExceedNum;
	unsigned long totalLowQualNum;
	unsigned long totalLowMeanNum;
	unsigned long totalSmallInsertNum;
	unsigned long totalPolyNNum;
	unsigned long singleReadNum;
	unsigned long singleBaseNum;

	//base distributions by read position(Raw)
	unsigned long base[MAX_LENGTH][5]; //ACGT
	unsigned long clean_base[MAX_LENGTH][5];

	//DistributionofQ20+/Q30+basesbyreadposition(Raw)
	unsigned long q20q30[MAX_LENGTH][2];  //Q20 Q30
	unsigned long clean_q20q30[MAX_LENGTH][2];

	//Basequalityvaluedistributionbyreadposition(Raw)
	unsigned long qual[MAX_LENGTH][MAX_QUALITY + 1];
	unsigned long clean_qual[MAX_LENGTH][MAX_QUALITY + 1];

	int maxQualityValue;  //记录fq文件中碱基的最大质量值.

	FqInfo()
	{
		bzero(this, sizeof(FqInfo));
		maxQualityValue = 41;
	}

	void add(const FqInfo &src)
	{
		rawTotalReadNum += src.rawTotalReadNum;
		cleanTotalReadNum += src.cleanTotalReadNum;
		rawTotalBaseNum += src.rawTotalBaseNum;
		cleanTotalBaseNum += src.cleanTotalBaseNum;
		rawBaseA += src.rawBaseA;
		cleanBaseA += src.cleanBaseA;
		rawBaseC += src.rawBaseC;
		cleanBaseC += src.cleanBaseC;
		rawBaseG += src.rawBaseG;
		cleanBaseG += src.cleanBaseG;
		rawBaseT += src.rawBaseT;
		cleanBaseT += src.cleanBaseT;
		rawBaseN += src.rawBaseN;
		cleanBaseN += src.cleanBaseN;

		rawQ20 += src.rawQ20;
		cleanQ20 += src.cleanQ20;
		rawQ30 += src.rawQ30;
		cleanQ30 += src.cleanQ30;

		adapterNum += src.adapterNum;
		nExceedNum += src.nExceedNum;
		lowQualNum += src.lowQualNum;
		lowMeanNum += src.lowMeanNum;
		smallInsertNum += src.smallInsertNum;
		polyNNum += src.polyNNum;
		duplicationNum += src.duplicationNum;

		totalDuplicationNum += src.totalDuplicationNum;
		totalAdapterNum += src.totalAdapterNum;
		totalNExceedNum += src.totalNExceedNum;
		totalLowQualNum += src.totalLowQualNum;
		totalLowMeanNum += src.totalLowMeanNum;
		totalSmallInsertNum += src.totalSmallInsertNum;
		totalPolyNNum += src.totalPolyNNum;
		singleReadNum += src.singleReadNum;
		singleBaseNum += src.singleBaseNum;

		for (unsigned int i = 0; i < rawReadLength; ++i)
		{
			base[i][0] += src.base[i][0];
			base[i][1] += src.base[i][1];
			base[i][2] += src.base[i][2];
			base[i][3] += src.base[i][3];
			base[i][4] += src.base[i][4];

			q20q30[i][0] += src.q20q30[i][0];
			q20q30[i][1] += src.q20q30[i][1];

			for (int j = 0; j <= MAX_QUALITY; ++j)
			{
				qual[i][j] += src.qual[i][j];
			}
		}

		for (unsigned int i = 0; i < cleanReadLength; ++i)
		{
			clean_base[i][0] += src.clean_base[i][0];
			clean_base[i][1] += src.clean_base[i][1];
			clean_base[i][2] += src.clean_base[i][2];
			clean_base[i][3] += src.clean_base[i][3];
			clean_base[i][4] += src.clean_base[i][4];

			clean_q20q30[i][0] += src.clean_q20q30[i][0];
			clean_q20q30[i][1] += src.clean_q20q30[i][1];

			for (int j = 0; j <= MAX_QUALITY; ++j)
			{
				clean_qual[i][j] += src.clean_qual[i][j];
			}
		}
	}

	void clear()
	{
		bzero(this, sizeof(FqInfo));
	}
};

/**
 * if the seq is only contain these 'ACGTacgt' characters
 * 		return true
 * else
 * 		return false
 */
bool isSequence(const string &seq);

/**
 * return the adapter's type
 * 1: adapter list
 * 2: adapter sequence
 * 3: error
 */
int adapterType(bool isPE, const string &adapter1, const string &adapter2);

/**
 * print all the fqInfo's member variables
 */
void printFqInfo(const string outDir, const string prefix, const FqInfo *fqInfo1, const FqInfo *fqInfo2);

/**
 * get the output file name for fq file
 */
string getOutputFileName(string filename, string prefix, string path);

void encode(string &seq, ReadSeq &base);

/**
 * get the sequence's reverse complementary strand
 */
string reverseComplementary(string &seq);

/**
 * upper the string
 */
void upper(char* str);


}  // namespace PreProcessTool

#endif /* PREPROCESSTOOL_3_LOCAL_COMMON_H_ */
