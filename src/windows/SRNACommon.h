/*
 * Common.h
 *
 * This header file provide the common definition and interface to
 * DNAProcessor and RNAProcessor Class.
 *
 *  Created on: 2012-6-3
 *      Author: shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */

#ifndef SRNAPROCESSTOOL_COMMON_H_
#define SRNAPROCESSTOOL_COMMON_H_
#include "CommonInclude.h"
#include "Logger.h"

namespace SRNAProcessTool {

//本工具所支持的最大读长
const static int MAX_LENGTH = 256;
//本工具所支持的最大质量值
const static int MAX_QUALITY = 42;

//本工具所支持的质量体系
enum QualitySystem
{
	ILLUMINA_ = 64,
	SANGER_ = 33,
};

/////////////////////////////////////
/////default output files name///////
/////////////////////////////////////
//log file name
const static string LOGOUT_NAME = "out.log";
//statistics files name
const static string SEQUENCING_QUALITY = "Basic_Statistics_of_Sequencing_Quality.txt";
const static string FILTERED_READS = "Statistics_of_Filtered_Reads.txt";
const static string BASE_DISTRIBUTIONS = "Base_distributions_by_read_position";
const static string DISTRIBUTION_OF_Q20_Q30 = "Distribution_of_Q20_Q30_bases_by_read_position";
const static string BASE_QUALITY_VALUE_DISTRIBUTION = "Base_quality_value_distribution_by_read_position";
const static string LENGTH_DISTRIBUTION = "Length_distribution.txt";
//clean fq file name's prefix
const static string CLEAN_FQ_PREFIX = "Clean_";



//代表FqInfo接头体中成员base二维数组的每行的下标
enum Base
{
	A_ = 0, C_, G_, T_, N_,
};

//FqInfo结构用于存储过滤的所有统计信息。
struct FqInfo
{
	//BasicStatisticsofSequencingQuality
	unsigned int readslength; //rawfq文件的read长度最大值（一般测序长度都一样）
	unsigned long rawTotalReads; //rawfq文件的read总数
	unsigned long cleanTotalReads; //cleanfq文件的read总数
	unsigned long rawTotalBases; //rawfq文件的碱基总数
	unsigned long cleanTotalBases; //cleanfq文件的碱基总数
	unsigned long rawBaseA; //rawfq文件所有read中碱基A的总数
	unsigned long cleanBaseA; //cleanfq文件所有read中碱基A的总数
	unsigned long rawBaseC; //rawfq文件所有read中碱基C的总数
	unsigned long cleanBaseC; //cleanfq文件所有read中碱基C的总数
	unsigned long rawBaseG; //rawfq文件所有read中碱基G的总数
	unsigned long cleanBaseG; //cleanfq文件所有read中碱基C的总数
	unsigned long rawBaseT; //rawfq文件所有read中碱基T的总数
	unsigned long cleanBaseT; //cleanfq文件所有read中碱基T的总数
	unsigned long rawBaseN; //rawfq文件所有read中碱基为N的总数
	unsigned long cleanBaseN; //cleanfq文件所有read中碱基为N的总数
	unsigned long rawQ20; //rawfq文件中碱基质量>=20的碱基总数
	unsigned long cleanQ20; //cleanfq文件中碱基质量>=20的碱基总数
	unsigned long rawQ30; //rawfq文件中碱基质量>=30的碱基总数
	unsigned long cleanQ30; //cleanfq文件中碱基质量>=30的碱基总数

	//StatisticsofFilteredReads
	unsigned long readsWithAdapter; //rawfq文件中带有adapter的read的总数
	unsigned long readsWithNrate; //rawfq文件中N值比例超过阀值的read的总数
	unsigned long readsWithLowQual;//rawfq文件中read中低质量值碱基比例大于阀值的//read的总数
	unsigned long readWithShortValidLength; //去掉两端接头后，长度小于18的read的总数
	unsigned long readWithPolyA; //带PolyA的read的总数
	unsigned long adapter3Null;
	unsigned long insertNull;
	unsigned long adapter5Pollute;
	int lengthStart;
	int lengthEnd;

	//base distributions by read position(Raw)
	unsigned long base[MAX_LENGTH][5]; //ACGTN，使用枚举类型

	//DistributionofQ20+/Q30+basesbyreadposition(Raw)
	unsigned long q20q30[MAX_LENGTH][2]; //Q20和Q30

	//Basequalityvaluedistributionbyreadposition(Raw)
	unsigned long qual[MAX_LENGTH][MAX_QUALITY]; //2-41
	unsigned long lengthDis[MAX_LENGTH];


	FqInfo()
	{
		readslength = 0;
		rawTotalReads = 0;
		cleanTotalReads = 0;
		rawTotalBases = 0;
		cleanTotalBases = 0;
		rawBaseA = 0;
		cleanBaseA = 0;
		rawBaseC = 0;
		cleanBaseC = 0;
		rawBaseG = 0;
		cleanBaseG = 0;
		rawBaseT = 0;
		cleanBaseT = 0;
		rawBaseN = 0;
		cleanBaseN = 0;
		rawQ20 = 0;
		cleanQ20 = 0;
		rawQ30 = 0;
		cleanQ30 = 0;

		readsWithAdapter = 0;
		readsWithNrate = 0;
		readsWithLowQual = 0;
		readWithShortValidLength = 0;
		readWithPolyA = 0;
		adapter3Null = 0;
		insertNull = 0;
		adapter5Pollute = 0;

		lengthStart = 18;
		lengthEnd = 44;

		memset(base, 0, sizeof(long) * 5 * MAX_LENGTH);
		memset(q20q30, 0, sizeof(long) * 2 * MAX_LENGTH);
		memset(qual, 0,	sizeof(long) * MAX_QUALITY * MAX_LENGTH);
		memset (lengthDis, 0, sizeof(long) * MAX_LENGTH);
	}

	/**
	 * add the src to this
	 */
	void add(FqInfo &src, int readLen)
	{
		rawTotalReads += src.rawTotalReads;
		cleanTotalReads += src.cleanTotalReads;
		rawTotalBases += src.rawTotalBases;
		cleanTotalBases += src.cleanTotalBases;
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

		readsWithAdapter += src.readsWithAdapter;
		readsWithNrate += src.readsWithNrate;
		readsWithLowQual += src.readsWithLowQual;
		readWithShortValidLength += src.readWithShortValidLength;
		readWithPolyA += src.readWithPolyA;

		adapter3Null += adapter3Null;
		insertNull += insertNull;
		adapter5Pollute += adapter5Pollute;

		for (int i = 0; i < readLen; ++i)
		{
			base[i][0] += src.base[i][0];
			base[i][1] += src.base[i][1];
			base[i][2] += src.base[i][2];
			base[i][3] += src.base[i][3];
			base[i][4] += src.base[i][4];

		}

		for (int i = 0; i < readLen; ++i)
		{
			q20q30[i][0] += src.q20q30[i][0];
			q20q30[i][1] += src.q20q30[i][1];
		}

		for (int i = 0; i < readLen; ++i)
			for (int j = 0; j < MAX_QUALITY; ++j)
			{
				qual[i][j] += src.qual[i][j];
			}

		for (int i = 0; i < readLen; i++)
		{
			lengthDis[i] += src.lengthDis[i];
		}

	}
};

/**
 * if the seq is only contain these 'ACGTacgt' characters
 * 		return true
 * else
 * 		return false
 */
bool isSequence(const std::string &seq);

/**
 * return the adapter's type
 * 1: adapter list
 * 2: adapter sequence
 * 3: error
 */
int adapterType(const std::string &adapter1, const std::string &adapter2);

/**
 * assignment the fqInfo's member variables to 0
 */
void clearFqInfo(FqInfo &fqInfo);

/**
 * get the output file name for fq file
 */
string getOutputFileName(string filename, string path);

/**
 * print all the fqInfo's member variables
 */
void printFqInfo(const string outDir, FqInfo &fqInfo);

}  // namespace SRNAProcessTool
#endif /* SRNAPROCESSTOOL_COMMON_H_ */
