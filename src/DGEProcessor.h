/*
 * DGEProcessor.h
 *
 *  Created on: 2012-6-3
 *      Author: Shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */

#ifndef DGEPROCESSTOOL_DGEPROCESSOR_H_
#define DGEPROCESSTOOL_DGEPROCESSOR_H_

#include "DGECommon.h"
#include "FqBuffer.h"
#include "Logger.h"
#include "Common.h"

namespace DGEProcessTool {

//线程任务函数的参数
typedef struct TaskParam
{
	PreProcessTool::Read *reads1;
	int left, right;
	FqInfo info1;
	int *result;
	TaskParam()
	{
		reads1 = NULL;
	}
} TaskParam;

//
typedef struct StatisInfo
{
	int a, c, g, t, n, ns; //ns 是指第截取tag的碱基N的个数//
	int q20, q30;
	StatisInfo()
	{
		a = 0;
		c = 0;
		g = 0;
		t = 0;
		n = 0;
		ns = 0;
		q20 = 0;
		q30 = 0;
	}
}StatisInfo;

//memory limit's unit is MB
const int MEM_UNIT = 1024 * 1024;

class DGEProcessor
{
public:
	DGEProcessor():
		outDir_("."),
		outPfx_("clean"),
		site_("CATG"),
		chr_("ACGTN"),
		tagLength_(17),
		tagStart_(17),
		tagEnd_(18),
		misMatch_(0),
		headTrim_(0),
		tailTrim_(0),
		cutOff_(1),
		readLen_(49),
		score_(64),
		threadNum_(4),
		standard_(1000000),
		nSeed_(1),
		number_(0),
		memLimit_(256 * MEM_UNIT),
		qualSystem_(ILLUMINA_),
		filterLowQual_(true),
		filterIndex_(false),
		outfq_(false),
		filterTile_(false)
	{};

	static void printUsage();
	static void printVersion();
	int processRNA(int argc, char **argv);

private:

	/**
	 * 初始化类的各个成员变量
	 */
	int processParams(int argc, char **argv);

	int processDGE();
	void taskDGE(TaskParam *param);//线程函数
	int statisticsDGE(PreProcessTool::Read &read, FqInfo &info);//统计函数,并判断是否过滤read


	/**
	 * called by statistics
	 */
	StatisInfo auxStatistics(PreProcessTool::Read &read, FqInfo &info);//统计fq的函数

	/*
	 * 处理3'接头
	 */
	void trimAdapter(char *sequence, int startPos);
	bool findAdapter(const char *sequence, const char *adapter);
	int smithWatermanAlign(string query, string target);

	void output(FqInfo *info1, FqInfo *info2);

//	void getTiles(string tiles, set<int> &tileSet);

private:
	map<string,int> nSequence_;
	map<string,int> adptSequence_;
	map<string,int> tagSequence_;
	map< unsigned long,vector<unsigned long> > sortClean_;
	vector< map<string,int>::iterator > cleanFa_;

	string fqFile1_;
	string fqFile2_;
	string adapter_;
	
	string outDir_;
	string outPfx_;
	string site_;
	string chr_;
	int tagLength_;
	int tagStart_;
	int tagEnd_;
	int misMatch_;
	int headTrim_;
	int tailTrim_;
	int cutOff_;
	int readLen_;
	int score_;
	int threadNum_;
	int standard_;
	int nSeed_;
	unsigned long number_;
	long memLimit_;
	QualitySystem qualSystem_;
	bool filterLowQual_;
	bool filterIndex_;
	bool outfq_;
	long tagPercent[7][2];
	float match_[5][5];

	bool filterTile_;
	set<string> tiles_;
};

} // namespace DGEProcessTool
#endif /* DGEPROCESSTOOL_DGEPROCESSOR_H_ */
