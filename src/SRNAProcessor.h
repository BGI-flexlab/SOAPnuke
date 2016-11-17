/*
 * RNAProcessor.h
 *
 *  Created on: 2012-6-3
 *      Author: Shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */

#ifndef SRNAPROCESSTOOL_RNAPROCESSOR_H_
#define SRNAPROCESSTOOL_RNAPROCESSOR_H_

#include "SRNACommon.h"
#include "FqBuffer.h"
#include "Logger.h"
#include "Common.h"
using namespace boost::threadpool;
namespace SRNAProcessTool {

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
		int a, c, g, t, n, ns; //ns 是指第1-30个碱基内的N个数//
		int q20, q30, lowQual1,lowQual2;//lowQual1是1-30碱基q<10的个数, lowQual2是1-30碱基q<13的个数.//
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
			lowQual1 = 0;
			lowQual2 = 0;
		}
	}StatisInfo;

	//memory limit's unit is MB
	const int MEM_UNIT = 1024 * 1024;

	class RNAProcessor
	{
		public:
			RNAProcessor():
				strict_(false),
				filterTmpNumber_(0),
				seedLength_(30),
				minInsertSize_(18),
				minSize_(100),
				adapter5_("GTTCAGAGTTCTACAGTCCGACGATC"),
				adapter3_("TCGTATGCCGTCTTCTGCTTG"),
				mrna_(false),
				fastq_(false),
				memLimit_(MEM_UNIT * 512),
				threadNum_(6),
				qualSystem_(ILLUMINA_),
				readLen_(49),
				outDir_("."),
				onlyStat_(false),
				outPfx_("clean"),
				nSeed_(0),
				nRead_(2),
				lowQualSeed1_(10),
				lowQualSeed1Num_(4),
				lowQualSeed2_(13),
				lowQualSeed2Num_(6),
				error_(0.01),
				removeBad_(false),
				trim_(true),
				filterPolyN_(0.7),
				filterPolyA_(0.7),
				filterLowQual_(true),
				continuousAlign_(6),
				alignRate_(0.8),
				miniAlignLength_(5),
				errorRate_(0.4),
				misMatch_(4),
				filterIndex_(false),
				outfq_(false),
				headTrim_(0),
				tailTrim_(0),
				cutOff_(0),
				score_(64),
				seqType_(0),
				pl_(6),
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

			/**
			 * 处理输出fq
			 */
			int processFQ();
			void taskFQ(TaskParam *param);//线程函数
			int statisticsFQ(PreProcessTool::Read &read, FqInfo &info);//统计函数,并判断是否过滤read


			/**
			 * 处理输出fa
			 */

			int processFA();
			void taskFA(TaskParam *param);//线程函数
			int statisticsFA(PreProcessTool::Read &read, FqInfo &info);//统计函数
			bool filterFA(const string &sequence, int &count, FqInfo &info);	//过滤函数
			//	void filterTags(FqInfo &info);
			/*	void sortCleanFa(vector<long> &arrayPair,long x,long y);*/
			/*	bool sortFuction(map<string,int>::iterator x, map<string,int>::iterator y);*/


			/**
			 * called by statistics
			 */
			StatisInfo auxStatistics(PreProcessTool::Read &read, FqInfo &info);//统计fq的函数

			/*
			 * 处理3'接头
			 */
			void trimAdapter(char *sequence, int startPos);
			int findAdapter(const char *sequence, const char *adapter);

			/*
			 * 处理5'接头
			 */
			bool hasAdapter(const char *sequence, const char *adapter);
			/*
			 * 输出统计文件
			 */
			void output(FqInfo *info1, FqInfo *info2);
			bool hasPolyN(char *sequence, int title, int num, ofstream &outFile, long polyNCount[][2], int wT);

			string mergeAndSortFilesByAscii(string &outpfx, vector<string> &files);
			string mergeAndSortFilesByNumber(string &outpfx, vector<string> &files);
			long splitFile(string &outpfx, string &inFile, vector<string> &outFiles);
			void mergeAndSortFilesByNumberTask(string inFile1, string inFile2, string outFile, int *result);
			void mergeAndSortFilesByAsciiTask(string inFile1, string inFile2, string outFile, int *result);
			int mergeTmpFilesAndPrint(FqInfo &globleInfo);
			void filterRawTags(map<string, int> &rawSequence, FqInfo &info);
			int printFiles(string &numberTmpFile, FqInfo &globleInfo, long length1);
			string intToString(int i);

			//void getTiles(string tiles, set<int> &tileSet);

		private:
			map<string,int> rawSequence_;
			map<string,int> shortSequence_;
			vector<string> tmpFile_;


			string fqFile1_;
			string fqFile2_;

			char * key_;

			bool strict_;
			int filterTmpNumber_;
			int seedLength_;
			int minInsertSize_;
			int minSize_;
			string adapter5_;
			string adapter3_;
			string outFqName_;
			bool mrna_;
			bool fastq_;
			long memLimit_;
			int threadNum_;
			QualitySystem qualSystem_;
			int readLen_;
			string outDir_;
			bool onlyStat_;
			string outPfx_;
			int nSeed_;
			int nRead_;
			int lowQualSeed1_;
			int lowQualSeed1Num_;
			int lowQualSeed2_;
			int lowQualSeed2Num_;
			float error_;
			bool removeBad_;
			bool trim_;
			float filterPolyN_;
			float filterPolyA_;
			bool filterLowQual_;
			int continuousAlign_;
			float alignRate_;
			int miniAlignLength_;
			float errorRate_;
			int misMatch_;
			bool filterIndex_;
			bool outfq_;
			int headTrim_;
			int tailTrim_;
			unsigned long cutOff_;
			long tagPercent[7];
			int score_;
			int seqType_;
			pool pl_;
			bool filterTile_;
			set<string> tiles_;
	};

} // namespace SRNAProcessTool
#endif /* SRNAPROCESSTOOL_RNAPROCESSOR_H_ */
