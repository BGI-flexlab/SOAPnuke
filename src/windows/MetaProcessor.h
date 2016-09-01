/* MetaProcessor.h
 *
 * Created on: 2012-12-19
 *     Author: Shuai Jiang
 *       Mail: jiangshuai@genomics.cn
 */

#ifndef PREPROCESSTOOL_3_LOCAL_METAPROCESSOR_H_
#define PREPROCESSTOOL_3_LOCAL_METAPROCESSOR_H_

#include "MetaCommon.h"
#include "Common.h"
namespace MetaPreProcessTool
{
	typedef struct TaskParam
	{
		PreProcessTool::Read *reads1;
		PreProcessTool::Read *reads2;
		int left, right;
		FqInfo info1, info2;

		TaskParam()
		{
			reads1 = NULL;
			reads2 = NULL;
		}
	} TaskParam;

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

	typedef struct StatisResult
	{
		bool hasAdpt;
		bool isLowQual;
		bool nExceed;
		bool isPolyN;

		StatisResult()
		{
			hasAdpt = false;
			isLowQual = false;
			nExceed = false;
			isPolyN = false;
		}
	}StatisResult;

	const int MEM_UNIT = 1024 * 1024;
	class MetaProcessor
	{
		public:
			MetaProcessor();
			static void printUsage();
			static void printVersion();
			int filterMeta(int argc, char** argv);
		private:
			bool hasAdapter(set<string> &readsName, const char *seqName);
			bool hasAdapter(const char *sequence, int readLen, const char *adapter, int adptLen);
			bool isQualityCutOff(const char*quality, int left, int right);
			bool isAlignLengthOK(char *line);
			bool statisticsPE(PreProcessTool::Read *read1, PreProcessTool::Read* read2, FqInfo *info1, FqInfo* info2, int index);
			int get3PLengthFromInfo(FqInfo &info);
			int getTrim3PLength(string fqFile, int threadNum);
			int getTrim5PPosition(const char* quality, int readLen);
			int getReadsNameFromFile(string filename, set<string> &readsName);
			int processParams(int argc, char** argv);
			int printMetaStats(FqInfo &info1, FqInfo &info2);
			StatisInfo auxStatistics(PreProcessTool::Read *read, int trimLeft, string adapter, int adptLen, set<string> &readsName, FqInfo &info, StatisResult &sr);
			void getTrim3PLengthTask(TaskParam *param);
			void statisticsTmp(PreProcessTool::Read *read, FqInfo &info);
			void task(TaskParam *param);
			void outputCleanDataTask(gzFile &file, PreProcessTool::Read *reads);
			void outputRawDataTask(gzFile &file, PreProcessTool::Read *reads, unsigned int size);
			void outputCleanDataTaskNew(gzFile &file, PreProcessTool::Read *reads, int size, int *result, int type);
			void outputSingleDataTask(gzFile &file, PreProcessTool::Read *reads1, PreProcessTool::Read *reads2, int size, int *result);
		private:
			bool isAdptList_;
			bool removeIndex_;
			bool trim_;

			int adapterLen1_;
			int adapterLen2_;
			int lengthThreshold_;
			int minAlignLength_;
			int matchNumber_;
			int nBaseNumber_;
			int qualityThreshold_;
			int trimLeft1_;
			int trimLeft2_;

			unsigned int PROCESS_THREAD_NUM;
			unsigned int *cleanDataIndexs_;
			int *singleResult_;
			long memLimit_;
			unsigned long cleanBaseNum_;

			float lowQualityRate_;
			float misMatchRate_;
			float polyN_;

			QualitySystem cleanQualSys_;
			QualitySystem qualSys_;

			boost::threadpool::pool pl_;

			boost::uint32_t doneNum_;
			boost::uint32_t size_;

			boost::mutex printSingleMutex_;
			boost::mutex sizeMutex_;
			gzFile outCleanFileS_;

			string adapter1_;
			string adapter2_;
			string fqFile1_;
			string fqFile2_;
			string outCleanPfx_;
			string outDir_;
			string outRawPfx_;

			set<string> readsName1_;
			set<string> readsName2_;

			bool filterTile_;
			set<int> tiles_;

			int seqType_;//测序fa name类型
			int polyAType_;//poly A过滤方案

	};
}//namespace
#endif
