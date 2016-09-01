/*
 *MetaProcessor.cpp
 *
 * Created on: 2012-12-15
 *     Author: Shuai Jiang
 *     Mail  : jiangshuai@genomics.cn
 */

#include "MetaProcessor.h"
#include "PeBuffer.h"
#include "FqBuffer.h"
#include "threadpool.hpp"
#include <boost/thread.hpp>
#include <boost/interprocess/detail/atomic.hpp>
using namespace boost;
using namespace boost::threadpool;
using namespace boost::interprocess::ipcdetail;

namespace MetaPreProcessTool
{
	void MetaProcessor::printVersion()
	{
		cout << "soapnuke filterMeta tools version 1.5.0\n";
		cout << "Author : jiangshuai\n";
		cout << "Mail   : jiangshuai@genomics.cn" << endl;
	}

	void MetaProcessor::printUsage()
	{
		cout << "Usage: FILTERMETA [OPTION] ... \n";
		cout << "  -f, --adapter1  : <s> 3'adapter sequence or adapter list file of fq1 file, (default: [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC])\n";
		cout << "  -r, --adapter2  : <s> 5'adapter sequence or adapter list file of fq2 file, (default: [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA])\n";
		cout << "  -1, --fq1       : <s> fq1 file, must set\n";
		cout << "  -2, --fq2       : <s> fq2 file, must set\n";
		cout << "  -c, --outCleanPfx    : <s> out clean fq files prefix name, must set\n";
		cout << "  --tile               : <s> tile number to ignore reads , such as [1101-1104,1205]\n";
		//		cout << "  -w, --outRawPfx      : <s> out raw fq files prefix name\n\n";

		cout << "  -o, --outDir         : <s> out directory, (default, [./]\n";
		cout << "  -Q, --qualitySystem  : <i> quality system of raw fq, 1:illumina, 2:sanger (default: [1])\n";
		cout << "  -S, --sanger         : <b> set clean data qualtiy system to sanger (default: illumina)\n";
		cout << "  -L, --alignLength    : <i> the min align length when match the adapter, (default: [15])\n";
		//		cout << "  -M, --misMatchRate  : <f> the max mismatch rate when match the adapter, !only for adapter sequence!, (default: [0.2])\n";
		//		cout << "  -N, --matchNumber   : <i> the min continuous match number when match the adapter, !only for adapter sequence!, (default: [10])\n";

		cout << "  -N, --nBaseNumber    : <i> filter reads contain N, (default, [3])\n";
		cout << "  -P, --polyN          : <f> filter polyN[A, T, G, C], 0 means do not filter, (default: [1.0])\n\n";
		cout << "  -U, --unTrim         : <b> don't trim reads (default, [off])\n";
		cout << "  when --untrim is off: the next 3 options are useful\n";
		cout << "    -q, --qualityThreshold : <i> quality threshold, (default: [20])\n";
		cout << "    -l, --lengthThreshold  : <i> length threshold, (default: [30])\n";
		cout << "    -R, --lowQualityRate   : <f> low quality(qual < -q) rate, (default: [0.5])\n\n";
		cout << "  -i, --index     : <b> remove index\n";
		cout << "  -C, --cut       : <f> the bases number you want to keep in all clean fq files\n";
		cout << "                        (unit:1024*1024*1024, 0 means not cut reads)\n";
		cout << "\n";
		cout << "\t-5, --seqType   : <i> Sequence fq name type, 0->old fastq name, 1->new fastq name[default: 0]\n";
		cout << "\t    old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
		cout << "\t    new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
		//	cout << "\t-6, --polyAType : <i> filter poly A type, 0->both two reads are poly a, 1->at least one reads is poly a, then filter, [default: 0]\n";

		cout << "  -a, --append    : <s> the log's output place : console or file (default: [console])\n";
		cout << "  -h, --help      : <b> help\n";
		cout << "  -v, --version   : <b> show version" << endl;
	}

	MetaProcessor::MetaProcessor():
		isAdptList_(false),
		removeIndex_(false),
		trim_(true),

		lengthThreshold_(30),
		minAlignLength_(15),
		matchNumber_(15),
		nBaseNumber_(3),
		qualityThreshold_(20),
		trimLeft1_(0),
		trimLeft2_(0),

		PROCESS_THREAD_NUM(4),
		cleanDataIndexs_(NULL),
		singleResult_(NULL),
		memLimit_(150 * MEM_UNIT),
		cleanBaseNum_(0),

		lowQualityRate_(0.5),
		misMatchRate_(0.2),
		polyN_(1.0),
		cleanQualSys_(ILLUMINA_),
		qualSys_(ILLUMINA_), 

		pl_(4),
		outCleanFileS_(NULL),
		adapter1_("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"),
		adapter2_("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"),
		outDir_("."),
		filterTile_(false),
		seqType_(0),
		polyAType_(0)
		{}

	int MetaProcessor::processParams(int argc, char** argv)
	{
		const char* shortOptions = "f:r:1:2:c:K:o:Q:SL:N:P:Uq:l:R:iC:a:5:6:hv";
		const struct option longOptions[] = 
		{
			{"adapter1",         1, NULL, 'f'},
			{"adapter2",         1, NULL, 'r'},
			{"fq1",              1, NULL, '1'},
			{"fq2",              1, NULL, '2'},
			{"outCleanPfx",      1, NULL, 'c'},
			{"tile",              1, NULL, 'K'},
			//			{"outRawPfx",        1, NULL, 'w'},
			{"outDir",           1, NULL, 'o'},
			{"qualitySystem",    1, NULL, 'Q'},
			{"sanger",           0, NULL, 'S'},
			{"alignLength",      1, NULL, 'L'},
			{"nBaseNumber",      1, NULL, 'N'},
			{"polyN",            1, NULL, 'P'},
			{"unTrim",           0, NULL, 'U'},
			{"qualityThreshold", 1, NULL, 'q'},
			{"lengthThreshold",  1, NULL, 'l'},
			{"lowQualityRate",   1, NULL, 'R'},
			{"index",            0, NULL, 'i'},
			{"cut",              1, NULL, 'C'},
			{"append",           1, NULL, 'a'},
			{ "seqType"  , 1, NULL, '5' },
			{ "polyAType" , 1, NULL, '6' },
			{"help",             0, NULL, 'h'},
			{"version",          0, NULL, 'v'},
		};

		if(argc==1)
		{
			printUsage();
			return 1;
		}

		string append;
		string tiles;
		int nextOpt;
		float num;
		while(-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions, NULL)))
		{
			switch(nextOpt)
			{
				case 'f':
					adapter1_.assign(optarg);
					break;
				case 'r':
					adapter2_.assign(optarg);
					break;
				case '1':
					fqFile1_.assign(optarg);
					break;
				case '2':
					fqFile2_.assign(optarg);
					break;
				case 'c':
					outCleanPfx_.assign(optarg);
					break;
				case 'K':
					filterTile_ = true;
					tiles.assign(optarg);
					PreProcessTool::getTiles(tiles, tiles_);
					break;
					//				case 'w':
					//					outRawPfx_.assign(optarg);
					//					break;
				case 'o':
					outDir_.assign(optarg);
					break;

				case 'Q':
					switch(optarg[0])
					{
						case '1':
							qualSys_ = ILLUMINA_;
							break;
						case '2':
							qualSys_ = SANGER_;
							break;
						default:
							cerr << "error quality system" <<endl;
							return 1;
					}
					break;
				case 'S':
					cleanQualSys_ = SANGER_;
					break;
				case 'L':
					minAlignLength_ = atoi(optarg);
					break;
				case 'N':
					nBaseNumber_ = atoi(optarg);
					break;
				case 'P':
					polyN_ = atof(optarg);
					break;
				case 'U':
					trim_ = false;
					break;
				case 'q':
					qualityThreshold_ = atoi(optarg);
					break;
				case 'l':
					lengthThreshold_ = atoi(optarg);
					break;
				case 'R':
					lowQualityRate_ = atof(optarg);
					break;
				case 'i':
					removeIndex_ = true;
					break;
				case 'C':
					num = atof(optarg);
					if(num < 1E-6)
					{
						cleanBaseNum_ = 0;
					}
					else
					{
						cleanBaseNum_ = (unsigned long)(num * 1024 * 1024 * 1024);
					}
					break;
				case 'a':
					append = optarg;
					break;
				case '5':
					seqType_ = atoi(optarg);
					break;
				case '6':
					polyAType_ = atoi(optarg);
					break;
				case 'h':
					printUsage();
					return 1;
				case 'v':
					printVersion();
					return 1;
				case '?':
					cout << "unkonwn option: -" << (char) optopt << endl;
					cout << "Print -h or --help for more information." << endl;
					return 1;
				default:
					cout << "Param : " << optarg << endl;
					return 1;
			}
		}

		if(argc != optind)
		{
			cerr << "options error, please check the options" << endl;
			return 1;
		}

		bool isPathNotExists = false;
		if(_access(outDir_.c_str(), 0) == -1)
		{
			isPathNotExists = true;
			int len = outDir_.size();
			char *path = (char *)malloc(len + 15);
			sprintf(path, "mkdir -p %s", outDir_.c_str());
			if (std::system(path) == -1)
			{
				cerr << "output directory " << outDir_ << " cannot create" << endl;
				return 1;
			}
		}

		if(fqFile1_.empty())
		{
			cerr << "fq1 file must be exists" << endl;
			return 1;
		}

		if(fqFile2_.empty())
		{
			cerr << "fq2 file must be exists" << endl;
			return 1;
		}

		if(outCleanPfx_.empty())
		{
			cerr << "outCleanPfx must be set" << endl;
			return 1;
		}
		if(outDir_.size() > 1 && outDir_[outDir_.size()-1] == '/')
		{
			outDir_.erase(outDir_.size()-1, 1);
		}

		string logoutPath = outDir_ + "/" + "out.log";
		if(!init_logger(append, logoutPath))
		{
			cerr << "Cannot Init Log:" << append << "-" << logoutPath << endl;
			return 1;
		}
		else
		{
			LOG(INFO, "Log Init Success");
		}

		if(isPathNotExists)
		{
			LOG(WARN, "output directory " << outDir_ << " does not exists, program will auto create");
			LOG(WARN, "output directory " << outDir_ << " has been created");
		}
		/*
		   if(!outRawPfx_.empty() && rawReadNum_ == 0)
		   {
		   char buf1[1024];
		   char buf2[1024];

		   string rawFqFile1 = outDir_ + "/" + outRawPfx_ + "raw.1.fq.gz";
		   string rawFqFile2 = outDir_ + "/" + outRawPfx_ + "raw.2.fq.gz";

		   if(fqFile1_.substr(fqFile1_.size()-2, 2) != "gz")
		   {
		   sprintf(buf1, "gzip -c %s > %s", fqFile1_.c_str(), rawFqFile1.c_str());
		   }
		   else
		   {
		   sprintf(buf1, "cp -uL %s  %s", fqFile1_.c_str(), rawFqFile1.c_str());
		   }

		   if(fqFile2_.substr(fqFile2_.size()-2, 2) != "gz")
		   {
		   sprintf(buf2, "gzip -c %s > %s", fqFile2_.c_str(), rawFqFile2.c_str());
		   }
		   else
		   {
		   sprintf(buf2, "cp -uL %s  %s", fqFile2_.c_str(), rawFqFile2.c_str());
		   }

		   LOG(INFO, buf1);
		   system(buf1);

		   LOG(INFO, buf2);
		   system(buf2);
		   }
		 */
		if(adapter1_.empty() || adapter2_.empty())
		{
			LOG(ERROR, "--adapter2 or --adapter2 should not be NULL");
			return 1;
		}
		else
		{
			int type = adapterType(true, adapter1_, adapter2_);
			if(type == 3)
			{
				LOG(ERROR, "adapter1 or adapter2 error");
				return 1;
			}
			else if(type == 2)
			{
				isAdptList_ = true;
			}
			else
			{
				isAdptList_ = false;
				for(int i=adapter1_.size()-1; i>=0; i--)
				{
					adapter1_[i] = toupper(adapter1_[i]);
				}
				for(int i=adapter2_.size()-1; i>=0; i--)
				{
					adapter2_[i] = toupper(adapter2_[i]);
				}
				adapterLen1_ = adapter1_.size();
				adapterLen2_ = adapter2_.size();
			}

		}
		return 0;		
	}

	int MetaProcessor::filterMeta(int argc, char** argv)
	{
		if(processParams(argc, argv) != 0)
		{
			return 1;
		}

		int threadNum = PROCESS_THREAD_NUM;
		if(trim_)
		{
			trimLeft1_ = getTrim3PLength(fqFile1_, threadNum);
			if(trimLeft1_ < 0)
			{
				LOG(ERROR, "we can't find trim number for 3 primer " << fqFile1_);
				return 1;
			}
			LOG(INFO,"fq1: trim 3p num: "<<trimLeft1_);

			trimLeft2_ = getTrim3PLength(fqFile2_, threadNum);
			if(trimLeft2_ < 0)
			{
				LOG(ERROR, "we can't find trim number for 3 primer " << fqFile2_);
				return 1;
			}
			LOG(INFO,"fq2: trim 3p num: "<<trimLeft2_);
		}
		FqInfo globleInfo1, globleInfo2;
		gzFile outRawFile1 = NULL, outRawFile2 = NULL, outCleanFile1 = NULL, outCleanFile2 = NULL;//, outCleanFileS = NULL;
		string outRawFileName1, outRawFileName2, outCleanFileName1, outCleanFileName2, outCleanFileNameS;
		/*
		   if(rawReadNum_ > 0)
		   {
		   if(outRawPfx_.empty())
		   {
		   outRawFileName1 = getOutputFileName(fqFile1_, "raw", outDir_);
		   outRawFileName2 = getOutputFileName(fqFile2_, "raw", outDir_);
		   }
		   else
		   {
		   outRawFileName1 = outDir_ + "/" + outRawPfx_ + ".raw.1.fq.gz";
		   outRawFileName2 = outDir_ + "/" + outRawPfx_ + ".raw.2.fq.gz";
		   }
		   outRawFile1 = gzopen(outRawFileName1.c_str(), "wb");
		   outRawFile2 = gzopen(outRawFileName2.c_str(), "wb");
		   if(!outRawFile1)
		   {
		   LOG(ERROR, "creat output file: " + outRawFileName1);
		   return 1;
		   }
		   if(!outRawFile2)
		   {
		   LOG(ERROR, "creat output file: " + outRawFileName2);
		   return 1;
		   }
		   }
		 */
		outCleanFileName1 = outDir_ + "/" + outCleanPfx_ + ".pair.1.fq.gz";
		outCleanFileName2 = outDir_ + "/" + outCleanPfx_ + ".pair.2.fq.gz";
		outCleanFileNameS = outDir_ + "/" + outCleanPfx_ + ".single.fq.gz";

		outCleanFile1 = gzopen(outCleanFileName1.c_str(), "wb");
		outCleanFile2 = gzopen(outCleanFileName2.c_str(), "wb");
		outCleanFileS_ = gzopen(outCleanFileNameS.c_str(), "wb");

		if(!outCleanFile1)
		{
			LOG(ERROR, "create output file: " + outCleanFileName1);
			return 1;
		}
		if(!outCleanFile2)
		{
			LOG(ERROR, "create output file: " + outCleanFileName2);
			return 1;
		}
		if(!outCleanFileS_)
		{
			LOG(ERROR, "create output file: " + outCleanFileNameS);
			return 1;
		}

		if(isAdptList_)
		{
			if(getReadsNameFromFile(adapter1_, readsName1_) != 0)
			{
				LOG(ERROR, "read file: " + adapter1_);
				return 1;
			}
			if(getReadsNameFromFile(adapter2_, readsName2_) != 0)
			{
				LOG(ERROR, "read file: " + adapter2_);
				return 1;
			}
			LOG(INFO, "read adapter.list file OK");
		}
		PreProcessTool::Read *reads1, *reads2;
		long capacity = static_cast<long>(memLimit_ / 2.5);
		PreProcessTool::PeBuffer buffer(fqFile1_.c_str(), fqFile2_.c_str(), capacity, PreProcessTool::PeBuffer::RB, filterTile_, tiles_);
		buffer.setSeqType(seqType_);
		TaskParam *params = new TaskParam[PROCESS_THREAD_NUM];

		int size=0;
		//		unsigned long sizeTmp=0;
		//		bool first = true;
		while(buffer.getReads()> 0)
		{
			size = buffer.getReadSize();
			if(size == 0)
				continue;
			//			sizeTmp += size;
			//			if(rawReadNum_>0 && sizeTmp>rawReadNum_)
			//			{
			//				size = size - (sizeTmp - rawReadNum_);
			//			}
			//			if(first)
			//			{
			//				first = false;
			//				cleanDataIndexs_ = new unsigned int[(unsigned int)(buffer.getInitReadSize()*1.1 + 1)];
			//			}
			singleResult_ = new int[size + 1];
			bzero(singleResult_, sizeof(int)*(size + 1));
			reads1 = buffer.getReadsOne();
			reads2 = buffer.getReadsTwo();
			/*
			   if(rawReadNum_ > 0)
			   {
			   pl_.schedule(boost::bind(&MetaProcessor::outputRawDataTask, this, outRawFile1, reads1, size));
			   pl_.schedule(boost::bind(&MetaProcessor::outputRawDataTask, this, outRawFile2, reads2, size));
			   pl_.wait();
			   }
			 */
			size_ = 0;
			doneNum_ = 0;
			int block = size / PROCESS_THREAD_NUM;
			int remain = size % PROCESS_THREAD_NUM;
			int index = 0;

			//			pl_.schedule(boost::bind(&MetaProcessor::outputCleanDataTask, this, outCleanFile1, reads1));
			//			pl_.schedule(boost::bind(&MetaProcessor::outputCleanDataTask, this, outCleanFile2, reads2));
			for(unsigned int i=0; i<PROCESS_THREAD_NUM; i++)
			{
				params[i].left = index;
				index += block;
				if(remain > 0)
				{
					index += 1;
					remain--;
				}
				params[i].right = index;
				params[i].reads1 = reads1;
				params[i].reads2 = reads2;
				bzero(&(params[i].info1), sizeof(FqInfo));
				bzero(&(params[i].info2), sizeof(FqInfo));
				pl_.schedule(boost::bind(&MetaProcessor::task, this, &params[i]));
			}
			pl_.wait();
			//			for (int i=0; i<size; i++)
			//			{
			//				if(singleResult_[i] == 1)
			//				{
			//					gzputs(outCleanFileS_, reads1[i].readName);
			//					gzputs(outCleanFileS_, "\n");
			//					gzputs(outCleanFileS_, reads1[i].baseSequence);
			//					gzputs(outCleanFileS_, "\n");
			//					gzputs(outCleanFileS_, reads1[i].optionalName);
			//					gzputs(outCleanFileS_, "\n");
			//					gzputs(outCleanFileS_, reads1[i].baseQuality);
			//					gzputs(outCleanFileS_, "\n");
			//				}
			//				else if(singleResult_[i] == 2)
			//				{
			//					gzputs(outCleanFileS_, reads2[i].readName);
			//					gzputs(outCleanFileS_, "\n");
			//					gzputs(outCleanFileS_, reads2[i].baseSequence);
			//					gzputs(outCleanFileS_, "\n");
			//					gzputs(outCleanFileS_, reads2[i].optionalName);
			//					gzputs(outCleanFileS_, "\n");
			//					gzputs(outCleanFileS_, reads2[i].baseQuality);
			//					gzputs(outCleanFileS_, "\n");
			//				}
			//			}
			pl_.schedule(boost::bind(&MetaProcessor::outputCleanDataTaskNew, this, outCleanFile1, reads1, size, singleResult_, 3));
			pl_.schedule(boost::bind(&MetaProcessor::outputCleanDataTaskNew, this, outCleanFile2, reads2, size, singleResult_, 3));
			pl_.schedule(boost::bind(&MetaProcessor::outputSingleDataTask, this, outCleanFileS_, reads1, reads2, size, singleResult_));
			pl_.wait();
			for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
			{
				if(globleInfo1.cleanReadLengthTmp < params[i].info1.cleanReadLengthTmp)
				{
					globleInfo1.cleanReadLengthTmp = params[i].info1.cleanReadLengthTmp;
				}
				if(globleInfo2.cleanReadLengthTmp < params[i].info2.cleanReadLengthTmp)
				{
					globleInfo2.cleanReadLengthTmp = params[i].info2.cleanReadLengthTmp;
				}
				if (params[i].info1.rawReadLength > globleInfo1.rawReadLength)
				{
					globleInfo1.rawReadLength = params[i].info1.rawReadLength;
				}
				if (params[i].info2.rawReadLength > globleInfo2.rawReadLength)
				{
					globleInfo2.rawReadLength = params[i].info2.rawReadLength;
				}
				if (params[i].info1.maxQualityValue > globleInfo1.maxQualityValue)
				{
					globleInfo1.maxQualityValue = params[i].info1.maxQualityValue;
				}
				if (params[i].info2.maxQualityValue > globleInfo2.maxQualityValue)
				{
					globleInfo2.maxQualityValue = params[i].info2.maxQualityValue;
				}
				globleInfo1.add(params[i].info1);
				globleInfo2.add(params[i].info2);
				//calculte the max quality value

			}
			if(singleResult_ != NULL)
			{
				delete []singleResult_;
				singleResult_ = NULL;
			}

			if(cleanBaseNum_ > 0 && globleInfo1.cleanTotalBaseNum + globleInfo2.cleanTotalBaseNum >= cleanBaseNum_)
				break;
		}
		if (globleInfo1.maxQualityValue > globleInfo2.maxQualityValue)
		{
			globleInfo2.maxQualityValue = globleInfo1.maxQualityValue;
		}
		else
		{
			globleInfo1.maxQualityValue = globleInfo2.maxQualityValue;
		}
		if(globleInfo1.cleanReadLengthTmp < globleInfo2.cleanReadLengthTmp)
		{
			globleInfo1.cleanReadLengthTmp = globleInfo2.cleanReadLengthTmp;
		}

		if (size == -1)
		{
			LOG(ERROR, "read fq1 or fq2 error");
			return 1;
		}

		if (size == -2)
		{
			LOG(WARN, "fq1 and fq2's read number not equal");
		}

		if (params != NULL)
		{
			LOG(INFO, "delete params");
			delete []params;
		}
		/*
		   if (rawReadNum_)
		   {
		   gzclose(outRawFile1);
		   gzclose(outRawFile2);
		   }
		 */
		gzclose(outCleanFile1);
		gzclose(outCleanFile2);
		gzclose(outCleanFileS_);
		printFqInfo(outDir_, "", &globleInfo1, &globleInfo2);
		if(printMetaStats(globleInfo1, globleInfo2) > 0)
		{
			return 1;
		}
		LOG(INFO, "PreProcess Finish");
		return 0;

	}

	void MetaProcessor::task(TaskParam *param)
	{
		PreProcessTool::Read *reads1 = param->reads1;
		PreProcessTool::Read *reads2 = param->reads2;
		int start = param->left;
		int end = param->right;
		bool isClean;
		if(start == end)
		{
			atomic_inc32(&doneNum_);
			return;
		}

		for(int i=start; i<end; i++)
		{
			isClean = statisticsPE(&reads1[i], &reads2[i], &(param->info1), &(param->info2), i);
			//			if(isClean)
			//			{
			//				mutex::scoped_lock(sizeMutex_);
			//				cleanDataIndexs_[atomic_read32(&size_)] = i;
			//				atomic_inc32(&size_);
			//			}
		}
		//		atomic_inc32(&doneNum_);
		return;
	}

	bool MetaProcessor::statisticsPE(PreProcessTool::Read *read1, PreProcessTool::Read* read2, FqInfo *info1, FqInfo* info2, int index)
	{
		StatisResult sr1, sr2;
		StatisInfo si1 = auxStatistics(read1, trimLeft1_, adapter1_, adapterLen1_, readsName1_, *info1, sr1);
		StatisInfo si2 = auxStatistics(read2, trimLeft2_, adapter2_, adapterLen2_, readsName2_, *info2, sr2);
		bool isClean1=false, isClean2=false;
		if(sr1.hasAdpt)
		{
			info1->adapterNum++;
		}
		else if(sr1.nExceed)
		{
			info1->nExceedNum++;
		}
		else if(sr1.isPolyN)
		{
			info1->polyNNum++;
		}
		else
		{
			if(trim_)
			{
				if(si1.readLen < lengthThreshold_)
				{
					info1->smallInsertNum++;
				}
				else if(sr1.isLowQual)
				{
					info1->lowQualNum++;
				}
				else
				{
					isClean1 = true;
				}
			}
			else
			{
				isClean1 = true;
			}
		}

		if(sr2.hasAdpt)
		{
			info2->adapterNum++;
		}
		else if(sr2.nExceed)
		{
			info2->nExceedNum++;
		}
		else if(sr2.isPolyN)
		{
			info2->polyNNum++;
		}
		else
		{
			if(trim_)
			{
				if(si2.readLen < lengthThreshold_)
				{
					info2->smallInsertNum++;
				}
				else if(sr2.isLowQual)
				{
					info2->lowQualNum++;
				}
				else
				{
					isClean2 = true;
				}
			}
			else
			{
				isClean2 = true;
			}
		}

		if(isClean1 && isClean2)
		{
			info1->cleanBaseA += si1.a;
			info1->cleanBaseC += si1.c;
			info1->cleanBaseG += si1.g;
			info1->cleanBaseT += si1.t;
			info1->cleanBaseN += si1.n;
			info1->cleanQ20 += si1.q20;
			info1->cleanQ30 += si1.q30;
			info1->cleanTotalReadNum++;
			info1->cleanTotalBaseNum += si1.readLen;

			if(info1->cleanReadLengthTmp < si1.readLen)
			{
				info1->cleanReadLengthTmp = si1.readLen;
			}
			//			calculateBaseDistribute(read1, *info1, si1.readLen);

			info2->cleanBaseA += si2.a;
			info2->cleanBaseC += si2.c;
			info2->cleanBaseG += si2.g;
			info2->cleanBaseT += si2.t;
			info2->cleanBaseN += si2.n;
			info2->cleanQ20 += si2.q20;
			info2->cleanQ30 += si2.q30;
			info2->cleanTotalReadNum++;
			info2->cleanTotalBaseNum += si2.readLen;
			if(info2->cleanReadLengthTmp < si2.readLen)
			{
				info2->cleanReadLengthTmp = si2.readLen;
			}
			//			calculateBaseDistribute(read2, *info2, si2.readLen);
			singleResult_[index] = 3;
			return true;
		}
		else
		{
			if(sr1.hasAdpt || sr2.hasAdpt)
			{
				info1->totalAdapterNum++;
			}
			else if(sr1.nExceed || sr2.nExceed)
			{
				info1->totalNExceedNum++;
			}
			else if(sr1.isPolyN || sr2.isPolyN)
			{
				info1->totalPolyNNum++;
			}
			else if(si1.readLen < lengthThreshold_ || si2.readLen < lengthThreshold_)
			{
				info1->totalSmallInsertNum++;
			}
			else
			{
				info1->totalLowQualNum++;
			}
			if(isClean1)
			{
				info1->singleReadNum++;
				info1->singleBaseNum += si1.readLen;
				singleResult_[index] = 1;
				//				{
				//					mutex::scoped_lock(printSingleMutex_);
				//					gzputs(outCleanFileS_, read1->readName);
				//					gzputs(outCleanFileS_, "\n");
				//					gzputs(outCleanFileS_, read1->baseSequence);
				//					gzputs(outCleanFileS_, "\n");
				//					gzputs(outCleanFileS_, read1->optionalName);
				//					gzputs(outCleanFileS_, "\n");
				//					gzputs(outCleanFileS_, read1->baseQuality);
				//					gzputs(outCleanFileS_, "\n");
				//				}
			}
			if(isClean2)
			{
				info2->singleReadNum++;
				info2->singleBaseNum += si2.readLen;
				singleResult_[index] = 2;
				//				{
				//					mutex::scoped_lock(printSingleMutex_);
				//					gzputs(outCleanFileS_, read2->readName);
				//					gzputs(outCleanFileS_, "\n");
				//					gzputs(outCleanFileS_, read2->baseSequence);
				//					gzputs(outCleanFileS_, "\n");
				//					gzputs(outCleanFileS_, read2->optionalName);
				//					gzputs(outCleanFileS_, "\n");
				//					gzputs(outCleanFileS_, read2->baseQuality);
				//					gzputs(outCleanFileS_, "\n");
				//				}
			}
			return false;
		}
		return false;
	}

	StatisInfo MetaProcessor::auxStatistics(PreProcessTool::Read *read, int trimLeft, string adapter, int adptLen, set<string> &readsName, FqInfo &info, StatisResult &sr)
	{
		int qual;
		int readLen = strlen(read->baseSequence);

		info.rawTotalBaseNum += readLen;
		info.rawTotalReadNum++;

		if(isAdptList_)
		{
			sr.hasAdpt = hasAdapter(readsName, read->readName);
		}
		else
		{
			sr.hasAdpt = hasAdapter(read->baseSequence, readLen, adapter.c_str(), adptLen);
		}

		int trimRight = readLen;
		if(trim_)
		{
			trimRight = getTrim5PPosition(read->baseQuality, readLen);
			//			sr.lengthCutOff = ((trimRight - trimLeft) < lengthThreshold_);
			//			sr.qualityCutOff = isQualityCutOff(read->baseQuality, trimLeft, trimRight);
		}
		StatisInfo si;
		int a = 0, g = 0, c = 0, t = 0, n = 0;
		int q20 = 0, q30 = 0;

		for (int i=0; i<readLen; ++i)
		{
			switch (read->baseSequence[i])
			{
				case 'A':
					a++;
					info.base[i][A_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.a++;
					}
					break;
				case 'a':
					read->baseSequence[i] = 'A';
					a++;
					info.base[i][A_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.a++;
					}
					break;			
				case 'C':
					c++;
					info.base[i][C_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.c++;
					}
					break;
				case 'c':
					read->baseSequence[i] = 'C';
					c++;
					info.base[i][C_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.c++;
					}
					break;
				case 'G':
					g++;
					info.base[i][G_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.g++;
					}
					break;
				case 'g':
					read->baseSequence[i] = 'G';
					g++;
					info.base[i][G_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.g++;
					}
					break;
				case 'T':
					t++;
					info.base[i][T_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.t++;
					}
					break;
				case 't':
					read->baseSequence[i] = 'T';
					t++;
					info.base[i][T_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.t++;
					}
					break;
				case 'N':
					n++;
					info.base[i][N_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.n++;
					}
					break;
				case 'n':
					read->baseSequence[i] = 'N';
					n++;
					info.base[i][N_]++;
					if (i>=trimLeft && i<trimRight)
					{
						si.n++;
					}
					break;

			}

			qual = read->baseQuality[i] - qualSys_;
			if (qual > MAX_QUALITY)
			{
				LOG(WARN, "some bases' quality larger than " << MAX_QUALITY << ", they have been set to " << MAX_QUALITY);
				qual = MAX_QUALITY;
			}

			if (qual < 0)
			{
				LOG(WARN, "some bases' quality smaller than 0, they have been set to 0");
				qual = 0;
			}

			if (qual > info.maxQualityValue)
			{
				info.maxQualityValue = qual;
			}


			read->baseQuality[i] = qual + cleanQualSys_;

			info.qual[i][qual]++;

			if (qual >= 20)
			{
				q20++;
				info.q20q30[i][0]++;
				if (qual >= 30)
				{
					q30++;
					info.q20q30[i][1]++;
				}
			}

			if (i>=trimLeft && i<trimRight)
			{
				if (qual < qualityThreshold_)
				{
					si.lowQual++;
				}
				if (qual >=20)
				{
					si.q20++;
					if (qual >= 30)
					{
						si.q30++;
					}
				}
			}
		}

		info.rawBaseA += a;
		info.rawBaseC += c;
		info.rawBaseG += g;
		info.rawBaseT += t;
		info.rawBaseN += n;
		info.rawQ20 += q20;
		info.rawQ30 += q30;

		if(info.rawReadLength < readLen)
		{
			info.rawReadLength = readLen;
		}

		si.readLen = trimRight - trimLeft;
		sr.nExceed = (n >= nBaseNumber_);
		sr.isLowQual = (si.lowQual > si.readLen * lowQualityRate_);
		if(polyN_ >= 1E-6)
		{
			sr.isPolyN= (1.0*a/readLen>=(polyN_-1E-6) || 1.0*c/readLen>=(polyN_-1E-6) || 1.0*g/readLen>=(polyN_-1E-6) || 1.0*t/readLen>=(polyN_-1E-6));
		}

		if(trim_ && trimRight>trimLeft)
		{
			read->baseSequence[trimRight] = '\0';
			read->baseSequence = read->baseSequence + trimLeft;
			read->baseQuality[trimRight] = '\0';
			read->baseQuality = read->baseQuality + trimLeft;
		}

		if(removeIndex_)
		{
			int sharpIndex = 0;
			int i=0;
			while (read->readName[i++] != '#')
				;
			sharpIndex = i;
			while (read->readName[i++] != '/')
				;
			strcpy(read->readName + sharpIndex, read->readName + i - 1);
		}

		return si;
	}

	bool MetaProcessor::isQualityCutOff(const char*quality, int left, int right)
	{
		if(right <= left)
		{
			return true;
		}
		int lowQual = 0;
		for(int i=left; i<right; i++)
		{
			if(quality[i] - qualSys_ < qualityThreshold_)
			{
				lowQual++;
			}
		}
		if(1.0*lowQual/(right-left)<lowQualityRate_)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	int MetaProcessor::getTrim5PPosition(const char* quality, int readLen)
	{
		for(int i=readLen-1; i>=0; i--)
		{
			if(quality[i] - qualSys_ >= qualityThreshold_)
			{
				return i+1;
			}
		}
		return 0;
	}
	int MetaProcessor::getTrim3PLength(string fqFile, int threadNum)
	{
		FqInfo tmpFqInfo;
		long capacity = static_cast<long>(memLimit_ / 2.5);
		PreProcessTool::FqBuffer buffer(fqFile.c_str(), capacity, PreProcessTool::FqBuffer::RB, filterTile_, tiles_);
		buffer.setSeqType(seqType_);
		TaskParam *params = new TaskParam[threadNum];
		PreProcessTool::Read *reads;
		int size = 0;
		//		unsigned long sizeTmp = 0;
		while((reads = buffer.getReads()))
		{
			size = buffer.getRealReadSize();
			if(size == 0)
				continue;
			if(size == -1)
			{
				LOG(ERROR, "read fq file: " << fqFile);
				return -1;
			}

			//			sizeTmp += size;
			//			if(rawReadNum_ > 0 && sizeTmp > rawReadNum_)
			//			{
			//				size = size - (sizeTmp - rawReadNum_);
			//			}

			int block = size / threadNum;
			int remain = size % threadNum;
			int index = 0;
			for( int i=0; i<threadNum; i++)
			{
				params[i].left = index;
				index += block;
				if(remain > 0)
				{
					index +=1;
					remain--;
				}
				params[i].right = index;
				params[i].reads1 = reads;
				bzero(&(params[i].info1), sizeof(FqInfo));
				pl_.schedule(boost::bind(&MetaProcessor::getTrim3PLengthTask, this, &params[i]));
			}
			pl_.wait();
			for(int i=0; i<threadNum; i++)
			{
				if(tmpFqInfo.rawReadLength < params[i].info1.rawReadLength)
				{
					tmpFqInfo.rawReadLength = params[i].info1.rawReadLength;
				}
				tmpFqInfo.add(params[i].info1);
			}
			//			if(rawReadNum_ > 0 && sizeTmp >=rawReadNum_)
			//				break;
		}
		if(params != NULL)
		{
			delete []params;
		}
		LOG(INFO, "raw readLen :" << tmpFqInfo.rawReadLength << " : " <<fqFile);
		int len = get3PLengthFromInfo(tmpFqInfo);
		if(len + lengthThreshold_ > tmpFqInfo.rawReadLength)
		{
			LOG(ERROR, "3 primer trim number of " << fqFile << "is too large: " << len);
			len = -2;
		}
		return len;
	}

	int MetaProcessor::get3PLengthFromInfo(FqInfo &info)
	{
		double averageBaseA = 1.0 * info.rawBaseA / info.rawReadLength;
		double averageBaseC = 1.0 * info.rawBaseC / info.rawReadLength;
		double averageBaseG = 1.0 * info.rawBaseG / info.rawReadLength;
		double averageBaseT = 1.0 * info.rawBaseT / info.rawReadLength;
		//		LOG(INFO, "average A: "<<averageBaseA << " C: "<<averageBaseC<< " G: "<<averageBaseG<< " T: "<<averageBaseT);

		double baseAsdSum=0.0, baseCsdSum=0.0, baseGsdSum=0.0, baseTsdSum=0.0;

		for(unsigned int i=0; i<info.rawReadLength; i++)
		{
			baseAsdSum += (info.base[i][A_] - averageBaseA) * (info.base[i][A_] - averageBaseA);
			baseCsdSum += (info.base[i][C_] - averageBaseC) * (info.base[i][C_] - averageBaseC);
			baseGsdSum += (info.base[i][G_] - averageBaseG) * (info.base[i][G_] - averageBaseG);
			baseTsdSum += (info.base[i][T_] - averageBaseT) * (info.base[i][T_] - averageBaseT);
			//			LOG(INFO, info.base[i][A_] << "-" << averageBaseA << " :^2 "<< (info.base[i][A_] - averageBaseA) * (info.base[i][A_] - averageBaseA));
			//			LOG(INFO, info.base[i][C_] << "-" << averageBaseC << " :^2 "<< (info.base[i][C_] - averageBaseC) * (info.base[i][C_] - averageBaseC));
			//			LOG(INFO, info.base[i][G_] << "-" << averageBaseG<< " :^2 "<< (info.base[i][G_] - averageBaseG) * (info.base[i][G_] - averageBaseG));
			//			LOG(INFO, info.base[i][T_] << "-" << averageBaseT << " :^2 "<< (info.base[i][T_] - averageBaseT) * (info.base[i][T_] - averageBaseT));
		}
		//		LOG(INFO, "sum A: " << baseAsdSum <<" sum C: " << baseCsdSum<<" sum G: " << baseGsdSum<<" sum T: " << baseTsdSum);
		double baseAstd = sqrt(baseAsdSum / (info.rawReadLength - 1));
		double baseCstd = sqrt(baseCsdSum / (info.rawReadLength - 1));
		double baseGstd = sqrt(baseGsdSum / (info.rawReadLength - 1));
		double baseTstd = sqrt(baseTsdSum / (info.rawReadLength - 1));
		//		LOG(INFO, "baseStd A: " << baseAstd <<" C: " << baseCstd <<" G: " << baseGstd <<" T: " << baseTstd);

		double minBaseA = averageBaseA - baseAstd * 2;
		double maxBaseA = averageBaseA + baseAstd * 2;
		double minBaseC = averageBaseC - baseCstd * 2;
		double maxBaseC = averageBaseC + baseCstd * 2;
		double minBaseG = averageBaseG - baseGstd * 2;
		double maxBaseG = averageBaseG + baseGstd * 2;
		double minBaseT = averageBaseT - baseTstd * 2;
		double maxBaseT = averageBaseT + baseTstd * 2;
		//		LOG(INFO, "A: " << minBaseA << "," << maxBaseA << " C :" << minBaseC << "," << maxBaseC << " G :" << minBaseG << "," << maxBaseG << " T :" << minBaseT << "," << maxBaseT);

		for(unsigned int i=0; i<info.rawReadLength; i++)
		{
			if(info.base[i][A_]<=maxBaseA
					&& info.base[i][A_]>=minBaseA
					&& info.base[i][C_]<=maxBaseC
					&& info.base[i][C_]>=minBaseC
					&& info.base[i][G_]<=maxBaseG
					&& info.base[i][G_]>=minBaseG 
					&& info.base[i][T_]<=maxBaseT
					&& info.base[i][T_]>=minBaseT)
			{
				return i;
			}
		}
		return -1;
	}

	void MetaProcessor::getTrim3PLengthTask(TaskParam *param)
	{
		PreProcessTool::Read *read = param->reads1;
		int start = param->left;
		int end = param->right;
		if( start == end)
			return;

		for(int i=start; i<end; i++)
		{
			statisticsTmp(&read[i], param->info1);
		}
		return;
	}

	void MetaProcessor::statisticsTmp(PreProcessTool::Read *read, FqInfo &info)
	{
		int readLen = strlen(read->baseSequence);
		if(info.rawReadLength < readLen)
			info.rawReadLength = readLen;
		info.rawTotalBaseNum += readLen;
		info.rawTotalReadNum++;
		int a=0, g=0, c=0, t=0, n=0;
		for(int i=0; i<readLen; i++)
		{
			switch(read->baseSequence[i])
			{
				case 'A':
					a++;
					info.base[i][A_]++;
					break;
				case 'a':
					a++;
					info.base[i][A_]++;
					break;
				case 'C':
					c++;
					info.base[i][C_]++;
					break;
				case 'c':
					c++;
					info.base[i][C_]++;
					break;
				case 'G':
					g++;
					info.base[i][G_]++;
					break;
				case 'g':
					g++;
					info.base[i][G_]++;
					break;
				case 'T':
					t++;
					info.base[i][T_]++;
					break;
				case 't':
					t++;
					info.base[i][T_]++;
					break;
				case 'N':
					n++;
					info.base[i][N_]++;
					break;
				case 'n':
					n++;
					info.base[i][N_]++;
					break;
				default:
					break;
			}
		}
		info.rawBaseA += a;
		info.rawBaseC += c;
		info.rawBaseG += g;
		info.rawBaseT += t;
		info.rawBaseN += n;
		return;
	}

	int MetaProcessor::getReadsNameFromFile(string filename, set<string> &readsName)
	{
		gzFile file = gzopen(filename.c_str(), "rb");
		if (file == NULL)
		{
			LOG(ERROR, "cannot open file: " + filename);
			return 1;
		}

		char buf[512];
		int terminate = 0;
		int len;
		gzgets(file, buf, 512);
		if(seqType_==0){
			while (gzgets(file, buf, 512) != NULL)
			{
				if(!isAlignLengthOK(buf))
				{
					continue;
				}
				terminate = 0;
				len = strlen(buf);
				while (buf[terminate++] != '/' && (terminate < len))
					;
				if (terminate == len)
				{
					LOG(ERROR, "[FORMAT ERROR] " << filename);
					return 1;
				}
				buf[terminate - 1] = '\0';
				readsName.insert(buf);
			}
		}else{
			while (gzgets(file, buf, 512) != NULL)
			{
				if(!isAlignLengthOK(buf))
				{
					continue;
				}
				terminate = 0;
				len = strlen(buf);
				while (buf[terminate++] != ' ' && (terminate < len))
					;
				if (terminate == len)
				{
					LOG(ERROR, "[FORMAT ERROR] " << filename);
					return 1;
				}
				buf[terminate - 1] = '\0';
				readsName.insert(buf);
			}
		}
		gzclose(file);
		return 0;
	}

	bool MetaProcessor::isAlignLengthOK(char *line)
	{
		int len = strlen(line);
		int tab = 0;
		char *adptLenPos = NULL;
		char *endPos = NULL;
		for(int i=0; i<len; i++)
		{
			if(line[i] =='\t')
			{
				tab++;
				if(tab==8)
				{
					adptLenPos = line + i + 1;
				}
				if(tab==9)
				{
					endPos = line + i;
					break;
				}
			}
		}
		if(adptLenPos ==NULL || endPos == NULL)
		{
			LOG(WARN, "can't find adptLen " << line);
			return false;
		}
		*endPos = '\0';
		int adptLen = atoi(adptLenPos);
		if(adptLen >= minAlignLength_)
		{
			return true;
		}
		return false;
	}

	bool MetaProcessor::hasAdapter(set<string> &readsName, const char *seqName)
	{
		int i = 0;
		if(seqType_==0){
			while (seqName[i++] != '/')
				;
		}else{
			while (seqName[i++] != ' ')
				;
		}

		return readsName.count(string(seqName + 1, i-2));
	}
	bool MetaProcessor::hasAdapter(const char *sequence, int readLen, const char *adapter, int adptLen)
	{
		bool find = false;
		int minMatchLen = minAlignLength_;
		int a1 = adptLen - minMatchLen;
		int r1 = 0;
		int len, mis;

		int right = readLen - minMatchLen;

		for (r1 = 0; r1 <= right;)
		{
			int len1 = adptLen - a1;
			int len2 = readLen - r1;
			len = (len1 < len2) ? len1 : len2;
			mis = 0;
			int map[MAX_LENGTH];
			map[0] = 0;
			for (int c = 0; c < len; ++c)
			{
				if (adapter[a1 + c] == sequence[r1 + c])
				{
					map[mis]++;
				}
				else
				{
					mis++;
					map[mis] = 0;
				}
			}
			int max_map = 0;
			for (int c = 0; c <= mis; ++c)
			{
				if (map[c] > max_map)
				{
					max_map = map[c];
				}
			}
			if ((mis < misMatchRate_*len) || (max_map >= matchNumber_))
			{
				find = true;
				break;
			}
			if (a1 > 0)
			{
				a1--;
			}
			else
			{
				r1++;
			}
		}

		return find;
	}
	void MetaProcessor::outputRawDataTask(gzFile &file,PreProcessTool::Read *reads, unsigned int size)
	{
		for (unsigned int i=0; i<size; ++i)
		{
			gzputs(file, reads[i].readName);
			gzputs(file, "\n");
			gzputs(file, reads[i].baseSequence);
			gzputs(file, "\n");
			gzputs(file, reads[i].optionalName);
			gzputs(file, "\n");
			gzputs(file, reads[i].baseQuality);
			gzputs(file, "\n");
		}
		return;
	}
	void MetaProcessor::outputCleanDataTaskNew(gzFile &file, PreProcessTool::Read *reads, int size, int *result, int type)
	{
		for(int i=0; i<size; i++)
		{
			if(result[i] == type)
			{
				gzputs(file, reads[i].readName);
				gzputs(file, "\n");
				gzputs(file, reads[i].baseSequence);
				gzputs(file, "\n");
				gzputs(file, reads[i].optionalName);
				gzputs(file, "\n");
				gzputs(file, reads[i].baseQuality);
				gzputs(file, "\n");
			}
		}
		return;
	}

	void MetaProcessor::outputSingleDataTask(gzFile &file, PreProcessTool::Read *reads1, PreProcessTool::Read *reads2, int size, int *result)
	{
		for(int i=0; i<size; i++)
		{
			if(result[i] == 1)
			{
				gzputs(file, reads1[i].readName);
				gzputs(file, "\n");
				gzputs(file, reads1[i].baseSequence);
				gzputs(file, "\n");
				gzputs(file, reads1[i].optionalName);
				gzputs(file, "\n");
				gzputs(file, reads1[i].baseQuality);
				gzputs(file, "\n");
			}
			else if(result[i] == 2)
			{
				gzputs(file, reads2[i].readName);
				gzputs(file, "\n");
				gzputs(file, reads2[i].baseSequence);
				gzputs(file, "\n");
				gzputs(file, reads2[i].optionalName);
				gzputs(file, "\n");
				gzputs(file, reads2[i].baseQuality);
				gzputs(file, "\n");
			}
		}
	}
	void MetaProcessor::outputCleanDataTask(gzFile &file, PreProcessTool::Read *reads)
	{
		unsigned int index = 0;
		while (index == atomic_read32(&size_) && (atomic_read32(&doneNum_) < PROCESS_THREAD_NUM))
		{
			this_thread::sleep(posix_time::seconds(2)); //sleep 2 second
		}

		while ((atomic_read32(&doneNum_) < PROCESS_THREAD_NUM) || index < atomic_read32(&size_))
		{
			this_thread::sleep(posix_time::seconds(1));
			while (index < atomic_read32(&size_))
			{
				gzputs(file, reads[cleanDataIndexs_[index]].readName);
				gzputs(file, "\n");
				gzputs(file, reads[cleanDataIndexs_[index]].baseSequence);
				gzputs(file, "\n");
				gzputs(file, reads[cleanDataIndexs_[index]].optionalName);
				gzputs(file, "\n");
				gzputs(file, reads[cleanDataIndexs_[index]].baseQuality);
				gzputs(file, "\n");

				index++;
			}
		}
		return;
	}
	int MetaProcessor::printMetaStats(FqInfo &info1, FqInfo &info2)
	{
		string mode = "fastx";
		string cleanStatFile = outDir_ + "/" + outCleanPfx_ + ".readtrimfilter." + mode + ".stats";
		string rawStatFile1 = outDir_ + "/" + outCleanPfx_ + ".pair.1.fq.gz.raw.reads.stats";
		string rawStatFile2 = outDir_ + "/" + outCleanPfx_ + ".pair.2.fq.gz.raw.reads.stats";
		ofstream cleanStat(cleanStatFile.c_str());
		ofstream rawStat1(rawStatFile1.c_str());
		ofstream rawStat2(rawStatFile2.c_str());
		if(!cleanStat)
		{
			LOG(ERROR,"can't create file: "<<cleanStatFile);
			return 1;
		}
		if(!rawStat1)
		{
			LOG(ERROR,"can't create file: "<<rawStatFile1);
			return 1;
		}
		if(!rawStat2)
		{
			LOG(ERROR,"can't create file: "<<rawStatFile2);
			return 1;
		}

		rawStat1 << "Reads\tBases\n";
		rawStat1 << info1.rawTotalReadNum << "\t" << info1.rawTotalBaseNum << "\n";
		rawStat2 << "Reads\tBases\n";
		rawStat2 << info2.rawTotalReadNum << "\t" << info2.rawTotalBaseNum << "\n";
		rawStat1.close();
		rawStat2.close();

		cleanStat << "Reads\tBases\tMax\tAvg\tKmer\tInserts\n";
		unsigned long reads = info1.cleanTotalReadNum + info2.cleanTotalReadNum + info1.singleReadNum + info2.singleReadNum;
		if(reads == 0)
		{
			cleanStat.close();
			return 2;
		}
		unsigned long bases = info1.cleanTotalBaseNum + info2.cleanTotalBaseNum + info1.singleBaseNum + info2.singleBaseNum;
		unsigned long inserts = info1.cleanTotalReadNum + info1.singleReadNum + info2.singleReadNum;
		int averageLen = int (1.0 * bases / reads + 0.5);
		int kmer = 0;
		if(averageLen%2==0)
		{
			kmer = averageLen / 2;
			if(kmer%2==0)
			{
				kmer += 1;
			}
			else
			{
				kmer += 2;
			}
		}
		else
		{
			kmer = (averageLen + 1) / 2;
			if(kmer%2==0)
			{
				kmer += 1;
			}
		}
		cleanStat << reads << "\t" << bases << "\t" << info1.cleanReadLengthTmp << "\t" << averageLen << "\t" << kmer << "\t" << inserts << "\n";
		//other
		//		reads = info1.cleanTotalReadNum + info2.cleanTotalReadNum;
		//		if(reads == 0)
		//		{
		//			cleanStat.close();
		//			return 2;
		//		}
		//		bases = info1.cleanTotalBaseNum + info2.cleanTotalBaseNum;
		//		inserts = info1.cleanTotalReadNum;
		//		averageLen = int (1.0 * bases / reads + 0.5);
		//		if(averageLen%2==0)
		//		{
		//			kmer = averageLen / 2;
		//			if(kmer%2==0)
		//			{
		//				kmer += 1;
		//			}
		//			else
		//			{
		//				kmer += 2;
		//			}
		//		}
		//		else
		//		{
		//			kmer = (averageLen + 1) / 2;
		//			if(kmer%2==0)
		//			{
		//				kmer += 1;
		//			}
		//		}
		//		cleanStat << reads << "\t" << bases << "\t" << info1.cleanReadLengthTmp << "\t" << averageLen << "\t" << kmer << "\t" << inserts << "\n";
		cleanStat.close();
		return 0;
	}
}//name space
