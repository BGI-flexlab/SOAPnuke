/*
 * DGEProcessor.cpp
 *
 *  Created on: 2012-6-3
 *      Author: Shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */
#include "DGEProcessor.h"


using namespace boost::threadpool;

namespace DGEProcessTool {

	void DGEProcessor::printVersion()
	{
		cout << "Preprocessing DGE version 1.3.0\n";
		cout << "Author:	chenyongsheng\n";
		cout << "Email :	chenyongsheng@genomics.cn\n";
	}

	void DGEProcessor::printUsage()
	{
		cout << "Useage: [OPTION]... FILE [FILE]\n";
		cout << "must arg:\n";
		cout << "\t-f, --fq         STR       fastq file\n";
		cout << "\t-A, --adapter    STR       adapter sequenc\n";
		cout << "\t--tile           STR       tile number to ignore reads , such as [1101-1104,1205]\n";

		cout << "\nusual args:\n";
		cout << "\t-F, --outfq                print out original fq file (default: off)\n";
		cout << "\t-o, --outDir     STR       out directory (default: current directory)\n";
		cout << "\t-x, --outPfx     STR       out file prefix or sample name (default: clean)\n";
		cout << "\t-s, --site       STR       Extra bases before clean tag (default: CATG)\n";
		cout << "\t-l, --tagLen     INT       Tag length (default: 17)\n";
		cout << "\t-r, --tagRange   STR       Tag length's range when find adapter (default: [17,18])\n";
		cout << "\t-M, --misMatch   INT       Max mismatch number when find adapter (dfault: 1)\n";
		cout << "\t-t, --trim       STR       trim some bp of the read's head and tail (default: [0,0])\n";
		cout << "\t-c, --cut        INT       copy number lower limit (default: 1)\n";
		cout << "\t-N, --number     FLOAT     reserve read number in each fq file (K reads(1024 reads), 0 means not cut reads)\n";

		cout << "\t-Q, --qualSys    INT       quality system, 1:illumina, 2:sanger (default: 1)\n";
		cout << "\t-i, --index                remove index of fq file\n";
		cout << "\t-G, --sanger               out put sanger quality score system fq. (defaul: off illumina)\n";

		cout << "\t-u, --unlowQ                do not filter low quality reads (N reads)(default: off)\n";
		cout << "\t-L, --readLen    INT       Max read length in fq file (default: 49)\n";

		cout << "\nhelp args:\n";
		cout << "\t-a, --append     STR       logger's appender: console or file (defualt: console)\n";
		cout << "\t-h, --help                 help\n";
		cout << "\t-v, --version              version information" << endl;
	}

	int DGEProcessor::processParams(int argc, char **argv)
	{
		const char *shortOptions = "f:FA:K:o:x:s:l:r:M:t:c:N:Q:iGuL:a:hv";
		const struct option longOptions[] =
		{
			{ "fq",       1, NULL, 'f' },
			{ "outfq",    0, NULL, 'F' },
			{ "adapter",  1, NULL, 'A' },
			{ "tile",     1, NULL, 'K' },
			{ "outDir",   1, NULL, 'o' },
			{ "outPfx",   1, NULL, 'x' },
			{ "site",     1, NULL, 's' },
			{ "tagLen",   1, NULL, 'l' },
			{ "tagRange", 1, NULL, 'r' },
			{ "misMatch", 1, NULL, 'M' },
			{ "trim",     1, NULL, 't' },
			{ "cut",      1, NULL, 'c' },
			{ "number",   1, NULL, 'N' },
			{ "qualSys",  1, NULL, 'Q' },
			{ "index",    0, NULL, 'i' },
			{ "sanger",   0, NULL, 'G' },
			{ "unlowQ",   0, NULL, 'u' },
			{ "readLen",  1, NULL, 'L' },
			{ "append",   1, NULL, 'a' },
			{ "help", 	  0, NULL, 'h' },
			{ "version",  0, NULL, 'v' },
		};

		string append;

		if (argc == 1)
		{
			cout << "Print -h or --help for more information." << endl;
			return 1;
		}

		int nextOpt;
		string trim;
		string tiles;
		string tagRange;
		while (-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions,
						NULL)))
		{
			switch (nextOpt)
			{
				case 'f':
					fqFile1_.assign(optarg);
					break;
				case 'F':
					outfq_ = true;
					break;
				case 'A':
					adapter_.assign(optarg);
					break;
				case 'K':
					filterTile_ = true;
					tiles.assign(optarg);
					PreProcessTool::getTiles(tiles, tiles_);
					break;
				case 'o':
					outDir_.assign(optarg);
					break;
				case 'x':
					outPfx_.assign(optarg);
					break;
				case 'u':
					filterLowQual_ = false;
					break;
				case 's':
					site_.assign(optarg);
					break;
				case 'l':
					tagLength_ = atoi(optarg);
					break;
				case 'M':
					misMatch_ = atoi(optarg);
					break;
				case 'Q':
					switch (optarg[0])
					{
						case '1':
							qualSystem_ = ILLUMINA_;
							break;
						case '2':
							qualSystem_ = SANGER_;
							break;
						default:
							cout << "error quality system" << endl;
							return 1;
					}
					break;
				case 'i':
					filterIndex_= true;
					break;
				case 'G':
					score_= 33;
					break;
				case 'r':
					tagRange.assign(optarg);
					if(tagRange.find(',') == string::npos)
					{
						cout << "-t/--tagRange options error" << endl;
						return 1;
					}
					tagStart_ = atoi(tagRange.substr(0, tagRange.find(',')).c_str());
					tagEnd_ = atoi(tagRange.substr(tagRange.find(',')+1).c_str());
					break;

				case 't':
					trim.assign(optarg);
					if (trim.find(',') == string::npos)
					{
						cout << "-t/--trim options error" << endl;
						return 1;
					}
					headTrim_ = atoi(trim.substr(0, trim.find(',')).c_str());
					tailTrim_ = atoi(trim.substr(trim.find(',') + 1).c_str());
					break;
				case 'L':
					readLen_ = atoi(optarg);
					break;
				case 'N':
					number_ = (unsigned long int)(atof(optarg) * 1024);
					break;
				case 'c':
					cutOff_ = atoi(optarg);
					break;
				case 'a':
					append = optarg;
					break;
				case 'h':
					printUsage();
					return 1;
					break;
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

		bool isPathNotExists = false; //指示输出目录是否存在
		if (access(outDir_.c_str(), F_OK) == -1)  //输出路径不存在
		{
			isPathNotExists = true;
			int len = outDir_.size();
			char *path = (char *)malloc(len + 15);
			sprintf(path, "mkdir -p %s", outDir_.c_str());
			if (system(path) == -1)
			{
				cerr << "output directory " << outDir_ << " cannot create" << endl;
				return 1;
			}
		}

		string logoutPath = outDir_ + "/" + LOGOUT_NAME;
		if (!init_logger(append, logoutPath))
		{
			cerr << "Cannot Init Log:" << append << "-" << logoutPath << endl;
			return 1;
		}
		else
		{
			LOG(INFO, "Log Init Success");
		}

		if (isPathNotExists)
		{
			LOG(WARN, "output directory " << outDir_ << " does not exists, program will auto create");
			LOG(WARN, "output directory " << outDir_ << " has been created");
		}

		if (fqFile1_.empty())
		{
			LOG(ERROR, "fq file must be exists");
			return 1;
		}

		if (adapter_.empty())
		{
			LOG(ERROR, "adapter options must exists");
			return 1;
		}

		return 0;
	}

	int DGEProcessor::processRNA(int argc, char **argv)
	{
		//process the command line params, and initial the member variables
		if (processParams(argc, argv) > 0)
		{
			return 1;
		}

		bool isSeq = isSequence(adapter_);
		if (!isSeq)
		{
			LOG(ERROR, "adapter must be sequence");
			return 1;
		}

		//to upper
		for (int i=0; adapter_[i]; ++i)
		{
			adapter_[i] = toupper(adapter_[i]);
		}

		if (tagStart_ < 0 || tagStart_ > readLen_ || tagEnd_ < 0 || tagEnd_ > readLen_ || tagStart_ > tagEnd_)
		{
			LOG(ERROR, "tagRange must for 0~" << readLen_ << " and ,first <= second number");
			return 1;
		}

		return processDGE();

		return 0;
	}
	int DGEProcessor::processDGE()
	{
		int chrLen = chr_.length();
		for (int i = 0; i < chrLen; i++)
		{
			for (int j = 0; j< chrLen; j++)
			{
				if (chr_[j] == 'N' || chr_[i] == 'N')
				{
					match_[i][j] = 1;
				}
				else if (chr_[i] == chr_[j])
				{
					match_[i][j] = 2;
				} 
				else
				{
					match_[i][j] = -4.5;
				}
			}
		}

		FqInfo globleInfo;
		if (readLen_ - headTrim_ - tailTrim_ < 0){
			LOG(ERROR, "trim["<<headTrim_<<","<<tailTrim_<<"] is large then readLen"<<readLen_);
			exit(1);
		}
		globleInfo.readslength = readLen_ - headTrim_ - tailTrim_;

		gzFile outFile = NULL; 
		string outFileName;
		outFileName = outDir_ + "/" + outPfx_ + ".fq.gz";

		if(outfq_){
			if(number_ == 0 && headTrim_ == 0 && tailTrim_ == 0 && !filterTile_ && !filterIndex_ && fqFile1_.substr(fqFile1_.size() - 2, 2) == "gz")
			{
				string sys = "cp -u " + fqFile1_ + " " + outFileName;
				system(sys.c_str());
				outfq_ = false;
			}
			else
			{
				outFile = gzopen(outFileName.c_str(), "wb");
				if (!outFile)
				{
					LOG(ERROR, "create output file: " + outFileName + " error");
					return 1;
				}
			}
		}

		long capacity = memLimit_;

		PreProcessTool::FqBuffer fqBuffer(fqFile1_.c_str(), capacity, PreProcessTool::FqBuffer::RB, filterTile_, tiles_);

		pool pl(threadNum_);

		TaskParam *params = new TaskParam[threadNum_];

		PreProcessTool::Read *reads;
		int size;
		unsigned long sizeTmp = 0;

		while ((reads = fqBuffer.getReads()))
		{
			size = fqBuffer.getRealReadSize();
			if (size == 0)
				continue;
			if (size == -1)
			{
				LOG(ERROR, "read fq file: " + fqFile1_ + " error");
				return 1;
			}

			sizeTmp += size;
			if (number_ > 0 && sizeTmp > number_)
			{
				size = size - (sizeTmp - number_);
			}
			int *isFilter = new int[size];
			bzero(isFilter, sizeof(int) * size);

			int block = size / threadNum_;
			int remain = size % threadNum_;
			int index = 0;

			for (int i=0; i<threadNum_; ++i)
			{
				params[i].left = index;
				index += block;
				if (remain > 0)
				{
					index += 1;
					remain--;
				}

				params[i].right = index;
				params[i].reads1 = reads;
				params[i].result = isFilter;
				clearFqInfo(params[i].info1);

				pl.schedule(boost::bind(&DGEProcessor::taskDGE, this, &params[i]));
			}

			pl.wait();

			for (int i=0; i<threadNum_; ++i)
			{
				globleInfo.add(params[i].info1, readLen_);
			}

			for (int i=0; i<size; ++i)
			{
				if(outfq_){

					gzputs(outFile, reads[i].readName);
					gzputs(outFile, "\n");
					gzputs(outFile, reads[i].baseSequence);
					gzputs(outFile, "\n");
					gzputs(outFile, reads[i].optionalName);
					gzputs(outFile, "\n");
					gzputs(outFile, reads[i].baseQuality);
					gzputs(outFile, "\n");
				}
				trimAdapter (reads[i].baseSequence, tagLength_);
				if(isFilter[i] == 0)
				{
					//				trimAdapter (reads[i].baseSequence, tagLength_);
					tagSequence_[string(reads[i].baseSequence)]++;
				}
				else if(isFilter[i] == 1)
				{
					adptSequence_[string(reads[i].baseSequence)]++;
				}
				else if(isFilter[i] == 2)
				{
					nSequence_[string(reads[i].baseSequence)]++;
				}

			}

			delete []isFilter;

			if(number_ > 0 && sizeTmp >= number_)
				break;
		}


		delete []params;
		if(outfq_){
			gzclose(outFile);
		}

		long jt = 0;
		for (map<string,int>::iterator it=tagSequence_.begin(); it != tagSequence_.end(); it++)
		{
			sortClean_[(*it).second].push_back(jt);
			cleanFa_.push_back(it);
			jt++;
		}
		for (int i=0; i<7 ;i++)
		{
			tagPercent[i][0] = 0;
			tagPercent[i][1] = 0;
		}

		outFileName = outDir_ + "/" + outPfx_ + "_TagCopyNumber.txt";
		ofstream faOut(outFileName.c_str());

		if (!faOut)
		{
			LOG(ERROR, "open output file "<<outFileName<<" error");
			exit(1);
		}
		//	faOut << "Tag\tCopy Num\n";
		faOut << "Tag\tCopyNumber\tTPM(Transcript Per Million)\n";
		faOut << fixed;
		for (map<unsigned long,vector<unsigned long> >::reverse_iterator it = sortClean_.rbegin(); it != sortClean_.rend(); it++)
		{
			for (unsigned long jt = 0; jt != (*it).second.size(); jt++)
			{	
				//			if ((*((cleanFa_[(*it).second[jt]]))).second > cutOff_)
				//			{
				//				faOut << site_ << (*((cleanFa_[(*it).second[jt]]))).first << "\t" << (*(cleanFa_[(*it).second[jt]])).second <<"\n";
				//			}

				if ((*((cleanFa_[(*it).second[jt]]))).second > 100)
				{
					tagPercent[6][0]++;
					tagPercent[6][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}

				if ((*((cleanFa_[(*it).second[jt]]))).second > 50)
				{
					tagPercent[5][0]++;
					tagPercent[5][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}

				if ((*((cleanFa_[(*it).second[jt]]))).second > 20)
				{
					tagPercent[4][0]++;
					tagPercent[4][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}

				if ((*((cleanFa_[(*it).second[jt]]))).second > 10)
				{
					tagPercent[3][0]++;
					tagPercent[3][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}

				if ((*((cleanFa_[(*it).second[jt]]))).second > 5)
				{
					tagPercent[2][0]++;
					tagPercent[2][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}

				if ((*((cleanFa_[(*it).second[jt]]))).second > cutOff_)
				{
					tagPercent[1][0]++;
					tagPercent[1][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}

				if ((*((cleanFa_[(*it).second[jt]]))).second <= cutOff_)
				{
					tagPercent[0][0]++;
					tagPercent[0][1] += (*((cleanFa_[(*it).second[jt]]))).second;
				}
			}
		}

		for (map<unsigned long,vector<unsigned long> >::iterator it = sortClean_.begin(); it != sortClean_.end(); it++)
		{
			for (unsigned long jt = 0; jt != (*it).second.size(); jt++)
			{
				if ((*((cleanFa_[(*it).second[jt]]))).second > cutOff_)
				{
					faOut << site_ << (*((cleanFa_[(*it).second[jt]]))).first << "\t" 
						<< (*(cleanFa_[(*it).second[jt]])).second << "\t" << setprecision(2) 
						<< 1.0 * standard_ * (*(cleanFa_[(*it).second[jt]])).second / tagPercent[1][1]
						<< "\n";
				}
			}
		}

		faOut.close();

		outFileName = outDir_ + "/Solexa_Tag_Library.xls";

		long nUniq = nSequence_.size(), adptUniq = adptSequence_.size();
		long totalUniq = nUniq + adptUniq + tagPercent[0][0] + tagPercent[1][0];
		faOut.open(outFileName.c_str());
		faOut << fixed;
		faOut << "\tDistinct Tag\t\tTotal Tag Number\t\n";

		faOut << "Raw Data         \t" << totalUniq << "\t100.00%\t" << globleInfo.rawTotalReads << "\t100.00%\n";

		faOut << "Tags Containing N\t" << nUniq << "\t" << setprecision(2) << 100.0 * nUniq / totalUniq << "%\t"
			<< globleInfo.readsWithLowQual << "\t" << setprecision(2) << 100.0 * globleInfo.readsWithLowQual / globleInfo.rawTotalReads << "%\n";

		faOut << "Adaptors         \t" << adptUniq << "\t" << setprecision(2) << 100.0 * adptUniq / totalUniq << "%\t"
			<< globleInfo.readsWithAdapter << "\t" << setprecision(2) << 100.0 * globleInfo.readsWithAdapter / globleInfo.rawTotalReads << "%\n";

		faOut << "Tag CopyNum <" << cutOff_+1 << "   \t" << tagPercent[0][0] << "\t" << setprecision(2) << 100.0 * tagPercent[0][0] / totalUniq << "%\t"
			<< tagPercent[0][1] << "\t" << setprecision(2) << 100.0 * tagPercent[0][1] / globleInfo.rawTotalReads << "%\n";

		faOut << "Clean Tag        \t" << tagPercent[1][0] << "\t" << setprecision(2) << 100.0 * tagPercent[1][0] / totalUniq << "%\t"
			<< tagPercent[1][1] << "\t" << setprecision(2) << 100.0 * tagPercent[1][1] / globleInfo.rawTotalReads << "%\n";

		faOut << "Tag CopyNum >=" << cutOff_+1 << "  \t" << tagPercent[1][0] << "\t" << setprecision(2) << 100.0 * tagPercent[1][0] / tagPercent[1][0] << "%\t"
			<< tagPercent[1][1] << "\t" << setprecision(2) << 100.0 * tagPercent[1][1] / tagPercent[1][1] << "%\n";

		faOut << "Tag CopyNum >5   \t" << tagPercent[2][0] << "\t" << setprecision(2) << 100.0 * tagPercent[2][0] / tagPercent[1][0] << "%\t"
			<< tagPercent[2][1] << "\t" << setprecision(2) << 100.0 * tagPercent[2][1] / tagPercent[1][1] << "%\n";

		faOut << "Tag CopyNum >10  \t" << tagPercent[3][0] << "\t" << setprecision(2) << 100.0 * tagPercent[3][0] / tagPercent[1][0] << "%\t"
			<< tagPercent[3][1] << "\t" << setprecision(2) << 100.0 * tagPercent[3][1] / tagPercent[1][1] << "%\n";

		faOut << "Tag CopyNum >20  \t" << tagPercent[4][0] << "\t" << setprecision(2) << 100.0 * tagPercent[4][0] / tagPercent[1][0] << "%\t"
			<< tagPercent[4][1] << "\t" << setprecision(2) << 100.0 * tagPercent[4][1] / tagPercent[1][1] << "%\n";

		faOut << "Tag CopyNum >50  \t" << tagPercent[5][0] << "\t" << setprecision(2) << 100.0 * tagPercent[5][0] / tagPercent[1][0] << "%\t"
			<< tagPercent[5][1] << "\t" << setprecision(2) << 100.0 * tagPercent[5][1] / tagPercent[1][1] << "%\n";

		faOut << "Tag CopyNum >100 \t" << tagPercent[6][0] << "\t" << setprecision(2) << 100.0 * tagPercent[6][0] / tagPercent[1][0] << "%\t"
			<< tagPercent[6][1] << "\t" << setprecision(2) << 100.0 * tagPercent[6][1] / tagPercent[1][1] << "%\n";
		faOut.close();

		output(&globleInfo, NULL);

		LOG(INFO, "PreProcess Finish");

		return 0;
	}

	void DGEProcessor::taskDGE(TaskParam *param)
	{
		PreProcessTool::Read *reads1 = param->reads1;

		int start = param->left;
		int end = param->right;

		if (start == end) //ask#*****************
		{
			return;
		}

		for (int i=start; i < end; ++i)
		{
			param->result[i] = statisticsDGE(reads1[i], param->info1);
		}
	}
	int DGEProcessor::statisticsDGE(PreProcessTool::Read &read, FqInfo &info)
	{
		StatisInfo si = auxStatistics(read, info);

		bool low = (filterLowQual_ && si.ns >= nSeed_ );
		if (low) //filter low quality
		{
			info.readsWithLowQual++;
			return 2;
		}

		if(findAdapter (read.baseSequence, adapter_.c_str()))
		{
			info.cleanTotalReads++;
			return 0;
		}
		else
		{
			info.readsWithAdapter++;
			return 1;
		}

		return 1;
	}






	void DGEProcessor::output(FqInfo *info1, FqInfo *info2)
	{
		if (info2 == NULL)
		{
			printFqInfo(outDir_.c_str(), *info1);
		}
		else
		{
			printFqInfo(outDir_.c_str(), *info1);
			printFqInfo(outDir_.c_str(), *info2);
		}

	}



	StatisInfo DGEProcessor::auxStatistics(PreProcessTool::Read &read, FqInfo &info)
	{
		if (headTrim_ > 0)
		{
			read.baseSequence = read.baseSequence + headTrim_;
			read.baseQuality = read.baseQuality + headTrim_;
		}
		int readLen = strlen(read.baseSequence);
		if (tailTrim_ > 0){
			int right = readLen - tailTrim_;
			read.baseSequence[right] = '\0';
			read.baseQuality[right] = '\0';
			readLen = strlen(read.baseSequence);
		}
		int qual;

		info.rawTotalReads++;
		info.rawTotalBases += readLen;

		StatisInfo si;
		for (int i = 0; i < readLen; ++i)
		{
			switch (read.baseSequence[i])
			{
				case 'A':
					si.a++;
					info.base[i][A_]++;
					break;
				case 'C':
					si.c++;
					info.base[i][C_]++;
					break;
				case 'G':
					si.g++;
					info.base[i][G_]++;
					break;
				case 'T':
					si.t++;
					info.base[i][T_]++;
					break;
				case 'N':
					si.n++;
					info.base[i][N_]++;
					if (i < tagLength_)
						si.ns++;
					break;
			}

			qual = read.baseQuality[i] - qualSystem_;
			info.qual[i][qual]++;
			read.baseQuality[i] = qual + score_;

			if (qual >= 20)
			{
				si.q20++;
				info.q20q30[i][0]++;
				if (qual >= 30)
				{
					si.q30++;
					info.q20q30[i][1]++;
				}
			}
		}

		info.rawBaseA += si.a;
		info.rawBaseC += si.c;
		info.rawBaseG += si.g;
		info.rawBaseT += si.t;
		info.rawBaseN += si.n;
		info.rawQ20 += si.q20;
		info.rawQ30 += si.q30;
		if (filterIndex_)
		{
			int sharpIndex = 0;
			int i = 0;
			while (read.readName[i++] != '#')
				;
			sharpIndex = i;
			while (read.readName[i++] != '/')
				;
			strcpy(read.readName + sharpIndex, read.readName + i - 1);
		}
		return si;
	}

	bool DGEProcessor::findAdapter(const char *sequence, const char *adapter)
	{
		bool find = false;
		int readLen = strlen(sequence);
		int adptLen = strlen(adapter);
		int r1;
		int mis;
		int len;
		for (r1 = tagStart_; r1 <= tagEnd_; r1++)
		{
			int len2 = readLen - r1;
			len = (adptLen < len2) ? adptLen : len2;
			mis = 0;
			for (int c = 0; c < len; c++)
			{
				if (adapter[c] != sequence[r1 + c])
				{
					mis++;
				}
				if (mis > misMatch_)
					break;
			}
			if (mis <= misMatch_)
				find = true;
		}

		if(!find)
		{
			string query = adapter;
			string target = sequence;
			int offset = smithWatermanAign(query, target) - 1;
			for (r1 = tagStart_; r1 <= tagEnd_; r1++)
			{
				if (offset == r1){
					find = true;
				}
			}
		}
		return find;
	}
	void DGEProcessor::trimAdapter(char *sequence, int startPos)
	{
		sequence[startPos]='\0';
	}

	int DGEProcessor::smithWatermanAign(string query, string target) { // Smith-Waterman algorithm
		// the match in query must be started on base 0
		int queryLen = query.length();
		int targetLen = target.length();
		int i, j;
		int open = 100;

		float **scoreMatrix, **directMatrix; // 0: up, 1: left, 2: northwest, 3: itself
		int maxcol, maxrow;
		float maxScore;
		scoreMatrix = (float **) malloc(sizeof(float *) * (targetLen + 1));
		for (i = 0; i <= targetLen; i++) {
			scoreMatrix[i] = (float *) malloc(sizeof(float) * (targetLen + 1));
		}
		directMatrix = (float **) malloc(sizeof(float *) * (targetLen + 1));
		for (i = 0; i <= targetLen; i++) {
			directMatrix[i] = (float *) malloc(sizeof(float) * (targetLen + 1));
		}

		for (i = 0; i <= targetLen; i++) {
			scoreMatrix[0][i] = 0;
			directMatrix[0][i] = 1;
		}
		for (i = 0; i <= queryLen; i++) {
			scoreMatrix[i][0] = i * open * (-1);
			directMatrix[i][0] = 0;
		}
		directMatrix[0][0] = 3;
		for (i = 1; i <= queryLen; i++) {
			for (j = 1; j <= targetLen; j++) {
				float a = scoreMatrix[i - 1][j - 1] + match_[chr_.find(toupper(query[i - 1]))][chr_.find(toupper(target[j - 1]))];
				float b = scoreMatrix[i - 1][j] - open;
				float c = scoreMatrix[i][j - 1] - open;
				float temp1;
				if (b >= c) {
					directMatrix[i][j] = 0;
					temp1 = b;
				} else {
					directMatrix[i][j] = 1;
					temp1 = c;
				}
				if (a >= temp1) {
					directMatrix[i][j] = 2;
					scoreMatrix[i][j] = a;
				} else {
					scoreMatrix[i][j] = temp1;
				}
				if (maxScore < scoreMatrix[i][j]) {
					maxScore = scoreMatrix[i][j];
					maxcol = j;
					maxrow = i;
				}
			}
		}
		// trace back start point
		int row = queryLen;
		int col = targetLen;

		// output
		float colmaxS = scoreMatrix[0][targetLen];
		row = 0;
		for (j = 1; j <= queryLen; j++) {
			if (scoreMatrix[j][targetLen] >= colmaxS) {
				row = j;
				colmaxS = scoreMatrix[j][targetLen];
			}
		}

		int queryS = -1, queryE, targetS = -1, targetE;
		targetE = col;
		queryE = row;
		while (directMatrix[row][col] != 3) {
			if (directMatrix[row][col] == 0) {
				row--;
			} else if (directMatrix[row][col] == 1) {
				col--;
			} else if (directMatrix[row][col] == 2) {
				targetS = col;
				queryS = row;
				col--;
				row--;
			}
		}

		// free memory
		for (i = 0; i <= targetLen; i++) {
			free(scoreMatrix[i]);
		}
		for (i = 0; i <= targetLen; i++) {
			free(directMatrix[i]);
		}
		free(scoreMatrix);
		free(directMatrix);

		//cout << queryS << "," << queryE << "," << targetS << "," << targetE;
		return targetS;
	}
}  // namespace DGEProcessTool
