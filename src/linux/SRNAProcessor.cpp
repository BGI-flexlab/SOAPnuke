/*
 * RNAProcessor.cpp
 *
 *  Created on: 2012-6-3
 *      Author: Shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */
#include "SRNAProcessor.h"
#include <sstream>
using namespace boost::threadpool;

namespace SRNAProcessTool {

	void RNAProcessor::printVersion()
	{
		cout << "Preprocessing RNA version 1.5.0\n";
		cout << "Author:	chenyongsheng\n";
		cout << "Email :	chenyongsheng@genomics.cn\n";
	}

	void RNAProcessor::printUsage()
	{
		cout << "Useage: [OPTION]... FILE [FILE]\n";
		cout << "must arg:\n";
		cout << "\t-f, --fq         <string> :  fastq file\n";

		cout << "\nusual args:\n";
		cout << "\t-m, --mrna       <switch> :  mrna filter(default: off)\n";
		cout << "\t\t-n, --polyN      <float>  :  remove polyN[A, T, G, C], 0 means do not filter, (default: 0.7)\n";
		cout << "\t-F, --outfq      <string> :  prefix of out orignal fq name, Eg. if set -F XXX, will print out XXX.fq.gz, otherwise will not print\n";
		cout << "\t-3, --adapter3   <string> :  3' adaptor sequence (default: TCGTATGCCGTCTTCTGCTTG)\n";
		cout << "\t-5, --adapter5   <string> :  5' adaptor sequence (default: GTTCAGAGTTCTACAGTCCGACGATC)\n";
		cout << "\t--tile           <string> :  tile number to ignore reads , such as [1101-1104,1205]\n";
		cout << "\t-o, --outDir     <string> :  out directory (default: current directory)\n";
		cout << "\t-x, --outPfx     <string> :  out file prefix (default: clean)\n";
		cout << "\t-s, --strict     <switch> :  filter low quality reads strictly (default: off)\n";
		cout << "\t-z, --minSize    <int>    :  small insert size (default: 18)\n";
		cout << "\t-p, --polyA      <float>  :  filter poly A, percent of A, 0 means do not filter, (default: 0.7)\n";
		cout << "\t-Q, --qualSys    <int>    :  quality system, 1:illumina, 2:sanger (default: 1)\n";
		cout << "\t-q, --fastq      <switch> :  out file type: on:fastq, off:fasta (default: off)\n";
		cout << "\t\t-i, --index      <switch> :  remove index\n";
		cout << "\t\t-G, --sanger     <switch> :  out put sanger quality score system fq. (defaul: off illumina)\n";
		cout << "\t-u, --untrim     <switch> :  do not trim 3' adapter (default: off)\n";
		cout << "\t-w, --unlowQ     <switch> :  do not filter low quality reads (default: off)\n";
		cout << "\t-L, --readLen    <int>    :  Max read length in fq file (default: 49)\n";
		cout << "\t-t, --trim       <string> :  trim some bp of the read's head and tail (default: [0,0])\n";
		cout << "\t-c, --cut        <float>  :  the read number you want to keep in each orignal fq file.\n"
			<< "\t                             Eg.: if set -c N, read number = N * 1024; default: N = 0, means reserve whole orignal fq;\n";
		cout << "\t-y, --seqType   : <i> Sequence fq name type, 0->old fastq name, 1->new fastq name HighSeq4000[default: 0]\n";
		cout << "\nhelp args:\n";
		cout << "\t-a, --append     <string> :  logger's appender: console or file (defualt: console)\n";
		cout << "\t-h, --help       <switch> :  help\n";
		cout << "\t-v, --version    <switch> :  version information\n";

		cout << "\nunusual arg:\n";
		cout << " find 5' adapter\n";
		cout << "\t-C, --continuous <int>    :  mini 5' adapter continuous alignment length (default: 6)\n";
		cout << "\t-A, --alignRate  <float>  :  mini alignment rate when find 5' adapter: alignment/tag (default: 0.8)\n";
		cout << " find 3' adapter\n";
		cout << "\t-l, --miniAlign  <int>    :  mini alignment length when find 3' adapter (default: 5)\n";
		cout << "\t-E, --errorRate  <float>  :  Max error rate when find 3' adapter (mismatch/match) (dfault: 0.4)\n";
		cout << "\t-M, --misMatch   <int>    :  Max mismatch number when find 3' adapter (dfault: 4)" << endl;
	}

	int RNAProcessor::processParams(int argc, char **argv)
	{
		const char *shortOptions = "f:mn:F:3:5:K:o:x:sz:p:Q:qy:GiuwL:t:c:a:hvC:A:l:E:M:";
		const struct option longOptions[] =
		{
			{ "fq",       1, NULL, 'f' },
			{ "mrna",     0, NULL, 'm' },
			{ "polyN",    1, NULL, 'n' },
			{ "outfq",    1, NULL, 'F' },
			{ "adapter3", 1, NULL, '3' },
			{ "adapter5", 1, NULL, '5' },
			{ "tile",     1, NULL, 'K' },
			{ "untrim",   0, NULL, 'u' },
			{ "unlowQ",   0, NULL, 'w' },
			{ "strict",   0, NULL, 's' },
			{ "polyA",    1, NULL, 'p' },
			{ "minSize",  1, NULL, 'z' },
			{ "outDir",   1, NULL, 'o' },
			{ "outPfx",   1, NULL, 'x' },
			{ "qualSys",  1, NULL, 'Q' },
			{ "fastq",    0, NULL, 'q' },
			{ "index",    0, NULL, 'i' },
			{ "seqType",  1, NULL, 'y' },
			{ "sanger",   0, NULL, 'G' },
			{ "readLen",  1, NULL, 'L' },
			{ "trim",     1, NULL, 't' },
			{ "cut",      1, NULL, 'c' },
			{ "continuous", 1, NULL, 'C' },
			{ "alignRate", 1, NULL, 'A' },
			{ "miniAlign", 1, NULL, 'l' },
			{ "errorRate", 1, NULL, 'E' },
			{ "misMatch",  1, NULL, 'M' },
			{ "append",    1, NULL, 'a' },
			{ "help", 	   0, NULL, 'h' },
			{ "version",   0, NULL, 'v' },
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
		while (-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions,
						NULL)))
		{
			switch (nextOpt)
			{
				case 'f':
					fqFile1_.assign(optarg);
					break;
				case 'm':
					mrna_ = true;
					break;
				case 'n':
					filterPolyN_ = atof(optarg);
					break;
				case 'F':
					outfq_ = true;
					outFqName_.assign(optarg);
					break;
				case '3':
					adapter3_.assign(optarg);
					break;
				case '5':
					adapter5_.assign(optarg);
					break;
				case 'K':
					filterTile_ = true;
					tiles.assign(optarg);
					PreProcessTool::getTiles(tiles, tiles_);
					break;
				case 'u':
					trim_ = false;
					break;
				case 'w':
					filterLowQual_ = false;
					break;
				case 's':
					strict_ = true;
					lowQualSeed1Num_ = 2;
					lowQualSeed2Num_ = 3;
					break;
				case 'p':
					filterPolyA_ = atof(optarg);
					break;
				case 'z':
					minInsertSize_ = atoi(optarg);
					break;
				case 'y':
					seqType_ = atoi(optarg);
					break;
				case 'o':
					outDir_.assign(optarg);
					break;
				case 'x':
					outPfx_.assign(optarg);
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
				case 'q':
					fastq_ = true;
					break;
				case 'i':
					filterIndex_= true;
					break;
				case 'G':
					score_= 33;
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
				case 'c':
					cutOff_ =(unsigned long int)(atof(optarg) * 1024);
					break;
				case 'C':
					continuousAlign_  = atoi(optarg);
					break;
				case 'A':
					alignRate_ = atof(optarg);
					break;
				case 'l':
					miniAlignLength_ = atoi(optarg);
					break;
				case 'E':
					errorRate_ = atof(optarg);
					break;
				case 'M':
					misMatch_ = atoi(optarg);
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
			LOG(ERROR, "fq1 file must be exists");
			return 1;
		}

		if (adapter3_.empty() || adapter5_.empty())
		{
			LOG(ERROR, "adapter3_ and adapter5_ options must exists");
			return 1;
		}

		return 0;
	}

	int RNAProcessor::processRNA(int argc, char **argv)
	{
		//process the command line params, and initial the member variables
		if (processParams(argc, argv) > 0)
		{
			return 1;
		}

		int adptType = adapterType(adapter3_, adapter5_);
		if (adptType != 1)
		{
			LOG(ERROR, "adapter3_ and adapter5_ must in the same type, (sequence)");
			return 1;
		}

		//to upper
		for (int i=0; adapter3_[i]; ++i)
		{
			adapter3_[i] = toupper(adapter3_[i]);
		}

		for (int i=0; adapter5_[i]; ++i)
		{
			adapter5_[i] = toupper(adapter5_[i]);
		}


		if (fastq_)
		{
			return processFQ();
		}
		else
		{
			return processFA();
		}

		return 0;
	}



	void RNAProcessor::output(FqInfo *info1, FqInfo *info2)
	{
		/*	ofstream out("Statistic_Info");
			if (!out)
			{
			LOG(ERROR, "open output file: Statistic_Info error");
			exit(1);
			}
		 */
		if (info2 == NULL)
		{
			//		out << "[SE]" << "\n";
			//		out << "[FQ1]" << "\n";
			//		info1.lengthStart = minInsertSize_;
			//		info1.lengthEnd = 44;
			printFqInfo(outDir_.c_str(), *info1);
		}
		else
		{
			//		out << "[PE]" << "\n";
			//		out << "[FQ1]" << "\n";
			printFqInfo(outDir_.c_str(), *info1);
			//		out << "[FQ2]" << "\n";
			printFqInfo(outDir_.c_str(), *info2);
		}

		//	out.close();
	}



	StatisInfo RNAProcessor::auxStatistics(PreProcessTool::Read &read, FqInfo &info)
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
					if (i < seedLength_)
						si.ns++;
					break;
			}

			qual = read.baseQuality[i] - qualSystem_;
			info.qual[i][qual]++;

			if (i < seedLength_)
			{
				if (qual < lowQualSeed2_)
				{
					si.lowQual2++;
					if (qual < lowQualSeed1_)
						si.lowQual1++;
				}
			}

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
			if(seqType_ == 0){
				int sharpIndex = 0;
				int i = 0;
				while (read.readName[i++] != '#')
					;
				sharpIndex = i;
				while (read.readName[i++] != '/')
					;
				strcpy(read.readName + sharpIndex, read.readName + i - 1);

			}else{
				int i = strlen(read.readName);
				while (read.readName[--i] !=':')
					;
				if(i>0)
					read.readName[i]='\0';
			}
		}


		return si;
	}

	int RNAProcessor::findAdapter(const char *sequence, const char *adapter)
	{
		int startPos = -1;
		int readLen = strlen(sequence);
		int adptLen = strlen(adapter);
		int a1 = 2;
		int r1 = 0;
		int len;
		int mis;
		bool flagType = false;
		int misTmp = 0;
		int totalMapTmp = 0;
		for (r1 = 0; r1 <= readLen - miniAlignLength_;)
		{
			int len1 = adptLen - a1;
			int len2 = readLen - r1;
			len = (len1 < len2) ? len1 : len2;
			mis = 0;
			int map[MAX_LENGTH];
			map[0] = 0;
			int totalMap = 0;
			for (int c = 0; c < len; c++)
			{
				if (sequence[r1 + c] =='N'){
					continue;
				}
				if (adapter[a1 + c] == sequence[r1 + c])
				{
					map[mis]++;
					totalMap++;
				}
				else
				{
					mis++;
					map[mis] = 0;
				}
			}
			int misAndMap = mis + totalMap;
			float rate =1.0 * mis / totalMap;
			if (mis <= misMatch_ && misAndMap >= miniAlignLength_ && rate <= errorRate_)
			{
				if(flagType)
				{
					if (mis <= misTmp && totalMap >= totalMapTmp)
					{
						startPos = r1;
						misTmp = mis;
						totalMapTmp = totalMap;
					}
				}
				else
				{
					startPos = r1;
					flagType = true;
					misTmp = mis;
					totalMapTmp = totalMap;
				}
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
		return startPos;
	}
	void RNAProcessor::trimAdapter(char *sequence, int startPos)
	{
		sequence[startPos]='\0';
	}
	bool RNAProcessor::hasAdapter(const char *sequence, const char *adapter)
	{
		bool find = false;
		int readLen = strlen(sequence);
		int adptLen = strlen(adapter);
		int a1 = adptLen - continuousAlign_;
		int r1 = 0;
		int len;
		int mis;
		int readLenSmall = (readLen - continuousAlign_ < 0) ? 0 : (readLen - continuousAlign_);
		for (r1 = 0; r1 <= readLenSmall;)
		{
			int len1 = adptLen - a1;
			int len2 = readLen - r1;
			len = (len1 < len2) ? len1 : len2;
			mis = 0;
			int map[MAX_LENGTH];
			map[0] = 0;
			int totalMap=0;
			for (int c = 0; c < len; c++)
			{
				if (adapter[a1 + c] == sequence[r1 + c])
				{
					map[mis]++;
					totalMap++;
				}
				else
				{
					mis++;
					map[mis] = 0;
				}
			}
			int max_map = 0;
			for (int c = 0; c <= mis; c++)
			{
				if (map[c] > max_map)
				{
					max_map = map[c];
				}
			}
			if (mis <= 4 && (max_map >= continuousAlign_ || readLen < 12) && (1.0 * totalMap / readLen >= alignRate_ || 1.0 * totalMap / adptLen >= alignRate_))
			{
				find = true;
				break;
			}
			/*		if (max_map >= 15){
					find = true;
					break;
					}*/
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

	string RNAProcessor::mergeAndSortFilesByAscii(string &outpfx, vector<string> &files)
	{
		if(files.size()<1)
			return "";
		int tmpNumber = 0;
		string outpfxnew=outpfx+".AsciiTmp";
		int *result = new int[10];
		while(files.size()>1)
		{
			int child=0;
			int fileNumber=files.size();
			int plNumber = fileNumber / 2;
			plNumber = plNumber > 6 ? 6 : plNumber;
			//			pool pl(plNumber);
			bzero(result, sizeof(int) * 10);
			string rm = "rm";
			while(fileNumber>=2 && child<plNumber)
			{
				string inFile1 = files[0];
				files.erase(files.begin());
				string inFile2 = files[0];
				files.erase(files.begin());
				rm+=" " + inFile1 + " " + inFile2;
				string outFile = outpfxnew + intToString(tmpNumber);
				tmpNumber++;
				files.push_back(outFile);
				pl_.schedule(boost::bind(&RNAProcessor::mergeAndSortFilesByAsciiTask, this, inFile1, inFile2, outFile, &result[child]));
				child++;
			}
			pl_.wait();
			system(rm.c_str());
			for(int i=0; i<plNumber; i++)
			{
				if(result[i]!=0)
					return "";
			}
			//			delete []result;
		}
		delete []result;
		return files[0];
	}

	string RNAProcessor::mergeAndSortFilesByNumber(string &outpfx, vector<string> &files)
	{
		if(files.size()<1)
			return "";
		int tmpNumber = 0;
		string outpfxnew=outpfx+".NumberTmp";
		int *result = new int[10];
		while(files.size()>1)
		{
			int child=0;
			int fileNumber=files.size();
			int plNumber = fileNumber / 2;
			plNumber = plNumber > 6 ? 6 : plNumber;
			//			pool pl(plNumber);
			//			int *result = new int[plNumber];
			bzero(result, sizeof(int) * 10);
			string rm = "rm";
			while(fileNumber>=2 && child<plNumber)
			{
				string inFile1 = files[0];
				files.erase(files.begin());
				string inFile2 = files[0];
				files.erase(files.begin());
				rm += " " + inFile1 + " " + inFile2;
				string outFile = outpfxnew + intToString(tmpNumber);
				tmpNumber++;
				files.push_back(outFile);
				pl_.schedule(boost::bind(&RNAProcessor::mergeAndSortFilesByNumberTask, this, inFile1, inFile2, outFile, &result[child]));
				child++;
			}
			pl_.wait();
			system(rm.c_str());
			for(int i=0; i<plNumber; i++)
			{
				if(result[i]!=0)
					return "";
			}
			//			delete []result;
		}
		delete []result;
		return files[0];
	}
	long RNAProcessor::splitFile(string &outpfx, string &inFile, vector<string> &outFiles)
	{
		ifstream in(inFile.c_str(), fstream::in);
		if(in==NULL)
			return -1;
		int lineNumber = 5000000;
		//		int lineNumber = 10000;
		string outpfxnew = outpfx + ".SplitTmp";
		multimap<int, string> tags;
		long uniqNumber=0;
		int tmpNumber=0;
		char *line = new char[256];
		char *pos1=NULL;
		int pos2=0;
		while(!in.eof())
		{
			in.getline(line, 256);
			if(strlen(line)<1)
				continue;
			char *tab = strchr(line, '\t');
			*tab = '\0';
			pos1 = line;
			pos2 = atoi(tab+1);
			tags.insert(make_pair(pos2, string(pos1)));
			uniqNumber++;
			if(uniqNumber % lineNumber==0)
			{
				string outFile = outpfxnew + intToString(tmpNumber);
				ofstream out(outFile.c_str(), fstream::out);
				for(multimap<int, string>::reverse_iterator it=tags.rbegin(); it!=tags.rend();it++)
				{
					out<<it->second<<"\t"<<it->first<<"\n";
				}
				out.flush();
				out.close();
				tags.clear();
				outFiles.push_back(outFile);
				tmpNumber++;
			}
		}
		if(uniqNumber % lineNumber > 0)
		{

			string outFile = outpfxnew + intToString(tmpNumber);
			ofstream out(outFile.c_str(), fstream::out);
			for(multimap<int, string>::reverse_iterator it=tags.rbegin(); it!=tags.rend();it++)
			{
				out<<it->second<<"\t"<<it->first<<"\n";
			}
			out.flush();
			out.close();
			tags.clear();
			outFiles.push_back(outFile);
		}
		string rm = "rm " + inFile;
		system(rm.c_str());
		return uniqNumber;
	}
	void RNAProcessor::mergeAndSortFilesByNumberTask(string inFile1, string inFile2, string outFile, int *result)
	{
		*result = 1;
		ifstream in1(inFile1.c_str(), fstream::in);
		ifstream in2(inFile2.c_str(), fstream::in);
		ofstream out(outFile.c_str(), ofstream::out);
		if(in1==NULL || in2==NULL || out==NULL)
		{
			*result = 1;
			return;
		}
		char *line1 = new char[256];
		char *line2 = new char[256];
		char *pos11=NULL;
		char *pos21=NULL;
		int pos12=0;
		int pos22=0;
		int tmp=0;
		while(true)
		{
			if(tmp==0)
			{
				in1.getline(line1, 256);
				in2.getline(line2, 256);

				if(strlen(line1)<1 || strlen(line2)<1)
					break;

				char *tab1=strchr(line1, '\t');
				char *tab2=strchr(line2, '\t');
				*tab1='\0';
				*tab2='\0';

				pos11=line1;
				pos21=line2;
				pos12=atoi(tab1+1);
				pos22=atoi(tab2+1);
			}
			else if(tmp<0)
			{
				in1.getline(line1, 256);

				if(strlen(line1)<1)
					break;

				char *tab1=strchr(line1, '\t');
				*tab1='\0';

				pos11=line1;
				pos12=atoi(tab1+1);
			}
			else
			{
				in2.getline(line2, 256);

				if(strlen(line2)<1)
					break;

				char *tab2=strchr(line2, '\t');
				*tab2='\0';

				pos21=line2;
				pos22=atoi(tab2+1);
			}

			tmp = pos22 - pos12;
			if(tmp>0)
			{
				out<<pos21<<"\t"<<pos22<<"\n";
			}
			else if(tmp<0)
			{
				out<<pos11<<"\t"<<pos12<<"\n";
			}
			else
			{
				tmp = strcmp(pos21, pos11);
				if(tmp<0)
				{
					out<<pos11<<"\t"<<pos12<<"\n";
				}
				else if(tmp>0)
				{
					out<<pos21<<"\t"<<pos22<<"\n";
				}
				else
				{
					*result = 2;
					return ;
				}
			}
		}
		if(tmp < 0)
		{
			out<<pos21<<"\t"<<pos22<<"\n";
		}
		else if(tmp>0)
		{
			out<<pos11<<"\t"<<pos12<<"\n";
		}
		else
		{
			if(strlen(line1)>0)
				out<<line1<<"\n";
			if(strlen(line2)>0)
				out<<line2<<"\n";
		}
		while(!in1.eof())
		{
			in1.getline(line1, 256);
			if(strlen(line1)>0)
				out<<line1<<"\n";
		}
		while(!in2.eof())
		{
			in2.getline(line2, 256);
			if(strlen(line2)>0)
				out<<line2<<"\n";
		}
		in1.close();
		in2.close();
		out.flush();
		out.close();
		//		string rm = "rm " + inFile1 + " " + inFile2;
		//		system(rm.c_str());
		*result=0;
		return;
	}
	void  RNAProcessor::mergeAndSortFilesByAsciiTask(string inFile1, string inFile2, string outFile, int *result)
	{
		*result = 1;
		ifstream in1(inFile1.c_str(), fstream::in);
		ifstream in2(inFile2.c_str(), fstream::in);
		ofstream out(outFile.c_str(), ofstream::out);
		if(in1==NULL || in2==NULL || out==NULL)
		{
			*result =1;
			return ;
		}
		char *line1 = new char[256];
		char *line2 = new char[256];
		char *pos11=NULL;
		char *pos12=NULL;
		char *pos21=NULL;
		char *pos22=NULL;
		int tmp=0;
		while(true)
		{
			if(tmp==0)
			{
				in1.getline(line1, 256);
				in2.getline(line2, 256);

				if(strlen(line1)<1 || strlen(line2)<1)
					break;

				char *tab1=strchr(line1, '\t');
				char *tab2=strchr(line2, '\t');
				*tab1='\0';
				*tab2='\0';

				pos11=line1;
				pos21=line2;
				pos12=tab1+1;
				pos22=tab2+1;
			}
			else if(tmp<0)
			{
				in1.getline(line1, 256);

				if(strlen(line1)<1)
					break;

				char *tab1=strchr(line1, '\t');
				*tab1='\0';

				pos11=line1;
				pos12=tab1+1;
			}
			else
			{
				in2.getline(line2, 256);

				if(strlen(line2)<1)
					break;

				char *tab2=strchr(line2, '\t');
				*tab2='\0';

				pos21=line2;
				pos22=tab2+1;
			}

			tmp = strcmp(pos11, pos21);

			if(tmp<0)
			{
				out<<pos11<<"\t"<<pos12<<"\n";
			}
			else if(tmp>0)
			{
				out<<pos21<<"\t"<<pos22<<"\n";
			}
			else
			{
				int number = atoi(pos12) + atoi(pos22);
				out<<pos11<<"\t"<<number<<"\n";
			}
		}

		if(tmp==0)
		{
			if(strlen(line1)>0)
				out<<line1<<"\n";
			if(strlen(line2)>0)
				out<<line2<<"\n";
		}
		else if(tmp<0)
		{
			out<<pos21<<"\t"<<pos22<<"\n";
		}
		else if(tmp>0)
		{
			out<<pos11<<"\t"<<pos12<<"\n";
		}
		while(!in1.eof())
		{
			in1.getline(line1, 256);
			if(strlen(line1)>0)
				out<<line1<<"\n";
		}
		while(!in2.eof())
		{
			in2.getline(line2, 256);
			if(strlen(line2)>0)
				out<<line2<<"\n";
		}
		in1.close();
		in2.close();
		out.flush();
		out.close();
		//		string rm = "rm " + inFile1 + " " + inFile2;
		//		system(rm.c_str());
		*result=0;
		return;
	}

	string RNAProcessor::intToString(int i)
	{
		stringstream ss;
		ss<<i;
		return ss.str();
	}
}  // namespace SRNAProcessTool
