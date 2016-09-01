/*
 * RNAProcessor.cpp
 *
 *  Created on: 2012-6-19
 *      Author: Shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */
#include "SRNAProcessor.h"


using namespace boost::threadpool;

namespace SRNAProcessTool {


	int RNAProcessor::processFA()
	{
		FqInfo globleInfo;
		if (readLen_ - headTrim_ - tailTrim_ < 0){
			LOG(ERROR, "trim["<<headTrim_<<","<<tailTrim_<<"] is large then readLen"<<readLen_);
			exit(1);
		}
		globleInfo.readslength = readLen_ - headTrim_ - tailTrim_;
		globleInfo.lengthStart = minInsertSize_;
		globleInfo.lengthEnd = 44;

		gzFile outRawFq;
		if(outfq_){
			string outFqName;
			if(outFqName_.empty()){
				outFqName = getOutputFileName(fqFile1_, outDir_);
			}
			else
			{
				outFqName = outDir_ + "/" + outFqName_ + ".fq.gz";
			}

			if(cutOff_ == 0 && headTrim_ == 0 && tailTrim_ == 0 && !filterTile_ &&!filterIndex_ && fqFile1_.substr(fqFile1_.size() - 2, 2) == "gz"){
				string sys = "cp -u " + fqFile1_ + " " + outFqName;
				system(sys.c_str());
				outfq_ = false;
			}
			else
			{
				//gzFile outRawFq;
				outRawFq = gzopen(outFqName.c_str(), "wb");
				if (!outRawFq)
				{
					LOG(ERROR, "create output file: " + outFqName + " error");
					return 1;
				}
			}
		}


		long capacity = memLimit_;

		PreProcessTool::FqBuffer fqBuffer(fqFile1_.c_str(), capacity, PreProcessTool::FqBuffer::RB, filterTile_, tiles_);

		//	pool pl(threadNum_);

		TaskParam *params = new TaskParam[threadNum_];

		PreProcessTool::Read *reads;

		int size;
		unsigned long sizeTmp = 0;
		key_ = new char[readLen_+256];

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
			if (cutOff_ > 0 && sizeTmp > cutOff_)
			{
				size = size - (sizeTmp - cutOff_);
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

				pl_.schedule(boost::bind(&RNAProcessor::taskFA, this, &params[i]));
			}

			pl_.wait();

			for (int i=0; i<threadNum_; ++i)
			{
				globleInfo.add(params[i].info1, readLen_);
			}

			for (int i=0; i<size; ++i)
			{
				if(outfq_){
					gzputs(outRawFq, reads[i].readName);
					gzputs(outRawFq, "\n");
					gzputs(outRawFq, reads[i].baseSequence);
					gzputs(outRawFq, "\n");
					gzputs(outRawFq, reads[i].optionalName);
					gzputs(outRawFq, "\n");
					gzputs(outRawFq, reads[i].baseQuality);
					gzputs(outRawFq, "\n");
				}
				if (isFilter[i] == 0)
				{
					rawSequence_[string(reads[i].baseSequence)]++;
				}
			}

			delete []isFilter;
			filterRawTags(rawSequence_, globleInfo);
			rawSequence_.clear();
			if(cutOff_ > 0 && sizeTmp >= cutOff_)
				break;
		}
		pl_.wait();
		delete []key_;
		delete []params;
		if(outfq_)
		{
			gzclose(outRawFq);
		}
//		fqBuffer.freeFqBuffer();
		LOG(INFO, "Begin to filter Tags");
		return mergeTmpFilesAndPrint(globleInfo);
	}

	int RNAProcessor::printFiles(string &numberTmpFile, FqInfo &globleInfo, long length1)
	{
		ifstream in(numberTmpFile.c_str(), fstream::in);
		if(!in)
		{
			LOG(ERROR, "can't open numberTmpFile");
			return 1;
		}
		int setwTmp = 0;
		do
		{
			setwTmp++;
			length1/=10;

		}while(length1>0);

		int iTmp=0;
		for (int i=0; i<7 ;i++)
		{
			tagPercent[i] = 0;
		}
		string outFileName = outDir_ + "/" + outPfx_ + ".fa";

		ofstream faOut(outFileName.c_str());

		if (!faOut)
		{
			LOG(ERROR, "open output file "<<outFileName<<" error");
			exit(1);
		}
		long polyNCount[4][2]={{0,0},{0,0},{0,0},{0,0}}; //mrna_;
		ofstream polyNTxt, polyNFa, polyNStat, cleanTxt, smallTxt;
		string polyNTxtName = outDir_ + "/" + "polyN.txt";
		string polyNFaName = outDir_ + "/" + "nopolyN.fa";
		string polyNStatName = outDir_ + "/" + "polyN.stat";
		string cleanTxtName = outDir_ + "/" + "clean.txt";
		string smalltxtName = outDir_ + "/" + "small.txt";
		if(mrna_)
		{
			polyNTxt.open(polyNTxtName.c_str());
			polyNFa.open(polyNFaName.c_str());
			polyNStat.open(polyNStatName.c_str());
		}//mrna_;
		else
		{
			cleanTxt.open(cleanTxtName.c_str());
			smallTxt.open(smalltxtName.c_str());
		}

		char *line = new char[256];
		char *pos1=NULL;
		int pos2=0;
		while(!in.eof())
		{
			in.getline(line,256);
			if(strlen(line)<1)
				continue;
			iTmp++;
			char *tab = strchr(line, '\t');
			*tab = '\0';
			pos1 = line;
			pos2 = atoi(tab+1);
			faOut << ">t" <<setw(setwTmp)<<setfill('0')<<iTmp<<"\t"<<pos2<<"\n"<<pos1<<"\n";
			if(mrna_)
			{
				bool isPolyN = hasPolyN(pos1, iTmp, pos2, polyNTxt, polyNCount, setwTmp);
				if(!isPolyN)
				{
					polyNFa << ">t" <<setw(setwTmp)<<setfill('0')<<iTmp<<"\t"<<pos2<<"\n"<<pos1<<"\n";
				}
			}
			else
			{
				cleanTxt<<strlen(pos1)<<"\t"<<pos2<<"\t"<<pos1<<"\n";
			}
		}
		delete []line;
		in.close();
		faOut.close();
		string rm = "rm " + numberTmpFile;
		system(rm.c_str());
		if(mrna_)
		{
			polyNStat <<fixed;
			polyNStat << "\t#unique\t%\t#total\t%\n";
			polyNStat << "Total\t"<<iTmp<<"\t"<<100<<"\t"<<globleInfo.cleanTotalReads<<"\t"<<100<<"\n";
			polyNStat << "polyA\t"<<polyNCount[0][0]<<"\t"<<setprecision(2)<<100.0*polyNCount[0][0]/iTmp<<"\t"
				<<polyNCount[0][1]<<"\t"<<setprecision(2)<<100.0*polyNCount[0][1]/globleInfo.cleanTotalReads<<"\n";
			polyNStat << "polyC\t"<<polyNCount[1][0]<<"\t"<<setprecision(2)<<100.0*polyNCount[1][0]/iTmp<<"\t"
				<<polyNCount[1][1]<<"\t"<<setprecision(2)<<100.0*polyNCount[1][1]/globleInfo.cleanTotalReads<<"\n";
			polyNStat << "polyG\t"<<polyNCount[2][0]<<"\t"<<setprecision(2)<<100.0*polyNCount[2][0]/iTmp<<"\t"
				<<polyNCount[2][1]<<"\t"<<setprecision(2)<<100.0*polyNCount[2][1]/globleInfo.cleanTotalReads<<"\n";
			polyNStat << "polyT\t"<<polyNCount[3][0]<<"\t"<<setprecision(2)<<100.0*polyNCount[3][0]/iTmp<<"\t"
				<<polyNCount[3][1]<<"\t"<<setprecision(2)<<100.0*polyNCount[3][1]/globleInfo.cleanTotalReads<<"\n";
			polyNStat.close();
			polyNTxt.close();
			polyNFa.close();
			globleInfo.lengthStart = minSize_;
		}
		else
		{
			cleanTxt.close();
			for(map<string,int>::iterator it=shortSequence_.begin(); it!=shortSequence_.end(); it++)
			{
				smallTxt<<it->first.size()<<"\t"<<it->second<<"\t"<<it->first<<"\n";
			}
			smallTxt.close();
			globleInfo.lengthStart = 10;
		}

		output(&globleInfo, NULL);
		LOG(INFO, "PreProcess Finish");
		return 0;
	}

	void RNAProcessor::taskFA(TaskParam *param)
	{
		PreProcessTool::Read *reads1 = param->reads1;

		int start = param->left;
		int end = param->right;

		if (start == end)
		{
			return;
		}

		for (int i=start; i < end; ++i) //ask#*****************
		{	
			param->result[i] = statisticsFA(reads1[i], param->info1);
		}
	}
	int RNAProcessor::statisticsFA(PreProcessTool::Read &read, FqInfo &info)
	{
		StatisInfo si = auxStatistics(read, info);
		if (filterLowQual_ && (si.n > nRead_ || si.ns > nSeed_ || si.lowQual2 > lowQualSeed2Num_ || si.lowQual1 > lowQualSeed1Num_))
		{
			info.readsWithLowQual++;
			return 1;
		}

		return 0;
	}

	/*bool RNAProcessor::filterFA(const std::basic_string<char, std::char_traits<char>, std::allocator<char> >&sequence, int & count, FqInfo &info)*/
	bool RNAProcessor::filterFA(const string & sequence, int & count, FqInfo &info)
	{
		bool hasAdpt = false;
		int ik = 0;
		for (ik = 0; sequence[ik]; ik++)
		{
			key_[ik] = sequence[ik];
		}
		key_[ik] = '\0';

		int pos = findAdapter (key_, adapter3_.c_str());
		if (pos < 0)
		{
			info.adapter3Null += count;
			hasAdpt = true;
		}
		else if (pos <= 2)
		{
			info.insertNull += count;
			hasAdpt = true;
		}
		if (!hasAdpt)
		{
			if (trim_)
			{
				trimAdapter (key_, pos);
			}
			hasAdpt = hasAdapter (key_, adapter5_.c_str());
			if (hasAdpt)
			{
				info.adapter5Pollute += count;
			}
		}

		if(hasAdpt) //filter adapter;
		{
			info.readsWithAdapter += count;
			return true;
		}

		int readLen = strlen(key_);
		if(minSize_ > readLen)
			minSize_ = readLen;

		if(mrna_)
		{
			info.lengthDis[readLen] += count;
		}
		else
		{
			if(readLen >= 10)
			{
				info.lengthDis[readLen] += count;
			}
		}

		if (readLen < minInsertSize_) //filter short insert
		{
			info.readWithShortValidLength += count;
			if(!mrna_)
			{
				shortSequence_[string(key_)] += count;
			}
			return true;
		}

		int ai = 0, ci = 0, gi = 0, ti = 0, ni = 0;
		for (int i = 0; i < readLen; i++)
		{
			switch (key_[i])
			{
				case 'A':
					ai++;
					break;
				case 'C':
					ci++;
					break;
				case 'G':
					gi++;
					break;
				case 'T':
					ti++;
					break;
				case 'N':
					ni++;
					break;
			}
		}

		if (!mrna_ && filterPolyA_ >= 1E-6 && (1.0 * ai / readLen) > filterPolyA_) //filtr polyA
		{
			info.readWithPolyA += count;
			return true;
		}

		info.cleanBaseA += ai * count;
		info.cleanBaseC += ci * count;
		info.cleanBaseG += gi * count;
		info.cleanBaseT += ti * count;
		info.cleanBaseN += ni * count;

		info.cleanTotalBases += readLen * count;
		info.cleanTotalReads += count;
		//if(!mrna_)
		//{
		//	info.lengthDis[readLen] += count;
		//}

		return false;	//clean 
	}

	void RNAProcessor::filterRawTags(map<string, int> &rawSequence, FqInfo &info)
	{
		if(rawSequence.size()<1)
			return;
		map<string, int> cleanTags;
		for (map<string, int>::iterator it=rawSequence.begin(); it!=rawSequence.end(); it++)
		{
			bool filterKey = filterFA(it->first, it->second, info);
			if(!filterKey)
			{
				cleanTags[string(key_)] += it->second;
			}
			rawSequence.erase(it);
		}
		string outFile = outDir_ + "/" + "Tmp.FilterTmp" + intToString(filterTmpNumber_);
		ofstream out(outFile.c_str());
		for(map<string, int>::iterator it=cleanTags.begin(); it!=cleanTags.end(); it++)
		{
			out<<it->first<<"\t"<<it->second<<"\n";
		}
		out.flush();
		out.close();
		cleanTags.clear();
		tmpFile_.push_back(outFile);
		filterTmpNumber_++;
		rawSequence.clear();
		return;
	}
	int RNAProcessor::mergeTmpFilesAndPrint(FqInfo &info)
	{
		LOG(INFO, "Begin to merge and sort FilterTmp files");
		string outpfx = outDir_ + "/Tmp";
		string filterTmpLast = mergeAndSortFilesByAscii(outpfx, tmpFile_);
		if(filterTmpLast=="")
		{
			LOG(ERROR, "Can't merge and sort FilterTmp files into AsciiTmp File, NO READS PASS FILTER, CHECK!!!!");
			return 1;
		}
		LOG(INFO, "Begin to splite AsciiTmp File to SplitTmp Files");
		vector<string> splitTmpFile;
		long uniqReads = splitFile(outpfx, filterTmpLast, splitTmpFile);
		if(uniqReads<=0)
		{
			LOG(ERROR, "Splite AsciiTmp File to SplitTmp Files Error");
			return 1;
		}
		LOG(INFO, "Begin to merge and Sort  SplitTmp Files into NumberTmp File");
		string numberTmpLast =  mergeAndSortFilesByNumber(outpfx, splitTmpFile);
		if(numberTmpLast=="")
		{
			LOG(ERROR, "Can't merge and sort SplitTmp Files into NumberTmp File");
			return 1;
		}
		LOG(INFO, "Begion to print clean.fa and statistics");
		int presult = printFiles(numberTmpLast, info, uniqReads);
//		string rm = "rm " + outpfx + ".NumberTmp*";
//		system(rm.c_str());
		if(presult !=0)
		{
			LOG(ERROR, "print clean.fa and statistics error");
			return 1;
		}
		return 0;
	}
	/*
	   void RNAProcessor::filterTags(FqInfo &info)
	   {
	   key_ = new char[readLen_+1];
	   for (map<string,int>::iterator it=rawSequence_.begin(); it != rawSequence_.end(); it ++)
	   {
	   bool filterKey = filterFA((*it).first, (*it).second, info);

	   if (!filterKey)
	   {	
	   cleanSequence_[string(key_)] += (*it).second;
	   }

	   rawSequence_.erase(it);
	   }
	   long jt = 0;
	   for (map<string,int>::iterator it=cleanSequence_.begin(); it != cleanSequence_.end(); it ++)
	   {
	   sortClean_[(*it).second].push_back(jt);
	   cleanFa_.push_back(it);
	//		cleanTmp_.push_back(jt);
	jt++;
	}
	delete []key_;
	}
	 */
	bool RNAProcessor::hasPolyN(char* sequence, int title, int num, ofstream &outFile, long polyNCount[][2], int wT)
	{
		int len = strlen(sequence);
		int a=0,t=0,g=0,c=0;
		for(int i=0; i<len; i++)
		{
			switch(sequence[i])
			{
				case 'A':
					a++;
					break;
				case 'a':
					a++;
					break;
				case 'T':
					t++;
					break;
				case 't':
					t++;
					break;
				case 'G':
					g++;
					break;
				case 'g':
					g++;
					break;
				case 'C':
					c++;
					break;
				case 'c':
					c++;
					break;
				default:
					break;
			}
		}
		float ap = 1.0 * a / len;
		float tp = 1.0 * t / len;
		float gp = 1.0 * g / len;
		float cp = 1.0 * c / len;
		if(filterPolyN_ > 10e-6)
		{
			if(ap>filterPolyN_)
			{
				outFile <<"t"<<setw(wT)<<setfill('0')<<title<<"\tpolyA\t"<<ap<<"\t"<<len<<"\t"<<num<<"\t"<<sequence<<"\n";
				polyNCount[0][0]++;
				polyNCount[0][1]+=num;
				return true;
			}
			else if(cp>filterPolyN_)
			{
				polyNCount[1][0]++;
				polyNCount[1][1]+=num;
				outFile <<"t"<<setw(wT)<<setfill('0')<<title<<"\tpolyC\t"<<cp<<"\t"<<len<<"\t"<<num<<"\t"<<sequence<<"\n";
				return true;
			}
			else if(gp>filterPolyN_)
			{
				outFile <<"t"<<setw(wT)<<setfill('0')<<title<<"\tpolyG\t"<<gp<<"\t"<<len<<"\t"<<num<<"\t"<<sequence<<"\n";
				polyNCount[2][0]++;
				polyNCount[2][1]+=num;
				return true;
			}
			else if(tp>filterPolyN_)
			{
				outFile <<"t"<<setw(wT)<<setfill('0')<<title<<"\tpolyT\t"<<tp<<"\t"<<len<<"\t"<<num<<"\t"<<sequence<<"\n";
				polyNCount[3][0]++;
				polyNCount[3][1]+=num;
				return true;
			}
		}
		return false;
	}
	/*
	   void RNAProcessor::sortCleanFa(vector<long> &arrayPair,long x,long y)
	   {
	   long xx=x;
	   long yy=y;
	   long k = arrayPair[x];
	   if(x > y)
	   return;
	   while(xx != yy)
	   {
	   while(xx < yy && (*(cleanFa_[arrayPair[yy]])).second <= (*cleanFa_[k]).second)
	   yy--;
	   arrayPair[xx] = arrayPair[yy];

	   while(xx < yy && (*(cleanFa_[arrayPair[xx]])).second >= (*cleanFa_[k]).second)
	   xx++;
	   arrayPair[yy] = arrayPair[xx];
	   }
	   arrayPair[xx] = k;

	   sortCleanFa(arrayPair, x, xx-1);
	   sortCleanFa(arrayPair, xx+1, y);
	   }
	 */

}  // namespace SRNAProcessTool
