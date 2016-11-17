/*
 * RNAProcessor.cpp
 *
 *  Created on: 2012-6-3
 *      Author: Shuai JIANG
 * 		Mail  : jiangshuai@genomics.cn
 */
#include "SRNAProcessor.h"


using namespace boost::threadpool;

namespace SRNAProcessTool {

	int RNAProcessor::processFQ()
	{
		FqInfo globleInfo;
		if (readLen_ - headTrim_ - tailTrim_ < 0){
			LOG(ERROR, "trim["<<headTrim_<<","<<tailTrim_<<"] is large then readLen"<<readLen_);
			exit(1);
		}
		globleInfo.readslength = readLen_ - headTrim_ - tailTrim_;
		globleInfo.lengthStart = 10;
		globleInfo.lengthEnd = 44;

		gzFile outFile;
		string outFileName;
		gzFile outRawFq = NULL;
		int headTrimTmp = headTrim_, tailTrimTmp = tailTrim_;
		bool filterIndexTmp = filterIndex_;
		int rightTmp;
		if(outfq_){
			string outFqName;
			if(outFqName_.empty()){
				outFqName = getOutputFileName(fqFile1_, outDir_);
			}
			else
			{
				outFqName = outDir_ + "/" + outFqName_ + ".fq.gz";
			}

			if(cutOff_ == 0 && headTrim_ == 0 && tailTrim_ == 0 && !filterTile_ && !filterIndex_ && fqFile1_.substr(fqFile1_.size() - 2, 2) == "gz"){
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
			headTrim_ = 0;
			tailTrim_ = 0;
			filterIndex_ = false;
		}

		//	outFileName = getOutputFileName(fqFile1_, outDir_);
		outFileName = outDir_ + "/" + outPfx_ + ".fq.gz";

		outFile = gzopen(outFileName.c_str(), "wb");
		if (!outFile)
		{
			LOG(ERROR, "create output file:¡¡" + outFileName + " error");
			return 1;
		}


		long capacity = memLimit_;

		PreProcessTool::FqBuffer fqBuffer(fqFile1_.c_str(), capacity, PreProcessTool::FqBuffer::RB, filterTile_, tiles_);
		fqBuffer.setSeqType(seqType_);
		//	pool pl(threadNum_);

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
			if (cutOff_ > 0 && sizeTmp > cutOff_)
			{
				size = size - (sizeTmp - cutOff_);
			}
			int *isFilter = new int[size];
			bzero(isFilter, sizeof(int) * size);

			int block = size / threadNum_;
			int remain = size % threadNum_;
			int index = 0;
			if(outfq_)
			{
				if(filterIndexTmp)
				{
					if(seqType_ == 0){
						for(int i=0; i<size; ++i)
						{
							int sharpIndex = 0;
							int j = 0;
							while (reads[i].readName[j++] != '#')
								;
							sharpIndex = j;
							while (reads[i].readName[j++] != '/')
								;
							strcpy(reads[i].readName + sharpIndex, reads[i].readName + j - 1);
						}

					}else{
						for(int i=0; i<size; ++i){
							int j = strlen(reads[i].readName);
							while (reads[i].readName[--j] !=':')
								;
							if(j>0)
								reads[i].readName[j]='\0';
						}
					}
				}
				if(headTrimTmp > 0)
				{
					for(int i=0; i<size; ++i)
					{
						reads[i].baseSequence = reads[i].baseSequence + headTrimTmp;
						reads[i].baseQuality = reads[i].baseQuality + headTrimTmp;
					}
				}
				if(tailTrimTmp > 0)
				{
					for(int i=0; i<size; ++i)
					{
						rightTmp = strlen(reads[i].baseSequence) - tailTrimTmp;
						reads[i].baseSequence[rightTmp] = '\0';
						reads[i].baseQuality[rightTmp] = '\0';
					}
				}
				for(int i=0; i<size; ++i)
				{
					gzputs(outRawFq, reads[i].readName);
					gzputs(outRawFq, "\n");
					gzputs(outRawFq, reads[i].baseSequence);
					gzputs(outRawFq, "\n");
					gzputs(outRawFq, reads[i].optionalName);
					gzputs(outRawFq, "\n");
					gzputs(outRawFq, reads[i].baseQuality);
					gzputs(outRawFq, "\n");
				}
			}

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
				//			bzero(isFilter, sizeof(int) * size);
				params[i].result = isFilter;

				//			params[i].type = type;

				clearFqInfo(params[i].info1);

				pl_.schedule(boost::bind(&RNAProcessor::taskFQ, this, &params[i]));
			}

			pl_.wait();

			for (int i=0; i<threadNum_; ++i)
			{
				globleInfo.add(params[i].info1, readLen_);
			}

			for (int i=0; i<size; ++i)
			{
				if (isFilter[i] == 0)
				{
					gzputs(outFile, reads[i].readName);
					gzputs(outFile, "\n");
					gzputs(outFile, reads[i].baseSequence);
					gzputs(outFile, "\n");
					gzputs(outFile, reads[i].optionalName);
					gzputs(outFile, "\n");
					gzputs(outFile, reads[i].baseQuality);
					gzputs(outFile, "\n");
				}
			}

			delete []isFilter;

			if(cutOff_ > 0 && sizeTmp >= cutOff_)
				break;
		}


		delete []params;
		if(outfq_)
		{
			gzclose(outRawFq);
		}

		gzclose(outFile);

		output(&globleInfo, NULL);

		LOG(INFO, "PreProcess Finish");

		return 0;
	}



	void RNAProcessor::taskFQ(TaskParam *param)
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
			param->result[i] = statisticsFQ(reads1[i], param->info1);
		}
	}




	int RNAProcessor::statisticsFQ(PreProcessTool::Read &read, FqInfo &info)
	{	
		StatisInfo si = auxStatistics(read, info);

		bool low = (filterLowQual_ && (si.n > nRead_ || si.ns > nSeed_ || si.lowQual2 > lowQualSeed2Num_ || si.lowQual1 > lowQualSeed1Num_));
		if (low) //filter low quality
		{
			info.readsWithLowQual++;
			return 1;
		}


		bool hasAdpt = false;
		int pos = findAdapter (read.baseSequence, adapter3_.c_str());
		if (pos < 0)
		{
			info.adapter3Null++;
			hasAdpt = true;
		}
		else if (pos <= 2)
		{
			info.insertNull++;
			hasAdpt = true;
		}

		if (!hasAdpt)
		{
			if (trim_)
			{
				trimAdapter (read.baseSequence, pos);
				trimAdapter (read.baseQuality, pos);
			}
			hasAdpt = hasAdapter (read.baseSequence, adapter5_.c_str());
			if (hasAdpt)
			{
				info.adapter5Pollute++;
			}
		}

		if (hasAdpt) //filter adapter
		{	
			info.readsWithAdapter++;
			return 1;
		}

		int readLen = strlen(read.baseSequence);
		if (readLen >= 10)
			info.lengthDis[readLen]++;

		if (readLen < minInsertSize_) //filter short insert
		{
			info.readWithShortValidLength++;
			return 1;
		}

		int ai = 0, ci = 0, gi = 0, ti = 0, ni = 0, q20i = 0, q30i = 0, quali;
		for (int i = 0; i < readLen; i++)
		{
			switch (read.baseSequence[i])
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

			quali = read.baseQuality[i] - qualSystem_;
			read.baseQuality[i] = quali + score_;
			if (quali >= 20)
			{
				q20i++;
				if (quali >= 30)
				{
					q30i++;
				}
			}
		}

		if (filterPolyA_ >= 1E-6 && (1.0 * ai / readLen) > filterPolyA_) //filPer polyA;
		{
			info.readWithPolyA++;
			return 1;
		}

		info.cleanBaseA += ai;
		info.cleanBaseC += ci;
		info.cleanBaseG += gi;
		info.cleanBaseT += ti;
		info.cleanBaseN += ni;
		info.cleanQ20 += q20i;
		info.cleanQ30 += q30i;

		info.cleanTotalBases += readLen;
		info.cleanTotalReads++;
		//	info.lengthDis[readLen]++;

		return 0;
	}

}  // namespace SRNAProcessTool
