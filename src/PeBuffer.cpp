/*
 * PeBuffer.cpp
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#include "PeBuffer.h"
#include "CommonInclude.h"
#include "Logger.h"

namespace PreProcessTool
{

	PeBuffer::PeBuffer(const char *fq1Filename, const char *fq2Filename,
			int capacity, MODE mode, bool filterTile, const set<string> &tiles) :
		capacity_(capacity), size1_(0), size2_(0), readSize_(0), realReadSize_(0), initReadSize_(
				0), lastIndex1_(0), lastIndex2_(0), lineNum1_(0), lineNum2_(0), seqType_(0)
	{
		if (fq1Filename == NULL || fq2Filename == NULL)
		{
			LOG(ERROR, "fqFile1 or fqFile2 is NULL");
			exit(1);
		}

		if (mode == RB)
		{
			file1_ = gzopen(fq1Filename, "rb");
			file2_ = gzopen(fq2Filename, "rb");
		}
		else
		{
			file1_ = gzopen(fq1Filename, "wb");
			file2_ = gzopen(fq2Filename, "wb");
		}

		buf1_ = new char[capacity_ + 1];
		buf2_ = new char[capacity_ + 1];

		filterTile_ = filterTile;
		tiles_ = tiles;
	}
	PeBuffer::PeBuffer(const char *fq1Filename, const char *fq2Filename,int capacity, MODE mode, bool filterTile, const set<string> &tiles, int seqType){
		PeBuffer(fq1Filename,fq2Filename,capacity,mode,filterTile,tiles);
		seqType_ = seqType;
	}
    
	void PeBuffer::setSeqType(int type){
		seqType_ = type;
	}
    
    void PeBuffer::setTileIsFov(bool b){
        tileIsFov_ = b;
    }

	Read* PeBuffer::getReadsOne()
	{
		return reads1_;
	}

	Read* PeBuffer::getReadsTwo()
	{
		return reads2_;
	}

	void PeBuffer::readTask(int file, int &result)
	{
		result = 0;
		if (file == 1)
		{
			int remain1 = size1_ - lastIndex1_;
			if (remain1 > 0)
			{
				memcpy(buf1_, buf1_ + lastIndex1_, remain1);
			}
			int receiveSize1;
			size1_ = remain1;

			while ((receiveSize1 = gzread(file1_, buf1_ + size1_,
							capacity_ - size1_)) > 0)
			{
				size1_ += receiveSize1;
				if (size1_ == capacity_)
				{
					break;
				}
			}

			if (receiveSize1 == -1) //error
			{
				result = -1;
			}
			else if (receiveSize1 == 0)
			{
				if (size1_ == 0)
					return;
				if (buf1_[size1_ - 1] != '\n')
				{
					buf1_[size1_++] = '\n';
				}
			}
		}
		else if (file == 2)
		{
			int remain2 = size2_ - lastIndex2_;
			if (remain2 > 0)
			{
				memcpy(buf2_, buf2_ + lastIndex2_, remain2);
			}
			int receiveSize2;
			size2_ = remain2;

			while ((receiveSize2 = gzread(file2_, buf2_ + size2_,
							capacity_ - size2_)) > 0)
			{
				size2_ += receiveSize2;
				if (size2_ == capacity_)
				{
					break;
				}
			}
			if (receiveSize2 == -1)
			{
				result = -1;
			}
			else if (receiveSize2 == 0)
			{
				if(size2_ == 0)
					return;
				if (buf2_[size2_ - 1] != '\n')
				{
					buf2_[size2_++] = '\n';
				}
			}
		}
		else
		{
			LOG(ERROR, "NO SUCH TYPE in PeBuffer::readTask");
			exit(1);
		}
	}

	int PeBuffer::getReads()
	{
		int result1 = 0, result2 = 0;
		readTask(1, result1);
		readTask(2, result2);

		//end
		if (result1 == 0 && size1_ == 0 && result2 == 0 && size2_ == 0)
		{
			return 0;
		}
		//one file end and another file not end,
		//so one file's read number is not equal to the other
		else if ((result1 == 0 && size1_ == 0 && (result2 != 0 || size2_ != 0))
				|| ((result1 != 0 || size1_ != 0) && result1 == 0 && size2_ == 0))
		{
			return -2;
		}
		else if (result1 == -1 || result2 == -1)
		{
			return -1;
		}

		//calculate the number of read in buffer
		if (initReadSize_ == 0)
		{
			int lines = 0;
			int tempSize1 = 0;
			for (int i = 0; i < size1_; ++i)
			{
				if (buf1_[i] == '\n')
				{
					lines++;
					if (lines == 4)
					{
						lines = 0;
						tempSize1++;
					}
				}
			}

			lines = 0;
			int tempSize2 = 0;
			for (int i = 0; i < size2_; ++i)
			{
				if (buf2_[i] == '\n')
				{
					lines++;
					if (lines == 4)
					{
						lines = 0;
						tempSize2++;
					}
				}
			}

			initReadSize_ = (tempSize1 > tempSize2 ? tempSize2 : tempSize1);

			reads1_ = new Read[initReadSize_];
			reads2_ = new Read[initReadSize_];
			nameIndex1_ = new char*[initReadSize_];
			nameIndex2_ = new char*[initReadSize_];
		}

		lineNum1_ += (readSize_ * 4);
		long tempNum = lineNum1_;

		int readNum1 = 0, readNumTmp1 = 0, index1 = 0, lines = 0;
		int pos[4], last1 = -1;

		while (index1 < size1_ && readNumTmp1 < initReadSize_)
		{
			if (buf1_[index1] == '\n' || buf1_[index1] == '\0')
			{
				pos[lines] = index1;
				lines++;
				tempNum++;
				if (lines == 4)
				{
					lines = 0;
					buf1_[pos[0]] = '\0';
					buf1_[pos[1]] = '\0';
					buf1_[pos[2]] = '\0';
					buf1_[pos[3]] = '\0';

					nameIndex1_[readNumTmp1++] = buf1_ + last1 + 1;
                    
                    if(filterTile_)
                    {
                        if(tileIsFov_)
                        {
                            if(isFilterFov(buf1_ + last1 + 1, seqType_))
                            {
                                last1 = pos[3];
                                index1++;
                                continue;
                            }
                        }else if(isFilterTile(buf1_ + last1 + 1, seqType_))
                        {
                            last1 = pos[3];
                            index1++;
                            continue;
                        }
                    }

					reads1_[readNum1].readName = buf1_ + last1 + 1;
					reads1_[readNum1].baseSequence = buf1_ + pos[0] + 1;
					reads1_[readNum1].optionalName = buf1_ + pos[1] + 1;
					reads1_[readNum1].baseQuality = buf1_ + pos[2] + 1;

					if (strlen(reads1_[readNum1].baseSequence)
							< strlen(reads1_[readNum1].baseQuality))
					{
						string temp = "" + tempNum;
						LOG(ERROR,
								"the length of base sequence and base quality are not equal in fq1, line number: " + temp);
						exit(1);
					}
					last1 = pos[3];
					readNum1++;
				}
			}
			index1++;
		}

		lineNum2_ += (readSize_ * 4);
		tempNum = lineNum2_;

		int readNum2 = 0, readNumTmp2 = 0, index2 = 0;
		int last2 = -1;
		lines = 0;
		while (index2 < size2_ && readNumTmp2 < initReadSize_)
		{
			if (buf2_[index2] == '\n' || buf2_[index2] == '\0')
			{
				pos[lines] = index2;
				lines++;
				tempNum++;
				if (lines == 4)
				{
					lines = 0;
					buf2_[pos[0]] = '\0';
					buf2_[pos[1]] = '\0';
					buf2_[pos[2]] = '\0';
					buf2_[pos[3]] = '\0';

					nameIndex2_[readNumTmp2++] = buf2_ + last2 + 1;
                    if(filterTile_)
                    {
                        if(tileIsFov_)
                        {
                            if(isFilterFov(buf2_ + last2 + 1, seqType_))
                            {
                                last2 = pos[3];
                                index2++;
                                continue;
                            }
                        }else if(isFilterTile(buf2_ + last2 + 1, seqType_))
                        {
                            last2 = pos[3];
                            index2++;
                            continue;
                        }
                    }
                    
					reads2_[readNum2].readName = buf2_ + last2 + 1;
					reads2_[readNum2].baseSequence = buf2_ + pos[0] + 1;
					reads2_[readNum2].optionalName = buf2_ + pos[1] + 1;
					reads2_[readNum2].baseQuality = buf2_ + pos[2] + 1;

					if (strlen(reads2_[readNum2].baseSequence)
							< strlen(reads2_[readNum2].baseQuality))
					{
						string temp = "" + tempNum;
						LOG(ERROR,
								"the length of base sequence and base quality are not equal in fq1, line number: " + temp);
						exit(1);
					}

					last2 = pos[3];

					readNum2++;
				}
			}
			index2++;
		}

		if (readNumTmp1 > readNumTmp2)
		{
			readSize_ = readNumTmp2;
			lastIndex1_ = nameIndex1_[readNumTmp2] - buf1_;
			lastIndex2_ = last2 + 1;
		}
		else if (readNumTmp1 < readNumTmp2)
		{
			readSize_ = readNumTmp1;
			lastIndex1_ = last1 + 1;
			lastIndex2_ = nameIndex2_[readNumTmp1] - buf2_;
		}
		else
		{
			readSize_ = readNumTmp1;
			lastIndex1_ = last1 + 1;
			lastIndex2_ = last2 + 1;
		}

		realReadSize_ = readNum1 < readNum2 ? readNum1 : readNum2;
		return readSize_;
	}

	char* PeBuffer::getBuf1()
	{
		return buf1_;
	}
	char* PeBuffer::getBuf2()
	{
		return buf2_;
	}
	int PeBuffer::getlastIndex1()
	{
		return lastIndex1_;
	}
	int PeBuffer::getlastIndex2()
	{
		return lastIndex2_;
	}

	int PeBuffer::getInitReadSize()
	{
		return initReadSize_;
	}
	int PeBuffer::getReadSize()
	{
		return realReadSize_;
	}

	PeBuffer::~PeBuffer()
	{
		if (buf1_ != NULL)
			delete[] buf1_;
		if (buf2_ != NULL)
			delete[] buf2_;
		if (reads1_ != NULL)
			delete[] reads1_;
		if (reads2_ != NULL)
			delete[] reads2_;
		if (file1_ != NULL)
			gzclose(file1_);
		if (file2_ != NULL)
			gzclose(file2_);
		if (nameIndex1_ != NULL)
			delete [] nameIndex1_;
		if (nameIndex2_ != NULL)
			delete [] nameIndex2_;
	}
	bool PeBuffer::isFilterTile(char *name, int seqType)
	{
		int num=0;
		int len = strlen(name);
		int i=0;
		if(seqType_==0){
			for(i=0; i<len; i++)
			{
				if(name[i]==':')
					num++;
				if(num>=2)
					break;
			}
			if(num<2 || i+4>=len)
				return false;
		}else{
			for(i=0; i<len; i++)
			{
				if(name[i]==':')
					num++;
				if(num>=4)
					break;
			}
			if(num<4 || i+4>=len)
				return false;
		}

		for(int j=0; j<5; j++)
		{
            if(name[i+j+1]>='0' && name[i+j+1] <= '9')
            {
                tile[j]=name[i+j+1];
            }else
            {
                tile[j] = '\0';
            }
		}
		tile[5]='\0';

        return tiles_.count(tile);
	}
    
    bool PeBuffer::isFilterFov(char *name, int seqType)
    {
        long len = strlen(name);
        int i=0;
        
        if(seqType==0){
            for(i=0; i<len; i++)
            {
                if(name[i]=='C')
                    if(i+8<len && name[i+4]=='R' && name[i+8]=='_')
                        break;
                
            }
        }else{
            LOG(ERROR, "Zebra-500 data(--fov), --seqType is 0");
        }
        
        for(int j=0; j<8; j++)
        {
            tile[j]=name[i+j];
        }
        tile[8]='\0';
        return tiles_.count(tile);
    }
} // namespace PreProcessTool

