/*
 * FqBuffer.cpp
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */
#include "FqBuffer.h"
#include "CommonInclude.h"
#include "Logger.h"

namespace PreProcessTool {

	FqBuffer::FqBuffer(const char *filename, int capacity, MODE mode, bool filterTile, const set<int> &tiles) :
		capacity_(capacity), size_(0), readSize_(0), realReadSize_(0), lastIndex_(0), lineNum_(0),seqType_(0)
	{
		if (filename == NULL)
		{
			LOG(ERROR, "No FileName");
			exit(1);
		}
		if (mode == RB)
			file_ = gzopen(filename, "rb");
		else
			file_ = gzopen(filename, "wb");

		if (file_ == NULL)
		{
			LOG(ERROR, "open file: " + string(filename) + " error");
			exit(1);
		}

		buf_ = new char[capacity_ + 1];

		filterTile_ = filterTile;
		tiles_ = tiles;
	}
	FqBuffer::FqBuffer(const char *filename, int capacity, MODE mode, bool filterTile, const set<int> &tiles, int seqType){
		FqBuffer(filename,capacity,mode,filterTile,tiles);
		seqType_ = seqType;
	}
	void FqBuffer::setSeqType(int type){
		seqType_ = type;
	}

	Read* FqBuffer::getReads()
	{
		int remain = size_ - lastIndex_;
		memcpy(buf_, buf_ + lastIndex_, remain);

		int receiveSize;
		size_ = remain;

		while ((receiveSize = gzread(file_, buf_ + size_, capacity_ - size_)) > 0)
		{
			size_ += receiveSize;
			if (size_ == capacity_)
			{
				break;
			}
		}

		if (receiveSize == -1) //error
		{
			readSize_ = -1;
			realReadSize_ = -1;
			return NULL;
		}
		else if (receiveSize == 0) //end of file
		{
			if (size_ == 0) //end of buffer
			{
				readSize_ = 0;
				realReadSize_ = 0;
				return NULL;
			}
			if (buf_[size_ - 1] != '\n')
			{
				buf_[size_++] = '\n';
			}
		}

		lineNum_ += (readSize_ * 4);

		if (readSize_ == 0)
		{
			int lines = 0;
			//calculate the number of read in buffer
			for (int i = 0; i < size_; ++i)
			{
				if (buf_[i] == '\n')
				{
					lines++;
					if (lines == 4)
					{
						lines = 0;
						readSize_++;
					}
				}
			}

			reads_ = new Read[readSize_];
		}

		long tempNum = lineNum_;
		int j = 0, i = 0, lines = 0;
		int pos[4], lastIndex = -1;
		realReadSize_ = 0;
		for (i = 0; i < readSize_; ++i)
		{
			while (j < size_)
			{
				if (buf_[j] == '\n')
				{
					pos[lines] = j;
					lines++;
					tempNum++;
					if (lines == 4)
					{
						lines = 0;
						buf_[pos[0]] = '\0';
						buf_[pos[1]] = '\0';
						buf_[pos[2]] = '\0';
						buf_[pos[3]] = '\0';

						if(filterTile_ && isFilterTile(buf_ + lastIndex + 1, seqType_))
						{
							lastIndex = pos[3];
							j++;
							break;
						}

						reads_[realReadSize_].readName = buf_ + lastIndex + 1;
						reads_[realReadSize_].baseSequence = buf_ + pos[0] + 1;
						reads_[realReadSize_].optionalName = buf_ + pos[1] + 1;
						reads_[realReadSize_].baseQuality = buf_ + pos[2] + 1;

						if (strlen(reads_[realReadSize_].baseSequence) > strlen(reads_[realReadSize_].baseQuality))
						{
							string temp = "" + tempNum;
							LOG(ERROR, "the length of base sequence and base quality are not equal, line number: " + temp);
							exit(1);
						}
						realReadSize_++;

						lastIndex = pos[3];
						j++;
						break;
					}
				}
				j++;
			}

			if (j >= size_)
			{
				break;
			}
		}

		if (j >= size_)
		{
			if (lastIndex == size_ - 1)
			{
				readSize_ = i + 1;
			}
			else
			{
				readSize_ = i;
			}
		}

		lastIndex_ = lastIndex + 1;

		return reads_;
	}

	char* FqBuffer::getBuf()
	{
		return buf_;
	}

	int FqBuffer::getLastIndex()
	{
		return lastIndex_;
	}

	int FqBuffer::getReadSize()
	{
		return readSize_;
	}

	int FqBuffer::getRealReadSize()
	{
		return realReadSize_;
	}

	FqBuffer::~FqBuffer()
	{
		if (buf_ != NULL)
			delete[] buf_;
		if (reads_ != NULL)
			delete[] reads_;
		if (file_ != NULL)
			gzclose(file_);
	}

	bool FqBuffer::isFilterTile(char *name, int seqType)
	{
		int num=0;
		int len = strlen(name);
		int i=0;
		int tileNum=0;

		if(seqType==0){
			for(i=0; i<len; i++)
			{
				if(name[i]==':')
					num++;
				if(num>=2)
					break;
			}
			if(num<2 || i+4>=len)
				return false;
			for(int j=0; j<5; j++)
			{
				tile[j]=name[i+j+1];
			}
			tile[5]='\0';
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
			for(int j=0; j<5; j++)
			{
				tile[j]=name[i+j+1];
			}
			tile[5]='\0';
		}

		tileNum = atoi(tile);
		return tiles_.count(tileNum);
	}
}  // namespace PreProcessTool

