/*
 * FqBuffer.h
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#ifndef PREPROCESSTOOL_3_LOCAL_FQBUFFER_H_
#define PREPROCESSTOOL_3_LOCAL_FQBUFFER_H_

#include "CommonInclude.h"
#include "Common.h"

namespace PreProcessTool {

class FqBuffer
{
public:
	enum MODE
	{
		RB, WB,
	};
public:
	FqBuffer(const char *filename, int capacity, MODE mode, bool filterTile, const set<int> &tiles);
	FqBuffer(const char *filename, int capacity, MODE mode, bool filterTile, const set<int> &tiles, int seqType);
	~FqBuffer();

	Read* getReads();

	int getReadSize();

	int getRealReadSize();

	char *getBuf();

	int getLastIndex();
	
	bool isFilterTile(char *name, int seqType);
	void setSeqType(int type);


private:
	bool filterTile_;
	gzFile file_;
	char *buf_;
	char tile[5];
	int capacity_;
	int size_;
	Read *reads_;
	int readSize_;
	int realReadSize_;

	int lastIndex_;

	long lineNum_;
	set<int> tiles_;
	int seqType_;
};

}  // namespace PreProcessTool


#endif /* PREPROCESSTOOL_3_LOCAL_FQBUFFER_H_ */
