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
	FqBuffer(gzFile fqStreaming, int capacity, bool filterTile, const set<string> &tiles);
	FqBuffer(const char *filename, int capacity, MODE mode, bool filterTile, const set<string> &tiles);
	FqBuffer(const char *filename, int capacity, MODE mode, bool filterTile, const set<string> &tiles, int seqType);
	~FqBuffer();

	Read* getStreamingReads();
	Read* getReads();

	int getReadSize();

	int getRealReadSize();

	char *getBuf();

	int getLastIndex();
	
	bool isFilterTile(char *name, int seqType);
    bool isFilterFov(char *name, int seqType);
    
	void setSeqType(int type);
    void setTileIsFov(bool b);


private:
    bool IS_STREAMING;
	bool filterTile_;
    bool tileIsFov_;
    char tile[9];
    set<string> tiles_;
    
	gzFile file_;
	char *buf_;
	int capacity_;
	int size_;
	Read *reads_;
	int readSize_;
	int realReadSize_;
	int lastIndex_;
	long lineNum_;
	int seqType_;
};

}  // namespace PreProcessTool


#endif /* PREPROCESSTOOL_3_LOCAL_FQBUFFER_H_ */
