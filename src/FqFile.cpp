#include "FqFile.h"
#include <string>
using namespace std;

namespace PreProcessTool {

bool FqFile::nextRead(StrRead& read1, StrRead& read2)
{
    string temp;
    int index;

    getline(ifs, temp); 
    if (ifs.eof())
    {
        return false;
    }
    index = sapaceIndex(temp);
    read1.readName = temp.substr(0, index);
    read2.readName = temp.substr(index+1);

    getline(ifs, temp);
    index = sapaceIndex(temp);
    read1.baseSequence = temp.substr(0, index);
    read2.baseSequence = temp.substr(index+1);

    getline(ifs, temp);
    index = sapaceIndex(temp);
    read1.optionalName = temp.substr(0, index);
    read2.optionalName = temp.substr(index+1);

    getline(ifs, temp);
    index = sapaceIndex(temp);
    read1.baseQuality = temp.substr(0, index);
    read2.baseQuality = temp.substr(index+1);

    return true;
}

bool FqFile::nextRead(StrRead& read)
{
    string temp;
    getline(ifs, temp);
    if (ifs.eof())
        return false;
    read.readName = temp;
    getline(ifs, temp);
    read.baseSequence = temp;
    getline(ifs, temp);
    read.optionalName = temp;
    getline(ifs, temp);
    read.baseQuality = temp;

    return true;
}

int FqFile::sapaceIndex(const string& str)
{
    int len = str.length();
    for (int i=0; i<len; ++i)
    {
        if (str[i] == '\t')
            return i;
    }

    return -1;
}

}
