#ifndef PREPROCESSTOOL_3_LOCAL_FQFILE_H_
#define PREPROCESSTOOL_3_LOCAL_FQFILE_H_ 

#include <iostream>
#include <fstream>
#include <string>

#include "Common.h"


namespace PreProcessTool {

class FqFile 
{
public:
    FqFile(){}
    /*FqFile(std::string filename):name(filename), ifs(filename)
    {
    }*/

    bool open(std::string filename)
    {
        name = filename;
        ifs.open(filename.c_str(), std::ifstream::in);
        if (ifs)
            return true;
        else
            return false;
    }

    ~FqFile()
    {
        ifs.close();
    }

    bool nextRead(StrRead& read1, StrRead& read2);
    bool nextRead(StrRead& read);

private:
    int sapaceIndex(const std::string& str);

private:
    std::string name;
    std::ifstream ifs;
};

}

#endif 
