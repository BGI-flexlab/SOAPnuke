/*
 * Main.cpp
 *
 *  Created on: 2012-6-17
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#include "CommonInclude.h"
#include "FilterProcessor.h"
#include "SRNAProcessor.h"
#include "DGEProcessor.h"
#include "MetaProcessor.h"

using namespace std;
using namespace PreProcessTool;
using namespace SRNAProcessTool;
using namespace DGEProcessTool;
using namespace MetaPreProcessTool;

#ifndef PACKAGEVERSION
#define PACKAGEVERSION "1.6.0"
#endif

int usage()
{
	cout << endl;
	cout << "Prpgram: soapnuke\n";
	cout << "Version: " << PACKAGEVERSION <<endl;
	cout << "Contact: YoungChan<chenyuxin@genomics.cn>\n";
	cout << "Command:\n";
	cout << "         filter        preprocessing sequences\n";
	cout << "         filtersRNA    preprocessing sRNA sequences\n";
	cout << "         filterDGE     preprocessing DGE sequences\n";
	cout << "         filterMeta    preprocessing Meta sequences\n";
	cout << endl;
	return 0;
}

int main(int argc,char *argv[])
{
	if(argc<2)
	{
		return usage();
	}
	for (int i=0; argv[1][i]; ++i)
	{
		argv[1][i] = toupper(argv[1][i]);
	}
	if(0==strcmp(argv[1],"FILTER"))
	{
		FilterProcessor filterProcessor;
		return filterProcessor.filter(argc-1, argv+1);
	}
	else if(0==strcmp(argv[1],"FILTERSRNA"))
	{
        RNAProcessor rnaProcessor;
        return rnaProcessor.processRNA(argc-1, argv+1);
	}
    else if (0 == strcmp(argv[1], "FILTERDGE"))
    {
        DGEProcessor dgeProcessor;
        return dgeProcessor.processRNA(argc-1, argv+1);
    }
	else if (0 == strcmp(argv[1], "FILTERMETA"))
	{
		MetaProcessor metaProcessor;
		return metaProcessor.filterMeta(argc-1, argv+1);
	}
	else
	{
		cerr<<"[main] unrecognized command "<<argv[1]<<endl;
		return 1;
	}
	return 0;

}


