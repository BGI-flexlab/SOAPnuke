#ifndef _PROCESS_ARGV_H
#define _PROCESS_ARGV_H

#include <string>
#include <vector>
//#include "module.h"
#include "global_parameter.h"
#include "gc.h"
using namespace::std;

#ifndef PACKAGEVERSION
#define PACKAGEVERSION "2.1"
#endif
#ifndef MINORVERSION
#define MINORVERSION "5"
#endif
void check_module(int argc,char* argv[]);
int global_parameter_initial(int argc,char* argv[],C_global_parameter& gp);
bool check_parameter(int argc,char* argv[],C_global_parameter& gp);
void printUsage(string c_module);
void printVersion();
void printModule();
void printHtsUsage();
void initFromConfigFile(C_global_parameter& gp,char* configFile);
#endif