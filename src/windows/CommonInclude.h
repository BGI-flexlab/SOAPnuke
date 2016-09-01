/*
 * CommonInclude.h
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#ifndef PREPROCESSTOOL_3_LOCAL_COMMONINCLUDE_H_
#define PREPROCESSTOOL_3_LOCAL_COMMONINCLUDE_H_

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <algorithm> 

#include <zlib.h>

#include <getopt.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>

#pragma comment(lib, "libeay32.lib")  
#pragma comment(lib, "ssleay32.lib")

#define bzero(a, b) memset(a, 0, b)

#include "threadpool.hpp"

using namespace std;

#endif /* PREPROCESSTOOL_3_LOCAL_COMMONINCLUDE_H_ */
