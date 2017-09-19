/*
 * DNAProcessor.cpp
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */

#include "FilterProcessor.h"

#include "Common.h"
#include "FqBuffer.h"
#include "PeBuffer.h"
#include "Logger.h"
#include "FqFile.h"

#include "threadpool.hpp"
#include <fstream>

#include <openssl/md5.h>
#include <boost/thread.hpp>
#include <boost/interprocess/detail/atomic.hpp>

#include <sys/stat.h>

using namespace boost;
using namespace boost::threadpool;
using namespace boost::interprocess::ipcdetail;

namespace PreProcessTool {

    void FilterProcessor::printVersion()
    {
        cerr << "soapnuke filter tools version 1.6.0\n";
        cerr << "Author:  chenhaosen\n";
        cerr << " Email:  chenhaosen@genomics.cn\n";
    }

    void FilterProcessor::printUsage()
    {
        cout << "Usage: filter [OPTION]... \n";
        cout << "\t-1, --fq1         FILE      fq1 file\n";
        cout << "\t-2, --fq2         FILE      fq2 file, used to pe\n";
        cout << "\t  --tile          STR       tile number to ignore reads, such as [1101-1104,1205]\n";
        cout << "\t  --fov           STR       fov number to ignore reads (only for zebra-platform data), such as [C001R003,C003R004]\n";
        cout << "\n";
        cout << "\t-f, --adapter1    STR       3' adapter sequence of fq1 file\n";
        cout << "\t-r, --adapter2    STR       5' adapter sequence of fq2 file (only for PE reads)\n";
//		cout << "\t    --cutAdaptor  INT[,INT] cut adaptor sequence, discard the read when the adaptor index of the read is less than INT, cut 3'-end or 5'-end (depend on -f/-r) [50,3]\n";
//      cout << "\t    --cutAdaptor  INT       cut adaptor sequence, discard the read when the adaptor index of the read is less than INT(depend on -f/-r) [50]\n";
        cout << "\t    --cutAdaptor  INT       cut adaptor sequence\n";
//		cout << "\t    --BaseNum     INT       the base number you want to keep in each clean fq file, (depend on --cutAdaptor)\n";
        cout << "\t    --misMatch    INT       the max mismatch number when match the adapter (depend on -f/-r)  [1]\n";
        cout << "\t    --matchRatio  FLOAT     adapter's shortest match ratio (depend on -f/-r)  [0.5]\n";
        cout << "\n";

        cout << "\t-l, --lowQual     INT       low quality threshold  [5]\n";
        cout << "\t-q, --qualRate    FLOAT     low quality rate  [0.5]\n";
        cout << "\t-n, --nRate       FLOAT     N rate threshold  [0.05]\n";
        cout << "\t-N, --maskLowQual INT       Turn those bases with low quality into N, set INT as the quality threshold  [-1]\n";
        cout << "\t-m, --mean        FLOAT     filter reads with low average quality, (<) \n";
        cout << "\t-p, --polyA       FLOAT     filter poly A, percent of A, 0 means do not filter   [0]\n";
//		cout << "\t-d, --rmdup                 remove PCR duplications\n";
//		cout << "\t-3, --dupRate               calculate PCR duplications rate only,but don't remove PCR duplication reads\n";
        cout << "\t-i, --index                 remove index\n";
        cout << "\t-c, --cut         FLOAT     the read number you want to keep in each clean fq file, (unit:1024*1024, 0 means not cut reads)\n";
        cout << "\t-t, --trim        INT,INT,INT,INT" << endl;
        cout <<	"\t                            trim some bp of the read's head and tail, they means: (read1's head and tail and read2's head and tail  [0,0,0,0]\n";
        cout << "\t  --trimBadTail   INT,INT   Trim from tail ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]\n";
        cout << "\t  --trimBadHead   INT,INT   Trim from head ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]\n";
        cout << "\n";

        cout << "\t-S, --small                 filter the small insert size\n";
        cout << "\t    --overlap     INT       minimun match length (depend on -S)  [10]\n";
        cout << "\t    --mis         FLOAT     the maximum miss match ratio (depend on -S)  [0.1]\n";
        cout << "\n";
        cout << "\t-e, --mem         INT       memory limit (MB default: [1024]MB), if rmdup was set, this memory limit has no effect\n";
        cout << "\t-T, --thread      INT       process thread number (default: [2])" << endl;
        cout << "\n";

        cout << "\t-Q, --qualSys     INT       quality system 1:illumina, 2:sanger  [1]\n";
        cout << "\t-G, --sanger                set clean data qualtiy system to sanger  [illumina]\n";
        cout << "\t-L, --read1Len    INT       read1 max length (default: all read1's length are equal, and auto acquire)\n";
        cout << "\t-I, --read2Len    INT       read2 max length (default: all read2's length are equal, and auto acquire)\n";
        cout << "\t    --minLen      INT       min length of read1&2, those reads longer than it will be filtered out\n";
        cout << "\n";

//      cout << "\t-a, --append      STR       the log's output place : console or file  [console]\n";
        cout << "\t-o, --outDir      STR       output directory, directory must exists  [.]\n";
        cout << "\t-C, --cleanFq1    STR       clean fq1 file name\n";
        cout << "\t-D, --cleanFq2    STR       clean fq2 file name\n";
        cout << "\t-R, --rawFq1      STR       raw fq1 file name\n";
        cout << "\t-W, --rawFq2      STR       raw fq2 file name\n";
        cout << "\n";
        cout << "\t-5, --seqType     INT       Sequence fq name type, 0->old fastq name, 1->new fastq name [0]\n";
        cout << "\t                                  old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
        cout << "\t                                  new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
        cout << "\t-6, --polyAType   INT       filter poly A type, 0->both two reads are poly a, 1->at least one reads is poly a, then filter  [0]\n";
        cout << "\t-7, --outType:    INT       Add /1, /2 at the end of fastq name, 0:not add, 1:add  [0]\n";
        cout << "\n";
        cout << "\t  --TtoU		               convert T to U when write data" << endl;
        cout << "\t  --UtoT		               convert U to T when read data" << endl;
        cout << "\n";
        cout << "\t-h, --help                  help" << endl;
        cout << "\t-v, --version               show version" << endl;
    }

    FilterProcessor::FilterProcessor() : PROCESS_THREAD_NUM(2), IS_STREAMING(false), filterTile_(false), tileIsFov_(false), misMatch_(1), matchRatio_(0.5), lowQual_(5),
                                         qualRate_(0.5), nRate_(0.05), polyA_(0), minMean_(0.0), filterIndex_(false),
                                         rmdup_(false), dupRateOnly_(false), cutReadNum_(0), headTrim_(0), tailTrim_(0), headTrim2_(0), tailTrim2_(0),
                                         memLimit_(700 * MEM_UNIT), qualSys_(ILLUMINA_), isFilterSmallInsertSize_(false), overlap_(10),
                                         mis_(0.1), readLen_(0), readLen2_(0), outDir_("."), onlyStat_(false), isPE_(true),minReadLength(50),cutAdaptor(false),cutBasesNumber(0),
                                         isAdptList_(true), isFull_(false), size_(0), cleanQualSys_(ILLUMINA_), filterAdapter_(true),seqType_(0),outType_(0),polyAType_(0)
    {
    }

    int FilterProcessor::processParams(int argc, char **argv)
    {

        trimBadHeadMaxLength = 0;
        trimBadHeadQuility = 0;
        trimBadTailMaxLength = 0;
        trimBadTailQuility = 0;
        maskLowQual = -1;
        TtoU = 0;
        UtoT = 0;
        cutAdaptorOri = 3;  //cut sequence from adaptor index

        const char *shortOptions = "f:r:1:2:K:M:A:l:T:q:n:m:p:d3in:N:t:e:c:SO:P:Q:L:I:Ga:o:C:D:R:W:5:6:7:Eb:x:y:z:hv";
        const struct option longOptions[] =
                {
                        { "adapter1", 1, NULL, 'f' },
                        { "adapter2", 1, NULL, 'r' },
                        { "fq1"     , 1, NULL, '1' },
                        { "fq2"     , 1, NULL, '2' },
                        { "tile"    , 1, NULL, 'K' },
                        { "fov"    , 1, NULL, 'F' },
                        { "misMatch", 1, NULL, 'M' },
                        { "matchRatio", 1, NULL, 'A' },
                        { "lowQual" , 1, NULL, 'l' },
                        { "qualRate", 1, NULL, 'q' },
                        { "nRate"   , 1, NULL, 'n' },
                        { "maskLowQual" , 1, NULL, 'N' },
                        { "mean"    , 1, NULL, 'm' },
                        { "polyA"   , 1, NULL, 'p' },
                        { "rmdup"   , 0, NULL, 'd' },
                        { "dupRate" , 0, NULL, '3' },
                        { "cutAdaptor" ,0, NULL, 'E'},
                        { "BaseNum" ,1, NULL, 'b'},
                        { "index"   , 0, NULL, 'i' },
                        { "cut"     , 1, NULL, 'c' },
                        { "trim"    , 1, NULL, 't' },
                        { "trimBadHead" , 1, NULL, 'x' },
                        { "trimBadTail" , 1, NULL, 'y' },

                        { "mem"     , 1, NULL, 'e' },
                        { "thread"  , 1, NULL, 'T' },
                        { "small"   , 0, NULL, 'S' },
                        { "overlap" , 1, NULL, 'O' },
                        { "mis"     , 1, NULL, 'P' },
                        { "qualSys" , 1, NULL, 'Q' },
                        { "read1Len", 1, NULL, 'L' },
                        { "read2Len", 1, NULL, 'I' },
                        { "minLen", 1, NULL, 'z' },
                        { "sanger"  , 0, NULL, 'G' },
                        { "append"  , 1, NULL, 'a' },
                        { "outDir"  , 1, NULL, 'o' },
                        { "cleanFq1", 1, NULL, 'C' },
                        { "cleanFq2", 1, NULL, 'D' },
                        { "rawFq1"  , 1, NULL, 'R' },
                        { "rawFq2"  , 1, NULL, 'W' },
                        { "seqType"  , 1, NULL, '5' },
                        { "polyAType" , 1, NULL, '6' },
                        { "outType" , 1, NULL, '7' },

                        { "TtoU" , 0, &TtoU, 1 },
                        { "UtoT" , 0, &UtoT, 1 },

                        { "help"    , 0, NULL, 'h' },
                        { "version" , 0, NULL, 'v' },
                        {     0,    0,    0,    0},
                };

        string append;

        if (argc == 1)
        {
            printUsage();
            return 1;
        }

        int nextOpt;
        string trim;
        string trimBadHeadStr;
        string trimBadTailStr;
        string cutAdaptorStr;
        size_t i;
        float num;
        string tiles;

        while (-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions, NULL)))
        {
            switch (nextOpt)
            {
                case 'f':
                    adapter1_.assign(optarg);
                    break;
                case 'r':
                    adapter2_.assign(optarg);
                    break;
                case '1':
                    fqFile1_.assign(optarg);
                    break;
                case '2':
                    fqFile2_.assign(optarg);
                    break;
                case 'K':
                    filterTile_ = true;
                    tiles.assign(optarg);
                    break;
                case 'F':
                    filterTile_ = true;
                    tileIsFov_ = true;
                    tiles.assign(optarg);
                    break;
                case 'M':
                    misMatch_ = atoi(optarg);
                    break;
                case 'A':
                    matchRatio_ = atof(optarg);
                    break;
                case 'l':
                    lowQual_ = atoi(optarg);
                    break;
                case 'q':
                    qualRate_ = atof(optarg);
                    break;
                case 'n':
                    nRate_ = atof(optarg);
                    break;
                case 'N':
                    maskLowQual = atoi(optarg);
                    break;
                case 'm':
                    minMean_ = atof(optarg);
                    break;
                case 'T':
                    PROCESS_THREAD_NUM = atoi(optarg);
                    break;
                case 'p':
                    polyA_ = atof(optarg);
                    break;
                case 'd':
                    rmdup_ = true;
                    break;
                case '3':
                    dupRateOnly_ = true;
                    break;
                case 'i':
                    filterIndex_ = true;
                    break;
                case 'c':
                    num = atof(optarg);
                    if (num < 1E-6)
                    {
                        cutReadNum_ = 0;
                    }
                    else
                    {
                        cutReadNum_ = (unsigned long)(num * 1024 * 1024);
                    }
                    break;
                case 't':
                    trim.assign(optarg);
                    i = trim.find_first_of(',');
                    if (i == string::npos)
                    {
                        cerr << "-t/--trim options error" << endl;
                        return 1;
                    }
                    headTrim_ = atoi(trim.substr(0, i).c_str());

                    trim = trim.substr(i+1);
                    i = trim.find_first_of(',');
                    if (i == string::npos)
                    {
                        cerr << "-t/--trim options error" << endl;
                        return 1;
                    }
                    tailTrim_ = atoi(trim.substr(0, i).c_str());

                    trim = trim.substr(i+1);
                    i = trim.find_first_of(',');
                    if (i == string::npos)
                    {
                        cerr << "-t/--trim options error" << endl;
                        return 1;
                    }
                    headTrim2_ = atoi(trim.substr(0, i).c_str());

                    trim = trim.substr(i+1);
                    tailTrim2_ = atoi(trim.c_str());
                    break;
                case 'e':
                    memLimit_ = atol(optarg) * MEM_UNIT;
                    break;
                case 'x':
                    trimBadHeadStr.assign(optarg);
                    i = trimBadHeadStr.find_first_of(',');
                    if (i == string::npos)
                    {
                        cout << "  --trimBadHead options error" << endl;
                        return 1;
                    }
                    trimBadHeadQuility  = atoi(trimBadHeadStr.substr(0, i).c_str());
                    trimBadHeadMaxLength = atoi(trimBadHeadStr.substr(i+1).c_str());
                    break;
                case 'y':
                    trimBadTailStr.assign(optarg);
                    i = trimBadTailStr.find_first_of(',');
                    if (i == string::npos)
                    {
                        cout << "  --trimBadTail options error" << endl;
                        return 1;
                    }
                    trimBadTailQuility = atoi(trimBadTailStr.substr(0, i).c_str());
                    trimBadTailMaxLength = atoi(trimBadTailStr.substr(i+1).c_str());
                    break;
                case 'S':
                    isFilterSmallInsertSize_ = true;
                    break;
                case 'O':
                    overlap_ = atoi(optarg);
                    break;
                case 'P':
                    mis_ = atof(optarg);
                    break;
                case 'Q':
                    switch (optarg[0])
                    {
                        case '1':
                            qualSys_ = ILLUMINA_;
                            break;
                        case '2':
                            qualSys_ = SANGER_;
                            break;
                        default:
                            cerr << "error quality system" << endl;
                            return 1;
                    }
                    break;
                case 'L':
                    readLen_ = atoi(optarg);
                    break;
                case 'I':
                    readLen2_ = atoi(optarg);
                    break;
                case 'G':
                    cleanQualSys_ = SANGER_;
                    break;
                case 'a':
                    append = optarg;
                    break;
                case 'o':
                    outDir_.assign(optarg);
                    if(outDir_[outDir_.length() - 1] == '\\')
                    {
                        outDir_ = outDir_.substr(0, outDir_.length() - 1);
                    }
                    break;
                case 'C':
                    cleanFq1_.assign(optarg);
                    break;
                case 'D':
                    cleanFq2_.assign(optarg);
                    break;
                case 'R':
                    rawFq1_.assign(optarg);
                    break;
                case 'W':
                    rawFq2_.assign(optarg);
                    break;
                case '5':
                    seqType_ = atoi(optarg);
                    break;
                case '6':
                    polyAType_ = atoi(optarg);
                    break;
                case '7':
                    outType_ = atoi(optarg);
                    break;
                case 's':
                    onlyStat_ = true;
                    break;
                case 'h':
                    printUsage();
                    return 1;
                case 'E':
                    cutAdaptor = true;
//                  cutAdaptorStr.assign(optarg);
//                  i = cutAdaptorStr.find_first_of(',');
//                  minReadLength = atoi(cutAdaptorStr.substr(0, i).c_str());
//                  if (i == string::npos)
//                      break;
//                  cutAdaptorOri = atoi(cutAdaptorStr.substr(i+1).c_str());
                    break;
                case 'z':
                    minReadLength = atoi(optarg);
                    break;
                case 'b':
                    cutBasesNumber = strtoul(optarg,NULL,10);
                    break;
                case 'v':
                    printVersion();
                    return 1;
                case 0:
                    break;
                case '?':
                    cerr << "unkonwn option: -" << (char) optopt << endl;
                    cerr << "Print -h or --help for more information." << endl;
                    return 1;
                default:
                    cerr << "Param : " << optarg << endl;
                    return 1;
            }
        }

        if(filterTile_)
        {
            if(tileIsFov_)
            {
                PreProcessTool::getFovs(tiles, tiles_);

            }else
            {
                PreProcessTool::getTiles(tiles, tiles_);
            }
        }

        if( cutAdaptor && cutBasesNumber != 0 ){
            cutReadNum_ = 0;
        }

        if(argc - optind == 1){
            streamingInput = strdup(argv[optind]);
            IS_STREAMING = true;
        }else if (argc != optind)
        {
            cerr << "options error, please check the options" << endl;
            return 1;
        }

        if(IS_STREAMING){
            if (!init_logger())
            {
                cerr << "Cannot Init Log."<< endl;
                return 1;
            }
            else
            {
                LOG(INFO, "Log Init Success");
            }

            if (rmdup_)
            {
                LOG(WARN, "-d,--rmdup parameter is unavailable in STREAMING mode!");
            }

            fqStreaming = strcmp(streamingInput, "-")? gzopen(streamingInput, "rb") : gzdopen(fileno(stdin), "r");

            char *buf = new char[1024];
            gzgets(fqStreaming, buf, 1024);
            if(strlen(buf) <= 10){
                LOG(ERROR, "Input is empty!");
                exit(1);
            }
            vector<string> s = split(buf,'\t');

            gzgets(fqStreaming, buf, 1024);
            vector<string> s2 = split(buf,'\t');

            if(s[2]!="2" && s2[2]!="2"){
                isPE_ = false;
            }

            if(isPE_){
                if(s[2]=="2"){
                    readLen_ = s2[3].size();
                    readLen2_ = s[3].size();
                }else{
                    readLen_ = s[3].size();
                    readLen2_ = s2[3].size();
                }
                LOG(INFO, "fq1 read length: " << readLen_);
                LOG(INFO, "fq2 read length: " << readLen2_);
            }else
            {
                readLen_ = s[3].size();
                LOG(INFO, "read length: " << readLen_);
            }
            delete[] buf;
            //gzclose(fqStreaming);
        }
        else{
            bool isPathNotExists = false; //指示输出目录是否存在
            if (access(outDir_.c_str(), F_OK) == -1)  //输出路径不存在
            {
                isPathNotExists = true;
                if (mkdir(outDir_.c_str(), 0755) != 0)
                {
                    cerr << "output directory " << outDir_ << " cannot create" << endl;
                    return 1;
                }
            }

            if (fqFile1_.empty())
            {
                cerr << "fq1 file must exists" << endl;
                return 1;
            }

            string logoutPath = outDir_ + "/" + LOG_FILE;
            if (!init_logger(append, logoutPath))
            {
                cerr << "Cannot Init Log:" << append << "-" << logoutPath << endl;
                return 1;
            }
            else
            {
                LOG(INFO, "Log Init Success");
            }

            if (isPathNotExists)
            {
                LOG(WARN, "output directory " << outDir_ << " does not exists, program will auto create");
                LOG(WARN, "output directory " << outDir_ << " has been created");
            }

            if (readLen_ == 0)
            {
                gzFile file = gzopen(fqFile1_.c_str(), "rb");
                if (file == NULL)
                {
                    LOG(ERROR, "No file: " << fqFile1_);
                    return 1;
                }
                char *buf = new char[512];
                gzgets(file, buf, 512);
                gzgets(file, buf, 512);
                readLen_ = strlen(buf) - 1;
                delete[] buf;
                gzclose(file);
            }

            LOG(INFO, "fq1 read length: " << readLen_);

            if (fqFile2_.empty())
            {
                isPE_ = false;
            }

            if (isPE_)
            {
                if (readLen2_ == 0)
                {
                    gzFile file = gzopen(fqFile2_.c_str(), "rb");
                    if (file == NULL)
                    {
                        LOG(ERROR, "No file: " << fqFile2_);
                        return 1;
                    }
                    char *buf = new char[512];
                    gzgets(file, buf, 512);
                    gzgets(file, buf, 512);
                    readLen2_ = strlen(buf) - 1;
                    delete[] buf;
                    gzclose(file);
                }

                LOG(INFO, "fq2 read length: " << readLen2_);
            }

            if( cutAdaptor && cutBasesNumber != 0 ){
                cutReadNum_ = 0;
            }

            //when rawFq1 or rawFq2 were set and not cut data
            if ((!rawFq1_.empty() || !rawFq2_.empty()) && cutReadNum_ == 0 && !filterTile_)
            {
                if (!rawFq1_.empty() && !fqFile1_.empty())
                {
                    string file = getOutputFileName(rawFq1_, "", outDir_);
                    if (fqFile1_.substr(fqFile1_.size() - 2, 2) != "gz")
                    {
                        gzLoad(fqFile1_.c_str(), file.c_str());
                    }
                    else
                    {
                        if ( CopyFile(fqFile1_.c_str(), file.c_str()) != 0 )
                        {
                            printf("Copy file error : cp %s %s \n", fqFile2_.c_str(), file.c_str());
                            return 2;
                        }
                    }
                }
                if (!rawFq2_.empty() && !fqFile2_.empty())
                {
                    string file = getOutputFileName(rawFq2_, "", outDir_);
                    if (fqFile2_.substr(fqFile2_.size()-2, 2) != "gz")
                    {
                        gzLoad(fqFile2_.c_str(), file.c_str());
                    }
                    else
                    {
                        if ( CopyFile(fqFile2_.c_str(), file.c_str()) != 0 )
                        {
                            printf("Copy file error : cp %s %s \n", fqFile2_.c_str(), file.c_str());
                            return 2;
                        }
                    }
                }
            }
        }

        if (readLen_ > MAX_LENGTH || readLen2_ > MAX_LENGTH)
        {
            LOG(ERROR, "the read's length must not larger than " << MAX_LENGTH);
            exit(1);
        }

        if (adapter1_.empty() && adapter2_.empty())
        {
            filterAdapter_ = false;
        }
        else
        {
            int type = adapterType(isPE_, adapter1_, adapter2_);
            if (type == 3)
            {
                cerr << "adapter1 or adapter2 type error" << endl;
                LOG(ERROR, "adapter1 or adapter2 type error");
                return 1;
            }
            else if (type == 2)
            {
                isAdptList_ = true;
                LOG(ERROR, "adapter1 or adapter2 type error.Not support adaptor list file.");
                return 2;
            }
            else
            {
                isAdptList_ = false;
                //convert all leters to capital
                for (int i=adapter1_.size(); i>=0; i--)
                {
                    adapter1_[i] = toupper(adapter1_[i]);
                }
                for (int i=adapter2_.size(); i>=0; i--)
                {
                    adapter2_[i] = toupper(adapter2_[i]);
                }
                adapterLen1_ = adapter1_.size();
                adapterLen2_ = adapter2_.size();
            }
        }

        return 0;
    }

    int FilterProcessor::filter(int argc, char **argv)
    {
        //process the command line params, and initial the member variables
        if (processParams(argc, argv) != 0)
        {
            return  1;
        }

        if(IS_STREAMING){
            FqInfo globleInfo1, globleInfo2;
            globleInfo1.rawReadLength = readLen_;
            globleInfo1.cleanReadLength = readLen_ - headTrim_ - tailTrim_;

            pool pl(PROCESS_THREAD_NUM + 1);
            if (isPE_)
            {
                globleInfo2.rawReadLength = readLen2_;
                globleInfo2.cleanReadLength = readLen2_ - headTrim2_ - tailTrim2_;

                Read *reads1, *reads2;
                long capacity = static_cast<long>(memLimit_ / 2.5);

                PeBuffer buffer(fqStreaming, capacity, filterTile_, tiles_);
                buffer.setSeqType(seqType_);
                buffer.setTileIsFov(tileIsFov_);
                TaskParam *params = new TaskParam[PROCESS_THREAD_NUM];

                int size = 0;

                bool first = true;

                while (buffer.getReads() > 0)
                {
                    size = buffer.getReadSize();
                    if(size == 0)
                        continue;
                    if (first)
                    {
                        first = false;
                        cleanDataIndexs_ = new unsigned int[buffer.getInitReadSize() + 1];
                    }

                    reads1 = buffer.getReadsOne();
                    reads2 = buffer.getReadsTwo();

                    size_ = 0;
                    doneNum_ = 0;

                    int block = size / PROCESS_THREAD_NUM;
                    int remain = size % PROCESS_THREAD_NUM;
                    int index = 0;

                    pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTaskStreamingPE, this, reads1, reads2));

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
                    {
                        params[i].left = index;
                        index += block;
                        if (remain > 0)
                        {
                            index += 1;
                            remain--;
                        }

                        params[i].right = index;
                        params[i].reads1 = reads1;
                        params[i].reads2 = reads2;

                        bzero(&(params[i].info1), sizeof(FqInfo));
                        bzero(&(params[i].info2), sizeof(FqInfo));

                        pl.schedule(boost::bind(&FilterProcessor::task, this, &params[i]));
                    }

                    pl.wait();

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
                    {
                        globleInfo1.add(params[i].info1);
                        globleInfo2.add(params[i].info2);
                        //calculte the max quality value
                        if (params[i].info1.maxQualityValue > globleInfo1.maxQualityValue)
                        {
                            globleInfo1.maxQualityValue = params[i].info1.maxQualityValue;
                        }
                        if (params[i].info2.maxQualityValue > globleInfo2.maxQualityValue)
                        {
                            globleInfo2.maxQualityValue = params[i].info2.maxQualityValue;
                        }
                    }

                    if (cutReadNum_ > 0 && cutReadNum_  < globleInfo1.cleanTotalReadNum)
                    {
                        break;
                    }

                    if(cutAdaptor && cutBasesNumber > 0 && cutBasesNumber < globleInfo1.cleanTotalBaseNum){
                        break;
                    }

                }

                if (globleInfo1.maxQualityValue > globleInfo2.maxQualityValue)
                {
                    globleInfo2.maxQualityValue = globleInfo1.maxQualityValue;
                }
                else
                {
                    globleInfo1.maxQualityValue = globleInfo2.maxQualityValue;
                }

                if (size == -1)
                {
                    LOG(ERROR, "read fq1 or fq2 error");
                    return 1;
                }

                if (size == -2)
                {
                    LOG(WARN, "fq1 and fq2's read number not equal");
                }

                if (params != NULL)
                {
                    LOG(INFO, "delete params");
                    delete []params;
                }
                printFqInfo(&globleInfo1, &globleInfo2);

            }else{
                long capacity = static_cast<long>(memLimit_ / 2.5);
                FqBuffer buffer(fqStreaming, capacity, filterTile_, tiles_);
                buffer.setSeqType(seqType_);
                buffer.setTileIsFov(tileIsFov_);

                TaskParam *params = new TaskParam[PROCESS_THREAD_NUM];

                Read *reads;
                int size;

                bool first = true;

                while ((reads = buffer.getReads()))
                {
                    size = buffer.getRealReadSize();
                    if(size==0)
                        continue;

                    if (size == -1)
                    {
                        LOG(ERROR, "read fq file: " + fqFile1_ + " error");
                        return 1;
                    }

                    if (first)
                    {
                        first = false;
                        cleanDataIndexs_ = new unsigned int[buffer.getReadSize() + 1];
                    }

                    size_ = 0;
                    doneNum_ = 0;

                    int block = size / PROCESS_THREAD_NUM;
                    int remain = size % PROCESS_THREAD_NUM;
                    int index = 0;

                    pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTaskStreamingSE, this, reads));

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
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

                        bzero(&(params[i].info1), sizeof(FqInfo));

                        pl.schedule(boost::bind(&FilterProcessor::task, this, &params[i]));
                    }

                    pl.wait();

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
                    {
                        globleInfo1.add(params[i].info1);
                        if (globleInfo1.maxQualityValue < params[i].info1.maxQualityValue)
                        {
                            globleInfo1.maxQualityValue = params[i].info1.maxQualityValue;
                        }
                    }


                    if (cutReadNum_ > 0 && cutReadNum_ < globleInfo1.cleanTotalReadNum)
                    {
                        break;
                    }

                    if(cutAdaptor && cutBasesNumber > 0 && cutBasesNumber < globleInfo1.cleanTotalBaseNum){
                        break;
                    }

                }


                if (params != NULL)
                {
                    LOG(INFO, "delete params");
                    delete []params;
                }

                printFqInfo(&globleInfo1, NULL);
            }
        }else{
            FqInfo globleInfo1, globleInfo2;
            //store clean read
            gzFile outCleanFile1 = NULL, outCleanFile2 = NULL;
            string outCleanFileName1, outCleanFileName2;
            gzFile outRawFile1 = NULL, outRawFile2 = NULL;
            string outRawFileName1, outRawFileName2;

            globleInfo1.rawReadLength = readLen_;
            globleInfo1.cleanReadLength = readLen_ - headTrim_ - tailTrim_;
            //get the output fq files' name
            if (!cleanFq1_.empty())  //有指定输出文件名
            {
                outCleanFileName1 = getOutputFileName(cleanFq1_, "", outDir_);
            }
            else
            {
                outCleanFileName1 = getOutputFileName(fqFile1_, CLEAN_FQ_PREFIX, outDir_);
            }

            outCleanFile1 = gzopen(outCleanFileName1.c_str(), "wb");

            if (outCleanFile1 == NULL)
            {
                LOG(ERROR, "create output file: " + outCleanFileName1 + " error");
                return 1;
            }

            if (isPE_)
            {
                globleInfo2.rawReadLength = readLen2_;
                globleInfo2.cleanReadLength = readLen2_ - headTrim2_ - tailTrim2_;

                if (!cleanFq2_.empty())
                {
                    outCleanFileName2 = getOutputFileName(cleanFq2_, "", outDir_);
                }
                else
                {
                    outCleanFileName2 = getOutputFileName(fqFile2_, CLEAN_FQ_PREFIX, outDir_);
                }

                outCleanFile2 = gzopen(outCleanFileName2.c_str(), "wb");
                if (!outCleanFile2)
                {
                    LOG(ERROR, "create output file: " + outCleanFileName2 + " error");
                    return 1;
                }
            }

            if (cutReadNum_ || filterTile_)
            {
                if (!rawFq1_.empty())
                {
                    outRawFileName1 = getOutputFileName(rawFq1_, "", outDir_);
                }
                else
                {
                    outRawFileName1 = getOutputFileName(fqFile1_, RAW_FQ_PREFIX, outDir_);
                }

                outRawFile1 = gzopen(outRawFileName1.c_str(), "wb");
                if (!outRawFile1)
                {
                    LOG(ERROR, "create output file: " + outRawFileName1);
                    return 1;
                }
                if (isPE_)
                {
                    if (!rawFq2_.empty())
                    {
                        outRawFileName2 = getOutputFileName(rawFq2_, "", outDir_);
                    }
                    else
                    {
                        outRawFileName2 = getOutputFileName(fqFile2_, RAW_FQ_PREFIX, outDir_);
                    }

                    outRawFile2 = gzopen(outRawFileName2.c_str(), "wb");
                    if (!outRawFile2)
                    {
                        LOG(ERROR, "create output file: " + outRawFileName2);
                        return 1;
                    }
                }
            }

            if (filterAdapter_ && isAdptList_)
            {
                if (getReadsNameFromFile(adapter1_, readsName1_) != 0)
                {
                    return 1;
                }
                if (isPE_)
                {
                    if (getReadsNameFromFile(adapter2_, readsName2_) != 0)
                    {
                        return 1;
                    }
                }
            }

            int threadNum;
            if (isPE_)
            {
                threadNum = PROCESS_THREAD_NUM + 2;
            }
            else
            {
                threadNum = PROCESS_THREAD_NUM + 1;
            }
            //
            pool pl(threadNum);
            if (isPE_)
            {
                Read *reads1, *reads2;

                long capacity = static_cast<long>(memLimit_ / 2.5);

                PeBuffer buffer(fqFile1_.c_str(), fqFile2_.c_str(), capacity, PeBuffer::RB, filterTile_, tiles_);
                buffer.setSeqType(seqType_);
                buffer.setTileIsFov(tileIsFov_);
                TaskParam *params = new TaskParam[PROCESS_THREAD_NUM];

                int size = 0;

                bool first = true;
                //循环一次，产生一个已排序的临时文件
                int fileNum = 0;

                while (buffer.getReads() > 0)
                {
                    size = buffer.getReadSize();
                    if(size == 0)
                        continue;
                    if (first)
                    {
                        first = false;
                        cleanDataIndexs_ = new unsigned int[buffer.getInitReadSize() + 1];
                    }

                    reads1 = buffer.getReadsOne();
                    reads2 = buffer.getReadsTwo();

                    if (cutReadNum_ || filterTile_) //output raw data
                    {
                        pl.schedule(boost::bind(&FilterProcessor::outputRawDataTask,
                                                this, outRawFile1, reads1, size)
                        );
                        pl.schedule(boost::bind(&FilterProcessor::outputRawDataTask,
                                                this, outRawFile2, reads2, size)
                        );

                        pl.wait();
                    }
                    size_ = 0;
                    doneNum_ = 0;

                    int block = size / PROCESS_THREAD_NUM;
                    int remain = size % PROCESS_THREAD_NUM;
                    int index = 0;
                    ofstream tempOFS;

                    if (!rmdup_)
                    {
                        if(outType_!=0){
                            pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTask2, this, outCleanFile1, reads1, "/1\n"));
                            pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTask2, this, outCleanFile2, reads2, "/2\n"));
                        }else{
                            pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTask, this, outCleanFile1, reads1));
                            pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTask, this, outCleanFile2, reads2));
                        }
                    }

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
                    {
                        params[i].left = index;
                        index += block;
                        if (remain > 0)
                        {
                            index += 1;
                            remain--;
                        }

                        params[i].right = index;
                        params[i].reads1 = reads1;
                        params[i].reads2 = reads2;

                        bzero(&(params[i].info1), sizeof(FqInfo));
                        bzero(&(params[i].info2), sizeof(FqInfo));

                        pl.schedule(boost::bind(&FilterProcessor::task, this, &params[i]));
                    }

                    pl.wait();

                    if (rmdup_)
                    {
                        char tempFile[1024];
                        sprintf(tempFile, "%s/%d.sort.temp", outDir_.c_str(), fileNum);
                        tempOFS.open(tempFile);
                        outputTempData(tempOFS, reads1, reads2);
                        if(dupRateOnly_){
                            ofstream tempOFS1;
                            char dupTempFile[1024];
                            sprintf(dupTempFile, "%s/%d.dup.temp", outDir_.c_str(), fileNum);
                            tempOFS1.open(dupTempFile);
                            outputDupData(tempOFS1, reads1, reads2,size);
                        }
                        duplications_.clear();
                    }

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
                    {
                        globleInfo1.add(params[i].info1);
                        globleInfo2.add(params[i].info2);
                        //calculte the max quality value
                        if (params[i].info1.maxQualityValue > globleInfo1.maxQualityValue)
                        {
                            globleInfo1.maxQualityValue = params[i].info1.maxQualityValue;
                        }
                        if (params[i].info2.maxQualityValue > globleInfo2.maxQualityValue)
                        {
                            globleInfo2.maxQualityValue = params[i].info2.maxQualityValue;
                        }
                    }

                    if (rmdup_ && cutReadNum_ > 0 && cutReadNum_ * 1.1 < globleInfo1.greyTotalReadNum)
                    {
                        if (rmdup_)
                        {
                            tempOFS.close();
                            fileNum++;
                        }
                        break;
                    }
                    else if (cutReadNum_ > 0 && cutReadNum_  < globleInfo1.cleanTotalReadNum)
                    {
                        break;
                    }

                    if(cutAdaptor && cutBasesNumber > 0 && cutBasesNumber < globleInfo1.cleanTotalBaseNum){
                        break;
                    }

                    if (rmdup_)
                    {
                        fileNum++;
                        tempOFS.close();
                    }
                }

                //merge all sorted files
                if (rmdup_)
                {
                    mergeSortedFiles(fileNum, &globleInfo1, &globleInfo2, outCleanFile1, outCleanFile2);
                    for (int i=0; i<fileNum; ++i)
                    {
                        char tempFile[1024];
                        sprintf(tempFile, "%s/%d.sort.temp", outDir_.c_str(), i);
                        remove(tempFile);
                        if(dupRateOnly_){
                            char dupTempFile[1024];
                            sprintf(dupTempFile, "%s/%d.dup.temp", outDir_.c_str(), i);
                            remove(dupTempFile);
                        }
                    }
                }

                if (globleInfo1.maxQualityValue > globleInfo2.maxQualityValue)
                {
                    globleInfo2.maxQualityValue = globleInfo1.maxQualityValue;
                }
                else
                {
                    globleInfo1.maxQualityValue = globleInfo2.maxQualityValue;
                }

                if (size == -1)
                {
                    LOG(ERROR, "read fq1 or fq2 error");
                    return 1;
                }

                if (size == -2)
                {
                    LOG(WARN, "fq1 and fq2's read number not equal");
                }

                if (params != NULL)
                {
                    LOG(INFO, "delete params");
                    delete []params;
                }

                if (cutReadNum_ || filterTile_)
                {
                    gzclose(outRawFile1);
                    gzclose(outRawFile2);
                }

                gzclose(outCleanFile1);
                gzclose(outCleanFile2);


            }
            else //SE
            {
                long capacity = static_cast<long>(memLimit_ / 2.5);
                FqBuffer buffer(fqFile1_.c_str(), capacity, FqBuffer::RB, filterTile_, tiles_);
                buffer.setSeqType(seqType_);
                buffer.setTileIsFov(tileIsFov_);

                TaskParam *params = new TaskParam[PROCESS_THREAD_NUM];

                Read *reads;
                int size;

                bool first = true;
                int fileNum = 0;

                while ((reads = buffer.getReads()))
                {
                    size = buffer.getRealReadSize();
                    if(size==0)
                        break;
                    //	continue;

                    if (size == -1)
                    {
                        LOG(ERROR, "read fq file: " + fqFile1_ + " error");
                        return 1;
                    }

                    if (first)
                    {
                        first = false;
                        cleanDataIndexs_ = new unsigned int[buffer.getReadSize() + 1];
                    }

                    size_ = 0;
                    doneNum_ = 0;

                    if (cutReadNum_ || filterTile_) //output raw data
                    {
                        outputRawDataTask(outRawFile1, reads, size);
                    }

                    int block = size / PROCESS_THREAD_NUM;
                    int remain = size % PROCESS_THREAD_NUM;
                    int index = 0;
                    ofstream tempOFS;

                    if (!rmdup_)
                    {
                        if(outType_!=0){
                            pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTask2, this, outCleanFile1, reads,"/1\n"));

                        }else{
                            pl.schedule(boost::bind(&FilterProcessor::outputCleanDataTask, this, outCleanFile1, reads));
                        }
                    }

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
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

                        bzero(&(params[i].info1), sizeof(FqInfo));

                        pl.schedule(boost::bind(&FilterProcessor::task, this, &params[i]));
                    }

                    pl.wait();

                    if (rmdup_)
                    {
                        char tempFile[1024];
                        sprintf(tempFile, "%s/%d.sort.temp", outDir_.c_str(), fileNum);
                        tempOFS.open(tempFile);
                        outputTempData(tempOFS, reads);
                        if(dupRateOnly_){
                            ofstream tempOFS1;
                            char dupTempFile[1024];
                            sprintf(dupTempFile, "%s/%d.dup.temp", outDir_.c_str(), fileNum);
                            tempOFS1.open(dupTempFile);
                            outputDupData(tempOFS1, reads, size);
                        }
                        duplications_.clear();
                    }

                    for (unsigned int i=0; i<PROCESS_THREAD_NUM; ++i)
                    {
                        globleInfo1.add(params[i].info1);
                        if (globleInfo1.maxQualityValue < params[i].info1.maxQualityValue)
                        {
                            globleInfo1.maxQualityValue = params[i].info1.maxQualityValue;
                        }
                    }

                    if (rmdup_ && cutReadNum_ > 0 && cutReadNum_*1.1 < globleInfo1.greyTotalReadNum)
                    {
                        if (rmdup_)
                        {
                            fileNum++;
                            tempOFS.close();
                        }
                        break;
                    }
                    else if (cutReadNum_ > 0 && cutReadNum_ < globleInfo1.cleanTotalReadNum)
                    {

                        break;
                    }

                    if(cutAdaptor && cutBasesNumber > 0 && cutBasesNumber < globleInfo1.cleanTotalBaseNum){
                        break;
                    }

                    if (rmdup_)
                    {
                        fileNum++;
                        tempOFS.close();
                    }
                }

                if (rmdup_)
                {
                    mergeSortedFiles(fileNum, &globleInfo1, outCleanFile1);
                    for (int i=0; i<fileNum; ++i)
                    {
                        char tempFile[1024];
                        sprintf(tempFile, "%s/%d.sort.temp", outDir_.c_str(), i);
                        remove(tempFile);
                        if(dupRateOnly_){
                            char dupTempFile[1024];
                            sprintf(dupTempFile, "%s/%d.dup.temp", outDir_.c_str(), i);
                            remove(dupTempFile);
                        }
                    }
                }

                if (params != NULL)
                {
                    LOG(INFO, "delete params");
                    delete []params;
                }

                if (cutReadNum_ || filterTile_)
                {
                    gzclose(outRawFile1);
                }
                gzclose(outCleanFile1);

                /*gzclose(filteredFile1);*/
            }

            //////////////////////////////////////////////////

            //output the statistic infomation
            if (!isPE_)
            {
                printFqInfo(outDir_, lanID_, &globleInfo1, NULL);
            }
            else
            {
                printFqInfo(outDir_, lanID_, &globleInfo1, &globleInfo2);
            }
        }
        LOG(INFO, "PreProcess Finish");

        return 0;
    }

    void FilterProcessor::outputRawDataTask(gzFile &file,Read *reads, unsigned int size)
    {
        for (unsigned int i=0; i<size; ++i)
        {
            gzputs(file, reads[i].readName);
            gzputs(file, "\n");
            gzputs(file, reads[i].baseSequence);
            gzputs(file, "\n");
            gzputs(file, reads[i].optionalName);
            gzputs(file, "\n");
            gzputs(file, reads[i].baseQuality);
            gzputs(file, "\n");
        }
    }

    void FilterProcessor::mergeSortedFiles(int num, FqInfo *info1, FqInfo *info2, gzFile &outFile1, gzFile &outFile2)
    {
        FqFile ifs[num];
        StrRead *reads1 = new StrRead[num];
        StrRead *reads2 = new StrRead[num];
        bool result = false;
        char tempFile[1024];
        map<string, int> heap;
        map<string, int>::iterator iter;
        pair<map<string,int>::iterator, bool> ret;
        string temp;
        string previous("");

        //initial
        for (int i=0; i<num; ++i)
        {
            sprintf(tempFile, "%s/%d.sort.temp", outDir_.c_str(), i);
            ifs[i].open(tempFile);
            result = ifs[i].nextRead(reads1[i], reads2[i]);
            if (!result)
            {
                continue;
            }
            temp = reads1[i].baseSequence + reads2[i].baseSequence;
            MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
            result = true;
            while (result)
            {
                ret = heap.insert(pair<string, int>(string((const char*)md5Seq_, MD5_DIGEST_LENGTH), i));
                if (ret.second)
                {
                    break;
                }
                else
                {
                    info1->duplicationNum++;
                    info2->duplicationNum++;
                    info1->totalDuplicationNum++;
                    if(dupRateOnly_){
                        outputCleanData(outFile1, reads1[i]);
                        outputCleanData(outFile2, reads2[i]);
                    }
                    result = ifs[i].nextRead(reads1[i], reads2[i]);
                    temp = reads1[i].baseSequence + reads2[i].baseSequence;
                    MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
                }
            }
        }

        int index;
        while (true)
        {
            if (heap.empty())
            {
                break;
            }

            iter = heap.begin();
            index = iter->second;
            if (iter->first != previous)
            {
                previous = iter->first;

                StatisInfo si1 = auxStatistics(&reads1[index]);
                StatisInfo si2 = auxStatistics(&reads2[index]);

                info1->cleanBaseA += si1.a;
                info1->cleanBaseC += si1.c;
                info1->cleanBaseG += si1.g;
                info1->cleanBaseT += si1.t;
                info1->cleanBaseN += si1.n;
                info1->cleanQ20 += si1.q20;
                info1->cleanQ30 += si1.q30;
                info1->cleanTotalReadNum++;
                info1->cleanTotalBaseNum += reads1[index].baseSequence.length();

                calculateBaseDistribute(&reads1[index], *info1, reads1[index].baseSequence.length());

                info2->cleanBaseA += si2.a;
                info2->cleanBaseC += si2.c;
                info2->cleanBaseG += si2.g;
                info2->cleanBaseT += si2.t;
                info2->cleanBaseN += si2.n;
                info2->cleanQ20 += si2.q20;
                info2->cleanQ30 += si2.q30;
                info2->cleanTotalReadNum++;
                info2->cleanTotalBaseNum += reads2[index].baseSequence.length();
                calculateBaseDistribute(&reads2[index], *info2, reads2[index].baseSequence.length());

                outputCleanData(outFile1, reads1[index]);
                outputCleanData(outFile2, reads2[index]);
            }
            else //duplication
            {
                info1->duplicationNum++;
                info2->duplicationNum++;
                info1->totalDuplicationNum++;

                if(dupRateOnly_){
                    outputCleanData(outFile1, reads1[index]);
                    outputCleanData(outFile2, reads2[index]);
                }
            }

            heap.erase(iter);
            result = ifs[index].nextRead(reads1[index], reads2[index]);
            while (result)
            {
                temp = reads1[index].baseSequence + reads2[index].baseSequence;
                MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
                ret = heap.insert(pair<string, int>(string((const char*)md5Seq_, MD5_DIGEST_LENGTH), index));
                if (!ret.second)
                {
                    info1->duplicationNum++;
                    info2->duplicationNum++;
                    info1->totalDuplicationNum++;

                    if(dupRateOnly_){
                        outputCleanData(outFile1, reads1[index]);
                        outputCleanData(outFile2, reads2[index]);
                    }
                    result = ifs[index].nextRead(reads1[index], reads2[index]);
                }
                else
                {
                    break;
                }
            }
        }

        if(dupRateOnly_){
            FqFile dup[num];
            for (int i=0; i<num; ++i)
            {
                char dupFile[1024];
                sprintf(dupFile, "%s/%d.dup.temp", outDir_.c_str(), i);
                dup[i].open(dupFile);
                while(true){
                    result = dup[i].nextRead(reads1[i], reads2[i]);
                    if (!result)
                    {
                        break;
                    }

                    outputCleanData(outFile1, reads1[i]);
                    outputCleanData(outFile2, reads2[i]);
                }
            }
        }

        delete[] reads1;
        delete[] reads2;
    }

    void FilterProcessor::mergeSortedFiles(int num, FqInfo *info, gzFile &outFile)
    {
        FqFile ifs[num];
        StrRead *reads = new StrRead[num];

        bool result = false;
        char tempFile[1024];
        map<string, int> heap;
        map<string, int>::iterator iter;
        pair<map<string,int>::iterator, bool> ret;
        string temp;
        string previous("");

        //initial
        for (int i=0; i<num; ++i)
        {
            sprintf(tempFile, "%s/%d.sort.temp", outDir_.c_str(), i);
            ifs[i].open(tempFile);
            result = ifs[i].nextRead(reads[i]);
            if (!result)
            {
                continue;
            }
            temp = reads[i].baseSequence;
            MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
            result = true;
            while (result)
            {
                ret = heap.insert(pair<string, int>(string((const char*)md5Seq_, MD5_DIGEST_LENGTH), i));
                if (ret.second)
                {
                    break;
                }
                else
                {
                    info->duplicationNum++;
                    if(dupRateOnly_){
                        outputCleanData(outFile, reads[i]);
                    }
                    result = ifs[i].nextRead(reads[i]);
                    temp = reads[i].baseSequence;
                    MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
                }
            }
        }

        int index;
        while (true)
        {
            if (heap.empty())
            {
                break;
            }
            iter = heap.begin();
            index = iter->second;
            if (iter->first != previous)
            {
                previous = iter->first;

                StatisInfo si = auxStatistics(&reads[index]);

                info->cleanBaseA += si.a;
                info->cleanBaseC += si.c;
                info->cleanBaseG += si.g;
                info->cleanBaseT += si.t;
                info->cleanBaseN += si.n;
                info->cleanQ20 += si.q20;
                info->cleanQ30 += si.q30;

                info->cleanTotalBaseNum += reads[index].baseSequence.length();
                info->cleanTotalReadNum++;

                calculateBaseDistribute(&reads[index], *info, reads[index].baseSequence.length());

                outputCleanData(outFile, reads[index]);
            }
            else
            {
                info->duplicationNum++;

                if(dupRateOnly_){
                    outputCleanData(outFile, reads[index]);
                }
            }
            heap.erase(iter);
            result = ifs[index].nextRead(reads[index]);
            while (result)
            {
                temp = reads[index].baseSequence;
                MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
                ret = heap.insert(pair<string, int>(string((const char*)md5Seq_, MD5_DIGEST_LENGTH), index));
                if (!ret.second)
                {
                    info->duplicationNum++;
                    if(dupRateOnly_){
                        outputCleanData(outFile, reads[index]);
                    }
                    result = ifs[index].nextRead(reads[index]);
                }
                else
                {
                    break;
                }
            }

        }

        if(dupRateOnly_){
            FqFile dup[num];
            for (int i=0; i<num; ++i)
            {
                char dupFile[1024];
                sprintf(dupFile, "%s/%d.dup.temp", outDir_.c_str(), i);
                dup[i].open(dupFile);
                while(true){
                    result = dup[i].nextRead(reads[i]);
                    if (!result)
                    {
                        break;
                    }
                    outputCleanData(outFile, reads[i]);
                }
            }
        }

        delete[] reads;
    }

    void FilterProcessor::outputCleanData(gzFile &file, const StrRead &read)
    {
        gzputs(file, read.readName.c_str());
        gzputs(file, "\n");
        gzputs(file, read.baseSequence.c_str());
        gzputs(file, "\n");
        gzputs(file, read.optionalName.c_str());
        gzputs(file, "\n");
        gzputs(file, read.baseQuality.c_str());
        gzputs(file, "\n");
    }

    void FilterProcessor::outputTempData(ofstream& file, Read* reads1, Read* reads2)
    {
        map<string ,int>::iterator iter;
        int index;
        for (iter = duplications_.begin(); iter != duplications_.end(); ++iter)
        {
            index = iter->second;
            file << reads1[index].readName << '\t'
                 << reads2[index].readName << '\n'
                 << reads1[index].baseSequence << '\t'
                 << reads2[index].baseSequence << '\n'
                 << reads1[index].optionalName << '\t'
                 << reads2[index].optionalName << '\n'
                 << reads1[index].baseQuality << '\t'
                 << reads2[index].baseQuality << '\n';
        }
    }

    void FilterProcessor::outputTempData(ofstream& file, Read* reads)
    {
        map<string ,int>::iterator iter;
        int index;
        for (iter = duplications_.begin(); iter != duplications_.end(); ++iter)
        {
            index = iter->second;
            file << reads[index].readName << '\n'
                 << reads[index].baseSequence << '\n'
                 << reads[index].optionalName << '\n'
                 << reads[index].baseQuality << '\n';
        }
    }

    void FilterProcessor::outputDupData(ofstream& file, Read* reads1, Read* reads2, unsigned int size)
    {
        bool* outIndex = new bool[size];
        map<string ,int>::iterator iter;
        unsigned int i;

        for (i=0; i < size; i++){
            outIndex[i] = false;
        }

        for (iter = duplications_.begin(); iter != duplications_.end(); ++iter)
        {
            i = iter->second;
            outIndex[i] = true;
        }

        for(i=0; i < atomic_read32(&size_); i++)
        {
            int index = cleanDataIndexs_[i];
            if(!outIndex[index])
            {
                file << reads1[index].readName << '\t'
                     << reads2[index].readName << '\n'
                     << reads1[index].baseSequence << '\t'
                     << reads2[index].baseSequence << '\n'
                     << reads1[index].optionalName << '\t'
                     << reads2[index].optionalName << '\n'
                     << reads1[index].baseQuality << '\t'
                     << reads2[index].baseQuality << '\n';
            }
        }
        delete[] outIndex;
    }

    void FilterProcessor::outputDupData(ofstream& file, Read* reads, unsigned int size)
    {
        bool* outIndex = new bool[size];
        map<string ,int>::iterator iter;
        unsigned int i;
        for(i=0; i < size; i++){
            outIndex[i] = false;
        }
        for (iter = duplications_.begin(); iter != duplications_.end(); ++iter)
        {
            i = iter->second;
            outIndex[i] = true;
        }

        for(i = 0; i < atomic_read32(&size_); i++)
        {
            int index = cleanDataIndexs_[i];
            if(!outIndex[index])
            {
                file << reads[index].readName << '\n'
                     << reads[index].baseSequence << '\n'
                     << reads[index].optionalName << '\n'
                     << reads[index].baseQuality << '\n';
            }
        }
        delete[] outIndex;
    }

    void FilterProcessor::outputCleanDataTaskStreamingPE(Read *reads, Read *reads2)
    {
        unsigned int index = 0;
        while (index == atomic_read32(&size_) && (atomic_read32(&doneNum_) < PROCESS_THREAD_NUM))
        {
            this_thread::sleep(posix_time::seconds(2)); //sleep 2 second
        }

        while ((atomic_read32(&doneNum_) < PROCESS_THREAD_NUM) || index < atomic_read32(&size_))
        {
            this_thread::sleep(posix_time::seconds(1));
            while (index < atomic_read32(&size_))
            {
                char *rName = reads[cleanDataIndexs_[index]].readName;
                cout << "@" << rName << "/1" << endl;
                cout << reads[cleanDataIndexs_[index]].baseSequence << endl;
                cout << "+" << endl;
                cout << reads[cleanDataIndexs_[index]].baseQuality << endl;

                cout << "@" << rName << "/2" << endl;
                cout << reads2[cleanDataIndexs_[index]].baseSequence << endl;
                cout << "+" << endl;
                cout << reads2[cleanDataIndexs_[index]].baseQuality << endl;

                index++;
            }
        }
    }

    void FilterProcessor::outputCleanDataTaskStreamingSE(Read *reads)
    {
        unsigned int index = 0;
        while (index == atomic_read32(&size_) && (atomic_read32(&doneNum_) < PROCESS_THREAD_NUM))
        {
            this_thread::sleep(posix_time::seconds(2)); //sleep 2 second
        }

        while ((atomic_read32(&doneNum_) < PROCESS_THREAD_NUM) || index < atomic_read32(&size_))
        {
            this_thread::sleep(posix_time::seconds(1));
            while (index < atomic_read32(&size_))
            {
                cout << "@" << reads[cleanDataIndexs_[index]].readName << "/1" << endl;
                cout << reads[cleanDataIndexs_[index]].baseSequence << endl;
                cout << "+" << endl;
                cout << reads[cleanDataIndexs_[index]].baseQuality << endl;
                index++;
            }
        }
    }

    void FilterProcessor::outputCleanDataTask(gzFile &file, Read *reads)
    {
        unsigned int index = 0;
        while (index == atomic_read32(&size_) && (atomic_read32(&doneNum_) < PROCESS_THREAD_NUM))
        {
            this_thread::sleep(posix_time::seconds(2)); //sleep 2 second
        }

        while ((atomic_read32(&doneNum_) < PROCESS_THREAD_NUM) || index < atomic_read32(&size_))
        {
            this_thread::sleep(posix_time::seconds(1));
            while (index < atomic_read32(&size_))
            {
                gzputs(file, reads[cleanDataIndexs_[index]].readName);
                gzputs(file, "\n");
                gzputs(file, reads[cleanDataIndexs_[index]].baseSequence);
                gzputs(file, "\n");
                gzputs(file, reads[cleanDataIndexs_[index]].optionalName);
                gzputs(file, "\n");
                gzputs(file, reads[cleanDataIndexs_[index]].baseQuality);
                gzputs(file, "\n");

                index++;
            }
        }
    }
    void FilterProcessor::outputCleanDataTask2(gzFile &file, Read *reads, const char* name)
    {
        unsigned int index = 0;
        while (index == atomic_read32(&size_) && (atomic_read32(&doneNum_) < PROCESS_THREAD_NUM))
        {
            this_thread::sleep(posix_time::seconds(2)); //sleep 2 second
        }

        while ((atomic_read32(&doneNum_) < PROCESS_THREAD_NUM) || index < atomic_read32(&size_))
        {
            this_thread::sleep(posix_time::seconds(1));
            while (index < atomic_read32(&size_))
            {
                char *rName = reads[cleanDataIndexs_[index]].readName;
                int ii = strlen(rName);
                while(ii>0 && rName[ii]!=' ')
                    ii --;
                if(ii > 0)
                    rName[ii]='\0';
                gzputs(file, rName);
                gzputs(file, name);
                gzputs(file, reads[cleanDataIndexs_[index]].baseSequence);
                gzputs(file, "\n");
                gzputs(file, reads[cleanDataIndexs_[index]].optionalName);
                gzputs(file, "\n");
                gzputs(file, reads[cleanDataIndexs_[index]].baseQuality);
                gzputs(file, "\n");

                index++;
            }
        }
    }

    void FilterProcessor::task(TaskParam *param)
    {
        Read *reads1 = param->reads1;
        Read *reads2 = param->reads2;

        int start = param->left;
        int end = param->right;

        bool isClean;

        if (start == end)
        {
            return;
        }

        if (reads2 != NULL) //PE
        {
            for (int i=start; i<end; ++i)
            {
                upper(reads1[i].baseSequence);
                upper(reads2[i].baseSequence);
                if(UtoT){
                    turnBase(reads1[i].baseSequence, 'U', 'T');
                    turnBase(reads2[i].baseSequence, 'U', 'T');
                }
                isClean = statisticsPE(reads1, reads2, i, &(param->info1), &(param->info2));
                if ((!rmdup_ && isClean) || (rmdup_ && dupRateOnly_ && isClean))
                {
                    boost::mutex::scoped_lock lock(sizeMutex_);
                    cleanDataIndexs_[atomic_read32(&size_)] = i;
                    atomic_inc32(&size_);
                }
            }
        }
        else //SE
        {
            for (int i=start; i<end; ++i)
            {
                upper(reads1[i].baseSequence);
                if(UtoT){
                    turnBase(reads1[i].baseSequence, 'U', 'T');
                }
                isClean = statisticsSE(reads1, i, &param->info1);
                if ((!rmdup_ && isClean) || (rmdup_ && dupRateOnly_ && isClean))
                {
                    boost::mutex::scoped_lock lock(sizeMutex_);
                    cleanDataIndexs_[atomic_read32(&size_)] = i;
                    atomic_inc32(&size_);
                }
            }
        }

        atomic_inc32(&doneNum_);
    }

    bool FilterProcessor::statisticsSE(Read *reads, int index, FqInfo *info)
    {
        StatisResult sr;
        StatisInfo si;

        Read* read = &reads[index];

        int tailTrimTemp1_ = tailTrim_ ;

        int index1_ = adaptorIndex(read,adapter1_,adapterLen1_,readsName1_,sr);

        if(index1_ != -1 && cutAdaptor){
            if(index1_ >= minReadLength){
                int cutLen1_ = strlen(read->baseSequence) - index1_;

                if(cutLen1_ > tailTrim_)
                    tailTrimTemp1_ = cutLen1_;

                info->totalCutAdaptorNum++;
            }
        }

        si = auxStatistics(read, headTrim_, tailTrimTemp1_, adapter1_, adapterLen1_, readsName1_, *info, sr);

        if(strlen(read->baseSequence) < minReadLength){
            return false;
        }

        if (!sr.hasAdpt) //no adapter
        {
            if (!sr.nExceed) // N rate not exceed
            {
                if (!sr.isLowQual) // not low quality
                {
                    if (sr.sumQuality >= minMean_ * si.readLen) //not low mean quality
                    {
                        if(!sr.isPolyA) //not polyA
                        {
                            if (rmdup_)
                            {
                                pair<map<string,int>::iterator,bool> ret;
                                {
                                    boost::mutex::scoped_lock lock(dupMutex_);
                                    MD5((const unsigned char *)(read->baseSequence), strlen(read->baseSequence), md5Seq_);
                                    ret = duplications_.insert(pair<string, int>(string((const char*)md5Seq_, MD5_DIGEST_LENGTH), index));
                                }

                                if (ret.second)
                                {
                                    info->greyTotalReadNum++;
                                    return true;
                                }
                                else
                                {
                                    info->duplicationNum++;
                                    if(dupRateOnly_)
                                        return true;
                                }
                            }
                            else
                            {
                                info->cleanBaseA += si.a;
                                info->cleanBaseC += si.c;
                                info->cleanBaseG += si.g;
                                info->cleanBaseT += si.t;
                                info->cleanBaseN += si.n;
                                info->cleanQ20 += si.q20;
                                info->cleanQ30 += si.q30;

                                info->cleanTotalBaseNum += si.readLen;
                                info->cleanTotalReadNum++;

                                calculateBaseDistribute(read, *info, si.readLen);

                                return true;
                            }
                        }
                        else //is polyA
                        {
                            info->polyANum++;
                        }
                    }
                    else
                    {
                        info->lowMeanNum++;
                    }
                }
                else // low quality
                {
                    info->lowQualNum++;
                }
            }
            else //n rate exceed
            {
                info->nExceedNum++;
            }
        }
        else //has adapter
        {
            info->adapterNum++;
        }

        return false;

    }

    bool FilterProcessor::statisticsPE(Read *reads1, Read *reads2, int index, FqInfo *info1, FqInfo *info2)
    {
        StatisResult sr1, sr2;
        Read* read1 = &reads1[index];
        Read* read2 = &reads2[index];

        int tailTrimTemp1_ = tailTrim_ ;
        int tailTrimTemp2_ = tailTrim2_;

//		int headTrimTemp1_ = headTrim_ ;
//		int headTrimTemp2_ = headTrim2_;

        int index1_ = adaptorIndex(read1,adapter1_,adapterLen1_,readsName1_,sr1);
        int index2_ = adaptorIndex(read2,adapter2_,adapterLen2_,readsName2_,sr2);
        int cutLen1_, cutLen2_;

        if(index1_ != -1 || index2_ != -1){
            int minLen;
            if(index1_ != -1 && index2_ != -1)
                minLen = index1_ > index2_ ? index1_ : index2_;
            else
                minLen = index1_ == -1 ? index2_ : index1_;

            if(minLen >= minReadLength){

                cutLen1_ = strlen(read1->baseSequence) - minLen;
                cutLen2_ = strlen(read2->baseSequence) - minLen;

                if(cutLen1_ > tailTrim_)
                    tailTrimTemp1_ = cutLen1_;
                if(cutLen2_ > tailTrim2_)
                    tailTrimTemp2_ = cutLen2_;

                info1->totalCutAdaptorNum++;
                info2->totalCutAdaptorNum++;
            }
        }

        //fq1
        StatisInfo si1 = auxStatistics(read1, headTrim_, tailTrimTemp1_, adapter1_, adapterLen1_, readsName1_, *info1, sr1);
        //fq2
        StatisInfo si2 = auxStatistics(read2, headTrim2_, tailTrimTemp2_, adapter2_, adapterLen2_, readsName2_, *info2, sr2);

        if(cutLen1_ < minReadLength || cutLen2_ < minReadLength){
            return false;
        }

        if (!sr1.hasAdpt && !sr2.hasAdpt) //no adapter
        {
            if (!sr1.nExceed && !sr2.nExceed) //n rate not exceed
            {
                if (!sr1.isLowQual && !sr2.isLowQual) //quality not too low
                {
                    if ((sr1.sumQuality >= minMean_ * si1.readLen) && (sr2.sumQuality >= minMean_ * si2.readLen)) //not low mean quality
                    {
                        if (!isFilterSmallInsertSize_
                            || !isSmallSize(read1->baseSequence, si1.readLen, read2->baseSequence, si2.readLen))
                        {
                            if ((polyAType_==0 && (!sr1.isPolyA || !sr2.isPolyA)) ||(polyAType_==1 && !sr1.isPolyA && !sr2.isPolyA)) //not polyA
                            {
                                if (rmdup_) // remove duplication
                                {
                                    pair<map<string,int>::iterator,bool> ret;
                                    {
                                        boost::mutex::scoped_lock lock(dupMutex_);
                                        string temp(read1->baseSequence);
                                        temp += read2->baseSequence;
                                        MD5((const unsigned char *)(temp.c_str()), temp.length(), md5Seq_);
                                        ret = duplications_.insert(pair<string, int>(string((const char*)md5Seq_, MD5_DIGEST_LENGTH), index));
                                    }

                                    if (ret.second) //not duplication
                                    {
                                        info1->greyTotalReadNum++;
                                        return true;
                                    }
                                    else //duplication
                                    {
                                        info1->duplicationNum++;
                                        info2->duplicationNum++;
                                        info1->totalDuplicationNum++;
                                        if(dupRateOnly_)
                                            return true;
                                    }
                                } // not remove duplication
                                else
                                {
                                    info1->cleanBaseA += si1.a;
                                    info1->cleanBaseC += si1.c;
                                    info1->cleanBaseG += si1.g;
                                    info1->cleanBaseT += si1.t;
                                    info1->cleanBaseN += si1.n;
                                    info1->cleanQ20 += si1.q20;
                                    info1->cleanQ30 += si1.q30;
                                    info1->cleanTotalReadNum++;
                                    info1->cleanTotalBaseNum += si1.readLen;
                                    calculateBaseDistribute(read1, *info1, si1.readLen);

                                    info2->cleanBaseA += si2.a;
                                    info2->cleanBaseC += si2.c;
                                    info2->cleanBaseG += si2.g;
                                    info2->cleanBaseT += si2.t;
                                    info2->cleanBaseN += si2.n;
                                    info2->cleanQ20 += si2.q20;
                                    info2->cleanQ30 += si2.q30;
                                    info2->cleanTotalReadNum++;
                                    info2->cleanTotalBaseNum += si2.readLen;
                                    calculateBaseDistribute(read2, *info2, si2.readLen);

                                    return true;
                                }
                            }
                            else //has polyA
                            {
                                if(polyAType_==0){
                                    info1->polyANum++;
                                    info2->polyANum++;
                                }else{
                                    if(sr1.isPolyA)
                                        info1->polyANum++;
                                    if(sr2.isPolyA)
                                        info2->polyANum++;
                                }
                                info1->totalPolyANum++;
                            }
                        }
                        else  //small insert size
                        {
                            info1->smallInsertNum++;
                            info2->smallInsertNum++;
                            info1->totalSmallInsertNum++;
                        }
                    }
                    else //low mean quality
                    {
                        info1->lowMeanNum++;
                        info2->lowMeanNum++;
                        info1->totalLowMeanNum++;
                    }
                }
                else //quality low
                {
                    info1->totalLowQualNum++;
                    if (sr1.isLowQual)
                    {
                        info1->lowQualNum++;
                    }
                    if (sr2.isLowQual)
                    {
                        info2->lowQualNum++;
                    }
                }

            }
            else //n rate exceed
            {
                info1->totalNExceedNum++;
                if (sr1.nExceed)
                {
                    info1->nExceedNum++;
                }
                if (sr2.nExceed)
                {
                    info2->nExceedNum++;
                }
            }
        }
        else // has adapter
        {
            info1->totalAdapterNum++;
            if (sr1.hasAdpt)
            {
                info1->adapterNum++;
            }
            if (sr2.hasAdpt)
            {
                info2->adapterNum++;
            }
        }

        return false;
    }

    StatisInfo FilterProcessor::auxStatistics(StrRead *read)
    {
        StatisInfo si;
        int readLen = read->baseSequence.length();
        int qual;
        for (int i=0; i<readLen; ++i)
        {
            switch (read->baseSequence[i])
            {
                case 'A':
                    si.a++;
                    break;
                case 'C':
                    si.c++;
                    break;
                case 'G':
                    si.g++;
                    break;
                case 'T':
                    si.t++;
                    break;
                case 'N':
                    si.n++;
                    break;
            }

            qual = read->baseQuality[i] - cleanQualSys_;

            if (qual <= lowQual_)
            {
                si.lowQual++;
            }
            else if (qual >=20)
            {
                si.q20++;
                if (qual >= 30)
                {
                    si.q30++;
                }
            }
        }

        return si;
    }

    int FilterProcessor::adaptorIndex(Read *read,string adapter,int adpLen,set<string> &readsName,StatisResult &sr)
    {
        int readLen = strlen(read->baseSequence);
        int index = -1;
        if(filterAdapter_)
        {
            if(isAdptList_)
            {
                sr.hasAdpt = hasAdapter(readsName,read->readName);
            }
            else
            {
                index = hasAdapter(read->baseSequence, readLen, adapter.c_str(), adpLen);
                if(index != -1)
                {
                    if(!cutAdaptor || (cutAdaptor && index < minReadLength))
                    {
                        sr.hasAdpt = true;
                    }
                }
            }
        }
        return index;
    }

    StatisInfo FilterProcessor::auxStatistics(Read *read, int headTrim, int tailTrim, string adapter, int adptLen, set<string> &readsName, FqInfo &info, StatisResult &sr)
    {
        int qual;

        int readLen = strlen(read->baseSequence);

        info.rawTotalBaseNum += readLen;
        info.rawTotalReadNum++;

        int badHeadTrim = 0;
        for(; badHeadTrim < trimBadHeadMaxLength; badHeadTrim++){
            if(read->baseQuality[badHeadTrim] - qualSys_ >= trimBadHeadQuility)
                break;
        }

        int badTailTrim = 0;
        for(; badTailTrim < trimBadTailMaxLength; badTailTrim++){
            if(read->baseQuality[readLen - 1 - badTailTrim] - qualSys_ >= trimBadTailQuility)
                break;
        }

        headTrim = headTrim > badHeadTrim ? headTrim : badHeadTrim;
        tailTrim = tailTrim > badTailTrim ? tailTrim : badTailTrim;

        /*if (filterAdapter_)
        {
            if (isAdptList_)
            {
                sr.hasAdpt = hasAdapter(readsName, read->readName);
            }
            else
            {
                sr.hasAdpt = hasAdapter(read->baseSequence, readLen, adapter.c_str(), adptLen);
            }
        }*/

        StatisInfo si;

        int a = 0, g = 0, c = 0, t = 0, n = 0;
        int q20 = 0, q30 = 0;

        int right = readLen - tailTrim;
        int sumQual = 0;

        for (int i=0; i<readLen; ++i)
        {
            switch (read->baseSequence[i])
            {
                case 'A':
                    a++;
                    info.base[i][A_]++;
                    if (i>=headTrim && i<right)
                    {
                        si.a++;
                    }
                    break;
                case 'C':
                    c++;
                    info.base[i][C_]++;
                    if (i>=headTrim && i<right)
                    {
                        si.c++;
                    }
                    break;
                case 'G':
                    g++;
                    info.base[i][G_]++;
                    if (i>=headTrim && i<right)
                    {
                        si.g++;
                    }
                    break;
                case 'T':
                    t++;
                    info.base[i][T_]++;
                    if (i>=headTrim && i<right)
                    {
                        si.t++;
                    }
                    break;
                case 'N':
                    n++;
                    info.base[i][N_]++;
                    if (i>=headTrim && i<right)
                    {
                        si.n++;
                    }
                    break;
            }

            qual = read->baseQuality[i] - qualSys_;
            if (qual > MAX_QUALITY)
            {
                LOG(WARN, "some bases' quality larger than " << MAX_QUALITY << ", they have been set to " << MAX_QUALITY);
                qual = MAX_QUALITY;
            }

            if (qual < 0)
            {
                LOG(WARN, "some bases' quality smaller than 0, they have been set to 0");
                qual = 0;
            }

            if (qual > info.maxQualityValue)
            {
                info.maxQualityValue = qual;
            }


            if (i>=headTrim && i<right)
            {
                sumQual += qual;
            }

            read->baseQuality[i] = qual + cleanQualSys_;

            info.qual[i][qual]++;

            if (qual >= 20)
            {
                q20++;
                info.q20q30[i][0]++;
                if (qual >= 30)
                {
                    q30++;
                    info.q20q30[i][1]++;
                }
            }

            if (i>=headTrim && i<right)
            {
                if (qual <= lowQual_)
                {
                    si.lowQual++;
                }
                else if (qual >=20)
                {
                    si.q20++;
                    if (qual >= 30)
                    {
                        si.q30++;
                    }
                }
            }
        }

        info.rawBaseA += a;
        info.rawBaseC += c;
        info.rawBaseG += g;
        info.rawBaseT += t;
        info.rawBaseN += n;
        info.rawQ20 += q20;
        info.rawQ30 += q30;

        //clean data read length
        readLen = readLen - headTrim - tailTrim;
        if(readLen < 0){
            readLen = 0;
            read->baseSequence[0] = '\0';
            read->baseQuality[0] = '\0';
        }else{
            //截断read的两端
            read->baseSequence[right] = '\0';
            read->baseSequence = read->baseSequence + headTrim;
            read->baseQuality[right] = '\0';
            read->baseQuality = read->baseQuality + headTrim;
        }

        si.readLen = readLen;

        sr.nExceed = (si.n >= readLen * nRate_);
        sr.isLowQual = (si.lowQual >= readLen * qualRate_);
        sr.sumQuality = sumQual;

        if (polyA_ >= 1E-6)
        {
            sr.isPolyA = (1.0 * si.a / readLen) >= (polyA_ - 1E-6);
        }

        if(maskLowQual >= 0){
            maskLowQualBase(read->baseQuality, read->baseSequence, maskLowQual + cleanQualSys_);
        }

        if(TtoU){
            turnBase(read->baseSequence, 'T', 'U');
        }

        //去掉index
        if (filterIndex_)
        {
            if(seqType_==0){
                int sharpIndex = 0;
                int i = 0;
                while (read->readName[i++] != '#' and read->readName[i] != '\0')
                    ;
                if (read->readName[i-1] != '\0'){
                    sharpIndex = i - 1;
                    while (read->readName[i++] != '/')
                        ;
                    strcpy(read->readName + sharpIndex, read->readName + i - 1);
                }
            }else{
                int i = strlen(read->readName);
                while (read->readName[--i] !=':')
                    ;
                if(i>0)
                    read->readName[i]='\0';
            }
        }

        return si;
    }

    bool FilterProcessor::hasAdapter(set<string> &readsName, const char *seqName)
    {
        int i = 0;
        if(seqType_==0){
            while (seqName[i++] != '/')
                ;
        }else{
            while (seqName[i++] != ' ')
                ;
        }

        return readsName.count(string(seqName + 1, i-2));
    }

    int FilterProcessor::hasAdapter(const char *sequence, int readLen, const char *adapter, int adptLen)
    {
        int find = -1;
        int minMatchLen = (int)ceil(adptLen * matchRatio_);
        int a1 = adptLen - minMatchLen;
        int r1 = 0;
        int len, mis;

        int right = readLen - minMatchLen;

        for (r1 = 0; r1 <= right;)
        {
            int len1 = adptLen - a1;
            int len2 = readLen - r1;
            len = (len1 < len2) ? len1 : len2;
            mis = 0;
            int map[MAX_LENGTH];
            map[0] = 0;
            for (int c = 0; c < len; ++c)
            {
                if (adapter[a1 + c] == sequence[r1 + c])
                {
                    map[mis]++;
                }
                else
                {
                    mis++;
                    map[mis] = 0;
                }
            }
            int max_map = 0;
            for (int c = 0; c <= mis; ++c)
            {
                if (map[c] > max_map)
                {
                    max_map = map[c];
                }
            }
            if ((mis <= misMatch_) || (max_map >= minMatchLen))
            {
                find = r1;
                break;
            }
            if (a1 > 0)
            {
                a1--;
            }
            else
            {
                r1++;
            }
        }

        return find;
    }

    int FilterProcessor::getReadsNameFromFile(string filename, set<string> &readsName)
    {
        gzFile file = gzopen(filename.c_str(), "rb");
        if (file == NULL)
        {
            LOG(ERROR, "cannot open file: " + filename);
            return 1;
        }

        char buf[512];
        int terminate = 0;
        int len;
        gzgets(file, buf, 512);
        if(seqType_ == 0){
            while (gzgets(file, buf, 512) != NULL)
            {
                terminate = 0;
                len = strlen(buf);
                while (buf[terminate++] != '/' && (terminate < len))
                    ;
                if (terminate == len)
                {
                    LOG(ERROR, "[FORMAT ERROR] " << filename);
                    return 1;
                }
                buf[terminate - 1] = '\0';
                readsName.insert(buf);
            }
        }else{
            while (gzgets(file, buf, 512) != NULL)
            {
                terminate = 0;
                len = strlen(buf);
                while (buf[terminate++] != ' ' && (terminate < len))
                    ;
                if (terminate == len)
                {
                    LOG(ERROR, "[FORMAT ERROR] " << filename);
                    return 1;
                }
                buf[terminate - 1] = '\0';
                readsName.insert(buf);
            }
        }

        return 0;
    }

    bool FilterProcessor::isSmallSize(const char *sequence1, int readLen1, const char *sequence2, int readLen2)
    {
        string temp(sequence2);
        string seq2 = reverseComplementary(temp);

        int mismatch, max_mismatch;
        int max_match_length = (readLen1 > readLen2) ? readLen2 : readLen1;
        for (int i = overlap_; i <= max_match_length; i++)
        {
            max_mismatch = static_cast<int>(mis_ * i);
            mismatch = 0;
            for (int j = 0; j < i; j++)
            {
                if (sequence1[readLen1 - i + j] == 'N' || seq2[j] == 'N')
                {
                    mismatch++;
                }
                else if (sequence1[readLen1 - i + j] != seq2[j])
                {
                    mismatch++;
                }
            }
            if (mismatch <= max_mismatch)
            {
                return true;
            }
        }

        return false;
    }

    void FilterProcessor::calculateBaseDistribute(Read* read, FqInfo &info, int readLen)
    {
        int qual;

        for (int i=0; i<readLen; ++i)
        {
            switch (read->baseSequence[i])
            {
                case 'A':
                    ++info.clean_base[i][A_];
                    break;
                case 'C':
                    ++info.clean_base[i][C_];
                    break;
                case 'G':
                    ++info.clean_base[i][G_];
                    break;
                case 'T':
                    ++info.clean_base[i][T_];
                    break;
                case 'N':
                    ++info.clean_base[i][N_];
                    break;
            }

            qual = read->baseQuality[i] - cleanQualSys_;
            ++info.clean_qual[i][qual];

            if (qual >= 20)
            {
                ++info.clean_q20q30[i][0];
                if (qual >= 30)
                {
                    ++info.clean_q20q30[i][1];
                }
            }
        }
    }

    void FilterProcessor::calculateBaseDistribute(StrRead* read, FqInfo &info, int readLen)
    {
        int qual;

        for (int i=0; i<readLen; ++i)
        {
            switch (read->baseSequence[i])
            {
                case 'A':
                    ++info.clean_base[i][A_];
                    break;
                case 'C':
                    ++info.clean_base[i][C_];
                    break;
                case 'G':
                    ++info.clean_base[i][G_];
                    break;
                case 'T':
                    ++info.clean_base[i][T_];
                    break;
                case 'N':
                    ++info.clean_base[i][N_];
                    break;
            }

            qual = read->baseQuality[i] - cleanQualSys_;

            ++info.clean_qual[i][qual];

            if (qual >= 20)
            {
                ++info.clean_q20q30[i][0];
                if (qual >= 30)
                {
                    ++info.clean_q20q30[i][1];
                }
            }
        }
    }
}  // namespace PreProcessTool
