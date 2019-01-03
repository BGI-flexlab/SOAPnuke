#include "process_argv.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <map>
#include "global_parameter.h"
//#include <string>
using namespace::std;
#define ADA_RATIO 0.4
map<string,vector<string>> wrong_paras;
void check_module(int argc,char* argv[]){
    if(argc<2){
        printModule();
    }else{
        set<string> modules;
        modules.insert("filter");
        modules.insert("filtersRNA");
        modules.insert("filterMeta");
        string module=argv[1];
        if(modules.find(module)==modules.end()){
            if(module=="-h" || module=="--help"){
                printModule();
            }else if(module=="-v" || module=="--version"){
                printVersion();
            }else{
                cerr<<"Error:no such module,type -h/--help for help"<<endl;
                exit(1);
            }
        }else{
            if(argc==2){
                printUsage(module);
            }
        }
    }
}

int global_parameter_initial(int argc,char* argv[],C_global_parameter& gp){
	//const char *shortOptions = "f:r:1:2:K:M:A:l:T:q:n:m:p:d3in:N:t:e:c:SO:P:Q:L:I:G:a:o:C:D:R:W:5:6:7:8:9:Eb:x:y:z:hv";
    string c_module(argv[1]);
    const char *shortOptions ="E:j1:2:R:W:C:D:o:5:8:Jaf:r:Z:z:c:d:k:Y:K:F:iQ:G:l:q:m:x:y:n:p:g:X:t:B:O:P:7e:T:3:4:w:M:A:9:S:s:U:u:b:0:hv";
    const struct option longOptions[] =
            {
                //common parameter
                    //input and output file
                    {"mode",1,NULL,'E'},
                    {"streaming",0,NULL,'j'},
                    { "fq1"     , 1, NULL, '1' },
                    { "fq2"     , 1, NULL, '2' },
                    { "trimFq1"  , 1, NULL, 'R' },
                    { "trimFq2"  , 1, NULL, 'W' },
                    { "cleanFq1", 1, NULL, 'C' },
                    { "cleanFq2", 1, NULL, 'D' },
                    { "outDir"  , 1, NULL, 'o' },
                    //input and output file type
                    { "seqType"  , 1, NULL, '5' },
                    {"outFileType",1,NULL,'8'},

                    {"ada_trim",0,NULL,'J'},
                    {"contam_trim",0,NULL,'a'},
                    { "adapter1", 1, NULL, 'f' },
                    { "adapter2", 1, NULL, 'r' },
                    { "contam1"   , 1, NULL, 'Z' },
                    { "contam2"   , 1, NULL, 'z' },
                    { "ctMatchR"  , 1,NULL,'Y'},
                    { "global_contams",1,NULL,'c'},
                    { "glob_cotm_mR",1,NULL,'d'},
                    { "glob_cotm_mM",1,NULL,'k'},
                    { "tile"    , 1, NULL, 'K' },
                    { "fov"    , 1, NULL, 'F' },
                    //index remove
                    { "index"   , 0, NULL, 'i' },
                    //base quality
                    { "qualSys" , 1, NULL, 'Q' },
                    { "outQualSys"  , 1, NULL, 'G' },
                    { "lowQual" , 1, NULL, 'l' },
                    { "qualRate", 1, NULL, 'q' },
                    { "mean"    , 1, NULL, 'm' },
                    { "trimBadHead" , 1, NULL, 'x' },
                    { "trimBadTail" , 1, NULL, 'y' },
                    //base content
                    { "nRate"   , 1, NULL, 'n' },
                    { "highA"   , 1, NULL, 'p' },
                    { "polyG_tail"   , 1, NULL, 'g' },
                    { "polyX"   , 1, NULL, 'X' },
                    { "trim"    , 1, NULL, 't' },
                    //  { "TtoU" , 0, &TtoU, 1 },
                  //  { "UtoT" , 0, &UtoT, 1 },
                    {"baseConvert",1,NULL,'B'},
                    //{ "small"   , 0, NULL, 'S' },
                    //PE reads
                    { "overlap" , 1, NULL, 'O' },
                    { "mis"     , 1, NULL, 'P' },
                    { "pe_info" , 0, NULL, '7' },
                    //computer resource
                    { "patch"     , 1, NULL, 'e' },
                    { "thread"  , 1, NULL, 'T' },
                    //{ "split_line",1,NULL, '6'},
                    //read length limit
                    {"maxReadLen",1,NULL,'3'},
                    {"minReadLen",1,NULL,'4'},
                    //reads number limit
                   // {"'totalReadsNum'"     , 1, NULL, 'c' },
                    {"output_clean",1,NULL,'w'},
                    {"adaMis",1,NULL,'M'},
                    {"adaMR",1,NULL,'A'},
                    {"adaEdge",1,NULL,'9'},

                    {"adaRCtg",1,NULL,'S'},
                    {"adaRAr",1,NULL,'s'},
                    {"adaRMa",1,NULL,'U'},
                    {"adaREr",1,NULL,'u'},
                    {"adaRMm",1,NULL,'b'},
                    
                 //   { "append"  , 1, NULL, 'a' },
                    { "log",1,NULL,'0'},
                    { "help"    , 0, NULL, 'h' },
                    { "version" , 0, NULL, 'v' },
                    {     0,    0,    0,    0}
            };

    int nextOpt;
    gp.module_name=argv[1];
    if(gp.module_name=="filtersRNA"){
        gp.min_read_length=18;
        gp.max_read_length=49;
    }
    gp.log="log";
    int error=0;
    while (-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions, NULL)))
    {
        switch (nextOpt)
        {
            case 'E':gp.mode.assign(optarg);break;
            case 'j':gp.is_streaming=true;break;
            case '1':gp.fq1_path.assign(optarg);break;
            case '2':gp.fq2_path.assign(optarg);break;
            case 'R':gp.trim_fq1.assign(optarg);break;
            case 'W':gp.trim_fq2.assign(optarg);break;
            case 'C':gp.clean_fq1.assign(optarg);break;
            case 'D':gp.clean_fq2.assign(optarg);break;
            case 'o':gp.output_dir.assign(optarg);break;
            case '5':gp.seq_type.assign(optarg);break;
            case '8':gp.output_file_type.assign(optarg);break;
            case 'J':gp.adapter_discard_or_trim="trim";break;
            case 'a':gp.contam_discard_or_trim="trim";break;
           	case 'f':gp.adapter1_seq.assign(optarg);break;
            case 'r':gp.adapter2_seq.assign(optarg);break;
            case 'Z':gp.contam1_seq.assign(optarg);break;
            case 'z':gp.contam2_seq.assign(optarg);break;
            case 'c':gp.global_contams.assign(optarg);break;
            case 'd':gp.g_mrs.assign(optarg);break;
            case 'k':gp.g_mms.assign(optarg);break;
            case 'Y':gp.ctMatchR.assign(optarg);break;
            case 'K':gp.tile.assign(optarg);break;
            case 'F':gp.fov.assign(optarg);break;
            case 'i':gp.index_remove = true;break;
            case 'Q':{
                gp.qualityPhred=atoi(optarg);
                if(gp.qualityPhred==1){
                    gp.qualityPhred=64;
                }else if(gp.qualityPhred==2){
                    gp.qualityPhred=33;
                }
                break;
            }
            case 'G':{
                gp.outputQualityPhred=atoi(optarg);
                if(gp.outputQualityPhred==1){
                    gp.outputQualityPhred=64;
                }else if(gp.outputQualityPhred==2){
                    gp.outputQualityPhred=33;
                }
                break;
            }
            case 'l':gp.lowQual=atoi(optarg);break;
            case 'q':gp.lowQualityBaseRatio=atof(optarg);break;
            case 'm':gp.meanQuality=atoi(optarg);break;
            case 'x':gp.trimBadHead.assign(optarg);break;
            case 'y':gp.trimBadTail.assign(optarg);break;
            case 'n':gp.n_ratio=atof(optarg);break;
            case 'p':gp.highA_ratio=atof(optarg);break;
            case 'g':gp.polyG_tail=atof(optarg);break;
            case 'X':gp.polyX_num=atof(optarg);break;
            case 't':gp.trim.assign(optarg);break;
            case 'B':gp.base_convert.assign(optarg);break;
            case 'O':gp.overlap_length=atoi(optarg);break;
            case 'P':gp.peMismatchRatio=atof(optarg);break;
            case '7':gp.whether_add_pe_info=true;break;
            case 'e':gp.patchSize=atoi(optarg);break;
            case 'T':gp.threads_num=atoi(optarg);break;
           // case '6':gp.split_line=atoi(optarg);break;
            case '3':gp.max_read_length=atoi(optarg);break;
            case '4':gp.min_read_length=atoi(optarg);break;
            //case 'c':gp.total_reads_num=atof(optarg)*1000*1000;break;
            case 'w':gp.output_clean=atoi(optarg);break;
            case 'M':gp.adaMis=atoi(optarg);wrong_paras["filtersRNA"].push_back("-M|--adaMis");break;
            case 'A':gp.adaMR=atof(optarg);wrong_paras["filtersRNA"].push_back("-A|adaMR");break;
            case '9':gp.adaEdge=atoi(optarg);wrong_paras["filtersRNA"].push_back("-9|--adaEdge");break;
            case 'S':gp.adaRCtg=atoi(optarg);wrong_paras["filter"].push_back("-S|--adaRCtg");break;
            case 's':gp.adaRAr=atof(optarg);wrong_paras["filter"].push_back("-s|--adaRAr");break;
            case 'U':gp.adaRMa=atoi(optarg);wrong_paras["filter"].push_back("-U|--adaRMa");break;
            case 'u':gp.adaREr=atof(optarg);wrong_paras["filter"].push_back("-u|--adaREr");break;
            case 'b':gp.adaRMm=atoi(optarg);wrong_paras["filter"].push_back("-b|--adaRMm");break;
           // case 'd':gp.rmdup = true;break;
            case '0':gp.log.assign(optarg);break;
            case 'v':printVersion();return 1;
            case 'h':printUsage(c_module);return 1;
            default:{
                exit(1);
            }
        }
    }
    if (argc != optind+1)
    {
        cerr << "Error:please check the options" << endl;
        exit(1);
    }
    if(gp.log.find("/")==string::npos){
        gp.log=gp.output_dir+"/"+gp.log;
    }
    
    if(gp.fq1_path.rfind(".gz")!=gp.fq1_path.size()-3){
        gp.mode="ssd";
    }
    if(gp.patchSize==0){
        gp.patchSize=gp.threads_num*10000/8;
    }
    /*
    int min_adapter_length=gp.adapter1_seq.size()>gp.adapter2_seq.size()?gp.adapter2_seq.size():gp.adapter1_seq.size();
    if(min_adapter_length>20){
        gp.adaEdge=min_adapter_length*ADA_RATIO;
    }
    */
}
bool check_parameter(int argc,char* argv[],C_global_parameter& gp){
    if(!gp.fq1_path.empty()){
        check_gz_file(gp.fq1_path);
    }else{
        cerr<<"Error:input fastq1 is required"<<endl;
        exit(1);
    }
    if(gp.output_dir.empty()){
        cerr<<"Error:output directory is required"<<endl;
        exit(1);
    }
    bool pe_data=false;
    if(!gp.fq2_path.empty()){
        pe_data=true;
        check_gz_file(gp.fq2_path);
    }
    if(gp.clean_fq1.empty()){
        cerr<<"Error:output clean fastq is required"<<endl;
        exit(1);
    }else{
        if(pe_data){
             if(gp.clean_fq2.empty()){
                cerr<<"Error:output clean fastq2 is required"<<endl;
                exit(1);
            }
        }
    }
    if(!pe_data && gp.module_name!="filtersRNA"){
        if(!gp.adapter2_seq.empty()){
            cerr<<"Error:no need adapter2"<<endl;
            exit(1);
        }
    }
    if(!pe_data){
        if(!gp.trim_fq2.empty() || !gp.clean_fq2.empty()){
            cerr<<"Error:input file is not pe data"<<endl;
            exit(1);
        }
    }else{
        if(!(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3 && gp.fq2_path.find(".gz")==gp.fq2_path.size()-3) && !(gp.fq1_path.find(".gz")!=gp.fq1_path.size()-3 && gp.fq2_path.find(".gz")!=gp.fq2_path.size()-3)){
            cerr<<"Error:the format of input fastq1 is inconsistent with fastq2"<<endl;
            exit(1);
        }
    }
    if(gp.seq_type!="0" && gp.seq_type!="1"){
        cerr<<"Error:seq_type value should be 0 or 1"<<endl;
        exit(1);
    }
    if(gp.output_file_type!="fastq" && gp.output_file_type!="fasta"){
        cerr<<"Error:output_file_type value should be fastq or fasta"<<endl;
        exit(1);
    }
    if(!gp.adapter_discard_or_trim.empty()){
        if(gp.adapter_discard_or_trim!="trim" && gp.adapter_discard_or_trim!="discard"){
            cerr<<"Error:adapter_discard_or_trim value should be trim or discard"<<endl;
            exit(1);
        }
    }
    if(!gp.tile.empty()){
        if(gp.tile.find("-")!=string::npos){
            vector<string> tmp_eles;
            line_split(gp.tile,'-',tmp_eles);
            string start,end;
            for(string::size_type ix=tmp_eles[0].size()-1;ix>=0;ix--){
                if(!isalpha(tmp_eles[0][ix])){
                    start.insert(start.begin(),tmp_eles[0][ix]);
                }
            }
            for(string::size_type ix=tmp_eles[1].size()-1;ix>=0;ix--){
                if(!isalpha(tmp_eles[1][ix])){
                    end.insert(end.begin(),tmp_eles[1][ix]);
                }
            }
            if(atoi(start.c_str())>atoi(end.c_str())){
                cerr<<"Error:tile value format error"<<endl;
                exit(1);
            }
        }
        for(string::size_type ix=0;ix!=gp.tile.size();ix++){
            if(!isalnum(gp.tile[ix]) && gp.tile[ix]!='-' && gp.tile[ix]!=','){
                cerr<<"Error:tile value format error"<<endl;
                exit(1);
            }
        }
    }
    if(gp.qualityPhred!=64 && gp.qualityPhred!=33){
        cerr<<"Error:qualityPhred value error"<<endl;
        exit(1);
    }
    if(gp.outputQualityPhred!=64 && gp.outputQualityPhred!=33){
        cerr<<"Error:outputQualityPhred value error"<<endl;
        exit(1);
    }
    if(gp.module_name=="filter" || gp.module_name=="filterMeta"){
        if(wrong_paras["filter"].size()>0){
            cerr<<"Error:these parameters should not appear in the module,"<<join_vector(wrong_paras["filter"],',')<<endl;
            exit(1);
        }
    }else{
        if(gp.module_name=="filtersRNA"){
            //adaMis,adaMR;
            if(wrong_paras["filtersRNA"].size()>0){
                cerr<<"Error:these parameters should not appear in the module,"<<join_vector(wrong_paras["filtersRNA"],',')<<endl;
                exit(1);
            }
        }
    }
    if(gp.output_clean!=0 && gp.output_clean<gp.patchSize){
        cerr<<"Error: output reads in each clean fastq file(-w) should be more than patch size(-e)"<<endl;
        exit(1);
    }
    if(!gp.trim.empty()){
        vector<string> tmp_eles;
        line_split(gp.trim,',',tmp_eles);
        if(pe_data){
            if(tmp_eles.size()!=4){
                cerr<<"Error:trim value format error"<<endl;
                exit(1);
            }
        }else{
            if(tmp_eles.size()!=2){
                cerr<<"Error:trim value format error"<<endl;
                exit(1);
            }
        }
    }
    if(!gp.trimBadHead.empty()){
        vector<string> tmp_eles;
        line_split(gp.trimBadHead,',',tmp_eles);
        if(pe_data){
            if(tmp_eles.size()!=2){
                cerr<<"Error:trimBadHead value format error"<<endl;
                exit(1);
            }
        }else{
            if(tmp_eles.size()!=1){
                cerr<<"Error:trimBadHead value format error"<<endl;
                exit(1);
            }
        }
        
    }
    if(!gp.trimBadTail.empty()){
        vector<string> tmp_eles;
        line_split(gp.trimBadTail,',',tmp_eles);
        if(pe_data){
            if(tmp_eles.size()!=2){
                cerr<<"Error:trimBadTail value format error"<<endl;
                exit(1);
            }
        }else{
            if(tmp_eles.size()!=1){
                cerr<<"Error:trimBadTail value format error"<<endl;
                exit(1);
            }
        }
    }
    if(!gp.base_convert.empty()){
        set<string> acgt_s;
        string tmp_s="ACGTacgt";
        string tmp_single;
        for(string::size_type ix=0;ix!=tmp_s.size();ix++){
            tmp_single.insert(tmp_single.end(),tmp_s[ix]);
            acgt_s.insert(tmp_single);
            tmp_single="";
        }
        if(gp.base_convert.find("TO")==string::npos && gp.base_convert.find("2")==string::npos){
            cerr<<"Error:base_convert value format error"<<endl;
            exit(1);
        }else{
            char sep=gp.base_convert.find("TO")!=string::npos?'O':'2';
            gp.base_convert.replace(gp.base_convert.find("TO"),2,"O");
            vector<string> tmp_eles;
            line_split(gp.base_convert,sep,tmp_eles);
            if(tmp_eles[0].size()!=1 || acgt_s.find(tmp_eles[0])==acgt_s.end() || tmp_eles[1].size()!=1 || acgt_s.find(tmp_eles[1])==acgt_s.end()){
                cerr<<"Error:base_convert value format error"<<endl;
                exit(1);
            }
        }
    }
    if(gp.threads_num>48){
        cerr<<"Error:threads number is limited to 48"<<endl;
        exit(1);
    }
    if(gp.patchSize>5000000){
        cerr<<"Error:patchSize cannot exceed 5M considering memory usage"<<endl;
        exit(1);
    }
}
void printModule(){
    cout << endl;
    cout << "Program: SOAPnuke\n";
    cout << "Version: " << PACKAGEVERSION <<endl;
    cout << "Contact: GongChun<gongchun@genomics.cn>  ChenYuXin<chenyuxin@genomics.cn>"<<endl;
    cout << "Command:\n";
    cout << "         filter        preprocessing sequences\n";
    cout << "         filtersRNA    preprocessing sRNA sequences\n";
    //cout << "         filterDGE     preprocessing DGE sequences\n";   //not include filterDGE in thie version
    cout << "         filterMeta    preprocessing Meta sequences\n";
    cout << endl;
    exit(1);
}
void printUsage(string c_module){
    cout << "Usage: "<<c_module<<" [OPTION]... \n";
    cout<<"\ncommon options\n";
    cout << "\t-E, --mode\t\tSTR\t\tif pigz software is available,you can assign ssd mode for accelerate(-E ssd). Default is non SSD\n";
    cout << "\t-1, --fq1\t\tFILE\t\tfq1 file(required)\n";
    cout << "\t-2, --fq2\t\tFILE\t\tfq2 file, used when pe\n";
    cout << "\t-C, --cleanFq1\t\tSTR\t\tclean fq1 file name(required,gz format)\n";
    cout << "\t-D, --cleanFq2\t\tSTR\t\tclean fq2 file name\n";
    cout << "\t-o, --outDir\t\tSTR\t\toutput directory, directory must exists\n";
    cout << "\t-8, --outFileType\tSTR\t\toutput file format: fastq or fasta[fastq]\n";
    cout << "\t-0, --log\t\tSTR\t\tlog file\n";
    cout << "\n";
    cout << "\t-f, --adapter1\t\tSTR\t\tadapter sequence of fq1 file (5' adapter when filtersRNA mode)\n";
    cout << "\t-r, --adapter2\t\tSTR\t\tadapter sequence of fq2 file (for PE reads or 3' adapter when filtersRNA)\n";
    cout << "\t-Z, --contam1\t\tSTR\t\tcontaminant sequence(s) for fq1 file, split by comma\n";
    cout << "\t-z, --contam2\t\tSTR\t\tcontaminant sequence(s) for fq2 file, split by comma\n";
    cout << "\t-Y, --ctMatchR\t\tFLOAT/STR\tcontam's shortest consistent matching ratio [0.2]\n";
    cout << "\t-c, --global_contams\tSTR\t\tglobal contaminant sequences which need to be detected, split by comma if more than 1\n";
    cout << "\t-d, --glob_cotm_mR\tSTR\t\tminimum match ratio in global contaminant sequences detection\n";
    cout << "\t-k, --glob_cotm_mM\tSTR\t\tmaximum mismatch number in global contaminant sequences detection\n";
    cout << "\t-5, --seqType\t\tINT\t\tSequence fq name type, 0->old fastq name, 1->new fastq name [0]\n";
    cout << "\t\t\t\t\t\t\t\told fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
    cout << "\t\t\t\t\t\t\t\tnew fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
    cout << "\t-R, --trimFq1\t\tSTR\t\ttrim fq1 file name(gz format)\n";
    cout << "\t-W, --trimFq2\t\tSTR\t\ttrim fq2 file name\n";
    cout << "\t-K, --tile\t\tSTR\t\ttile number to ignore reads, such as [1101-1104,1205]\n";
    cout << "\t-F, --fov\t\tSTR\t\tfov number to ignore reads (only for zebra-platform data), such as [C001R003,C003R004]\n";
    cout << "\tAdapter related:\n";
    cout << "\t-J, --ada_trim\t\t\t\ttrim read when find adapter[discard]\n";
    cout << "\t-a, --contam_trim\t\t\ttrim read when find contam[discard]\n";
    if(c_module=="filter" || c_module=="filtermeta"){
        //cout << "filter and filtermeta module adapter related parameter\n";
        cout << "\t-M, --adaMis\t\tINT\t\tthe max mismatch number when match the adapter (depend on -f/-r)  [1]\n";
        cout << "\t-A, --adaMR\t\tFLOAT\t\tadapter's shortest match ratio (depend on -f/-r)  [0.5]\n";
        cout << "\t-9, --adaEdge\t\tINT\t\tthe min length for segmental alignment [6]\n";
    }else if(c_module=="filtersRNA"){
        //cout << "filtersRNA module adapter related parameter\n";
    cout << "  find 5' adapter\n";
    cout << "\t-S, --adaRCtg\tINT\t\tmini 5' adapter continuous alignment length (default: 6)\n";
    cout << "\t-s, --adaRAr\tFLOAT\t\tmini alignment rate when find 5' adapter: alignment/tag (default: 0.8)\n";
    cout << "  find 3' adapter\n";
    cout << "\t-U, --adaRMa\tINT\t\tmini alignment length when find 3' adapter (default: 5)\n";
    cout << "\t-u, --adaREr\tFLOAT\t\tMax error rate when find 3' adapter (mismatch/match) (dfault: 0.4)\n";
    cout << "\t-b, --adaRMm\tINT\t\tMax mismatch number when find 3' adapter (dfault: 4)" << endl;
    }
    cout<<endl;
   
    cout << "\t-l, --lowQual\t\tINT\t\tlow quality threshold  [5]\n";
    cout << "\t-q, --qualRate\t\tFLOAT\t\tlow quality rate  [0.5]\n";
    cout << "\t-n, --nRate\t\tFLOAT\t\tN rate threshold  [0.05]\n";
    //cout << "\t-N, --maskLowQual INT       Turn those bases with low quality into N, set INT as the quality threshold  [-1]\n";
    cout << "\t-m, --mean\t\tFLOAT\t\tfilter reads with low average quality\n";
    cout << "\t-p, --highA\t\tFLOAT\t\tfilter reads if ratio of A in a read exceed [FLOAT]\n";
    cout << "\t-g, --polyG_tail\tFLOAT\t\tfilter reads if found polyG in tail [INT]\n";
    cout << "\t-X, --polyX\t\tINT\t\tfilter reads if a read contains polyX [INT]\n";
    //cout << "\t-d, --rmdup                 remove PCR duplications\n";
    //cout << "\t-3, --dupRate               keep PCR duplicated reads and calculate duplications rate\n";
    cout << "\t-i, --index\t\t\t\tremove index\n";
    //cout << "\t-c, --'totalReadsNum'\tFLOAT\t\tthe read number you want to keep in each clean fq file, (unit:1000*1000, 0 means not cut reads)\n";
    cout << "\t-t, --trim\t\tINT,INT,INT,INT\n";
    cout << "\t    \t\t\t\t\ttrim some bp of the read's head and tail, they means: (PE type:read1's head and tail and read2's head and tail  [0,0,0,0]; SE type:read head and tail [0,0])\n";
    
    cout << "\t-x, --trimBadHead\tINT,INT\t\tTrim from head ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]\n";
    cout << "\t-y, --trimBadTail\tINT,INT\t\tTrim from tail ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]\n";
    cout << "\n";

    cout << "\t-O, --overlap\t\tINT\t\tfilter the small insert size.Not filter until the value exceed 1[-1]\n";
    cout << "\t-P, --mis\t\tFLOAT\t\tthe maximum mismatch ratio when find overlap between PE reads(depend on -O)[0.1]\n";
    cout << "\n";
    //cout << "\t-6, --split_line\tINT\t\tsplit raw fastq by <split_line>, default 10M reads per file (if ssd mode is open)\n";
    cout << "\t-e, --patch\t\tINT\t\treads number of a patch processed[400000]\n";
    cout << "\t-T, --thread\t\tINT\t\tprocess thread number[6]" << endl;
    cout << "\n";
    cout << "\t-Q, --qualSys\t\tINT\t\tquality system 1:64, 2:33[default:2]\n";
    cout << "\t-G, --outQualSys\tINT\t\tout quality system 1:64, 2:33[default:2]\n";
    cout << "\t-3, --maxReadLen\tINT\t\tread max length,default 49 for filtersRNA\n";
    cout << "\t-4, --minReadLen\tINT\t\tread min length,default 18 for filtersRNA,30 for other modules\n";
    cout << "\t-w, --output_clean\tINT\t\tmax reads number in each output clean fastq file\n";
//      cout << "\t-a, --append      STR       the log's output place : console or file  [console]\n";
    
    
    cout << "\n";
    
    cout << "\t-7, --pe_info\t\t\t\tAdd /1, /2 at the end of fastq name.[default:not add]\n";
    cout << "\t-B, --baseConvert\tSTR\t\tconvert base when write data,example:TtoU\n";
    cout << "\n";
    cout << "\t-h, --help\t\t\t\thelp" << endl;
    cout << "\t-v, --version\t\t\t\tshow version" << endl;
    if(c_module=="filter"){
        cout << "\tExample:  ./SOAPnuke filter -l 10 -q 0.1 -n 0.01 -M 2 --adaMR 0.5 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG  -1 test.r1.fq.gz -2 test.r2.fq.gz -C clean_1.fq.gz -D clean_2.fq.gz  -o result -T 16"<<endl;
    }
    exit(1);
}
void printVersion(){
    cerr << "SOAPnuke filter tools version "<<PACKAGEVERSION<<"\n";
    exit(1);
}
