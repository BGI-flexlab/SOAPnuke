#include "process_argv.h"
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <map>
#include <fstream>
#include "global_parameter.h"
#include "sys/sysinfo.h"
// #include <string>
using namespace ::std;
#define ADA_RATIO 0.4

map<string, vector<string>> wrong_paras;

void check_module(int argc, char *argv[])
{
	if (argc < 2)
	{
		printModule();
	}
	else
	{
		set<string> modules;
		modules.insert("filter");
		modules.insert("filtersRNA");
		modules.insert("filterMeta");
#ifdef _PROCESSHTS
		modules.insert("filterHts");
#endif
		modules.insert("filterStLFR");
		string module = argv[1];
		if (modules.find(module) == modules.end())
		{
			if (module == "-h" || module == "--help")
			{
				printModule();
			}
			else if (module == "-v" || module == "--version")
			{
				printVersion();
			}
			else
			{
				if (module == "filterHts")
				{
					cout
						<< "Error:filterHts module not installed, please re-make after set USEHTS true in Makefile, details see Install part in Readme"
						<< endl;
				}
				cerr << "Error:no such module,type -h/--help for help" << endl;
				exit(1);
			}
		}
		else
		{
			if (argc == 2)
			{
				if (module == "filterHts")
				{
					printHtsUsage();
				}
				else
				{
					printUsage(module);
				}
			}
		}
	}
}

int global_parameter_initial(int argc, char *argv[], C_global_parameter &gp)
{
	// const char *shortOptions = "f:r:1:2:K:M:A:l:T:q:n:m:p:d3in:N:t:e:c:SO:P:Q:L:I:G:a:o:C:D:R:W:5:6:7:8:9:Eb:x:y:z:hv";
	string c_module(argv[1]);
	//    const char *shortOptions ="j1:2:R:W:C:D:o:5:8:E:Jaf:r:Z:z:c:d:k:Y:K:F:iQ:G:l:q:m:x:y:n:p:g:X:t:B:O:P:7e:T:3:4:L:w:M:A:9:S:s:U:u:b:0:hv";
	const char *shortOptions = "j1:2:C:D:o:c:E:Jf:r:l:q:m:x:y:n:p:g:X:t:T:3:4:L:w:hv";
	// move some uncommonly used parameters to config file
	// include RW58aZzYcdkKFiBOP7eMA9SsUub0QG
	const struct option longOptions[] = {
		// common parameter
		// input and output file
		//{"mode",1,NULL,'E'},
		{"streaming", 0, NULL, 'j'},
		{"fq1", 1, NULL, '1'},
		{"fq2", 1, NULL, '2'},

		//                    { "trimFq1"  , 1, NULL, 'R' },//uncommonly used parameters
		//                    { "trimFq2"  , 1, NULL, 'W' },//uncommonly used parameters

		{"cleanFq1", 1, NULL, 'C'},
		{"cleanFq2", 1, NULL, 'D'},
		{"outDir", 1, NULL, 'o'},
		{"configFile", 1, NULL, 'c'},

		// input and output file type
		//                    { "seqType"  , 1, NULL, '5' },//uncommonly used parameters
		//                    {"outFileType",1,NULL,'8'},//uncommonly used parameters

		// reference for cram input
		{"ref", 1, NULL, 'E'},
		{"ada_trim", 0, NULL, 'J'},
		//                    {"contam_trim",0,NULL,'a'},//uncommonly used parameters
		{"adapter1", 1, NULL, 'f'},
		{"adapter2", 1, NULL, 'r'},
		//                    { "contam1"   , 1, NULL, 'Z' },//uncommonly used parameters
		//                    { "contam2"   , 1, NULL, 'z' },//uncommonly used parameters

		//                    { "ctMatchR"  , 1,NULL,'Y'},//uncommonly used parameters
		//
		//                    { "global_contams",1,NULL,'c'},//uncommonly used parameters
		//                    //uncommonly used parameters
		//                    { "glob_cotm_mR",1,NULL,'d'},//uncommonly used parameters
		//                    { "glob_cotm_mM",1,NULL,'k'},//uncommonly used parameters
		//                    { "tile"    , 1, NULL, 'K' },//uncommonly used parameters
		//                    { "fov"    , 1, NULL, 'F' },//uncommonly used parameters
		//                    //index remove
		//                    { "index"   , 0, NULL, 'i' },//uncommonly used parameters
		//
		//                    //base quality
		//                    { "qualSys" , 1, NULL, 'Q' },//uncommonly used parameters
		//                    { "outQualSys"  , 1, NULL, 'G' },//uncommonly used parameters
		{"lowQual", 1, NULL, 'l'},
		{"qualRate", 1, NULL, 'q'},
		{"mean", 1, NULL, 'm'},
		{"trimBadHead", 1, NULL, 'x'},
		{"trimBadTail", 1, NULL, 'y'},
		// base content
		{"nRate", 1, NULL, 'n'},
		{"highA", 1, NULL, 'p'},
		{"polyG_tail", 1, NULL, 'g'},
		{"polyX", 1, NULL, 'X'},
		{"trim", 1, NULL, 't'},
		//  { "TtoU" , 0, &TtoU, 1 },
		//  { "UtoT" , 0, &UtoT, 1 },

		//                    {"baseConvert",1,NULL,'B'},//uncommonly used parameters
		//                    //{ "small"   , 0, NULL, 'S' },
		//                    //PE reads
		//                    { "overlap" , 1, NULL, 'O' },//uncommonly used parameters
		//                    { "mis"     , 1, NULL, 'P' },//uncommonly used parameters
		//                    { "pe_info" , 0, NULL, '7' },//uncommonly used parameters
		//                    //computer resource
		//                    { "patch"     , 1, NULL, 'e' },//uncommonly used parameters

		{"thread", 1, NULL, 'T'},
		//{ "split_line",1,NULL, '6'},
		// read length limit
		//                    {"maxReadLen",1,NULL,'3'},//uncommonly used parameters
		{"minReadLen", 1, NULL, '4'},
		// reads number limit
		// {"'totalReadsNum'"     , 1, NULL, 'c' },
		//                    {"totalReadsNum",1,NULL,'L'},
		{"output_clean", 1, NULL, 'w'},

		//                    {"adaMis",1,NULL,'M'},//uncommonly used parameters
		//                    {"adaMR",1,NULL,'A'},//uncommonly used parameters
		//                    {"adaEdge",1,NULL,'9'},//uncommonly used parameters
		//
		//                    {"adaRCtg",1,NULL,'S'},//uncommonly used parameters
		//                    {"adaRAr",1,NULL,'s'},//uncommonly used parameters
		//                    {"adaRMa",1,NULL,'U'},//uncommonly used parameters
		//                    {"adaREr",1,NULL,'u'},//uncommonly used parameters
		//                    {"adaRMm",1,NULL,'b'},//uncommonly used parameters

		//   { "append"  , 1, NULL, 'a' },
		//                    { "log",1,NULL,'0'},//uncommonly used parameters
		{"help", 0, NULL, 'h'},
		{"version", 0, NULL, 'v'},
		{0, 0, 0, 0}};

	int nextOpt;
	gp.module_name = argv[1];
	if (gp.module_name == "filtersRNA")
	{
		gp.min_read_length = 18;
		gp.max_read_length = 49;
	}
	gp.log = "log";
	//    int error=0;
	while (-1 != (nextOpt = getopt_long(argc, argv, shortOptions, longOptions, NULL)))
	{
		switch (nextOpt)
		{
		case 'E':
			gp.reference.assign(optarg);
			break;
		case 'j':
			gp.is_streaming = true;
			break;
		case '1':
		{
			gp.fq1_path.assign(optarg);
			if (gp.fq1_path.rfind(".gz") == gp.fq1_path.size() - 3)
			{
				gp.inputGzformat = true;
			}
			else
			{
				gp.inputGzformat = false;
			}
			break;
		}
		case '2':
			gp.fq2_path.assign(optarg);
			break;
			//            case 'R':{
			//                gp.trim_fq1.assign(optarg);
			//                if(gp.trim_fq1.rfind(".gz")==gp.trim_fq1.size()-3){
			//                    gp.trimOutGzformat=true;
			//                }else{
			//                    gp.trimOutGzformat=false;
			//                }
			//                break;
			//            }
			//            case 'W':gp.trim_fq2.assign(optarg);break;
		case 'C':
		{
			gp.clean_fq1.assign(optarg);
			if (gp.clean_fq1.rfind(".gz") == gp.clean_fq1.size() - 3)
			{
				gp.cleanOutGzFormat = true;
			}
			else
			{
				gp.cleanOutGzFormat = false;
			}
			break;
		}
		case 'D':
			gp.clean_fq2.assign(optarg);
			break;
		case 'o':
			gp.output_dir.assign(optarg);
			break;
			//            case '5':gp.seq_type.assign(optarg);break;
			//            case '8':gp.output_file_type.assign(optarg);break;
		case 'J':
			gp.adapter_discard_or_trim = "trim";
			break;
			//            case 'a':gp.contam_discard_or_trim="trim";break;
		case 'f':
		{
			string validSeq = "ACGTacgtNn";
			ifstream ifAda1(optarg);
			if (!ifAda1)
			{
				gp.adapter1_seq.assign(optarg);
				gp.ada1s.push_back(gp.adapter1_seq);
				for (int i = 0; i < gp.adapter1_seq.size(); i++)
				{
					if (validSeq.find(gp.adapter1_seq[i]) == string::npos)
					{
						cerr << "Error:invalid character found in adapter:" << gp.adapter1_seq[i]
							 << ". Only ACGTacgtNn are supported" << endl;
						exit(1);
					}
				}
				break;
			}
			else
			{
				cout << "input adapter1 list file:" << optarg << endl;
				string adaLine;
				while (getline(ifAda1, adaLine))
				{
					gp.ada1s.push_back(adaLine);
				}
			}
			ifAda1.close();
			break;
		}
		case 'r':
		{
			//                gp.adapter2_seq.assign(optarg);break;
			string validSeq = "ACGTacgtNn";
			ifstream ifAda2(optarg);
			if (!ifAda2)
			{
				gp.adapter2_seq.assign(optarg);
				gp.ada2s.push_back(gp.adapter2_seq);
				for (int i = 0; i < gp.adapter2_seq.size(); i++)
				{
					if (validSeq.find(gp.adapter2_seq[i]) == string::npos)
					{
						cerr << "Error:invalid character found in adapter:" << gp.adapter2_seq[i]
							 << ". Only ACGTacgtNn are supported" << endl;
						exit(1);
					}
				}
				break;
			}
			else
			{
				cout << "input adapter2 list file:" << optarg << endl;
				string adaLine;
				while (getline(ifAda2, adaLine))
				{
					gp.ada2s.push_back(adaLine);
				}
			}
			ifAda2.close();
			break;
		}
			//            case 'Z':gp.contam1_seq.assign(optarg);break;
		case 'z':
			gp.contam2_seq.assign(optarg);
			break;
			//            case 'c':gp.global_contams.assign(optarg);break;
		case 'c':
			initFromConfigFile(gp, optarg);
			break; //-c means config file path ,no longer for global contams
				   //            case 'd':gp.g_mrs.assign(optarg);break;
				   //            case 'k':gp.g_mms.assign(optarg);break;
				   //            case 'Y':gp.ctMatchR.assign(optarg);break;
				   //            case 'K':gp.tile.assign(optarg);break;
				   //            case 'F':gp.fov.assign(optarg);break;
				   //            case 'i':gp.index_remove = true;break;
				   //            case 'Q':{
				   //                gp.qualityPhred=atoi(optarg);
				   //                if(gp.qualityPhred==1){
				   //                    gp.qualityPhred=64;
				   //                }else if(gp.qualityPhred==2){
				   //                    gp.qualityPhred=33;
				   //                }
				   //                break;
				   //            }
				   //            case 'G':{
				   //                gp.outputQualityPhred=atoi(optarg);
				   //                if(gp.outputQualityPhred==1){
				   //                    gp.outputQualityPhred=64;
				   //                }else if(gp.outputQualityPhred==2){
				   //                    gp.outputQualityPhred=33;
				   //                }
				   //                break;
				   //            }
		case 'l':
			gp.lowQual = atoi(optarg);
			break;
		case 'q':
			gp.lowQualityBaseRatio = atof(optarg);
			break;
		case 'm':
			gp.meanQuality = atoi(optarg);
			break;
		case 'x':
			gp.trimBadHead.assign(optarg);
			break;
		case 'y':
			gp.trimBadTail.assign(optarg);
			break;
		case 'n':
			gp.n_ratio = atof(optarg);
			break;
		case 'p':
			gp.highA_ratio = atof(optarg);
			break;
		case 'g':
			gp.polyG_tail = atof(optarg);
			break;
		case 'X':
			gp.polyX_num = atof(optarg);
			break;
		case 't':
			gp.trim.assign(optarg);
			break;
			//            case 'B':gp.base_convert.assign(optarg);break;
			//            case 'O':gp.overlap_length=atoi(optarg);break;
			//            case 'P':gp.peMismatchRatio=atof(optarg);break;
			//            case '7':gp.whether_add_pe_info=true;break;
			//            case 'e':gp.patchSize=atoi(optarg);break;
		case 'T':
			gp.threads_num = atoi(optarg);
			break;
			// case '6':gp.split_line=atoi(optarg);break;
			//            case '3':gp.max_read_length=atoi(optarg);break;
		case '4':
			gp.min_read_length = atoi(optarg);
			break;
			//            case 'L':{
			//                string tmp_str;
			//                tmp_str.assign(optarg);
			//                if(tmp_str.find("head")==string::npos){
			//                    gp.total_reads_num_random=true;
			//                    //gp.catWhenrunning=false;
			//                    for(int i=0;i!=tmp_str.size();i++){
			//                        if(!isdigit(tmp_str[i]) && tmp_str[i]!='.'){
			//                            cerr<<"Error:-L value should be a positive integer or float"<<endl;
			//                            exit(1);
			//                        }
			//                    }
			//                }else{
			//                    gp.total_reads_num_random=false;
			//                    tmp_str.erase(tmp_str.find("head"),4);
			//                    if(tmp_str.find(".")!=string::npos){
			//                        cerr<<"Error:-L value should be a integer when with head suffix"<<endl;
			//                        exit(1);
			//                    }else{
			//                        for(int i=0;i!=tmp_str.size();i++){
			//                            if(!isdigit(tmp_str[i])){
			//                                cerr<<"Error:-L value should be an integer when with head suffix"<<endl;
			//                                exit(1);
			//                            }
			//                        }
			//                    }
			//                }
			//
			//                float tmp_val=atof(optarg);
			//                if(tmp_val==0){
			//                    cerr<<"Error:-L value should be a positive integer or float"<<endl;
			//                    exit(1);
			//                }
			//                gp.total_reads_num=tmp_val;
			//                if(tmp_val<1){
			//                    gp.f_total_reads_ratio=tmp_val;
			//                }else{
			//                    istringstream is_str(tmp_str);
			//                    is_str>>gp.l_total_reads_num;
			//                }
			//                if(gp.f_total_reads_ratio>0 && gp.l_total_reads_num>0){
			//                    cerr<<"Error:reads number and ratio should not be both assigned at the same time"<<endl;
			//                    exit(1);
			//                }
			//                break;
			//            }
		case 'w':
		{
			string paraCheck(optarg);
			for (int i = 0; i != paraCheck.size(); i++)
			{
				if (!isdigit(paraCheck[i]))
				{
					cerr << "Error:-w value should be a positive integer" << endl;
					exit(1);
				}
			}
			gp.cleanOutSplit = atoi(paraCheck.c_str());
			if (gp.cleanOutSplit == 0)
			{
				cerr << "Error:-w value should be a positive integer" << endl;
				exit(1);
			}
			break;
		}
			//            case 'M':{
			//                wrong_paras["filtersRNA"].emplace_back("-M|--adaMis");
			//                string tmp_str;
			//                tmp_str.assign(optarg);
			//                if(tmp_str.find(",")==string::npos){
			//                    gp.adaMis=atoi(optarg);
			//                    gp.adaMis2=gp.adaMis;
			//                }else{
			//                    vector<string> values;
			//                    line_split(tmp_str,',',values);
			//                    if(values.size()<2){
			//                        cerr<<"Error:expected two values in -M parameter"<<endl;
			//                        exit(1);
			//                    }
			//                    gp.adaMis=atoi(values[0].c_str());
			//                    gp.adaMis2=atoi(values[1].c_str());
			//                }
			//                break;
			//            }
			//            case 'A':{
			//                wrong_paras["filtersRNA"].emplace_back("-A|adaMR");
			//                string tmp_str;
			//                tmp_str.assign(optarg);
			//                if(tmp_str.find(",")==string::npos){
			//                    gp.adaMR=atof(optarg);
			//                    gp.adaMR2=gp.adaMR;
			//                }else{
			//                    vector<string> values;
			//                    line_split(tmp_str,',',values);
			//                    if(values.size()<2){
			//                        cerr<<"Error:expected two values in -A parameter"<<endl;
			//                        exit(1);
			//                    }
			//                    gp.adaMR=atof(values[0].c_str());
			//                    gp.adaMR2=atof(values[1].c_str());
			//                }
			//                break;
			//            }
			//            case '9':{
			//                wrong_paras["filtersRNA"].emplace_back("-9|--adaEdge");
			//                string tmp_str;
			//                tmp_str.assign(optarg);
			//                if(tmp_str.find(",")==string::npos){
			//                    gp.adaEdge=atoi(optarg);
			//                    gp.adaEdge2=gp.adaEdge;
			//                }else{
			//                    vector<string> values;
			//                    line_split(tmp_str,',',values);
			//                    if(values.size()<2){
			//                        cerr<<"Error:expected two values in -9 parameter"<<endl;
			//                        exit(1);
			//                    }
			//                    gp.adaEdge=atoi(values[0].c_str());
			//                    gp.adaEdge2=atoi(values[1].c_str());
			//                }
			//                break;
			//            }
			//            case 'S':gp.adaRCtg=atoi(optarg);wrong_paras["filter"].emplace_back("-S|--adaRCtg");break;
			//            case 's':gp.adaRAr=atof(optarg);wrong_paras["filter"].emplace_back("-s|--adaRAr");break;
			//            case 'U':gp.adaRMa=atoi(optarg);wrong_paras["filter"].emplace_back("-U|--adaRMa");break;
			//            case 'u':gp.adaREr=atof(optarg);wrong_paras["filter"].emplace_back("-u|--adaREr");break;
			//            case 'b':gp.adaRMm=atoi(optarg);wrong_paras["filter"].emplace_back("-b|--adaRMm");break;
			//           // case 'd':gp.rmdup = true;break;
			//            case '0':gp.log.assign(optarg);break;
		case 'v':
			printVersion();
			return 1;
		case 'h':
			printUsage(c_module);
			return 1;
		default:
		{
			exit(1);
		}
		}
	}
	if (argc != optind + 1)
	{
		cerr << "Error:please check the options" << endl;
		exit(1);
	}
	if (!gp.rmdup || gp.cleanOutSplit <= 0)
	{
	}
	else
	{
		cerr << "Warning:generating split files(-w was set) would become slower when rmdup function was on" << endl;
	}
	if (gp.log.find("/") == string::npos)
	{
		gp.log = gp.output_dir + "/" + gp.log;
	}
	if (gp.fq1_path.rfind(".gz") != gp.fq1_path.size() - 3)
	{
		gp.mode = "ssd";
	}
	if (gp.patchSize == 0)
	{
		gp.patchSize = gp.threads_num * 20000 / 8;
	}
	/*
	int min_adapter_length=gp.adapter1_seq.size()>gp.adapter2_seq.size()?gp.adapter2_seq.size():gp.adapter1_seq.size();
	if(min_adapter_length>20){
		gp.adaEdge=min_adapter_length*ADA_RATIO;
	}
	*/
	return 0;
}

bool check_parameter(int argc, char *argv[], C_global_parameter &gp)
{
	// 对于filterHts模块不再做检查
	bool pe_data = false;
	if (gp.module_name != "filterHts")
	{
		if (!gp.fq1_path.empty())
		{
			if (file_exist_and_not_empty(gp.fq1_path) == 0)
			{
				cerr << "Error:input fastq1 is required" << endl;
				exit(1);
			}
		}
		else
		{
			cerr << "Error:input fastq1 is required" << endl;
			exit(1);
		}
		if (gp.output_dir.empty())
		{
			cerr << "Error:output directory is required" << endl;
			exit(1);
		}
		if (!gp.fq2_path.empty())
		{
			pe_data = true;
			if (file_exist_and_not_empty(gp.fq2_path) == 0)
			{
				cerr << "Error:input fastq2 is required" << endl;
				exit(1);
			}
			if (gp.fq1_path == gp.fq2_path)
			{
				cerr << "Error:input fq1 and fq2 are the same,please check the parameters" << endl;
				exit(1);
			}
		}
		if (gp.clean_fq1.empty())
		{
			cerr << "Error:output clean fastq is required" << endl;
			exit(1);
		}
		else
		{
			if (pe_data)
			{
				if (gp.clean_fq2.empty())
				{
					cerr << "Error:output clean fastq2 is required" << endl;
					exit(1);
				}
				if (!(
						gp.clean_fq1.rfind(".gz") == gp.clean_fq1.size() - 3 && gp.clean_fq2.rfind(".gz") == gp.clean_fq2.size() - 3) &&
					!(
						gp.clean_fq1.rfind(".gz") != gp.clean_fq1.size() - 3 && gp.clean_fq2.rfind(".gz") != gp.clean_fq2.size() - 3))
				{
					cerr << "Error:the format of clean fastq1 is inconsistent with fastq2" << endl;
					exit(1);
				}
				if (gp.cleanOutSplit > 0 || gp.total_reads_num > 0)
				{
					if (gp.clean_fq1.rfind(".gz") != gp.clean_fq1.size() - 3 && gp.clean_fq2.rfind(".gz") != gp.clean_fq2.size() - 3)
					{
						cerr << "Error:the clean out fastq should be non-gz format when clean output reads are limited"
							 << endl;
						exit(1);
					}
				}
			}
		}
		if (!pe_data && gp.module_name != "filtersRNA")
		{
			if (!gp.adapter2_seq.empty())
			{
				cerr << "Error:no need adapter2" << endl;
				exit(1);
			}
		}
		if (!pe_data)
		{
			if (!gp.trim_fq2.empty() || !gp.clean_fq2.empty())
			{
				cerr << "Error:input file is not pe data" << endl;
				exit(1);
			}
		}
		else
		{
			if (!(
					gp.fq1_path.rfind(".gz") == gp.fq1_path.size() - 3 && gp.fq2_path.rfind(".gz") == gp.fq2_path.size() - 3) &&
				!(
					gp.fq1_path.rfind(".gz") != gp.fq1_path.size() - 3 && gp.fq2_path.rfind(".gz") != gp.fq2_path.size() - 3))
			{
				cerr << "Error:the format of input fastq1 is inconsistent with fastq2" << endl;
				exit(1);
			}
		}

		if (gp.seq_type != "0" && gp.seq_type != "1")
		{
			cerr << "Error:seq_type value should be 0 or 1" << endl;
			exit(1);
		}
		if (gp.output_file_type != "fastq" && gp.output_file_type != "fasta")
		{
			cerr << "Error:output_file_type value should be fastq or fasta" << endl;
			exit(1);
		}
	}
	if (gp.module_name == "filterStLFR")
	{
		// check barcode list file
		if (gp.barcodeListPath == "")
		{
			cerr << "Error:barcode list not assigned" << endl;
			exit(1);
		}
		else
		{
			ifstream bl(gp.barcodeListPath);
			if (!bl)
			{
				cerr << "Error:cannot open such file," << gp.barcodeListPath << endl;
				exit(1);
			}
			bl.close();
		}
		// check barcode region format
		if (gp.barcodeRegionStr.find("_") == string::npos)
		{
			cerr << "Error:barcode region format error, it should be set as 101_10,117_10,133_10" << endl;
			exit(1);
		}
		else
		{
			vector<string> eles;
			line_split(gp.barcodeRegionStr, ',', eles);
			if (eles.size() != 3)
			{
				cerr << "Error:barcode region format error, it should be set as 101_10,117_10,133_10" << endl;
				exit(1);
			}
			for (vector<string>::iterator ix = eles.begin(); ix != eles.end(); ix++)
			{
				vector<string> eles2;
				line_split(*ix, '_', eles2);
				if (eles2.size() != 2)
				{
					cerr << "Error:barcode region format error, it should be set as 101_10,117_10,133_10" << endl;
					exit(1);
				}
			}
		}
	}
	if (!gp.adapter_discard_or_trim.empty())
	{
		if (gp.adapter_discard_or_trim != "trim" && gp.adapter_discard_or_trim != "discard")
		{
			cerr << "Error:adapter_discard_or_trim value should be trim or discard" << endl;
			exit(1);
		}
	}
	if (!gp.tile.empty())
	{
		if (gp.tile.find("-") != string::npos)
		{
			vector<string> tmp_eles;
			line_split(gp.tile, '-', tmp_eles);
			string start, end;
			for (string::size_type ix = tmp_eles[0].size() - 1; ix >= 0; ix--)
			{
				if (!isalpha(tmp_eles[0][ix]))
				{
					start.insert(start.begin(), tmp_eles[0][ix]);
				}
			}
			for (string::size_type ix = tmp_eles[1].size() - 1; ix >= 0; ix--)
			{
				if (!isalpha(tmp_eles[1][ix]))
				{
					end.insert(end.begin(), tmp_eles[1][ix]);
				}
			}
			if (atoi(start.c_str()) > atoi(end.c_str()))
			{
				cerr << "Error:tile value format error" << endl;
				exit(1);
			}
		}
		for (string::size_type ix = 0; ix != gp.tile.size(); ix++)
		{
			if (!isalnum(gp.tile[ix]) && gp.tile[ix] != '-' && gp.tile[ix] != ',')
			{
				cerr << "Error:tile value format error" << endl;
				exit(1);
			}
		}
	}
	if (gp.qualityPhred != 64 && gp.qualityPhred != 33)
	{
		cerr << "Error:qualityPhred value error" << endl;
		exit(1);
	}
	if (gp.outputQualityPhred != 64 && gp.outputQualityPhred != 33)
	{
		cerr << "Error:outputQualityPhred value error" << endl;
		exit(1);
	}
	if (gp.module_name == "filter" || gp.module_name == "filterMeta")
	{
		if (wrong_paras["filter"].size() > 0)
		{
			cerr << "Error:these parameters should not appear in the module," << join_vector(wrong_paras["filter"], ',')
				 << endl;
			exit(1);
		}
	}
	else
	{
		if (gp.module_name == "filtersRNA")
		{
			// adaMis,adaMR;
			if (wrong_paras["filtersRNA"].size() > 0)
			{
				cerr << "Error:these parameters should not appear in the module,"
					 << join_vector(wrong_paras["filtersRNA"], ',') << endl;
				exit(1);
			}
		}
	}
	if (gp.cleanOutSplit != 0 && gp.cleanOutSplit < gp.patchSize)
	{
		cerr << "Error: output reads in each clean fastq file(-w) should be more than patch size(-e)" << endl;
		exit(1);
	}
	if (gp.module_name != "filterHts")
	{
		if (!gp.trim.empty())
		{
			vector<string> tmp_eles;
			line_split(gp.trim, ',', tmp_eles);
			if (pe_data)
			{
				if (tmp_eles.size() != 4)
				{
					cerr << "Error:trim value format error" << endl;
					exit(1);
				}
			}
			else
			{
				if (tmp_eles.size() != 2)
				{
					cerr << "Error:trim value format error" << endl;
					exit(1);
				}
			}
			for (int i = 0; i < gp.trim.size(); i++)
			{
				if (!isdigit(gp.trim[i]) && gp.trim[i] != ',')
				{
					cerr << "Error:trim value format error:" << gp.trim << endl;
					cerr << "e.g.: -t 10 2 10 2" << endl;
					exit(1);
				}
			}
		}
		if (!gp.trimBadHead.empty())
		{
			vector<string> tmp_eles;
			line_split(gp.trimBadHead, ',', tmp_eles);
			if (pe_data)
			{
				if (tmp_eles.size() != 2)
				{
					cerr << "Error:trimBadHead value format error" << endl;
					exit(1);
				}
			}
			else
			{
				if (tmp_eles.size() != 1)
				{
					cerr << "Error:trimBadHead value format error" << endl;
					exit(1);
				}
			}
		}
		if (!gp.trimBadTail.empty())
		{
			vector<string> tmp_eles;
			line_split(gp.trimBadTail, ',', tmp_eles);
			if (pe_data)
			{
				if (tmp_eles.size() != 2)
				{
					cerr << "Error:trimBadTail value format error" << endl;
					exit(1);
				}
			}
			else
			{
				if (tmp_eles.size() != 1)
				{
					cerr << "Error:trimBadTail value format error" << endl;
					exit(1);
				}
			}
		}
	}
	if (!gp.base_convert.empty())
	{
		set<string> acgt_s;
		string tmp_s = "ACGTacgt";
		string tmp_single;
		for (string::size_type ix = 0; ix != tmp_s.size(); ix++)
		{
			tmp_single.insert(tmp_single.end(), tmp_s[ix]);
			acgt_s.insert(tmp_single);
			tmp_single = "";
		}
		if (gp.base_convert.find("TO") == string::npos && gp.base_convert.find("2") == string::npos)
		{
			cerr << "Error:base_convert value format error" << endl;
			exit(1);
		}
		else
		{
			string firstC = gp.base_convert.substr(0, 1);
			string lastC = gp.base_convert.substr(gp.base_convert.size() - 1, 1);
			if (acgt_s.find(firstC) == acgt_s.end() || acgt_s.find(lastC) == acgt_s.end())
			{
				cerr << "Error:base_convert value format error" << endl;
				exit(1);
			}
		}
	}
	if (gp.cleanOutSplit > 0 && gp.total_reads_num > 0)
	{
		cerr << "Error:-w and -L cannot be both assigned" << endl;
		exit(1);
	}
	if (gp.cleanOutSplit > 0)
	{
		gp.clean_file_reads = gp.cleanOutSplit;
	}
	else if (gp.l_total_reads_num > 0)
	{
		gp.clean_file_reads = gp.l_total_reads_num;
	}
	if (gp.threads_num > get_nprocs())
	{
		gp.threads_num = get_nprocs();
		cerr << "Warning:threads number exceeds the system cpu number" << endl;
		// exit(1);
	}
	if (gp.patchSize > 5000000)
	{
		cerr << "Error:patchSize cannot exceed 5M considering memory usage" << endl;
		exit(1);
	}
	return 0;
}

void printModule()
{
	cout << endl;
	cout << "Program: SOAPnuke\n";
	cout << "Version: " << PACKAGEVERSION << "." << MINORVERSION << endl;
	cout << "Contact: GongChun<gongchun@genomics.cn>  ChenYuXin<chenyuxin@genomics.cn>" << endl;
	cout << "Command:\n";
	cout << "         filter        preprocessing normal Fastq files\n";
#ifdef _PROCESSHTS
	cout << "         filterHts     preprocessing BAM/CRAM files\n";
#endif
	cout << "         filterStLFR   preprocessing stLFR Fastq files\n";
	cout << "         filtersRNA    preprocessing sRNA Fastq files\n";
	// cout << "         filterDGE     preprocessing DGE sequences\n";   //not include filterDGE in thie version
	cout << "         filterMeta    preprocessing Meta Fastq files\n";
	cout << endl;
	exit(1);
}

void printHtsUsage()
{
	cout << "Usage: "
		 << "filterHts"
		 << " [OPTION]... \n";
	cout << "\ncommonly used parameters\n";
	cout << "\t-E, --ref\t\tFILE\t\treference file(required when process cram format)\n";
	cout << "\t-1, \t\t\tFILE\t\tinput bam/cram file(required)\n";
	cout << "\t-2, \t\t\tFILE\t\toutput bam/cram file(required)\n";
	cout << "\t-o, --outDir\tSTR\t\toutput directory, directory must exists\n";
	cout << "\n";
	cout << "\t-f, --adapter1\t\tSTR\t\tadapter sequence or list file of read1\n";
	cout << "\t-r, --adapter2\t\tSTR\t\tadapter sequence or list file of read2 (if PE)\n";
	cout << "\t-J, --ada_trim\t\t\t\ttrim read when find adapter[default:discard]\n";
	cout << "\t-T, --thread\t\tINT\t\tthreads number used in process[6]\n";
	cout
		<< "\t-c, --configFile\tSTR\t\tconfig file which include uncommonly used parameters. Each line contains a parameter, for value needed parameter: adaMis=2, for bool parameter:contam_trim, which means change mode from discard to trim\n";
	cout << "\t-l, --lowQual\t\tINT\t\tlow quality threshold  [default:5]\n";
	cout << "\t-q, --qualRate\t\tFLOAT\t\tlow quality rate  [default:0.5]\n";
	cout << "\t-n, --nRate\t\tFLOAT\t\tN rate threshold  [default:0.05]\n";
	cout << "\t-m, --mean\t\tFLOAT\t\tfilter reads with low average quality\n";
	cout << "\t-p, --highA\t\tFLOAT\t\tfilter reads if ratio of A in a read exceed [FLOAT]\n";
	cout << "\t-g, --polyG_tail\tFLOAT\t\tfilter reads if found polyG in tail [INT]\n";
	cout << "\t-X, --polyX\t\tINT\t\tfilter reads if a read contains polyX [INT]\n";
	cout << "\t-i, --index\t\t\t\tremove index\n";
	cout
		<< "\t-L, --totalReadsNum\tINT/FLOAT\tnumber/fraction of reads you want to keep in the output clean fq file(cannot be assigned when -w is given).\n";
	cout
		<< "\t    \t\t\t\t\tIt will extract reads randomly through the total clean fq file by default, you also can get the head reads\n";
	cout << "\t    \t\t\t\t\tfor save time by add head suffix to the integer(e.g. -L 10000000head)\n";
	cout << "\t-4, --minReadLen\tINT\t\tread min length,default 18 for filtersRNA,30 for other modules\n";
	cout << "\t-h, --help\t\t\t\thelp" << endl;
	cout << "\t-v, --version\t\t\t\tshow version" << endl;
	cout
		<< "\tExample:  ./SOAPnuke filterHts -l 10 -q 0.1 -n 0.01 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG  --ref chr21.fa -1 input.bam -2 output.cram  -o result"
		<< endl;
	cout << "\n";
	cout << "uncommonly used parameters, which should be set in config file" << endl;
	//    cout << "\tcontam1\t\tSTR\t\tcontaminant sequence(s) for read1, split by comma\n";
	//    cout << "\tcontam2\t\tSTR\t\tcontaminant sequence(s) for read2, split by comma\n";
	cout << "\tctMatchR\tFLOAT/STR\tcontam's shortest consistent matching ratio [default:0.2]\n";
	//    cout << "\tglobal_contams\tSTR\t\tglobal contaminant sequences which need to be detected, split by comma if more than 1\n";
	//    cout << "\tglob_cotm_mR\tSTR\t\tminimum match ratio in global contaminant sequences detection\n";
	//    cout << "\tglob_cotm_mM\tSTR\t\tmaximum mismatch number in global contaminant sequences detection\n";
	cout << "\tseqType\t\tINT\t\tSequence fq name type, 0->old fastq name, 1->new fastq name [0]\n";
	cout << "\t\t\t\t\t\t\t\told fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
	cout << "\t\t\t\t\t\t\t\tnew fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
	cout << "\ttile\t\tSTR\t\ttile number to ignore reads, such as [1101-1104,1205]\n";
	cout << "\tfov\t\tSTR\t\tfov number to ignore reads (only for zebra-platform data), such as [C001R003,C003R004]\n";
	cout << "\tAdapter related:\n";

	//    cout << "\tcontam_trim\t\t\ttrim read when find contam[default:discard]\n";
	cout << "\tadaMis\t\tINT,[INT]\tthe max mismatch number when match the adapter (depend on -f/-r).If different values are required for fq1 and fq2, you can set two values seperated by comma,e.g. 1,2 (same as -A/-9)[1]\n";
	cout << "\tadaMR\t\tFLOAT,[FLOAT]\tadapter's shortest match ratio (depend on -f/-r)  [0.5]\n";
	cout << "\tadaEdge\t\tINT,[INT]\tthe min length for segmental alignment [6]\n";
	cout << endl;
	cout << "\toverlap\t\tINT\t\tfilter the small insert size.Not filter until the value exceed 1[-1]\n";
	cout << "\tmis\t\tFLOAT\t\tthe maximum mismatch ratio when find overlap between PE reads(depend on -O)[0.1]\n";
	cout << "\n";
	cout << "\tpatch\t\tINT\t\treads number of a patch processed[400000]\n";
	cout << "\n";
	cout << "\tqualSys\t\tINT\t\tquality system 1:64, 2:33[default:2]\n";
	cout << "\toutQualSys\tINT\t\tout quality system 1:64, 2:33[default:2]\n";
	cout << "\tmaxReadLen\tINT\t\tread max length,default 49 for filtersRNA\n";

	cout << "\n";

	cout << "\tpe_info\t\t\t\tAdd /1, /2 at the end of fastq name.[default:not add]\n";
	cout << "\tbaseConvert\tSTR\t\tconvert base when write data,example:TtoU, means convert base T to base U in the output\n";
	cout << "\tlog\t\tSTR\t\tlog file\n";
	cout << "\n";

	exit(1);
}

void printUsage(string c_module)
{
	cout << "Usage: " << c_module << " [OPTION]... \n";
	cout << "\ncommonly used parameters\n";
	// cout << "\t-E, --mode\t\tSTR\t\tif pigz software is available,you can assign ssd mode for accelerate(-E ssd). Default is non SSD\n";
	cout << "\t-1, --fq1\t\tFILE\t\tfq1 file(required), .gz or normal text format are both supported(required)\n";
	cout << "\t-2, --fq2\t\tFILE\t\tfq2 file(used when process PE data), format should be same as fq1 file, both are gz or both are normal text\n";
	cout << "\t-C, --cleanFq1\t\tSTR\t\treads which passed QC from fq1 file would output to this file\n";
	cout << "\t-D, --cleanFq2\t\tSTR\t\treads which passed QC from fq2 file would output to this file\n";
	cout << "\t-o, --outDir\t\tSTR\t\tOutput directory. Processed fq files and statistical results would be output to here\n";
	cout << "\t-f, --adapter1\t\tSTR\t\tadapter sequence or list file of read1\n";
	cout << "\t-r, --adapter2\t\tSTR\t\tadapter sequence or list file of read2 (if PE)\n";
	cout << "\t-J, --ada_trim\t\t\t\ttrim read when find adapter[default:discard]\n";
	cout << "\t-T, --thread\t\tINT\t\tthreads number used in process[6]\n";
	cout << "\t-c, --configFile\tSTR\t\tconfig file which include uncommonly used parameters. Each line contains a parameter, for value needed parameter: adaMis=2, for bool parameter:contam_trim, which means change mode from discard to trim\n";
	cout << "\t-l, --lowQual\t\tINT\t\tlow quality threshold  [default:5]\n";
	cout << "\t-q, --qualRate\t\tFLOAT\t\tlow quality rate  [default:0.5]\n";
	cout << "\t-n, --nRate\t\tFLOAT\t\tN rate threshold  [default:0.05]\n";
	cout << "\t-m, --mean\t\tFLOAT\t\tfilter reads with low average quality\n";
	cout << "\t-p, --highA\t\tFLOAT\t\tfilter reads if ratio of A in a read exceed [FLOAT]\n";
	cout << "\t-g, --polyG_tail\tFLOAT\t\tfilter reads if found polyG in tail [INT]\n";
	cout << "\t-X, --polyX\t\tINT\t\tfilter reads if a read contains polyX [INT]\n";
	cout << "\t-4, --minReadLen\tINT\t\tread min length,default 18 for filtersRNA,30 for other modules\n";
	cout << "\t-w, --cleanOutSplit\tINT\t\tmax reads number in each output clean fastq file\n";
	cout << "\t-h, --help\t\t\t\thelp" << endl;
	cout << "\t-v, --version\t\t\t\tshow version" << endl;
	if (c_module == "filter")
	{
		cout
			<< "\tExample:  ./SOAPnuke filter -l 10 -q 0.1 -n 0.01  -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG  -1 test.r1.fq.gz -2 test.r2.fq.gz -C clean_1.fq.gz -D clean_2.fq.gz  -o result -T 8"
			<< endl;
	}
	cout << "\n";
	cout << "uncommonly used parameters, which should be set in config file" << endl;
	//    cout << "\tcontam1\t\tSTR\t\tcontaminant sequence(s) for fq1 file, split by comma\n";
	//    cout << "\tcontam2\t\tSTR\t\tcontaminant sequence(s) for fq2 file, split by comma. The two parameters is similar to -f/-r, but support a list of sequences\n";
	cout << "\tAdapter related:\n";
	cout << "\tctMatchR\tFLOAT/STR\tcontam's shortest consistent matching ratio [default:0.2]\n";
	//    cout << "\tglobal_contams\tSTR\t\tglobal contaminant sequences which need to be detected, split by comma if more than 1\n";
	//    cout << "\tglob_cotm_mR\tSTR\t\tminimum match ratio in global contaminant sequences detection\n";
	//    cout << "\tglob_cotm_mM\tSTR\t\tmaximum mismatch number in global contaminant sequences detection\n";
	cout << "\tseqType\t\tINT\t\tSequence fq name type, 0->old fastq name, 1->new fastq name [0]\n";
	cout << "\t\t\t\t\t\t\t\told fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
	cout << "\t\t\t\t\t\t\t\tnew fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
	cout << "\ttrimFq1\t\tSTR\t\ttrim fq1 file name(gz format)[optional]\n";
	cout
		<< "\ttrimFq2\t\tSTR\t\ttrim fq2 file name[optional]. If trim related parameters were set on, these output files would include the total reads which only do trimming. For example, if read A failed QC after trimming, it will still output to -R/-W, but not to -C/-D\n";
	cout << "\ttile\t\tSTR\t\ttile number to ignore reads, such as [1101-1104,1205]\n";
	cout << "\tfov\t\tSTR\t\tfov number to ignore reads (only for zebra-platform data), such as [C001R003,C003R004]\n";
	//    cout << "\tcontam_trim\t\t\ttrim read when find contam[default:discard]\n";
	if (c_module == "filter" || c_module == "filtermeta")
	{
		// cout << "filter and filtermeta module adapter related parameter\n";
		cout
			<< "\tadaMis\t\tINT,[INT]\tthe max mismatch number when match the adapter (depend on -f/-r).If different values are required for fq1 and fq2, you can set two values seperated by comma,e.g. 1,2 (same as -A/-9)[1]\n";
		cout << "\tadaMR\t\tFLOAT,[FLOAT]\tadapter's shortest match ratio (depend on -f/-r)  [0.5]\n";
		cout << "\tadaEdge\t\tINT,[INT]\tthe min length for segmental alignment [6]\n";
	}
	else if (c_module == "filtersRNA")
	{
		// cout << "filtersRNA module adapter related parameter\n";
		cout << "  find 5' adapter\n";
		cout << "\tadaRCtg\tINT\t\tmini 5' adapter continuous alignment length (default: 6)\n";
		cout << "\tadaRAr\tFLOAT\t\tmini alignment rate when find 5' adapter: alignment/tag (default: 0.8)\n";
		cout << "  find 3' adapter\n";
		cout << "\tadaRMa\tINT\t\tmini alignment length when find 3' adapter (default: 5)\n";
		cout << "\tadaREr\tFLOAT\t\tMax error rate when find 3' adapter (mismatch/match) (default: 0.4)\n";
		cout << "\tadaRMm\tINT\t\tMax mismatch number when find 3' adapter (default: 4)" << endl;
	}
	else if (c_module == "filterStLFR")
	{
		/*if(para=="barcodeListPath"){
			gp.barcodeListPath.assign(value.c_str());
		}
		if(para=="barcodeRegionStr"){
			gp.barcodeRegionStr.assign(value.c_str());
		}
		if(para=="notCutNoLFR"){
			gp.notCutNoLFR=true;
		}*/
		cout << "\tbarcodeListPath\tSTR\t\tbarcode list of two columns:sequence and barcodeID\n";
		cout
			<< "\tbarcodeRegionStr\tSTR\t\tbarcode regions, such as: 101_10,117_10,145_10 or 101_10,117_10,133_10\n";
		cout << "\tnotCutNoLFR\t\t\t\tdo not cut sequence when fail found barcode\n";
		cout << "\tinputAsList\t\t\t\tinput file list not a file\n";
		cout << "\ttenX\t\t\t\t\toutput tenX format\n";
	}
	cout << "\toutFileType\tSTR\t\toutput file format: fastq or fasta[fastq]\n";

	cout << endl;
	if (c_module == "filterStLFR" || c_module == "filter")
	{
		cout
			<< "\trmdup\t\t\tremove duplicate reads. The function contains a certain false positive rate which means a small number of non-duplicate reads would be marked as duplication, and limited in large reads number\n";
		//        cout<<"\tapproximateReadsNum\tSTR/LONG\tapproximate reads number. We suggest you set this parameter if you know. e.g. 100m"<<endl;
		//        cout<<"\tmemSizeUsedInRmdup\tSTR/LONG\tmaximum memory size used in rmdup. e.g. 1g\n";
		//        cout<<"\texpectedFalsePositive\tFLOAT\texpected false positive rate which means maybe a <expectedFalsePositive> of reads would be marked as duplication which actually not. e.g. 1e-7"<<endl;
	}
	cout << "\tindex\t\t\t\tremove index\n";
	cout
		<< "\ttotalReadsNum\tINT/FLOAT\tnumber/fraction of reads you want to keep in the output clean fq file(cannot be assigned when -w is given).\n";
	cout
		<< "\t    \t\t\t\t\tIt will extract reads randomly through the total clean fq file by default, you also can get the head reads\n";
	cout << "\t    \t\t\t\t\tfor save time by add head suffix to the integer(e.g. -L 10000000head)\n";
	cout << "\ttrim\t\tINT,INT,INT,INT\n";
	cout
		<< "\t    \t\t\t\t\ttrim some bp of the read's head and tail, they means: (PE type:read1's head and tail and read2's head and tail  [0,0,0,0]; SE type:read head and tail [0,0])\n";

	cout
		<< "\ttrimBadHead\tINT,INT\t\tTrim from head ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]\n";
	cout
		<< "\ttrimBadTail\tINT,INT\t\tTrim from tail ends until meeting high-quality base or reach the length threshold, set (quality threshold,MaxLengthForTrim)  [0,0]\n";
	cout << "\n";

	cout << "\toverlap\t\tINT\t\tfilter the small insert size.Not filter until the value exceed 1[-1]\n";
	cout << "\tmis\t\tFLOAT\t\tthe maximum mismatch ratio when find overlap between PE reads(depend on -O)[0.1]\n";
	cout << "\n";
	// cout << "\tsplit_line\tINT\t\tsplit raw fastq by <split_line>, default 10M reads per file (if ssd mode is open)\n";
	cout << "\tpatch\t\tINT\t\treads number of a patch processed[400000]\n";
	cout << "\n";
	cout << "\tqualSys\t\tINT\t\tquality system 1:64, 2:33[default:2]\n";
	cout << "\toutQualSys\tINT\t\tout quality system 1:64, 2:33[default:2]\n";
	cout << "\tmaxReadLen\tINT\t\tread max length,default 49 for filtersRNA\n";

	//    cout << "\tcleanOutSplit\tINT\t\tmax reads number in each output clean fastq file\n";
	//      cout << "\tappend      STR       the log's output place : console or file  [console]\n";

	cout << "\n";

	cout << "\tpe_info\t\t\t\tAdd /1, /2 at the end of fastq name.[default:not add]\n";
	cout
		<< "\tbaseConvert\tSTR\t\tconvert base when write data,example:TtoU, means convert base T to base U in the output\n";
	cout << "\tmaxBaseQuality\tINT\t\tmaximum base quality[42]\n";
	cout << "\tlog\t\tSTR\t\tlog file\n";
	cout << "\n";

	exit(1);
}

void printVersion()
{
	cerr << "SOAPnuke filter tools version " << PACKAGEVERSION << "." << MINORVERSION << "\n";
	exit(1);
}

void initFromConfigFile(C_global_parameter &gp, char *configFile)
{
	set<string> legalParas, boolParas;
	legalParas.insert("trimFq1");
	legalParas.insert("trimFq2");
	legalParas.insert("seqType");
	legalParas.insert("outFileType");
	legalParas.insert("contam_trim");
	legalParas.insert("contam1");
	legalParas.insert("contam2");
	legalParas.insert("ctMatchR");
	legalParas.insert("global_contams");
	legalParas.insert("glob_cotm_mR");
	legalParas.insert("glob_cotm_mM");
	legalParas.insert("tile");
	legalParas.insert("fov");
	legalParas.insert("index");
	legalParas.insert("qualSys");
	legalParas.insert("outQualSys");
	legalParas.insert("baseConvert");
	legalParas.insert("maxBaseQuality");
	legalParas.insert("overlap");
	legalParas.insert("mis");
	legalParas.insert("pe_info");
	legalParas.insert("patch");
	legalParas.insert("maxReadLen");
	legalParas.insert("adaMis");
	legalParas.insert("adaMR");
	legalParas.insert("adaEdge");
	legalParas.insert("adaRCtg");
	legalParas.insert("adaRAr");
	legalParas.insert("adaRMa");
	legalParas.insert("adaREr");
	legalParas.insert("adaRMm");
	legalParas.insert("log");
	legalParas.insert("totalReadsNum");
	legalParas.insert("cleanOutSplit");
	legalParas.insert("trim");
	legalParas.insert("trimBadHead");
	legalParas.insert("trimBadTail");
	// stLFR
	legalParas.insert("barcodeListPath");
	legalParas.insert("barcodeRegionStr");
	legalParas.insert("notCutNoLFR");
	legalParas.insert("inputAsList");
	legalParas.insert("tenX");
	legalParas.insert("rmdup");
	//    legalParas.insert("approximateReadsNum");
	//    legalParas.insert("memSizeUsedInRmdup");
	//    legalParas.insert("expectedFalsePositive");

	boolParas.insert("index");
	boolParas.insert("pe_info");
	boolParas.insert("contam_trim");
	// stLFR
	boolParas.insert("notCutNoLFR");
	boolParas.insert("inputAsList");
	boolParas.insert("tenX");

	boolParas.insert("rmdup");
	string configFilePath(configFile);
	ifstream ifConfig;
	ifConfig.open(configFilePath.c_str());
	if (!ifConfig)
	{
		cerr << "Error:cannot open such file," << configFile << endl;
		exit(1);
	}
	string line;
	while (getline(ifConfig, line))
	{
		string para, value;
		if (line.find("#") == 0)
		{
			continue;
		}
		if (line.find("=") != string::npos)
		{
			vector<string> eles;
			line_split(line, '=', eles);
			if (eles.size() != 2)
			{
				cerr << "Error:unrecgonized format parameter," << line << endl;
				exit(1);
			}
			para = eles[0];
			value = eles[1];
			chomp_space(para, "all");
			chomp_space(value, "all");
		}
		else
		{
			para = line;
			if (boolParas.find(para) == boolParas.end())
			{
				cerr << "Error:this parameter should set a value," << para << endl;
				exit(1);
			}
		}
		if (legalParas.find(para) == legalParas.end())
		{
			cerr << "Error:no such parameter," << para << endl;
			exit(1);
		}
		if (para == "trimFq1")
		{
			gp.trim_fq1.assign(value);
			if (gp.trim_fq1.rfind(".gz") == gp.trim_fq1.size() - 3)
			{
				gp.trimOutGzformat = true;
			}
			else
			{
				gp.trimOutGzformat = false;
			}
		}
		if (para == "trimFq2")
		{
			gp.trim_fq2.assign(value);
		}
		if (para == "seqType")
		{
			gp.seq_type.assign(value);
		}
		if (para == "outFileType")
		{
			gp.output_file_type.assign(value);
		}
		if (para == "contam_trim")
		{
			gp.contam_discard_or_trim = "trim";
		}
		if (para == "contam1")
		{
			gp.contam1_seq.assign(value);
		}
		if (para == "contam2")
		{
			gp.contam2_seq.assign(value);
		}
		if (para == "ctMatchR")
		{
			gp.ctMatchR.assign(value.c_str());
		}
		if (para == "global_contams")
		{
			gp.global_contams.assign(value.c_str());
		}
		if (para == "glob_cotm_mR")
		{
			gp.g_mrs.assign(value.c_str());
		}
		if (para == "glob_cotm_mM")
		{
			gp.g_mms.assign(value.c_str());
		}
		if (para == "tile")
		{
			gp.tile.assign(value.c_str());
		}
		if (para == "fov")
		{
			gp.fov.assign(value.c_str());
		}
		if (para == "index")
		{
			gp.index_remove = true;
		}
		if (para == "qualSys")
		{
			gp.qualityPhred = atoi(value.c_str());
			if (gp.qualityPhred == 1)
			{
				gp.qualityPhred = 64;
			}
			else if (gp.qualityPhred == 2)
			{
				gp.qualityPhred = 33;
			}
		}
		if (para == "outQualSys")
		{
			gp.outputQualityPhred = atoi(value.c_str());
			if (gp.outputQualityPhred == 1)
			{
				gp.outputQualityPhred = 64;
			}
			else if (gp.outputQualityPhred == 2)
			{
				gp.outputQualityPhred = 33;
			}
		}
		if (para == "baseConvert")
		{
			gp.base_convert.assign(value.c_str());
		}
		if (para == "maxBaseQuality")
		{
			gp.maxBaseQuality = atoi(value.c_str());
		}
		if (para == "overlap")
		{
			gp.overlap_length = atoi(value.c_str());
		}
		if (para == "mis")
		{
			gp.peMismatchRatio = atof(value.c_str());
		}
		if (para == "pe_info")
		{
			gp.whether_add_pe_info = true;
		}
		if (para == "patch")
		{
			gp.patchSize = atoi(value.c_str());
		}
		if (para == "maxReadLen")
		{
			gp.max_read_length = atoi(value.c_str());
		}
		if (para == "adaMis")
		{
			wrong_paras["filtersRNA"].emplace_back("-M|--adaMis");
			string tmp_str;
			tmp_str.assign(value.c_str());
			if (tmp_str.find(",") == string::npos)
			{
				gp.adaMis = atoi(value.c_str());
				gp.adaMis2 = gp.adaMis;
			}
			else
			{
				vector<string> values;
				line_split(tmp_str, ',', values);
				if (values.size() < 2)
				{
					cerr << "Error:expected two values in -M parameter" << endl;
					exit(1);
				}
				gp.adaMis = atoi(values[0].c_str());
				gp.adaMis2 = atoi(values[1].c_str());
			}
		}
		if (para == "adaMR")
		{
			wrong_paras["filtersRNA"].emplace_back("-A|adaMR");
			string tmp_str;
			tmp_str.assign(value.c_str());
			if (tmp_str.find(",") == string::npos)
			{
				gp.adaMR = atof(value.c_str());
				gp.adaMR2 = gp.adaMR;
			}
			else
			{
				vector<string> values;
				line_split(tmp_str, ',', values);
				if (values.size() < 2)
				{
					cerr << "Error:expected two values in -A parameter" << endl;
					exit(1);
				}
				gp.adaMR = atof(values[0].c_str());
				gp.adaMR2 = atof(values[1].c_str());
			}
		}
		if (para == "adaEdge")
		{
			wrong_paras["filtersRNA"].emplace_back("-9|--adaEdge");
			string tmp_str;
			tmp_str.assign(value.c_str());
			if (tmp_str.find(",") == string::npos)
			{
				gp.adaEdge = atoi(value.c_str());
				gp.adaEdge2 = gp.adaEdge;
			}
			else
			{
				vector<string> values;
				line_split(tmp_str, ',', values);
				if (values.size() < 2)
				{
					cerr << "Error:expected two values in -9 parameter" << endl;
					exit(1);
				}
				gp.adaEdge = atoi(values[0].c_str());
				gp.adaEdge2 = atoi(values[1].c_str());
			}
		}
		if (para == "adaRCtg")
		{
			gp.adaRCtg = atoi(value.c_str());
			wrong_paras["filter"].emplace_back("-S|--adaRCtg");
		}
		if (para == "adaRAr")
		{
			gp.adaRAr = atof(value.c_str());
			wrong_paras["filter"].emplace_back("-s|--adaRAr");
		}
		if (para == "adaRMa")
		{
			gp.adaRMa = atoi(value.c_str());
			wrong_paras["filter"].emplace_back("-U|--adaRMa");
		}
		if (para == "adaREr")
		{
			gp.adaREr = atof(value.c_str());
			wrong_paras["filter"].emplace_back("-u|--adaREr");
		}
		if (para == "adaRMm")
		{
			gp.adaRMm = atoi(value.c_str());
			wrong_paras["filter"].emplace_back("-b|--adaRMm");
		}
		if (para == "log")
		{
			gp.log.assign(value.c_str());
		}
		if (para == "totalReadsNum")
		{
			string tmp_str;
			tmp_str.assign(value.c_str());
			if (tmp_str.find("head") == string::npos)
			{
				gp.total_reads_num_random = true;
				// gp.catWhenrunning=false;
				for (int i = 0; i != tmp_str.size(); i++)
				{
					if (!isdigit(tmp_str[i]) && tmp_str[i] != '.')
					{
						cerr << "Error:-L value should be a positive integer or float" << endl;
						exit(1);
					}
				}
			}
			else
			{
				gp.total_reads_num_random = false;
				tmp_str.erase(tmp_str.find("head"), 4);
				if (tmp_str.find(".") != string::npos)
				{
					cerr << "Error:-L value should be a integer when with head suffix" << endl;
					exit(1);
				}
				else
				{
					for (int i = 0; i != tmp_str.size(); i++)
					{
						if (!isdigit(tmp_str[i]))
						{
							cerr << "Error:-L value should be an integer when with head suffix" << endl;
							exit(1);
						}
					}
				}
			}

			float tmp_val = atof(value.c_str());
			if (tmp_val == 0)
			{
				cerr << "Error:-L value should be a positive integer or float" << endl;
				exit(1);
			}
			gp.total_reads_num = tmp_val;
			if (tmp_val < 1)
			{
				gp.f_total_reads_ratio = tmp_val;
			}
			else
			{
				istringstream is_str(tmp_str);
				is_str >> gp.l_total_reads_num;
			}
			if (gp.f_total_reads_ratio > 0 && gp.l_total_reads_num > 0)
			{
				cerr << "Error:reads number and ratio should not be both assigned at the same time" << endl;
				exit(1);
			}
		}
		if (para == "cleanOutSplit")
		{
			string paraCheck = value;
			for (int i = 0; i != paraCheck.size(); i++)
			{
				if (!isdigit(paraCheck[i]))
				{
					cerr << "Error:-w value should be a positive integer" << endl;
					exit(1);
				}
			}
			gp.cleanOutSplit = atoi(paraCheck.c_str());
			if (gp.cleanOutSplit == 0)
			{
				cerr << "Error:-w value should be a positive integer" << endl;
				exit(1);
			}
		}
		if (para == "trim")
		{
			gp.trim = value;
		}
		if (para == "trimBadHead")
		{
			gp.trimBadHead = value;
		}
		if (para == "trimBadTail")
		{
			gp.trimBadTail = value;
		}
		/*legalParas.insert("cleanOutSplit");
		legalParas.insert("trim");
		legalParas.insert("trimBadHead");
		legalParas.insert("trimBadTail");
		 */
		// stLFR
		/*string barcodeListPath;
	string barcodeRegionStr;
	bool notCutNoLFR;
		 */
		if (para == "barcodeListPath")
		{
			gp.barcodeListPath.assign(value.c_str());
		}
		if (para == "barcodeRegionStr")
		{
			gp.barcodeRegionStr.assign(value.c_str());
		}
		if (para == "notCutNoLFR")
		{
			gp.notCutNoLFR = true;
		}
		if (para == "inputAsList")
		{
			gp.inputAsList = true;
		}
		if (para == "tenX")
		{
			gp.tenX = true;
		}

		if (para == "rmdup")
		{
			gp.rmdup = true;
		}

		//        if(para=="approximateReadsNum"){
		//            if(isdigit(value[value.size()-1])){
		//                gp.approximateReadsNum=atol(value.c_str());
		//            }else if(value[value.size()-1]=='M' || value[value.size()-1]=='m'){
		//                string realValue=value.substr(0,value.size()-1);
		//                gp.approximateReadsNum=atol(realValue.c_str())*1024*1024;
		//            }else if(value[value.size()-1]=='G' || value[value.size()-1]=='g') {
		//                string realValue = value.substr(0, value.size() - 1);
		//                gp.approximateReadsNum = atol(realValue.c_str()) * 1024 * 1024 * 1024;
		//            }else{
		//                cerr<<"Error:cannot recognize value,"<<value<<endl;
		//                exit(1);
		//            }
		//        }
		//        if(para=="memSizeUsedInRmdup"){
		//            if(isdigit(value[value.size()-1])){
		//                gp.memSizeUsedInRmdup=atol(value.c_str());
		//            }else if(value[value.size()-1]=='M' || value[value.size()-1]=='m'){
		//                string realValue=value.substr(0,value.size()-1);
		//                gp.memSizeUsedInRmdup=atol(realValue.c_str())*1024*1024;
		//            }else if(value[value.size()-1]=='G' || value[value.size()-1]=='g') {
		//                string realValue = value.substr(0, value.size() - 1);
		//                gp.memSizeUsedInRmdup = atol(realValue.c_str()) * 1024 * 1024 * 1024;
		//            }else{
		//                cerr<<"Error:cannot recognize value,"<<value<<endl;
		//                exit(1);
		//            }
		//        }
		//        if(para=="expectedFalsePositive"){
		//            gp.expectedFalsePositive=stof(value.c_str());
		//            if(gp.expectedFalsePositive==0 || gp.expectedFalsePositive==1){
		//                cout<<"Warning:expected FP you set is "<<gp.expectedFalsePositive<<", please check the parameter"<<endl;
		//            }
		//        }
	}
}
