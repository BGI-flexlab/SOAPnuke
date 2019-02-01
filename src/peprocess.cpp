#include <iostream>
#include <string>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <thread>
#include <mutex>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <time.h>
#include <dirent.h>
#include <math.h>
#include "peprocess.h"
#include "process_argv.h"
#include "zlib.h"
#include "gc.h"
#include "sequence.h"
using namespace::std;
#define READBUF 500
#define random(x) (rand()%x)
mutex pe_cal_m;
mutex read_m,stat_m,write_m;
unsigned long long buffer(1024*1024*2);
string new_fq1_path,new_fq2_path;
C_filter_stat local_fs[max_thread];
C_fastq_file_stat local_raw_stat1[max_thread],local_raw_stat2[max_thread],local_trim_stat1[max_thread],local_trim_stat2[max_thread],local_clean_stat1[max_thread],local_clean_stat2[max_thread];
queue<vector<C_fastq> > thread_read_buffer1[max_thread],thread_read_buffer2[max_thread],thread_write_buffer1[max_thread],thread_write_buffer2[max_thread];
string sticky_reads1[max_thread],sticky_reads2[max_thread];
gzFile gz_fq1,gz_fq2;
map<int,int> pe1_out,pe2_out;
int queue_buffer=50;
int exceed_output1(0);
int exceed_output2(0);
peProcess::peProcess(C_global_parameter m_gp){ //initialize
	gp=m_gp;
	gv=C_global_variable();
	used_threads_num=0;
	processed_reads=0;
	srand((unsigned)time(NULL));
	ostringstream tmpstring;
	tmpstring<<rand()%100;
	random_num=tmpstring.str();
	file_end=false;
	limit_end=0;
}
void peProcess::print_stat(){	//print statistic information to the file
	string filter_out=gp.output_dir+"/Statistics_of_Filtered_Reads.txt";
	string general_out=gp.output_dir+"/Basic_Statistics_of_Sequencing_Quality.txt";
	string bs1_out=gp.output_dir+"/Base_distributions_by_read_position_1.txt";
	string bs2_out=gp.output_dir+"/Base_distributions_by_read_position_2.txt";
	string qs1_out=gp.output_dir+"/Base_quality_value_distribution_by_read_position_1.txt";
	string qs2_out=gp.output_dir+"/Base_quality_value_distribution_by_read_position_2.txt";
	string q20_out1=gp.output_dir+"/Distribution_of_Q20_Q30_bases_by_read_position_1.txt";
	string q20_out2=gp.output_dir+"/Distribution_of_Q20_Q30_bases_by_read_position_2.txt";
	string trim_stat1=gp.output_dir+"/Statistics_of_Trimming_Position_of_Reads_1.txt";
	string trim_stat2=gp.output_dir+"/Statistics_of_Trimming_Position_of_Reads_2.txt";
	ofstream of_filter_stat(filter_out.c_str());
	ofstream of_general_stat(general_out.c_str());
	ofstream of_readPos_base_stat1(bs1_out.c_str());
	ofstream of_readPos_base_stat2(bs2_out.c_str());
	ofstream of_readPos_qual_stat1(qs1_out.c_str());
	ofstream of_readPos_qual_stat2(qs2_out.c_str());
	ofstream of_q2030_stat1(q20_out1.c_str());
	ofstream of_q2030_stat2(q20_out2.c_str());
	ofstream of_trim_stat1(trim_stat1.c_str());
	ofstream of_trim_stat2(trim_stat2.c_str());
	if(!of_filter_stat){
		cerr<<"Error:cannot open such file,"<<filter_out<<endl;
		exit(1);
	}
	if(!of_general_stat){
		cerr<<"Error:cannot open such file,"<<general_out<<endl;
		exit(1);
	}
	if(!of_readPos_base_stat1 || !of_readPos_base_stat2){
		cerr<<"Error:cannot open such file,Base_distributions_by_read_position*.txt"<<endl;
		exit(1);
	}
	if(!of_readPos_qual_stat1 || !of_readPos_qual_stat2){
		cerr<<"Error:cannot open such file,Base_quality_value_distribution_by_read_position*.txt"<<endl;
		exit(1);
	}
	if(!of_q2030_stat1 || !of_q2030_stat2){
		cerr<<"Error:cannot open such file,Distribution_of_Q20_Q30_bases_by_read_position*.txt"<<endl;
		exit(1);
	}
	of_filter_stat<<"Item\t\t\t\tTotal\tPercentage\tfastq1\tfastq2\toverlap"<<endl;
	vector<string> filter_items;
	filter_items.push_back("Reads limited to output number");
	filter_items.push_back("Reads with filtered tile");
	filter_items.push_back("Reads with filtered fov");
	filter_items.push_back("Reads too short");
	filter_items.push_back("Reads too long");
	filter_items.push_back("Reads with global contam sequence");
	filter_items.push_back("Reads with contam sequence");
	filter_items.push_back("Reads with n rate exceed");
	filter_items.push_back("Reads with highA");
	filter_items.push_back("Reads with polyX");
	filter_items.push_back("Reads with low quality");
	filter_items.push_back("Reads with low mean quality");
	filter_items.push_back("Reads with small insert size");
	filter_items.push_back("Reads with adapter");

	map<string,unsigned long long> filter_number,filter_pe1,filter_pe2,filter_overlap;
	filter_number["Reads with global contam sequence"]=gv.fs.include_global_contam_seq_num;
	filter_number["Reads with contam sequence"]=gv.fs.include_contam_seq_num;
	filter_number["Reads too short"]=gv.fs.short_len_num;
	filter_number["Reads with adapter"]=gv.fs.include_adapter_seq_num;
	filter_number["Reads with low quality"]=gv.fs.low_qual_base_ratio_num;
	filter_number["Reads with low mean quality"]=gv.fs.mean_quality_num;
	filter_number["Reads with n rate exceed"]=gv.fs.n_ratio_num;
	filter_number["Reads with small insert size"]=gv.fs.over_lapped_num;
	filter_number["Reads with highA"]=gv.fs.highA_num;
	filter_number["Reads with polyX"]=gv.fs.polyX_num;
	filter_number["Reads with filtered tile"]=gv.fs.tile_num;
	filter_number["Reads with filtered fov"]=gv.fs.fov_num;
	filter_number["Reads too long"]=gv.fs.long_len_num;
	//filter_number["Reads limited to output number"]=gv.fs.output_reads_num;

	filter_pe1["Reads with contam sequence"]=gv.fs.include_contam_seq_num1;
	filter_pe1["Reads with global contam sequence"]=gv.fs.include_global_contam_seq_num1;
	filter_pe1["Reads too short"]=gv.fs.short_len_num1;
	filter_pe1["Reads with adapter"]=gv.fs.include_adapter_seq_num1;
	filter_pe1["Reads with low quality"]=gv.fs.low_qual_base_ratio_num1;
	filter_pe1["Reads with low mean quality"]=gv.fs.mean_quality_num1;
	filter_pe1["Reads with n rate exceed"]=gv.fs.n_ratio_num1;
	filter_pe1["Reads with small insert size"]=gv.fs.over_lapped_num;
	filter_pe1["Reads with highA"]=gv.fs.highA_num1;
	filter_pe1["Reads with polyX"]=gv.fs.polyX_num1;
	filter_pe1["Reads with filtered tile"]=gv.fs.tile_num;
	filter_pe1["Reads with filtered fov"]=gv.fs.fov_num;
	filter_pe1["Reads too long"]=gv.fs.long_len_num1;
	//filter_pe1["Reads limited to output number"]=gv.fs.output_reads_num;

	filter_pe2["Reads with contam sequence"]=gv.fs.include_contam_seq_num2;
	filter_pe2["Reads with global contam sequence"]=gv.fs.include_global_contam_seq_num2;
	filter_pe2["Reads too short"]=gv.fs.short_len_num2;
	filter_pe2["Reads with adapter"]=gv.fs.include_adapter_seq_num2;
	filter_pe2["Reads with low quality"]=gv.fs.low_qual_base_ratio_num2;
	filter_pe2["Reads with low mean quality"]=gv.fs.mean_quality_num2;
	filter_pe2["Reads with n rate exceed"]=gv.fs.n_ratio_num2;
	filter_pe2["Reads with small insert size"]=gv.fs.over_lapped_num;
	filter_pe2["Reads with highA"]=gv.fs.highA_num2;
	filter_pe2["Reads with polyX"]=gv.fs.polyX_num2;
	filter_pe2["Reads with filtered tile"]=gv.fs.tile_num;
	filter_pe2["Reads with filtered fov"]=gv.fs.fov_num;
	filter_pe2["Reads too long"]=gv.fs.long_len_num2;
	//filter_pe2["Reads limited to output number"]=gv.fs.output_reads_num;

	filter_overlap["Reads with contam sequence"]=gv.fs.include_contam_seq_num_overlap;
	filter_overlap["Reads with global contam sequence"]=gv.fs.include_global_contam_seq_num_overlap;
	filter_overlap["Reads too short"]=gv.fs.short_len_num_overlap;
	filter_overlap["Reads with adapter"]=gv.fs.include_adapter_seq_num_overlap;
	filter_overlap["Reads with low quality"]=gv.fs.low_qual_base_ratio_num_overlap;
	filter_overlap["Reads with low mean quality"]=gv.fs.mean_quality_num_overlap;
	filter_overlap["Reads with n rate exceed"]=gv.fs.n_ratio_num_overlap;
	filter_overlap["Reads with small insert size"]=gv.fs.over_lapped_num;
	filter_overlap["Reads with highA"]=gv.fs.highA_num_overlap;
	filter_overlap["Reads with polyX"]=gv.fs.polyX_num_overlap;
	filter_overlap["Reads with filtered tile"]=gv.fs.tile_num;
	filter_overlap["Reads with filtered fov"]=gv.fs.fov_num;
	filter_overlap["Reads too long"]=gv.fs.long_len_num_overlap;
	//filter_overlap["Reads limited to output number"]=gv.fs.output_reads_num;
	unsigned long long total_filter_fq1_num=0;
	for(map<string,unsigned long long>::iterator ix=filter_number.begin();ix!=filter_number.end();ix++){
		total_filter_fq1_num+=ix->second;
	}
	of_filter_stat<<setiosflags(ios::fixed);
	of_filter_stat<<"Total filtered read pair number\t"<<total_filter_fq1_num<<"\t100.00%\t\t"<<total_filter_fq1_num<<"\t"<<total_filter_fq1_num<<"\t"<<total_filter_fq1_num<<endl;
	for(vector<string>::iterator ix=filter_items.begin();ix!=filter_items.end();ix++){
		if(filter_number[*ix]>0){
			of_filter_stat<<*ix<<"\t"<<filter_number[*ix]<<"\t";
			of_filter_stat<<setprecision(2)<<100*(float)filter_number[*ix]/total_filter_fq1_num<<"%\t";
			of_filter_stat<<filter_pe1[*ix]<<"\t"<<filter_pe2[*ix]<<"\t"<<filter_overlap[*ix]<<endl;
		}
	}
	/*
	int total_filter_fq1_num=gv.fs.output_reads_num+gv.fs.short_len_num+gv.fs.include_contam_seq_num+gv.fs.include_adapter_seq_num+gv.fs.n_ratio_num+gv.fs.highA_num+gv.fs.tile_num+gv.fs.fov_num+gv.fs.low_qual_base_ratio_num+gv.fs.mean_quality_num+gv.fs.short_len_num+gv.fs.over_lapped_num;
	
	of_filter_stat<<"Reads too short\t\t\t"<<gv.fs.short_len_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.short_len_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.short_len_num1<<"\t"<<gv.fs.short_len_num2<<"\t"<<gv.fs.short_len_num_overlap<<endl;
	of_filter_stat<<"Reads with contam sequence\t"<<gv.fs.include_contam_seq_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.include_contam_seq_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.include_contam_seq_num1<<"\t"<<gv.fs.include_contam_seq_num2<<"\t"<<gv.fs.include_contam_seq_num_overlap<<endl;
	of_filter_stat<<"Reads with adapter\t\t"<<gv.fs.include_adapter_seq_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.include_adapter_seq_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.include_adapter_seq_num1<<"\t"<<gv.fs.include_adapter_seq_num2<<"\t"<<gv.fs.include_adapter_seq_num_overlap<<endl;
	of_filter_stat<<"Reads with low quality\t\t"<<gv.fs.low_qual_base_ratio_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.low_qual_base_ratio_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.low_qual_base_ratio_num1<<"\t"<<gv.fs.low_qual_base_ratio_num2<<"\t"<<gv.fs.low_qual_base_ratio_num_overlap<<endl;
	of_filter_stat<<"Reads with low mean quality\t"<<gv.fs.mean_quality_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.mean_quality_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.mean_quality_num1<<"\t"<<gv.fs.mean_quality_num2<<"\t"<<gv.fs.mean_quality_num_overlap<<endl;
	of_filter_stat<<"Reads with n rate exceed\t"<<gv.fs.n_ratio_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.n_ratio_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.n_ratio_num1<<"\t"<<gv.fs.n_ratio_num2<<"\t"<<gv.fs.n_ratio_num_overlap<<endl;
	of_filter_stat<<"Reads with small insert size\t"<<gv.fs.over_lapped_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.over_lapped_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.over_lapped_num<<"\t"<<gv.fs.over_lapped_num<<"\t"<<gv.fs.over_lapped_num<<endl;
	of_filter_stat<<"Reads with highA\t\t"<<gv.fs.highA_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.highA_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.highA_num1<<"\t"<<gv.fs.highA_num2<<"\t"<<gv.fs.highA_num_overlap<<endl;
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
	}
	of_filter_stat.close();
	*/
	of_general_stat<<"Item\traw reads(fq1)\tclean reads(fq1)\traw reads(fq2)\tclean reads(fq2)"<<endl;
	float raw1_rl(0),raw2_rl(0),clean1_rl(0),clean2_rl(0);
	char filter_r1_ratio[100],filter_r2_ratio[100];
	char raw_r1[7][100];
	char clean_r1[7][100];
	char raw_r2[7][100];
	char clean_r2[7][100];
	/*for(int i=0;i!=7;i++){
		raw_r1[i]={0};
		clean_r1[i]={0};
		raw_r2[i]={0};
		clean_r2[i]={0};
	}*/
	if(gv.raw1_stat.gs.reads_number!=0){
		raw1_rl=1.0*gv.raw1_stat.gs.base_number/gv.raw1_stat.gs.reads_number;
		sprintf(filter_r1_ratio,"%.2f",100*(float)total_filter_fq1_num/gv.raw1_stat.gs.reads_number);
		sprintf(raw_r1[0],"%.2f",100*(float)gv.raw1_stat.gs.a_number/gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[1],"%.2f",100*(float)gv.raw1_stat.gs.c_number/gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[2],"%.2f",100*(float)gv.raw1_stat.gs.g_number/gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[3],"%.2f",100*(float)gv.raw1_stat.gs.t_number/gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[4],"%.2f",100*(float)gv.raw1_stat.gs.n_number/gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[5],"%.2f",100*(float)gv.raw1_stat.gs.q20_num/gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[6],"%.2f",100*(float)gv.raw1_stat.gs.q30_num/gv.raw1_stat.gs.base_number);
	}
	if(gv.clean1_stat.gs.reads_number!=0){
		clean1_rl=1.0*gv.clean1_stat.gs.base_number/gv.clean1_stat.gs.reads_number;
		sprintf(clean_r1[0],"%.2f",100*(float)gv.clean1_stat.gs.a_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[1],"%.2f",100*(float)gv.clean1_stat.gs.c_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[2],"%.2f",100*(float)gv.clean1_stat.gs.g_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[3],"%.2f",100*(float)gv.clean1_stat.gs.t_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[4],"%.2f",100*(float)gv.clean1_stat.gs.n_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[5],"%.2f",100*(float)gv.clean1_stat.gs.q20_num/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[6],"%.2f",100*(float)gv.clean1_stat.gs.q30_num/gv.clean1_stat.gs.base_number);
	}
	if(gv.raw2_stat.gs.reads_number!=0){
		raw2_rl=1.0*gv.raw2_stat.gs.base_number/gv.raw2_stat.gs.reads_number;
		sprintf(filter_r2_ratio,"%.2f",100*(float)total_filter_fq1_num/gv.raw2_stat.gs.reads_number);
		sprintf(raw_r2[0],"%.2f",100*(float)gv.raw2_stat.gs.a_number/gv.raw2_stat.gs.base_number);
		sprintf(raw_r2[1],"%.2f",100*(float)gv.raw2_stat.gs.c_number/gv.raw2_stat.gs.base_number);
		sprintf(raw_r2[2],"%.2f",100*(float)gv.raw2_stat.gs.g_number/gv.raw2_stat.gs.base_number);
		sprintf(raw_r2[3],"%.2f",100*(float)gv.raw2_stat.gs.t_number/gv.raw2_stat.gs.base_number);
		sprintf(raw_r2[4],"%.2f",100*(float)gv.raw2_stat.gs.n_number/gv.raw2_stat.gs.base_number);
		sprintf(raw_r2[5],"%.2f",100*(float)gv.raw2_stat.gs.q20_num/gv.raw2_stat.gs.base_number);
		sprintf(raw_r2[6],"%.2f",100*(float)gv.raw2_stat.gs.q30_num/gv.raw2_stat.gs.base_number);
	}
	if(gv.clean2_stat.gs.reads_number!=0){
		clean2_rl=1.0*gv.clean2_stat.gs.base_number/gv.clean2_stat.gs.reads_number;
		sprintf(clean_r2[0],"%.2f",100*(float)gv.clean2_stat.gs.a_number/gv.clean2_stat.gs.base_number);
		sprintf(clean_r2[1],"%.2f",100*(float)gv.clean2_stat.gs.c_number/gv.clean2_stat.gs.base_number);
		sprintf(clean_r2[2],"%.2f",100*(float)gv.clean2_stat.gs.g_number/gv.clean2_stat.gs.base_number);
		sprintf(clean_r2[3],"%.2f",100*(float)gv.clean2_stat.gs.t_number/gv.clean2_stat.gs.base_number);
		sprintf(clean_r2[4],"%.2f",100*(float)gv.clean2_stat.gs.n_number/gv.clean2_stat.gs.base_number);
		sprintf(clean_r2[5],"%.2f",100*(float)gv.clean2_stat.gs.q20_num/gv.clean2_stat.gs.base_number);
		sprintf(clean_r2[6],"%.2f",100*(float)gv.clean2_stat.gs.q30_num/gv.clean2_stat.gs.base_number);
	}
	of_general_stat<<setiosflags(ios::fixed)<<setprecision(1)<<"Read length\t"<<raw1_rl<<"\t"<<clean1_rl<<"\t"<<raw2_rl<<"\t"<<clean2_rl<<endl;
	of_general_stat<<"Total number of reads\t"<<setprecision(15)<<gv.raw1_stat.gs.reads_number<<" (100.00%)\t"<<gv.clean1_stat.gs.reads_number<<" (100.00%)\t"<<gv.raw2_stat.gs.reads_number<<" (100.00%)\t"<<gv.clean2_stat.gs.reads_number<<" (100.00%)"<<endl;
	of_general_stat<<"Number of filtered reads\t"<<total_filter_fq1_num<<" ("<<filter_r1_ratio<<"%)\t-\t"<<total_filter_fq1_num<<" ("<<filter_r2_ratio<<"%)\t-"<<endl;
	of_general_stat<<"Total number of bases\t"<<setprecision(15)<<gv.raw1_stat.gs.base_number<<" (100.00%)\t"<<gv.clean1_stat.gs.base_number<<" (100.00%)\t"<<gv.raw2_stat.gs.base_number<<" (100.00%)\t"<<gv.clean2_stat.gs.base_number<<" (100.00%)"<<endl;
	unsigned long long filter_base1=total_filter_fq1_num*gv.raw1_stat.gs.read_length;
	unsigned long long filter_base2=total_filter_fq1_num*gv.raw1_stat.gs.read_length;
	of_general_stat<<"Number of filtered bases\t"<<setprecision(15)<<filter_base1<<" ("<<filter_r1_ratio<<"%)\t-\t"<<filter_base2<<" ("<<filter_r2_ratio<<"%)\t-"<<endl;
	of_general_stat<<"Number of base A\t"<<setprecision(15)<<gv.raw1_stat.gs.a_number<<" ("<<raw_r1[0]<<"%)\t"<<gv.clean1_stat.gs.a_number<<" ("<<clean_r1[0]<<"%)\t"<<gv.raw2_stat.gs.a_number<<" ("<<raw_r2[0]<<"%)\t"<<gv.clean2_stat.gs.a_number<<" ("<<clean_r2[0]<<"%)"<<endl;
	of_general_stat<<"Number of base C\t"<<setprecision(15)<<gv.raw1_stat.gs.c_number<<" ("<<raw_r1[1]<<"%)\t"<<gv.clean1_stat.gs.c_number<<" ("<<clean_r1[1]<<"%)\t"<<gv.raw2_stat.gs.c_number<<" ("<<raw_r2[1]<<"%)\t"<<gv.clean2_stat.gs.c_number<<" ("<<clean_r2[1]<<"%)"<<endl;
	of_general_stat<<"Number of base G\t"<<setprecision(15)<<gv.raw1_stat.gs.g_number<<" ("<<raw_r1[2]<<"%)\t"<<gv.clean1_stat.gs.g_number<<" ("<<clean_r1[2]<<"%)\t"<<gv.raw2_stat.gs.g_number<<" ("<<raw_r2[2]<<"%)\t"<<gv.clean2_stat.gs.g_number<<" ("<<clean_r2[2]<<"%)"<<endl;
	of_general_stat<<"Number of base T\t"<<setprecision(15)<<gv.raw1_stat.gs.t_number<<" ("<<raw_r1[3]<<"%)\t"<<gv.clean1_stat.gs.t_number<<" ("<<clean_r1[3]<<"%)\t"<<gv.raw2_stat.gs.t_number<<" ("<<raw_r2[3]<<"%)\t"<<gv.clean2_stat.gs.t_number<<" ("<<clean_r2[3]<<"%)"<<endl;
	of_general_stat<<"Number of base N\t"<<setprecision(15)<<gv.raw1_stat.gs.n_number<<" ("<<raw_r1[4]<<"%)\t"<<gv.clean1_stat.gs.n_number<<" ("<<clean_r1[4]<<"%)\t"<<gv.raw2_stat.gs.n_number<<" ("<<raw_r2[4]<<"%)\t"<<gv.clean2_stat.gs.n_number<<" ("<<clean_r2[4]<<"%)"<<endl;
	of_general_stat<<"Q20 number\t"<<setprecision(15)<<gv.raw1_stat.gs.q20_num<<" ("<<raw_r1[5]<<"%)\t"<<gv.clean1_stat.gs.q20_num<<" ("<<clean_r1[5]<<"%)\t"<<gv.raw2_stat.gs.q20_num<<" ("<<raw_r2[5]<<"%)\t"<<gv.clean2_stat.gs.q20_num<<" ("<<clean_r2[5]<<"%)"<<endl;
	/*
	float clean1_q20_ratio(0),clean2_q20_ratio(0);
	if(gv.clean1_stat.gs.base_number!=0)
		clean1_q20_ratio=(float)gv.clean1_stat.gs.q20_num/gv.clean1_stat.gs.base_number;
	if(gv.clean2_stat.gs.base_number!=0)
		clean2_q20_ratio=(float)gv.clean2_stat.gs.q20_num/gv.clean2_stat.gs.base_number;
	of_general_stat<<"Q20 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q20_num/gv.raw1_stat.gs.base_number<<"\t"<<clean1_q20_ratio<<"\t"<<(float)gv.raw2_stat.gs.q20_num/gv.raw2_stat.gs.base_number<<"\t"<<clean2_q20_ratio<<endl;
	*/
	of_general_stat<<"Q30 number\t"<<setprecision(15)<<gv.raw1_stat.gs.q30_num<<" ("<<raw_r1[6]<<"%)\t"<<gv.clean1_stat.gs.q30_num<<" ("<<clean_r1[6]<<"%)\t"<<gv.raw2_stat.gs.q30_num<<" ("<<raw_r2[6]<<"%)\t"<<gv.clean2_stat.gs.q30_num<<" ("<<clean_r2[6]<<"%)"<<endl;
	/*
	float clean1_q30_ratio(0),clean2_q30_ratio(0);
	if(gv.clean1_stat.gs.base_number!=0)
		clean1_q30_ratio=(float)gv.clean1_stat.gs.q30_num/gv.clean1_stat.gs.base_number;
	if(gv.clean2_stat.gs.base_number!=0)
		clean2_q30_ratio=(float)gv.clean2_stat.gs.q30_num/gv.clean2_stat.gs.base_number;
	//of_general_stat<<"Q30 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q30_num/gv.raw1_stat.gs.base_number<<"\t"<<clean1_q30_ratio<<"\t"<<(float)gv.raw2_stat.gs.q30_num/gv.raw2_stat.gs.base_number<<"\t"<<clean2_q30_ratio<<endl;
	*/
	of_general_stat.close();
	of_readPos_base_stat1<<"Pos\tA\tC\tG\tT\tN\tclean A\tclean C\tclean G\tclean T\tclean N"<<endl;
	of_readPos_base_stat2<<"Pos\tA\tC\tG\tT\tN\tclean A\tclean C\tclean G\tclean T\tclean N"<<endl;
	for(int i=0;i<gv.raw1_stat.gs.read_length;i++){
		of_readPos_base_stat1<<i+1<<"\t";
		of_readPos_base_stat2<<i+1<<"\t";
		string base_set="ACGTN";
		float raw1_cur_pos_total_base=0;
		float clean1_cur_pos_total_base=0;
		float raw2_cur_pos_total_base=0;
		float clean2_cur_pos_total_base=0;
		for(int j=0;j!=base_set.size();j++){
			raw1_cur_pos_total_base+=gv.raw1_stat.bs.position_acgt_content[i][j];
			raw2_cur_pos_total_base+=gv.raw2_stat.bs.position_acgt_content[i][j];
			clean1_cur_pos_total_base+=gv.clean1_stat.bs.position_acgt_content[i][j];
			clean2_cur_pos_total_base+=gv.clean2_stat.bs.position_acgt_content[i][j];
		}
		for(int j=0;j!=base_set.size();j++){
			of_readPos_base_stat1<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw1_stat.bs.position_acgt_content[i][j]/raw1_cur_pos_total_base<<"%\t";
		}
		for(int j=0;j!=base_set.size();j++){
			of_readPos_base_stat1<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean1_stat.bs.position_acgt_content[i][j]/clean1_cur_pos_total_base<<"%";
			if(base_set[j]!='N'){
				of_readPos_base_stat1<<"\t";
			}else{
				of_readPos_base_stat1<<endl;
			}
		}
		for(int j=0;j!=base_set.size();j++){
			of_readPos_base_stat2<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw2_stat.bs.position_acgt_content[i][j]/raw2_cur_pos_total_base<<"%\t";
		}
		for(int j=0;j!=base_set.size();j++){
			of_readPos_base_stat2<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean2_stat.bs.position_acgt_content[i][j]/clean2_cur_pos_total_base<<"%";
			if(base_set[j]!='N'){
				of_readPos_base_stat2<<"\t";
			}else{
				of_readPos_base_stat2<<endl;
			}
		}
	}
	of_readPos_base_stat1.close();
	of_readPos_base_stat2.close();

	of_readPos_qual_stat1<<"#raw fastq1 quality distribution"<<endl;
	of_readPos_qual_stat2<<"#raw fastq2 quality distribution"<<endl;
	of_readPos_qual_stat1<<"Pos\t";
	of_readPos_qual_stat2<<"Pos\t";
	int max_qual=0;
	for(int i=0;i<gv.raw1_stat.gs.read_length;i++){
		for(int j=1;j<=MAX_QUAL;j++){
			if(gv.raw1_stat.qs.position_qual[i][j]>0)
				max_qual=j;
		}
	}
	for(int i=0;i<=max_qual;i++){
		of_readPos_qual_stat1<<"Q"<<i<<"\t";
		of_readPos_qual_stat2<<"Q"<<i<<"\t";
	}
	of_readPos_qual_stat1<<"Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile"<<endl;
	of_readPos_qual_stat2<<"Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile"<<endl;
	
	float raw1_q20[gv.raw1_stat.gs.read_length],raw1_q30[gv.raw1_stat.gs.read_length];
	float raw2_q20[gv.raw1_stat.gs.read_length],raw2_q30[gv.raw1_stat.gs.read_length];
	float clean1_q20[gv.raw1_stat.gs.read_max_length],clean1_q30[gv.raw1_stat.gs.read_max_length];
	float clean2_q20[gv.raw1_stat.gs.read_max_length],clean2_q30[gv.raw1_stat.gs.read_max_length];
	for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
		of_readPos_qual_stat1<<i+1<<"\t";
		of_readPos_qual_stat2<<i+1<<"\t";
		unsigned long long raw1_q20_num(0),raw1_q30_num(0),raw1_total(0);
		unsigned long long raw2_q20_num(0),raw2_q30_num(0),raw2_total(0);
		for(int j=0;j<=max_qual;j++){
			if(j>=20){
				raw1_q20_num+=gv.raw1_stat.qs.position_qual[i][j];
				raw2_q20_num+=gv.raw2_stat.qs.position_qual[i][j];
			}
			if(j>=30){
				raw1_q30_num+=gv.raw1_stat.qs.position_qual[i][j];
				raw2_q30_num+=gv.raw2_stat.qs.position_qual[i][j];
			}
			raw1_total+=gv.raw1_stat.qs.position_qual[i][j];
			of_readPos_qual_stat1<<setiosflags(ios::fixed);
			of_readPos_qual_stat1<<setprecision(0)<<gv.raw1_stat.qs.position_qual[i][j]<<"\t";
			raw2_total+=gv.raw2_stat.qs.position_qual[i][j];
			of_readPos_qual_stat2<<setiosflags(ios::fixed);
			of_readPos_qual_stat2<<setprecision(0)<<gv.raw2_stat.qs.position_qual[i][j]<<"\t";
		}
		raw1_q20[i]=(float)raw1_q20_num/raw1_total;
		raw1_q30[i]=(float)raw1_q30_num/raw1_total;
		raw2_q20[i]=(float)raw2_q20_num/raw2_total;
		raw2_q30[i]=(float)raw2_q30_num/raw2_total;
		quartile_result raw1_quar=cal_quar_from_array(gv.raw1_stat.qs.position_qual[i],max_qual+1);
		quartile_result raw2_quar=cal_quar_from_array(gv.raw2_stat.qs.position_qual[i],max_qual+1);//lower_quar,upper_quar
		of_readPos_qual_stat1<<setiosflags(ios::fixed)<<setprecision(2)<<raw1_quar.mean<<"\t";
		of_readPos_qual_stat1<<setprecision(0)<<raw1_quar.median<<"\t"<<raw1_quar.lower_quar<<"\t"<<raw1_quar.upper_quar<<"\t"<<raw1_quar.first10_quar<<"\t"<<raw1_quar.last10_quar<<endl;
		of_readPos_qual_stat2<<setiosflags(ios::fixed)<<setprecision(2)<<raw2_quar.mean<<"\t";
		of_readPos_qual_stat2<<setprecision(0)<<raw2_quar.median<<"\t"<<raw2_quar.lower_quar<<"\t"<<raw2_quar.upper_quar<<"\t"<<raw2_quar.first10_quar<<"\t"<<raw2_quar.last10_quar<<endl;
	}
	of_readPos_qual_stat1<<"#clean fastq1 quality distribution"<<endl;
	of_readPos_qual_stat2<<"#clean fastq2 quality distribution"<<endl;
	of_readPos_qual_stat1<<"Pos\t";
	of_readPos_qual_stat2<<"Pos\t";
	for(int i=0;i<=max_qual;i++){
		of_readPos_qual_stat1<<"Q"<<i<<"\t";
		of_readPos_qual_stat2<<"Q"<<i<<"\t";
	}
	of_readPos_qual_stat1<<"Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile"<<endl;
	of_readPos_qual_stat2<<"Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile"<<endl;


	of_q2030_stat1<<"Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\tPercentage of Clean Q20+\tPercentage of Clean Q30+"<<endl;
	of_q2030_stat2<<"Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\tPercentage of Clean Q20+\tPercentage of Clean Q30+"<<endl;
	for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
		of_readPos_qual_stat1<<i+1<<"\t";
		of_readPos_qual_stat2<<i+1<<"\t";
		unsigned long long clean1_q20_num(0),clean1_q30_num(0),clean1_total(0);
		unsigned long long clean2_q20_num(0),clean2_q30_num(0),clean2_total(0);
		for(int j=0;j<=max_qual;j++){
			if(j>=20){
				clean1_q20_num+=gv.clean1_stat.qs.position_qual[i][j];
				clean2_q20_num+=gv.clean2_stat.qs.position_qual[i][j];
			}
			if(j>=30){
				clean1_q30_num+=gv.clean1_stat.qs.position_qual[i][j];
				clean2_q30_num+=gv.clean2_stat.qs.position_qual[i][j];
			}
			of_readPos_qual_stat1<<setiosflags(ios::fixed);
			clean1_total+=gv.clean1_stat.qs.position_qual[i][j];
			of_readPos_qual_stat1<<setprecision(0)<<gv.clean1_stat.qs.position_qual[i][j]<<"\t";
			clean2_total+=gv.clean2_stat.qs.position_qual[i][j];
			of_readPos_qual_stat2<<setiosflags(ios::fixed);
			of_readPos_qual_stat2<<setprecision(0)<<gv.clean2_stat.qs.position_qual[i][j]<<"\t";
		}

		clean1_q20[i]=(float)clean1_q20_num/clean1_total;
		clean1_q30[i]=(float)clean1_q30_num/clean1_total;
		clean2_q20[i]=(float)clean2_q20_num/clean2_total;
		clean2_q30[i]=(float)clean2_q30_num/clean2_total;
		quartile_result clean1_quar=cal_quar_from_array(gv.clean1_stat.qs.position_qual[i],max_qual+1);
		quartile_result clean2_quar=cal_quar_from_array(gv.clean2_stat.qs.position_qual[i],max_qual+1);//lower_quar,upper_quar
		of_readPos_qual_stat1<<setiosflags(ios::fixed)<<setprecision(2)<<clean1_quar.mean<<"\t";
		of_readPos_qual_stat1<<setprecision(0)<<clean1_quar.median<<"\t"<<clean1_quar.lower_quar<<"\t"<<clean1_quar.upper_quar<<"\t"<<clean1_quar.first10_quar<<"\t"<<clean1_quar.last10_quar<<endl;
		of_readPos_qual_stat2<<setiosflags(ios::fixed)<<setprecision(2)<<clean2_quar.mean<<"\t";
		of_readPos_qual_stat2<<setprecision(0)<<clean2_quar.median<<"\t"<<clean2_quar.lower_quar<<"\t"<<clean2_quar.upper_quar<<"\t"<<clean2_quar.first10_quar<<"\t"<<clean2_quar.last10_quar<<endl;
		of_q2030_stat1<<i+1<<setiosflags(ios::fixed)<<setprecision(2)<<"\t"<<100*raw1_q20[i]<<"%\t"<<100*raw1_q30[i]<<"%\t"<<100*clean1_q20[i]<<"%\t"<<100*clean1_q30[i]<<"%"<<endl;
		of_q2030_stat2<<i+1<<setiosflags(ios::fixed)<<setprecision(2)<<"\t"<<100*raw2_q20[i]<<"%\t"<<100*raw2_q30[i]<<"%\t"<<100*clean2_q20[i]<<"%\t"<<100*clean2_q30[i]<<"%"<<endl;
	}
	of_readPos_qual_stat1.close();
	of_readPos_qual_stat2.close();
	of_q2030_stat1.close();
	of_q2030_stat2.close();
	of_trim_stat1<<"Pos\tHeadLowQual\tHeadFixLen\tTailAdapter\tTailLowQual\tTailFixLen\tCleanHeadLowQual\tCleanHeadFixLen\tCleanTailAdapter\tCleanTailLowQual\tCleanTailFixLen"<<endl;
	of_trim_stat2<<"Pos\tHeadLowQual\tHeadFixLen\tTailAdapter\tTailLowQual\tTailFixLen\tCleanHeadLowQual\tCleanHeadFixLen\tCleanTailAdapter\tCleanTailLowQual\tCleanTailFixLen"<<endl;
	long long head_total1(0),tail_total1(0),head_total2(0),tail_total2(0);
	long long head_total_clean1(0),tail_total_clean1(0),head_total_clean2(0),tail_total_clean2(0);
	for(int i=0;i<gv.raw1_stat.gs.read_length;i++){
		head_total1+=gv.raw1_stat.ts.ht[i]+gv.raw1_stat.ts.hlq[i];
		tail_total1+=gv.raw1_stat.ts.ta[i]+gv.raw1_stat.ts.tlq[i]+gv.raw1_stat.ts.tt[i];
		head_total2+=gv.raw2_stat.ts.ht[i]+gv.raw2_stat.ts.hlq[i];
		tail_total2+=gv.raw2_stat.ts.ta[i]+gv.raw2_stat.ts.tlq[i]+gv.raw2_stat.ts.tt[i];
		head_total_clean1+=gv.clean1_stat.ts.ht[i]+gv.clean1_stat.ts.hlq[i];
		tail_total_clean1+=gv.clean1_stat.ts.ta[i]+gv.clean1_stat.ts.tlq[i]+gv.clean1_stat.ts.tt[i];
		head_total_clean2+=gv.clean2_stat.ts.ht[i]+gv.clean2_stat.ts.hlq[i];
		tail_total_clean2+=gv.clean2_stat.ts.ta[i]+gv.clean2_stat.ts.tlq[i]+gv.clean2_stat.ts.tt[i];
		//of_trim_stat1<<i+1<<"\t"<<gv.trim1_stat.ts.hlq[i]<<"\t"<<gv.trim1_stat.ts.ht[i]<<"\t"<<gv.trim1_stat.ts.ta[i]
	}
	//cout<<head_total_clean1<<"\t"<<tail_total_clean1<<"\t"<<head_total_clean2<<"\t"<<tail_total_clean2<<endl;
	for(int i=1;i<=gv.raw1_stat.gs.read_length;i++){
		of_trim_stat1<<i<<"\t";
		if(head_total1>0){
			of_trim_stat1<<gv.raw1_stat.ts.hlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw1_stat.ts.hlq[i]/head_total1<<"%\t";
			of_trim_stat1<<gv.raw1_stat.ts.ht[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw1_stat.ts.ht[i]/head_total1<<"%\t";
		}else{
			of_trim_stat1<<gv.raw1_stat.ts.hlq[i]<<"\t0.00%\t";
			of_trim_stat1<<gv.raw1_stat.ts.ht[i]<<"\t0.00%\t";
		}
		if(tail_total1>0){
			of_trim_stat1<<gv.raw1_stat.ts.ta[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw1_stat.ts.ta[i]/tail_total1<<"%\t";
			of_trim_stat1<<gv.raw1_stat.ts.tlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw1_stat.ts.tlq[i]/tail_total1<<"%\t";
			of_trim_stat1<<gv.raw1_stat.ts.tt[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw1_stat.ts.tt[i]/tail_total1<<"%\t";
		}else{
			of_trim_stat1<<gv.raw1_stat.ts.ta[i]<<"\t0.00%\t";
			of_trim_stat1<<gv.raw1_stat.ts.tlq[i]<<"\t0.00%\t";
			of_trim_stat1<<gv.raw1_stat.ts.tlq[i]<<"\t0.00%\t";
		}
		if(head_total_clean1>0){
			of_trim_stat1<<gv.clean1_stat.ts.hlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean1_stat.ts.hlq[i]/head_total_clean1<<"%\t";
			of_trim_stat1<<gv.clean1_stat.ts.ht[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean1_stat.ts.ht[i]/head_total_clean1<<"%\t";
		}else{
			of_trim_stat1<<gv.clean1_stat.ts.hlq[i]<<"\t0.00%\t";
			of_trim_stat1<<gv.clean1_stat.ts.ht[i]<<"\t0.00%\t";
		}
		if(tail_total_clean1>0){
			of_trim_stat1<<gv.clean1_stat.ts.ta[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean1_stat.ts.ta[i]/tail_total_clean1<<"%\t";
			of_trim_stat1<<gv.clean1_stat.ts.tlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean1_stat.ts.tlq[i]/tail_total_clean1<<"%\t";
			of_trim_stat1<<gv.clean1_stat.ts.tt[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean1_stat.ts.tt[i]/tail_total_clean1<<"%"<<endl;
		}else{
			of_trim_stat1<<gv.clean1_stat.ts.ta[i]<<"\t0.00%\t";
			of_trim_stat1<<gv.clean1_stat.ts.tlq[i]<<"\t0.00%\t";
			of_trim_stat1<<gv.clean1_stat.ts.tlq[i]<<"\t0.00%"<<endl;
		}
		
		of_trim_stat2<<i<<"\t";
		if(head_total2>0){
			of_trim_stat2<<gv.raw2_stat.ts.hlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw2_stat.ts.hlq[i]/head_total2<<"%\t";
			of_trim_stat2<<gv.raw2_stat.ts.ht[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw2_stat.ts.ht[i]/head_total2<<"%\t";
		}else{
			of_trim_stat2<<gv.raw2_stat.ts.hlq[i]<<"\t0.00%\t";
			of_trim_stat2<<gv.raw2_stat.ts.ht[i]<<"\t0.00%\t";
		}
		if(tail_total2>0){
			of_trim_stat2<<gv.raw2_stat.ts.ta[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw2_stat.ts.ta[i]/tail_total2<<"%\t";
			of_trim_stat2<<gv.raw2_stat.ts.tlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw2_stat.ts.tlq[i]/tail_total2<<"%\t";
			of_trim_stat2<<gv.raw2_stat.ts.tt[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.raw2_stat.ts.tt[i]/tail_total2<<"%\t";
		}else{
			of_trim_stat2<<gv.raw2_stat.ts.ta[i]<<"\t0.00%\t";
			of_trim_stat2<<gv.raw2_stat.ts.tlq[i]<<"\t0.00%\t";
			of_trim_stat2<<gv.raw2_stat.ts.tlq[i]<<"\t0.00%\t";
		}
		if(head_total_clean2>0){
			of_trim_stat2<<gv.clean2_stat.ts.hlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean2_stat.ts.hlq[i]/head_total_clean2<<"%\t";
			of_trim_stat2<<gv.clean2_stat.ts.ht[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean2_stat.ts.ht[i]/head_total_clean2<<"%\t";
		}else{
			of_trim_stat2<<gv.clean2_stat.ts.hlq[i]<<"\t0.00%\t";
			of_trim_stat2<<gv.clean2_stat.ts.ht[i]<<"\t0.00%\t";
		}
		if(tail_total_clean2>0){
			of_trim_stat2<<gv.clean2_stat.ts.ta[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean2_stat.ts.ta[i]/tail_total_clean2<<"%\t";
			of_trim_stat2<<gv.clean2_stat.ts.tlq[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean2_stat.ts.tlq[i]/tail_total_clean2<<"%\t";
			of_trim_stat2<<gv.clean2_stat.ts.tt[i]<<"\t"<<setiosflags(ios::fixed)<<setprecision(2)<<100*(float)gv.clean2_stat.ts.tt[i]/tail_total_clean2<<"%"<<endl;
		}else{
			of_trim_stat2<<gv.clean2_stat.ts.ta[i]<<"\t0.00%\t";
			of_trim_stat2<<gv.clean2_stat.ts.tlq[i]<<"\t0.00%\t";
			of_trim_stat2<<gv.clean2_stat.ts.tlq[i]<<"\t0.00%"<<endl;
		}
	}
	of_trim_stat1.close();
	of_trim_stat2.close();
}
void peProcess::update_stat(C_fastq_file_stat& fq1s_stat,C_fastq_file_stat& fq2s_stat,C_filter_stat& fs_stat,string type){	//update statistic information from each thread
	if(type=="raw"){
		if(gv.raw1_stat.gs.read_length==0)
			gv.raw1_stat.gs.read_length=fq1s_stat.gs.read_length;	//generate stat
		gv.raw1_stat.gs.reads_number+=fq1s_stat.gs.reads_number;
		gv.raw1_stat.gs.base_number+=fq1s_stat.gs.base_number;
		gv.raw1_stat.gs.a_number+=fq1s_stat.gs.a_number;
		gv.raw1_stat.gs.c_number+=fq1s_stat.gs.c_number;
		gv.raw1_stat.gs.g_number+=fq1s_stat.gs.g_number;
		gv.raw1_stat.gs.t_number+=fq1s_stat.gs.t_number;
		gv.raw1_stat.gs.n_number+=fq1s_stat.gs.n_number;
		gv.raw1_stat.gs.q20_num+=fq1s_stat.gs.q20_num;
		gv.raw1_stat.gs.q30_num+=fq1s_stat.gs.q30_num;
		
		if(gv.raw2_stat.gs.read_length==0)
			gv.raw2_stat.gs.read_length=fq2s_stat.gs.read_length;
		gv.raw2_stat.gs.reads_number+=fq2s_stat.gs.reads_number;
		gv.raw2_stat.gs.base_number+=fq2s_stat.gs.base_number;
		gv.raw2_stat.gs.a_number+=fq2s_stat.gs.a_number;
		gv.raw2_stat.gs.c_number+=fq2s_stat.gs.c_number;
		gv.raw2_stat.gs.g_number+=fq2s_stat.gs.g_number;
		gv.raw2_stat.gs.t_number+=fq2s_stat.gs.t_number;
		gv.raw2_stat.gs.n_number+=fq2s_stat.gs.n_number;
		gv.raw2_stat.gs.q20_num+=fq2s_stat.gs.q20_num;
		gv.raw2_stat.gs.q30_num+=fq2s_stat.gs.q30_num;

		string base_set="ACGTN";	//base content and quality along read position stat
		int max_qual=0;
		for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
			for(int j=0;j!=base_set.size();j++){
				gv.raw1_stat.bs.position_acgt_content[i][j]+=fq1s_stat.bs.position_acgt_content[i][j];
				gv.raw2_stat.bs.position_acgt_content[i][j]+=fq2s_stat.bs.position_acgt_content[i][j];
			}
		}
		for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
			gv.raw1_stat.ts.ht[i]+=fq1s_stat.ts.ht[i];
			gv.raw1_stat.ts.hlq[i]+=fq1s_stat.ts.hlq[i];
			gv.raw1_stat.ts.tt[i]+=fq1s_stat.ts.tt[i];
			gv.raw1_stat.ts.tlq[i]+=fq1s_stat.ts.tlq[i];
			gv.raw1_stat.ts.ta[i]+=fq1s_stat.ts.ta[i];
			gv.raw2_stat.ts.ht[i]+=fq2s_stat.ts.ht[i];
			gv.raw2_stat.ts.hlq[i]+=fq2s_stat.ts.hlq[i];
			gv.raw2_stat.ts.tt[i]+=fq2s_stat.ts.tt[i];
			gv.raw2_stat.ts.tlq[i]+=fq2s_stat.ts.tlq[i];
			gv.raw2_stat.ts.ta[i]+=fq2s_stat.ts.ta[i];
		}
		for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
			for(int j=1;j<=MAX_QUAL;j++){
				if(fq1s_stat.qs.position_qual[i][j]>0)
					max_qual=j;
			}
		}
		
		for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
			for(int j=0;j<=max_qual;j++){
				gv.raw1_stat.qs.position_qual[i][j]+=fq1s_stat.qs.position_qual[i][j];
				gv.raw2_stat.qs.position_qual[i][j]+=fq2s_stat.qs.position_qual[i][j];
			}
		}
		//gv.fs.output_reads_num+=fs_stat.output_reads_num;	//filter stat
		gv.fs.in_adapter_list_num+=fs_stat.in_adapter_list_num;
		gv.fs.include_adapter_seq_num+=fs_stat.include_adapter_seq_num;
		gv.fs.include_adapter_seq_num1+=fs_stat.include_adapter_seq_num1;
		gv.fs.include_adapter_seq_num2+=fs_stat.include_adapter_seq_num2;
		gv.fs.include_adapter_seq_num_overlap+=fs_stat.include_adapter_seq_num_overlap;

		gv.fs.n_ratio_num+=fs_stat.n_ratio_num;
		gv.fs.n_ratio_num1+=fs_stat.n_ratio_num1;
		gv.fs.n_ratio_num2+=fs_stat.n_ratio_num2;
		gv.fs.n_ratio_num_overlap+=fs_stat.n_ratio_num_overlap;

		gv.fs.highA_num+=fs_stat.highA_num;
		gv.fs.highA_num1+=fs_stat.highA_num1;
		gv.fs.highA_num2+=fs_stat.highA_num2;
		gv.fs.highA_num_overlap+=fs_stat.highA_num_overlap;

		gv.fs.polyX_num+=fs_stat.polyX_num;
		gv.fs.polyX_num1+=fs_stat.polyX_num1;
		gv.fs.polyX_num2+=fs_stat.polyX_num2;
		gv.fs.polyX_num_overlap+=fs_stat.polyX_num_overlap;

		gv.fs.tile_num+=fs_stat.tile_num;
		gv.fs.fov_num+=fs_stat.fov_num;
		gv.fs.over_lapped_num+=fs_stat.over_lapped_num;

		gv.fs.low_qual_base_ratio_num+=fs_stat.low_qual_base_ratio_num;
		gv.fs.low_qual_base_ratio_num1+=fs_stat.low_qual_base_ratio_num1;
		gv.fs.low_qual_base_ratio_num2+=fs_stat.low_qual_base_ratio_num2;
		gv.fs.low_qual_base_ratio_num_overlap+=fs_stat.low_qual_base_ratio_num_overlap;

		gv.fs.mean_quality_num+=fs_stat.mean_quality_num;
		gv.fs.mean_quality_num1+=fs_stat.mean_quality_num1;
		gv.fs.mean_quality_num2+=fs_stat.mean_quality_num2;
		gv.fs.mean_quality_num_overlap+=fs_stat.mean_quality_num_overlap;

		gv.fs.short_len_num+=fs_stat.short_len_num;
		gv.fs.short_len_num1+=fs_stat.short_len_num1;
		gv.fs.short_len_num2+=fs_stat.short_len_num2;
		gv.fs.short_len_num_overlap+=fs_stat.short_len_num_overlap;

		gv.fs.long_len_num+=fs_stat.long_len_num;
		gv.fs.long_len_num1+=fs_stat.long_len_num1;
		gv.fs.long_len_num2+=fs_stat.long_len_num2;
		gv.fs.long_len_num_overlap+=fs_stat.long_len_num_overlap;

		gv.fs.include_contam_seq_num+=fs_stat.include_contam_seq_num;
		gv.fs.include_contam_seq_num1+=fs_stat.include_contam_seq_num1;
		gv.fs.include_contam_seq_num2+=fs_stat.include_contam_seq_num2;
		gv.fs.include_contam_seq_num_overlap+=fs_stat.include_contam_seq_num_overlap;

		gv.fs.include_global_contam_seq_num+=fs_stat.include_global_contam_seq_num;
		gv.fs.include_global_contam_seq_num1+=fs_stat.include_global_contam_seq_num1;
		gv.fs.include_global_contam_seq_num2+=fs_stat.include_global_contam_seq_num2;
		gv.fs.include_global_contam_seq_num_overlap+=fs_stat.include_global_contam_seq_num_overlap;

		
	}else if(type=="trim"){
			//generate stat
		gv.trim1_stat.gs.reads_number+=fq1s_stat.gs.reads_number;
		gv.trim1_stat.gs.base_number+=fq1s_stat.gs.base_number;
		if(gv.trim1_stat.gs.reads_number==0){
			gv.trim1_stat.gs.read_length=fq1s_stat.gs.read_length;
		}else{
			gv.trim1_stat.gs.read_length=gv.trim1_stat.gs.base_number/gv.trim1_stat.gs.reads_number;
		}
		if(gv.trim1_stat.gs.read_max_length<fq1s_stat.gs.read_length){
			gv.trim1_stat.gs.read_max_length=fq1s_stat.gs.read_length;
		}
		gv.trim1_stat.gs.a_number+=fq1s_stat.gs.a_number;
		gv.trim1_stat.gs.c_number+=fq1s_stat.gs.c_number;
		gv.trim1_stat.gs.g_number+=fq1s_stat.gs.g_number;
		gv.trim1_stat.gs.t_number+=fq1s_stat.gs.t_number;
		gv.trim1_stat.gs.n_number+=fq1s_stat.gs.n_number;
		gv.trim1_stat.gs.q20_num+=fq1s_stat.gs.q20_num;
		gv.trim1_stat.gs.q30_num+=fq1s_stat.gs.q30_num;
		
		
		//gv.trim2_stat.gs.read_length=fq2s_stat.gs.read_length;
		gv.trim2_stat.gs.reads_number+=fq2s_stat.gs.reads_number;
		gv.trim2_stat.gs.base_number+=fq2s_stat.gs.base_number;
		if(gv.trim2_stat.gs.reads_number==0){
			gv.trim2_stat.gs.read_length=fq2s_stat.gs.read_length;
		}else{
			gv.trim2_stat.gs.read_length=gv.trim2_stat.gs.base_number/gv.trim2_stat.gs.reads_number;
		}
		if(gv.trim2_stat.gs.read_max_length<fq2s_stat.gs.read_length){
			gv.trim2_stat.gs.read_max_length=fq2s_stat.gs.read_length;
		}
		gv.trim2_stat.gs.a_number+=fq2s_stat.gs.a_number;
		gv.trim2_stat.gs.c_number+=fq2s_stat.gs.c_number;
		gv.trim2_stat.gs.g_number+=fq2s_stat.gs.g_number;
		gv.trim2_stat.gs.t_number+=fq2s_stat.gs.t_number;
		gv.trim2_stat.gs.n_number+=fq2s_stat.gs.n_number;
		gv.trim2_stat.gs.q20_num+=fq2s_stat.gs.q20_num;
		gv.trim2_stat.gs.q30_num+=fq2s_stat.gs.q30_num;
		int max_qual(0);
		string base_set="ACGTN";	//base content and quality along read position stat
		for(int i=0;i!=gv.trim1_stat.gs.read_max_length;i++){
			for(int j=0;j!=base_set.size();j++){
				gv.trim1_stat.bs.position_acgt_content[i][j]+=fq1s_stat.bs.position_acgt_content[i][j];
				gv.trim2_stat.bs.position_acgt_content[i][j]+=fq2s_stat.bs.position_acgt_content[i][j];
			}
		}
		for(int i=0;i!=gv.trim1_stat.gs.read_max_length;i++){
			for(int j=1;j<=MAX_QUAL;j++){
				if(fq1s_stat.qs.position_qual[i][j]>0)
					max_qual=j;
			}
		}
		for(int i=0;i!=gv.trim1_stat.gs.read_max_length;i++){
			for(int j=0;j<=max_qual;j++){
				gv.trim1_stat.qs.position_qual[i][j]+=fq1s_stat.qs.position_qual[i][j];
				gv.trim2_stat.qs.position_qual[i][j]+=fq2s_stat.qs.position_qual[i][j];
			}
		}
	}else if(type=="clean"){
		//gp.output_reads_num+=fq1s_stat.gs.reads_number;
		//cout<<gp.output_reads_num<<"\t"<<gp.total_reads_num<<endl;
		//if(gv.clean1_stat.gs.read_length==0)
		//	gv.clean1_stat.gs.read_length=fq1s_stat.gs.read_length;	//generate stat
		gv.clean1_stat.gs.reads_number+=fq1s_stat.gs.reads_number;
		gv.clean1_stat.gs.base_number+=fq1s_stat.gs.base_number;
		if(gv.clean1_stat.gs.reads_number==0){
			gv.clean1_stat.gs.read_length=fq1s_stat.gs.read_length;
		}else{
			gv.clean1_stat.gs.read_length=gv.clean1_stat.gs.base_number/gv.clean1_stat.gs.reads_number;
		}
		if(gv.clean1_stat.gs.read_max_length<fq1s_stat.gs.read_length){
			gv.clean1_stat.gs.read_max_length=fq1s_stat.gs.read_length;
		}
		gv.clean1_stat.gs.a_number+=fq1s_stat.gs.a_number;
		gv.clean1_stat.gs.c_number+=fq1s_stat.gs.c_number;
		gv.clean1_stat.gs.g_number+=fq1s_stat.gs.g_number;
		gv.clean1_stat.gs.t_number+=fq1s_stat.gs.t_number;
		gv.clean1_stat.gs.n_number+=fq1s_stat.gs.n_number;
		gv.clean1_stat.gs.q20_num+=fq1s_stat.gs.q20_num;
		gv.clean1_stat.gs.q30_num+=fq1s_stat.gs.q30_num;
		
		//if(gv.clean2_stat.gs.read_length==0)
		//	gv.clean2_stat.gs.read_length=fq2s_stat.gs.read_length;
		gv.clean2_stat.gs.reads_number+=fq2s_stat.gs.reads_number;
		gv.clean2_stat.gs.base_number+=fq2s_stat.gs.base_number;
		if(gv.clean2_stat.gs.reads_number==0){
			gv.clean2_stat.gs.read_length=fq2s_stat.gs.read_length;
		}else{
			gv.clean2_stat.gs.read_length=gv.clean2_stat.gs.base_number/gv.clean2_stat.gs.reads_number;
		}
		if(gv.clean2_stat.gs.read_max_length<gv.clean2_stat.gs.read_length){
			gv.clean2_stat.gs.read_max_length=gv.clean2_stat.gs.read_length;
		}
		gv.clean2_stat.gs.a_number+=fq2s_stat.gs.a_number;
		gv.clean2_stat.gs.c_number+=fq2s_stat.gs.c_number;
		gv.clean2_stat.gs.g_number+=fq2s_stat.gs.g_number;
		gv.clean2_stat.gs.t_number+=fq2s_stat.gs.t_number;
		gv.clean2_stat.gs.n_number+=fq2s_stat.gs.n_number;
		gv.clean2_stat.gs.q20_num+=fq2s_stat.gs.q20_num;
		gv.clean2_stat.gs.q30_num+=fq2s_stat.gs.q30_num;
		string base_set="ACGTN";	//base content and quality along read position stat
		int max_qual(0);
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
			for(int j=0;j!=base_set.size();j++){
				gv.clean1_stat.bs.position_acgt_content[i][j]+=fq1s_stat.bs.position_acgt_content[i][j];
				gv.clean2_stat.bs.position_acgt_content[i][j]+=fq2s_stat.bs.position_acgt_content[i][j];
			}
		}
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
			gv.clean1_stat.ts.ht[i]+=fq1s_stat.ts.ht[i];
			gv.clean1_stat.ts.hlq[i]+=fq1s_stat.ts.hlq[i];
			gv.clean1_stat.ts.tt[i]+=fq1s_stat.ts.tt[i];
			gv.clean1_stat.ts.tlq[i]+=fq1s_stat.ts.tlq[i];
			gv.clean1_stat.ts.ta[i]+=fq1s_stat.ts.ta[i];
			gv.clean2_stat.ts.ht[i]+=fq2s_stat.ts.ht[i];
			gv.clean2_stat.ts.hlq[i]+=fq2s_stat.ts.hlq[i];
			gv.clean2_stat.ts.tt[i]+=fq2s_stat.ts.tt[i];
			gv.clean2_stat.ts.tlq[i]+=fq2s_stat.ts.tlq[i];
			gv.clean2_stat.ts.ta[i]+=fq2s_stat.ts.ta[i];
		}
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
			for(int j=1;j<=MAX_QUAL;j++){
				if(fq1s_stat.qs.position_qual[i][j]>0)
					max_qual=j;
			}
		}
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
			for(int j=0;j<=max_qual;j++){
				gv.clean1_stat.qs.position_qual[i][j]+=fq1s_stat.qs.position_qual[i][j];
				gv.clean2_stat.qs.position_qual[i][j]+=fq2s_stat.qs.position_qual[i][j];
			}
		}
	}else{
		cerr<<"Error:code error"<<endl;
		exit(1);
	}
	

}
void* peProcess::stat_pe_fqs(PEstatOption opt){	//statistic the pair-ends fastq
	opt.stat1->gs.reads_number+=opt.fq1s->size();
	opt.stat2->gs.reads_number+=opt.fq2s->size();
	for(vector<C_fastq>::iterator ix=opt.fq1s->begin();ix!=opt.fq1s->end();ix++){
		if((*ix).head_hdcut>0 || (*ix).head_lqcut>0){
			if((*ix).head_hdcut>=(*ix).head_lqcut){
				opt.stat1->ts.ht[(*ix).head_hdcut]++;
			}else{
				opt.stat1->ts.hlq[(*ix).head_lqcut]++;
			}
		}
		if((*ix).tail_hdcut>0 || (*ix).tail_lqcut>0 || (*ix).adacut_pos>0){
			if((*ix).tail_hdcut>=(*ix).tail_lqcut){
				if((*ix).tail_hdcut>=(*ix).adacut_pos){
					opt.stat1->ts.tt[(*ix).raw_length-(*ix).tail_hdcut+1]++;
				}else{

					opt.stat1->ts.ta[(*ix).raw_length-(*ix).adacut_pos+1]++;
				}
			}else{
				if((*ix).tail_lqcut>=(*ix).adacut_pos){
					opt.stat1->ts.tlq[(*ix).raw_length-(*ix).tail_lqcut+1]++;
				}else{
					opt.stat1->ts.ta[(*ix).raw_length-(*ix).adacut_pos+1]++;
				}
			}
		}
		for(string::size_type i=0;i!=(*ix).sequence.size();i++){
			switch(((*ix).sequence)[i]){
				case 'a':
				case 'A':opt.stat1->bs.position_acgt_content[i][0]++;opt.stat1->gs.a_number++;break;
				case 'c':
				case 'C':opt.stat1->bs.position_acgt_content[i][1]++;opt.stat1->gs.c_number++;break;
				case 'g':
				case 'G':opt.stat1->bs.position_acgt_content[i][2]++;opt.stat1->gs.g_number++;break;
				case 't':
				case 'T':opt.stat1->bs.position_acgt_content[i][3]++;opt.stat1->gs.t_number++;break;
				case 'n':
				case 'N':opt.stat1->bs.position_acgt_content[i][4]++;opt.stat1->gs.n_number++;break;
				default:{
					cerr<<"Error:unrecognized sequence,"<<((*ix).sequence)<<endl;
					exit(1);
				}
			}
		}
		for(string::size_type i=0;i!=(*ix).qual_seq.size();i++){	//process quality sequence
			int base_quality=((*ix).qual_seq)[i]-gp.qualityPhred;
			opt.stat1->qs.position_qual[i][base_quality]++;
			if(base_quality>=20)
				opt.stat1->gs.q20_num++;
			if(base_quality>=30)
				opt.stat1->gs.q30_num++;
		}
		opt.stat1->gs.read_length=(*ix).sequence.size();
		opt.stat1->gs.base_number+=opt.stat1->gs.read_length;
	}
	
	for(vector<C_fastq>::iterator ix=opt.fq2s->begin();ix!=opt.fq2s->end();ix++){
		if((*ix).head_hdcut>0 || (*ix).head_lqcut>0){
			if((*ix).head_hdcut>=(*ix).head_lqcut){
				opt.stat2->ts.ht[(*ix).head_hdcut]++;
			}else{
				opt.stat2->ts.hlq[(*ix).head_lqcut]++;
			}
		}
		if((*ix).tail_hdcut>0 || (*ix).tail_lqcut>0 || (*ix).adacut_pos>0){
			if((*ix).tail_hdcut>=(*ix).tail_lqcut){
				if((*ix).tail_hdcut>=(*ix).adacut_pos){
					opt.stat2->ts.tt[(*ix).sequence.size()-(*ix).tail_hdcut+1]++;
				}else{
					opt.stat2->ts.ta[(*ix).sequence.size()-(*ix).adacut_pos+1]++;
				}
			}else{
				if((*ix).tail_lqcut>=(*ix).adacut_pos){
					opt.stat2->ts.tlq[(*ix).sequence.size()-(*ix).tail_lqcut+1]++;
				}else{
					opt.stat2->ts.ta[(*ix).sequence.size()-(*ix).adacut_pos+1]++;
				}
			}
		}
		for(string::size_type i=0;i!=(*ix).sequence.size();i++){
			switch(((*ix).sequence)[i]){
				case 'a':
				case 'A':opt.stat2->bs.position_acgt_content[i][0]++;opt.stat2->gs.a_number++;break;
				case 'c':
				case 'C':opt.stat2->bs.position_acgt_content[i][1]++;opt.stat2->gs.c_number++;break;
				case 'g':
				case 'G':opt.stat2->bs.position_acgt_content[i][2]++;opt.stat2->gs.g_number++;break;
				case 't':
				case 'T':opt.stat2->bs.position_acgt_content[i][3]++;opt.stat2->gs.t_number++;break;
				case 'n':
				case 'N':opt.stat2->bs.position_acgt_content[i][4]++;opt.stat2->gs.n_number++;break;
				default:{
					cerr<<"Error:unrecognized sequence,"<<((*ix).sequence)<<endl;
					exit(1);
				}
			}
		}
		for(string::size_type i=0;i!=(*ix).qual_seq.size();i++){	//process quality sequence
			int base_quality=((*ix).qual_seq)[i]-gp.qualityPhred;
			if(base_quality>MAX_QUAL){
				cerr<<"Error:quality is too high,please check the quality system parameter or fastq file"<<endl;
				exit(1);
			}
			if(base_quality<MIN_QUAL){
				cerr<<"Error:quality is too low,please check the quality system parameter or fastq file"<<endl;
				exit(1);
			}
			opt.stat2->qs.position_qual[i][base_quality]++;
			if(base_quality>=20)
				opt.stat2->gs.q20_num++;
			if(base_quality>=30)
				opt.stat2->gs.q30_num++;
		}
		opt.stat2->gs.read_length=(*ix).sequence.size();
		opt.stat2->gs.base_number+=opt.stat2->gs.read_length;
	}
	
}
void peProcess::filter_pe_fqs(PEcalOption opt){
	//C_reads_trim_stat_2 cut_pos;
	vector<C_fastq>::iterator i2=opt.fq2s->begin();
	for(vector<C_fastq>::iterator i=opt.fq1s->begin();i!=opt.fq1s->end();i++){
		C_pe_fastq_filter pe_fastq_filter=C_pe_fastq_filter(*i,*i2,gp);
		/*int head_hdcut,head_lqcut,tail_hdcut,tail_lqcut,adacut_pos;
	int contam_pos;
	int global_contam_pos;
	int raw_length;*/
		pe_fastq_filter.pe_trim(gp);
		if(gp.adapter_discard_or_trim=="trim" || gp.contam_discard_or_trim=="trim" || !gp.trim.empty() || !gp.trimBadHead.empty() || !gp.trimBadTail.empty()){
			(*i).head_hdcut=pe_fastq_filter.fq1.head_hdcut;
			(*i).head_lqcut=pe_fastq_filter.fq1.head_lqcut;
			(*i).tail_hdcut=pe_fastq_filter.fq1.tail_hdcut;
			(*i).tail_lqcut=pe_fastq_filter.fq1.tail_lqcut;
			(*i).adacut_pos=pe_fastq_filter.fq1.adacut_pos;
			//(*i).contam_pos=pe_fastq_filter.fq1.contam_pos;
			//(*i).global_contam_pos=pe_fastq_filter.fq1.global_contam_pos;
			//(*i).raw_length=pe_fastq_filter.fq1.raw_length;
			(*i2).head_hdcut=pe_fastq_filter.fq2.head_hdcut;
			(*i2).head_lqcut=pe_fastq_filter.fq2.head_lqcut;
			(*i2).tail_hdcut=pe_fastq_filter.fq2.tail_hdcut;
			(*i2).tail_lqcut=pe_fastq_filter.fq2.tail_lqcut;
			(*i2).adacut_pos=pe_fastq_filter.fq2.adacut_pos;
			//(*i2).contam_pos=pe_fastq_filter.fq2.contam_pos;
			//(*i2).global_contam_pos=pe_fastq_filter.fq2.global_contam_pos;
			//(*i2).raw_length=pe_fastq_filter.fq2.raw_length;
		}
		if(!gp.trim_fq1.empty()){
			preOutput(1,pe_fastq_filter.fq1);
			preOutput(2,pe_fastq_filter.fq2);
			opt.trim_result1->push_back(pe_fastq_filter.fq1);
			opt.trim_result2->push_back(pe_fastq_filter.fq2);
		}
		if(pe_fastq_filter.pe_discard(opt.local_fs,gp)!=1){
			if(!gp.clean_fq1.empty()){
				preOutput(1,pe_fastq_filter.fq1);
				preOutput(2,pe_fastq_filter.fq2);
				opt.clean_result1->push_back(pe_fastq_filter.fq1);
				opt.clean_result2->push_back(pe_fastq_filter.fq2);
			}
		}
		i2++;
	}
	//return cut_pos;
}

void  peProcess::preOutput(int type,C_fastq& a){	//modify the sequences before output if necessary
	if(gp.whether_add_pe_info){
		if(type==1){
			a.seq_id=a.seq_id+"/1";
		}else{
			a.seq_id=a.seq_id+"/2";
		}
	}
	if(!gp.base_convert.empty()){
		gp.base_convert=gp.base_convert.replace(gp.base_convert.find("TO"),2,"");
		gp.base_convert=gp.base_convert.replace(gp.base_convert.find("2"),1,"");
		if(gp.base_convert.size()!=2){
			cerr<<"Error:base_conver value format error"<<endl;
			exit(1);
		}
		for(string::size_type ix=0;ix!=a.sequence.size();ix++){
			if(toupper(a.sequence[ix])==toupper(gp.base_convert[0]))
				a.sequence[ix]=gp.base_convert[1];
		}
	}
}
void peProcess::peWrite_split(vector<C_fastq>& pe1,vector<C_fastq>& pe2){
	/*if(gp.threads_num==1){
		output_split_fastqs("1",pe1);
		output_split_fastqs("2",pe2);
	}else{
		thread write1(bind(&peProcess::output_split_fastqs,this,"1",pe1));
		thread write2(bind(&peProcess::output_split_fastqs,this,"2",pe2));
		write1.join();
		write2.join();
	}
	*/
	output_split_fastqs("1",pe1);
	output_split_fastqs("2",pe2);
}
void peProcess::peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,string type,gzFile out1,gzFile out2){	//output the sequences to  files
	/*if(gp.threads_num==1){
		if(gp.mode=="nonssd" && gp.output_clean>0 && type=="clean"){
			output_split_fastqs("1",pe1);
			output_split_fastqs("2",pe2);
		}else{
			output_fastqs("1",pe1,out1);
			output_fastqs("2",pe2,out2);
		}
		
	}else{
		if(gp.mode=="nonssd" && gp.output_clean>0 && type=="clean"){
			thread write1(bind(&peProcess::output_split_fastqs,this,"1",pe1));
			thread write2(bind(&peProcess::output_split_fastqs,this,"2",pe2));
			write1.join();
			write2.join();
		}else{
			thread write1(bind(&peProcess::output_fastqs,this,"1",pe1,out1));
			thread write2(bind(&peProcess::output_fastqs,this,"2",pe2,out2));
			write1.join();
			write2.join();
		}
	}
	*/
	output_fastqs("1",pe1,out1);
	output_fastqs("2",pe2,out2);
}
void peProcess::C_fastq_init(C_fastq& a,C_fastq& b){
	a.seq_id="";
	a.sequence="";
	a.qual_seq="";
	b.seq_id="";
	b.sequence="";
	b.qual_seq="";
	a.adapter_seq=gp.adapter1_seq;
	b.adapter_seq=gp.adapter2_seq;
	a.contam_seq=gp.contam1_seq;
	b.contam_seq=gp.contam2_seq;
	//a.global_contams=gp.global_contams;
	//b.global_contams=gp.global_contams;
	a.head_hdcut=-1;
	a.head_lqcut=-1;
	a.tail_hdcut=-1;
	a.tail_lqcut=-1;
	a.adacut_pos=-1;
	a.contam_pos=-1;
	a.global_contam_pos=-1;
	b.head_hdcut=-1;
	b.head_lqcut=-1;
	b.tail_hdcut=-1;
	b.tail_lqcut=-1;
	b.adacut_pos=-1;
	b.contam_pos=-1;
	b.global_contam_pos=-1;
	a.raw_length=0;
	b.raw_length=0;
	if(!gp.trim.empty()){
		vector<string> tmp_eles=get_pe_hard_trim(gp.trim);
		a.head_trim_len=tmp_eles[0];
		a.tail_trim_len=tmp_eles[1];
		b.head_trim_len=tmp_eles[2];
		b.tail_trim_len=tmp_eles[3];
	}

}
int peProcess::read(vector<C_fastq>& pe1,vector<C_fastq>& pe2,ifstream& infile1,ifstream& infile2){	//read file from ASCII text,used in ssd mode
	
	string buf1,buf2;
	int file1_line_num(0),file2_line_num(0);
	C_fastq fastq1,fastq2;
	C_fastq_init(fastq1,fastq2);
	for(int i=0;i<gp.patchSize*4;i++){
		if(getline(infile1,buf1)){
			file1_line_num++;
			//string s_line(buf1);
			//s_line.erase(s_line.size()-1);
			if(file1_line_num%4==1){
				fastq1.seq_id=buf1;
			}
			if(file1_line_num%4==2){
				fastq1.sequence=buf1;
			}
			if(file1_line_num%4==0){
				fastq1.qual_seq=buf1;
				//fq1s.push_back(fastq1);
			}
		}
		
		if(getline(infile2,buf2)){
			file2_line_num++;
			//string s_line(buf2);
			if(file2_line_num%4==1){
				fastq2.seq_id=buf2;
			}
			if(file2_line_num%4==2){
				fastq2.sequence=buf2;
			}
			if(file2_line_num%4==0){
				fastq2.qual_seq=buf2;
				pe1.push_back(fastq1);
				pe2.push_back(fastq2);
			}
		}else{
			return -1;
		}
	}
	return 0;
}
int peProcess::read_gz(vector<C_fastq>& pe1,vector<C_fastq>& pe2){	//read file from gz format,used in non-ssd mode
	char buf1[READBUF],buf2[READBUF];
	C_fastq fastq1,fastq2;
	C_fastq_init(fastq1,fastq2);
	int file1_line_num(0),file2_line_num(0);
	for(int i=0;i<gp.patchSize*4;i++){
		if(gzgets(gzfp1,buf1,READBUF)!=NULL){
			file1_line_num++;
			if(file1_line_num%4==1){
				fastq1.seq_id.assign(buf1);
				fastq1.seq_id.erase(fastq1.seq_id.size()-1);
			}
			if(file1_line_num%4==2){
				fastq1.sequence.assign(buf1);
				fastq1.sequence.erase(fastq1.sequence.size()-1);
			}
			if(file1_line_num%4==0){
				fastq1.qual_seq.assign(buf1);
				fastq1.qual_seq.erase(fastq1.qual_seq.size()-1);
			}
		}
		
		if(gzgets(gzfp2,buf2,READBUF)!=NULL){
			file2_line_num++;
			if(file2_line_num%4==1){
				fastq2.seq_id.assign(buf2);
				fastq2.seq_id.erase(fastq2.seq_id.size()-1);
			}
			if(file2_line_num%4==2){
				fastq2.sequence.assign(buf2);
				fastq2.sequence.erase(fastq2.sequence.size()-1);
			}
			if(file2_line_num%4==0){
				fastq2.qual_seq.assign(buf2);
				fastq2.qual_seq.erase(fastq2.qual_seq.size()-1);
				pe1.push_back(fastq1);
				pe2.push_back(fastq2);
			}
		}else{
			return -1;
		}
	}
	return 0;
}
void peProcess::create_thread_trimoutputFile(int index){
	if(!gp.trim_fq1.empty()){	//create output trim files handle
		ostringstream trim_outfile1,trim_outfile2;
		trim_outfile1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".trim.r1.fq.gz";
		trim_outfile2<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".trim.r2.fq.gz";
		gz_trim_out1[index]=gzopen(trim_outfile1.str().c_str(),"wb");
		gzsetparams(gz_trim_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1[index],1024*1024*8);
		gz_trim_out2[index]=gzopen(trim_outfile2.str().c_str(),"wb");
		gzsetparams(gz_trim_out2[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out2[index],1024*1024*8);
	}
}
void peProcess::create_thread_cleanoutputFile(int index){
	if(gp.output_clean<=0){
		if(!gp.clean_fq1.empty()){	//create output clean files handle
			ostringstream outfile1,outfile2;
			outfile1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r1.fq.gz";
			outfile2<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r2.fq.gz";
			gz_clean_out1[index]=gzopen(outfile1.str().c_str(),"wb");
			if(!gz_clean_out1[index]){
				cerr<<"Error:cannot write to the file,"<<outfile1.str()<<endl;
				exit(1);
			}
			gzsetparams(gz_clean_out1[index], 2, Z_DEFAULT_STRATEGY);
			gzbuffer(gz_clean_out1[index],1024*1024*8);
			gz_clean_out2[index]=gzopen(outfile2.str().c_str(),"wb");
			gzsetparams(gz_clean_out2[index], 2, Z_DEFAULT_STRATEGY);
			gzbuffer(gz_clean_out2[index],1024*1024*8);
		}
	}
}
/*
void peProcess::add_raw_trim(C_fastq_file_stat& a,C_fastq_file_stat& a2,C_reads_trim_stat& b,C_reads_trim_stat& b2){
	for(int i=1;i<=a.gs.read_length;i++){
		a.ts.hlq[i]+=b.hlq[i];
		a.ts.ht[i]+=b.ht[i];
		a.ts.ta[i]+=b.ta[i];
		a.ts.tlq[i]+=b.tlq[i];
		a.ts.tt[i]+=b.tt[i];
		a2.ts.hlq[i]+=b2.hlq[i];
		a2.ts.ht[i]+=b2.ht[i];
		a2.ts.ta[i]+=b2.ta[i];
		a2.ts.tlq[i]+=b2.tlq[i];
		a2.ts.tt[i]+=b2.tt[i];
	}
}
*/
void peProcess::thread_process_reads(int index,vector<C_fastq> &fq1s,vector<C_fastq> &fq2s){
	vector<C_fastq> trim_result1,trim_result2,clean_result1,clean_result2;
	
	PEcalOption opt2;
	opt2.local_fs=&local_fs[index];
	opt2.fq1s=&fq1s;
	opt2.fq2s=&fq2s;
	opt2.trim_result1=&trim_result1;
	opt2.trim_result2=&trim_result2;
	opt2.clean_result1=&clean_result1;
	opt2.clean_result2=&clean_result2;
	filter_pe_fqs(opt2);		//filter raw fastqs by the given parameters
	PEstatOption opt_raw;
	opt_raw.fq1s=&fq1s;
	opt_raw.stat1=&local_raw_stat1[index];
	opt_raw.fq2s=&fq2s;
	opt_raw.stat2=&local_raw_stat2[index];
	stat_pe_fqs(opt_raw);		//statistic raw fastqs
	//add_raw_trim(local_raw_stat1[index],local_raw_stat2[index],raw_cut.stat1,raw_cut.stat2);
	fq1s.clear();
	fq2s.clear();
	PEstatOption opt_trim,opt_clean;
	if(!gp.trim_fq1.empty()){	//trim means only trim but not discard.
		opt_trim.fq1s=&trim_result1;
		opt_trim.stat1=&local_trim_stat1[index];
		opt_trim.fq2s=&trim_result2;
		opt_trim.stat2=&local_trim_stat2[index];
		stat_pe_fqs(opt_trim);	//statistic trim fastqs
	}

	//write_m.lock();
	if(!gp.trim_fq1.empty()){
		peWrite(trim_result1,trim_result2,"trim",gz_trim_out1[index],gz_trim_out2[index]);	//output trim files
		trim_result1.clear();
		trim_result2.clear();
	}
	if(!gp.clean_fq1.empty()){
		opt_clean.stat1=&local_clean_stat1[index];
		opt_clean.stat2=&local_clean_stat2[index];
		if(gp.output_clean>0 || (gp.total_reads_num_random==false && gp.l_total_reads_num>0)){
			write_m.lock();
			if(limit_end>0){
				write_m.unlock();
				clean_result1.clear();
				clean_result2.clear();
				return;
			}
			peWrite_split(clean_result1,clean_result2);
			if(limit_end==1){
				int to_remove=gp.have_output1-gp.clean_file_reads;
				if(to_remove>=clean_result1.size()){
					write_m.unlock();
					clean_result1.clear();
					clean_result2.clear();
					return;
				}else{
					clean_result1.erase(clean_result1.end()-to_remove,clean_result1.end());
					clean_result2.erase(clean_result2.end()-to_remove,clean_result2.end());
				}
			}
			write_m.unlock();
			//statistic clean fastqs
		}else{
			peWrite(clean_result1,clean_result2,"clean",gz_clean_out1[index],gz_clean_out2[index]);//output clean files
		}
		opt_clean.fq1s=&clean_result1;
		opt_clean.fq2s=&clean_result2;
		stat_pe_fqs(opt_clean);	
		/*thread_write_m[index].lock();
		thread pewrite_t(bind(&peProcess::peWrite,this,clean_result1,clean_result2,"clean",gz_clean_out1[index],gz_clean_out2[index]));
		thread_write_m[index].unlock();*/
		if(gp.is_streaming){
			write_m.lock();
			C_global_variable tmp_gv;
			tmp_gv.fs=*(opt2.local_fs);
			tmp_gv.raw1_stat=*(opt_raw.stat1);
			tmp_gv.raw2_stat=*(opt_raw.stat2);
			//tmp_gv.trim1_stat=local_trim_stat1[index];
			//tmp_gv.trim2_stat=local_trim_stat2[index];
			tmp_gv.clean1_stat=*(opt_clean.stat1);
			tmp_gv.clean2_stat=*(opt_clean.stat2);
			peStreaming_stat(tmp_gv);
			write_m.unlock();
		}
		
		clean_result1.clear();
		clean_result2.clear();
	}
}

void* peProcess::sub_thread_nonssd(int index){	//sub thread process in non-ssd mode 
	of_log<<get_local_time()<<"\tthread "<<index<<" start"<<endl;
	/*
	local_fs[index]=C_filter_stat();
	local_raw_stat1[index]=C_fastq_file_stat();
	local_raw_stat2[index]=C_fastq_file_stat();
	if(!gp.trim_fq1.empty()){
		local_trim_stat1[index]=C_fastq_file_stat();
		local_trim_stat2[index]=C_fastq_file_stat();
	}
	if(!gp.clean_fq1.empty()){
		local_clean_stat1[index]=C_fastq_file_stat();
		local_clean_stat2[index]=C_fastq_file_stat();
	}
	*/
	vector<C_fastq> fq1s,fq2s;
	vector<C_fastq> trim_result1,trim_result2,clean_result1,clean_result2;
	int done(0);
	while(1){
		read_m.lock();
		if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
			if(read_gz(fq1s,fq2s)==-1){		//read fastqs from raw files
				done=1;
			}
		}else{
			if(read(fq1s,fq2s,nongzfp1,nongzfp2)==-1){		
				done=1;
			}
		}
		processed_reads+=fq1s.size();
		if(processed_reads%(gp.patchSize*gp.threads_num)==0){
			of_log<<get_local_time()<<"\tprocessed reads number:"<<processed_reads<<endl;
		}
		read_m.unlock();
		PEcalOption opt2;
		opt2.local_fs=&local_fs[index];
		opt2.fq1s=&fq1s;
		opt2.fq2s=&fq2s;
		opt2.trim_result1=&trim_result1;
		opt2.trim_result2=&trim_result2;
		opt2.clean_result1=&clean_result1;
		opt2.clean_result2=&clean_result2;
		filter_pe_fqs(opt2);		//filter raw fastqs by the given parameters
		PEstatOption opt_raw;
		opt_raw.fq1s=&fq1s;
		opt_raw.stat1=&local_raw_stat1[index];
		opt_raw.fq2s=&fq2s;
		opt_raw.stat2=&local_raw_stat2[index];
		stat_pe_fqs(opt_raw);		//statistic raw fastqs
		fq1s.clear();
		fq2s.clear();
		PEstatOption opt_trim,opt_clean;
		if(!gp.trim_fq1.empty()){	//trim means only trim but not discard.
			opt_trim.fq1s=&trim_result1;
			opt_trim.stat1=&local_trim_stat1[index];
			opt_trim.fq2s=&trim_result2;
			opt_trim.stat2=&local_trim_stat2[index];
			stat_pe_fqs(opt_trim);	//statistic trim fastqs
		}
		if(!gp.clean_fq1.empty()){
			opt_clean.fq1s=&clean_result1;
			opt_clean.stat1=&local_clean_stat1[index];
			opt_clean.fq2s=&clean_result2;
			opt_clean.stat2=&local_clean_stat2[index];
			stat_pe_fqs(opt_clean);	//statistic clean fastqs
		}
		write_m.lock();
		if(!gp.trim_fq1.empty()){
			peWrite(trim_result1,trim_result2,"trim",gz_trim_out1_nonssd,gz_trim_out2_nonssd);	//output trim files
		}
		if(!gp.clean_fq1.empty()){
			peWrite(clean_result1,clean_result2,"clean",gz_clean_out1_nonssd,gz_clean_out2_nonssd);	//output clean files
			if(gp.is_streaming){
				C_global_variable tmp_gv;
				tmp_gv.fs=*(opt2.local_fs);
				tmp_gv.raw1_stat=*(opt_raw.stat1);
				tmp_gv.raw2_stat=*(opt_raw.stat2);
				//tmp_gv.trim1_stat=local_trim_stat1[index];
				//tmp_gv.trim2_stat=local_trim_stat2[index];
				tmp_gv.clean1_stat=*(opt_clean.stat1);
				tmp_gv.clean2_stat=*(opt_clean.stat2);
				peStreaming_stat(tmp_gv);
			}
		}
		write_m.unlock();
		if(!gp.trim_fq1.empty()){
			trim_result1.clear();
			trim_result2.clear();
		}
		if(!gp.clean_fq1.empty()){
			clean_result1.clear();
			clean_result2.clear();
		}
		if(done==1){
			break;
		}
		if(gp.is_streaming){
			/*
			local_fs[index]=C_filter_stat();
			local_raw_stat1[index]=C_fastq_file_stat();
			local_raw_stat2[index]=C_fastq_file_stat();
			local_clean_stat1[index]=C_fastq_file_stat();
			local_clean_stat2[index]=C_fastq_file_stat();
			*/
		}
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}
void* peProcess::sub_thread(int index){	//sub thread in ssd mode
	off_t start_pos=t_start_pos[index];
	off_t end_pos=t_end_pos[index];
	string head,tail;
	int flag(0);
	int self_fq1fd=open(new_fq1_path.c_str(),O_RDONLY);
	int self_fq2fd=open(new_fq2_path.c_str(),O_RDONLY);
	if(self_fq1fd==-1){
		cerr<<"Error:cannot open the file,"<<new_fq1_path<<endl;
		exit(1);
	}
	if(self_fq2fd==-1){
		cerr<<"Error:cannot open the file,"<<new_fq2_path<<endl;
		exit(1);
	}
	int min_len=1024*4*10;
	if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
		create_thread_trimoutputFile(index);
		create_thread_cleanoutputFile(index);
	}else if(!gp.trim_fq1.empty()){
		create_thread_trimoutputFile(index);
	}
	
	
	//char *fq1_buf,*fq2_buf;
	unsigned long long mmap_start(start_pos);
	unsigned long long copysz(buffer);
	int tmp_iter(0);
	char *buf1,*buf2;
	C_fastq fastq1,fastq2;
	C_fastq_init(fastq1,fastq2);
	vector<C_fastq> fq1s,fq2s;
	vector<C_fastq> trim_result1,trim_result2,clean_result1,clean_result2;
	while(1){
		if(mmap_start>=end_pos)
			break;
		if(end_pos - mmap_start<buffer){
			copysz=end_pos - mmap_start;
		}else{
			if(end_pos - mmap_start -buffer <buffer/4){
				copysz=end_pos - mmap_start;
			}else{
				copysz=buffer;
			}
		}
		if((buf1=(char*)mmap(0,copysz,PROT_READ,MAP_SHARED,self_fq1fd,mmap_start))==MAP_FAILED){
			cerr<<"Error:mmap error"<<endl;
			exit(1);
		}
		if((buf2=(char*)mmap(0,copysz,PROT_READ,MAP_SHARED,self_fq2fd,mmap_start))==MAP_FAILED){
			cerr<<"Error:mmap error"<<endl;
			exit(1);
		}
		string tmp_head1,tmp_tail1;
		string tmp_head2,tmp_tail2;
		int head_idx=0;
		for(;head_idx<copysz;head_idx++){
			if(head_idx==0 || head_idx==1){
				tmp_head1.insert(tmp_head1.end(),buf1[head_idx]);
				tmp_head2.insert(tmp_head2.end(),buf2[head_idx]);
			}else{
				if(buf1[head_idx]=='@' && buf1[head_idx-1]=='\n' && buf1[head_idx-2]!='+'){
					break;
				}else{
					tmp_head1.insert(tmp_head1.end(),buf1[head_idx]);
					tmp_head2.insert(tmp_head2.end(),buf2[head_idx]);
				}
			}
		}
		int tail_idx=copysz-1;
		for(;tail_idx>1;tail_idx--){
			if(tail_idx==copysz-1){
				tmp_tail1.insert(tmp_tail1.begin(),buf1[tail_idx]);
				tmp_tail2.insert(tmp_tail2.begin(),buf2[tail_idx]);
			}else{
				if(buf1[tail_idx]=='\n' && buf1[tail_idx+1]=='@' && buf1[tail_idx-1]!='+'){
					break;
				}else{
					tmp_tail1.insert(tmp_tail1.begin(),buf1[tail_idx]);
					tmp_tail2.insert(tmp_tail2.begin(),buf2[tail_idx]);
				}
			}
		}
		sticky_reads1[index]+=tmp_head1;
		sticky_reads1[index]+=tmp_tail1;
		sticky_reads2[index]+=tmp_head2;
		sticky_reads2[index]+=tmp_tail2;
		int line_num(0);
		
		
		for(int i=head_idx;i<=tail_idx;i++){
			if(buf1[i]=='\n'){
				if(line_num%4==3){
					fq1s.push_back(fastq1);
					fq2s.push_back(fastq2);
					fastq1.seq_id.clear();
					fastq1.sequence.clear();
					fastq1.qual_seq.clear();
					fastq2.seq_id.clear();
					fastq2.sequence.clear();
					fastq2.qual_seq.clear();
				}
				line_num++;
				continue;
			}
			if(line_num%4==0){
				fastq1.seq_id.insert(fastq1.seq_id.end(),buf1[i]);
				fastq2.seq_id.insert(fastq2.seq_id.end(),buf2[i]);
			}
			if(line_num%4==1){
				fastq1.sequence.insert(fastq1.sequence.end(),buf1[i]);
				fastq2.sequence.insert(fastq2.sequence.end(),buf2[i]);
			}
			if(line_num%4==3){
				fastq1.qual_seq.insert(fastq1.qual_seq.end(),buf1[i]);
				fastq2.qual_seq.insert(fastq2.qual_seq.end(),buf2[i]);
			}
		}
		thread_process_reads(index,fq1s,fq2s);
		if(limit_end>0){
			break;
		}
		mmap_start+=copysz;
		munmap(buf1,copysz);
		munmap(buf2,copysz);
		tmp_iter++;
	}
	close(self_fq1fd);
	close(self_fq2fd);
	self_fq1fd=open(new_fq1_path.c_str(),O_RDONLY);
	self_fq2fd=open(new_fq2_path.c_str(),O_RDONLY);
	if(self_fq1fd==-1){
		cerr<<"Error:cannot open the file,"<<new_fq1_path<<endl;
		exit(1);
	}
	if(self_fq2fd==-1){
		cerr<<"Error:cannot open the file,"<<new_fq2_path<<endl;
		exit(1);
	}
	close(self_fq1fd);
	close(self_fq2fd);
	if(!gp.trim_fq1.empty()){
		if(index!=0){
			gzclose(gz_trim_out1[index]);
			gzclose(gz_trim_out2[index]);
		}
		
	}
	if(!gp.clean_fq1.empty()){
		if(index!=0){
			if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
				gzclose(gz_clean_out1[index]);
				gzclose(gz_clean_out2[index]);
			}
		}
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}

void peProcess::run_pigz(int type){	//split raw files with pigz and "split" command
	ostringstream cmd1;
	int pigz_thread=gp.threads_num>2?gp.threads_num/2:1;
	if(type==1){
		cmd1<<"pigz -c -d -p "<<pigz_thread<<" "<<gp.fq1_path<<" > "<<gp.output_dir<<"/raw.r1.fq";
	}else if(type==2){
		cmd1<<"pigz -c -d -p "<<pigz_thread<<" "<<gp.fq2_path<<" > "<<gp.output_dir<<"/raw.r2.fq";
	}
	if(system(cmd1.str().c_str())==-1){
		cerr<<"Error:run pigz error"<<endl;
		exit(1);
	}
}
void peProcess::merge_stat(){
	for(int i=0;i!=gp.threads_num;i++){
		update_stat(local_raw_stat1[i],local_raw_stat2[i],local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(local_trim_stat1[i],local_trim_stat2[i],local_fs[i],"trim");
		}
		if(!gp.clean_fq1.empty()){
			update_stat(local_clean_stat1[i],local_clean_stat2[i],local_fs[i],"clean");
		}
	}
}
void peProcess::merge_stat(int index){
	for(int i=0;i<=gp.threads_num;i++){
		update_stat(local_raw_stat1[i],local_raw_stat2[i],local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(local_trim_stat1[i],local_trim_stat2[i],local_fs[i],"trim");
		}
	}
	for(int i=0;i<=index;i++){
		if(!gp.clean_fq1.empty()){
			update_stat(local_clean_stat1[i],local_clean_stat2[i],local_fs[i],"clean");
		}
	}
}
void peProcess::merge_stat_nonssd(){
	for(int i=0;i!=gp.threads_num;i++){
		update_stat(local_raw_stat1[i],local_raw_stat2[i],local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(local_trim_stat1[i],local_trim_stat2[i],local_fs[i],"trim");
		}
		if(!gp.clean_fq1.empty()){
			update_stat(local_clean_stat1[i],local_clean_stat2[i],local_fs[i],"clean");
		}
	}
}
void peProcess::run_cmd(string cmd){
	if(system(cmd.c_str())==-1){
		cerr<<"Error:when running "<<cmd<<endl;
		exit(1);
	}
}
void peProcess::merge_trim_data(){	//cat all output files generated by multi-threads to a single large file
	if(!gp.trim_fq1.empty()){
		string trim_file1,trim_file2;
		trim_file1=gp.output_dir+"/"+gp.trim_fq1;
		trim_file2=gp.output_dir+"/"+gp.trim_fq2;
		if(check_gz_empty(trim_file1)==1){
			string rm_cmd="rm "+trim_file1;
			system(rm_cmd.c_str());
		}
		if(check_gz_empty(trim_file2)==1){
			string rm_cmd="rm "+trim_file2;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream trim_out_fq1_tmp,trim_out_fq2_tmp;
			trim_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r1.fq.gz";
			trim_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r2.fq.gz";
			
			if(access(trim_out_fq1_tmp.str().c_str(),0)!=-1 && access(trim_out_fq2_tmp.str().c_str(),0)!=-1){

				ostringstream cat_cmd1,cat_cmd2;
				cat_cmd1<<"cat "<<trim_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq1<<";rm "<<trim_out_fq1_tmp.str();
				cat_cmd2<<"cat "<<trim_out_fq2_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq2<<";rm "<<trim_out_fq2_tmp.str();
				if(gp.threads_num>1){
					thread cat_t1(bind(&peProcess::run_cmd,this,cat_cmd1.str()));
					thread cat_t2(bind(&peProcess::run_cmd,this,cat_cmd2.str()));
					cat_t1.join();
					cat_t2.join();
				}else{
					run_cmd(cat_cmd1.str());
					run_cmd(cat_cmd2.str());
				}
			}
			
		}
	}
	
}
void peProcess::merge_clean_data(){	//cat all output files generated by multi-threads to a single large file
	if(!gp.clean_fq1.empty()){
		string clean_file1,clean_file2;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		clean_file2=gp.output_dir+"/"+gp.clean_fq2;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
		if(check_gz_empty(clean_file2)==1){
			string rm_cmd="rm "+clean_file2;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp,clean_out_fq2_tmp;
			clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			clean_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r2.fq.gz";
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1 && access(clean_out_fq2_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1,cat_cmd2;
				cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
				cat_cmd2<<"cat "<<clean_out_fq2_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq2<<";rm "<<clean_out_fq2_tmp.str();
				if(gp.threads_num>1){
					thread cat_t1(bind(&peProcess::run_cmd,this,cat_cmd1.str()));
					thread cat_t2(bind(&peProcess::run_cmd,this,cat_cmd2.str()));
					cat_t1.join();
					cat_t2.join();
				}else{
					run_cmd(cat_cmd1.str());
					run_cmd(cat_cmd2.str());
				}
			}
		}
	}
}
void peProcess::merge_clean_data(int index){	//cat all output files generated by multi-threads to a single large file
	if(!gp.clean_fq1.empty()){
		string clean_file1,clean_file2;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		clean_file2=gp.output_dir+"/"+gp.clean_fq2;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
		if(check_gz_empty(clean_file2)==1){
			string rm_cmd="rm "+clean_file2;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp,clean_out_fq2_tmp;
			if(i==index){
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/last.r1.fq.gz";
				clean_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/last.r2.fq.gz";
			}else{
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
				clean_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r2.fq.gz";
			}
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1 && access(clean_out_fq2_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1,cat_cmd2;
				if(i<=index){
					cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
					cat_cmd2<<"cat "<<clean_out_fq2_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq2<<";rm "<<clean_out_fq2_tmp.str();
					if(i==index){
						cat_cmd1<<";rm "<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
						cat_cmd2<<";rm "<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r2.fq.gz";
					}
				}else{
					cat_cmd1<<"rm "<<clean_out_fq1_tmp.str();
					cat_cmd2<<"rm "<<clean_out_fq2_tmp.str();
				}
				if(gp.threads_num>1){
					thread cat_t1(bind(&peProcess::run_cmd,this,cat_cmd1.str()));
					thread cat_t2(bind(&peProcess::run_cmd,this,cat_cmd2.str()));
					cat_t1.join();
					cat_t2.join();
				}else{
					run_cmd(cat_cmd1.str());
					run_cmd(cat_cmd2.str());
				}
			}
		}
	}
}
void peProcess::merge_data(int index){	//cat all output files to a single large file in limit output mode
	if(!gp.trim_fq1.empty()){
		string trim_file1,trim_file2;
		trim_file1=gp.output_dir+"/"+gp.trim_fq1;
		trim_file2=gp.output_dir+"/"+gp.trim_fq2;
		if(check_gz_empty(trim_file1)==1){
			string rm_cmd="rm "+trim_file1;
			system(rm_cmd.c_str());
		}
		if(check_gz_empty(trim_file2)==1){
			string rm_cmd="rm "+trim_file2;
			system(rm_cmd.c_str());
		}
	}
	if(!gp.clean_fq1.empty()){
		string clean_file1,clean_file2;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		clean_file2=gp.output_dir+"/"+gp.clean_fq2;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
		if(check_gz_empty(clean_file2)==1){
			string rm_cmd="rm "+clean_file2;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream trim_out_fq1_tmp,trim_out_fq2_tmp;
			trim_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r1.fq.gz";
			trim_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r2.fq.gz";
			
			if(access(trim_out_fq1_tmp.str().c_str(),0)!=-1 && access(trim_out_fq2_tmp.str().c_str(),0)!=-1){

				ostringstream cat_cmd1,cat_cmd2;
				cat_cmd1<<"cat "<<trim_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq1<<";rm "<<trim_out_fq1_tmp.str();
				cat_cmd2<<"cat "<<trim_out_fq2_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq2<<";rm "<<trim_out_fq2_tmp.str();
				if(gp.threads_num>1){
					thread cat_t1(bind(&peProcess::run_cmd,this,cat_cmd1.str()));
					thread cat_t2(bind(&peProcess::run_cmd,this,cat_cmd2.str()));
					cat_t1.join();
					cat_t2.join();
				}else{
					run_cmd(cat_cmd1.str());
					run_cmd(cat_cmd2.str());
				}
			}
			
		}
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp,clean_out_fq2_tmp;
			if(i==index){
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/last.r1.fq.gz";
				clean_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/last.r2.fq.gz";
			}else{
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
				clean_out_fq2_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r2.fq.gz";
			}
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1 && access(clean_out_fq2_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1,cat_cmd2;
				if(i<=index){
					cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
					cat_cmd2<<"cat "<<clean_out_fq2_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq2<<";rm "<<clean_out_fq2_tmp.str();
					if(i==index){
						cat_cmd1<<";rm "<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
						cat_cmd2<<";rm "<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r2.fq.gz";
					}
				}else{
					cat_cmd1<<"rm "<<clean_out_fq1_tmp.str();
					cat_cmd2<<"rm "<<clean_out_fq2_tmp.str();
				}
				if(gp.threads_num>1){
					thread cat_t1(bind(&peProcess::run_cmd,this,cat_cmd1.str()));
					thread cat_t2(bind(&peProcess::run_cmd,this,cat_cmd2.str()));
					cat_t1.join();
					cat_t2.join();
				}else{
					run_cmd(cat_cmd1.str());
					run_cmd(cat_cmd2.str());
				}
			}
		}
	}
	
}
void peProcess::create_thread_read(int index){
	multi_gzfq1[index]=gzopen((gp.fq1_path).c_str(),"rb");
	if(!multi_gzfq1[index]){
		cerr<<"Error:cannot open the file,"<<gp.fq1_path<<endl;
		exit(1);
	}
	gzsetparams(multi_gzfq1[index], 2, Z_DEFAULT_STRATEGY);
	gzbuffer(multi_gzfq1[index],2048*2048);
	multi_gzfq2[index]=gzopen((gp.fq2_path).c_str(),"rb");
	if(!multi_gzfq2[index]){
		cerr<<"Error:cannot open the file,"<<gp.fq2_path<<endl;
		exit(1);
	}
	gzsetparams(multi_gzfq2[index], 2, Z_DEFAULT_STRATEGY);
	gzbuffer(multi_gzfq2[index],2048*2048);
}
void* peProcess::sub_thread_nonssd_realMultiThreads(int index){
	of_log<<get_local_time()<<"\tthread "<<index<<" start"<<endl;
	if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
		//create_thread_outputFile(index);
		create_thread_trimoutputFile(index);
		create_thread_cleanoutputFile(index);
	}else if(!gp.trim_fq1.empty()){
		create_thread_trimoutputFile(index);
	}
	create_thread_read(index);
	
	char buf1[READBUF],buf2[READBUF];
	C_fastq fastq1,fastq2;
	C_fastq_init(fastq1,fastq2);
	unsigned long long file1_line_num(0),file2_line_num(0);
	unsigned long long block_line_num1(0),block_line_num2(0);
	int thread_read_block=4*gp.patchSize;
	vector<C_fastq> fq1s,fq2s;
	while(1){
		if(gzgets(multi_gzfq1[index],buf1,READBUF)!=NULL){
			if((file1_line_num/thread_read_block)%gp.threads_num==index){
				block_line_num1++;
				if(block_line_num1%4==1){
					fastq1.seq_id.assign(buf1);
					fastq1.seq_id.erase(fastq1.seq_id.size()-1);
				}
				if(block_line_num1%4==2){
					fastq1.sequence.assign(buf1);
					fastq1.sequence.erase(fastq1.sequence.size()-1);
				}
				if(block_line_num1%4==0){
					fastq1.qual_seq.assign(buf1);
					fastq1.qual_seq.erase(fastq1.qual_seq.size()-1);
				}
			}
			file1_line_num++;
		}
		if(gzgets(multi_gzfq2[index],buf2,READBUF)!=NULL){
			if((file2_line_num/thread_read_block)%gp.threads_num==index){
				block_line_num2++;
				if(block_line_num2%4==1){
					fastq2.seq_id.assign(buf2);
					fastq2.seq_id.erase(fastq2.seq_id.size()-1);
				}else if(block_line_num2%4==2){
					fastq2.sequence.assign(buf2);
					fastq2.sequence.erase(fastq2.sequence.size()-1);
				}else if(block_line_num2%4==0){
					fastq2.qual_seq.assign(buf2);
					fastq2.qual_seq.erase(fastq2.qual_seq.size()-1);
					fq1s.push_back(fastq1);
					fq2s.push_back(fastq2);
					if(fq1s.size()==gp.patchSize){
						if(index==0)
							of_log<<get_local_time()<<" processed_reads:\t"<<file1_line_num/4<<endl;
						thread_process_reads(index,fq1s,fq2s);
						if(limit_end>0){
							break;
						}
					}
				}
			}
			file2_line_num++;
		}else{
			if(fq1s.size()>0){
				thread_process_reads(index,fq1s,fq2s);
				if(limit_end){
					break;
				}
			}
			gzclose(multi_gzfq1[index]);
			gzclose(multi_gzfq2[index]);
			break;
		}
	}
	if(limit_end>0){
		gzclose(multi_gzfq1[index]);
		gzclose(multi_gzfq2[index]);
	}
	create_thread_read(index);
	gzclose(multi_gzfq1[index]);
	gzclose(multi_gzfq2[index]);
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1[index]);
		gzclose(gz_trim_out2[index]);
	}
	if(!gp.clean_fq1.empty()){
		if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
			gzclose(gz_clean_out1[index]);
			gzclose(gz_clean_out2[index]);
		}
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}
void* peProcess::monitor_read_thread(){
	//queue<vector<C_fastq> > thread_read_buffer[max_thread],thread_write_buffer[max_thread];
	//mutex thread_read_m[max_thread],thread_write_m[max_thread];
	vector<C_fastq> fq1s,fq2s;
	int done(0);
	while(1){
		if(file_end){
			break;
		}else{
			int ready(0);
			for(int i=0;i!=gp.threads_num;i++){
				if(thread_read_m[i].try_lock()){
					cout<<"monitor: sub thread "<<i<<" size "<<thread_read_buffer1[i].size()<<endl;
					if(thread_read_buffer1[i].size()<queue_buffer){
						if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
							if(read_gz(fq1s,fq2s)==-1){		//read fastqs from raw files
								done=1;
							}
						}else{
							if(read(fq1s,fq2s,nongzfp1,nongzfp2)==-1){		
								done=1;
							}
						}
						thread_read_buffer1[i].push(fq1s);
						thread_read_buffer2[i].push(fq2s);
						fq1s.clear();
						fq2s.clear();
						if(done==1){
							file_end=true;
							thread_read_m[i].unlock();
							break;
						}
						cout<<"done:"<<done<<" after monitor: sub thread "<<i<<" size "<<thread_read_buffer1[i].size()<<endl;
					}else{
						ready++;
					}
					thread_read_m[i].unlock();
					if(ready==gp.threads_num){
						usleep(1000000);
					}
				}
			}
		}
	}
}
void peProcess::process_nonssd(){
	string mkdir_str="mkdir -p "+gp.output_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir fail"<<endl;
		exit(1);
	}
	of_log.open(gp.log.c_str());
	if(!of_log){
		cerr<<"Error:cannot open such file,"<<gp.log<<endl;
		exit(1);
	}
	of_log<<get_local_time()<<"\tAnalysis start!"<<endl;
	if(gp.output_clean<=0){
		if(!(gp.l_total_reads_num>0 && gp.total_reads_num_random==false)){
			make_tmpDir();
		}else if(!gp.trim_fq1.empty()){
			make_tmpDir();
		}
	}else if(!gp.trim_fq1.empty()){
		make_tmpDir();
	}
	thread t_array[gp.threads_num];
	//thread read_monitor(bind(&peProcess::monitor_read_thread,this));
	//sleep(10);
	for(int i=0;i<gp.threads_num;i++){
		//t_array[i]=thread(bind(&peProcess::sub_thread_nonssd_multiOut,this,i));
		t_array[i]=thread(bind(&peProcess::sub_thread_nonssd_realMultiThreads,this,i));
	}
	for(int i=0;i<gp.threads_num;i++){
		t_array[i].join();
	}
	//gzclose(gz_fq1);
	//gzclose(gz_fq2);
	//read_monitor.join();
	if(gp.total_reads_num_random==true && gp.total_reads_num>0){
		run_extract_random();
		remove_tmpDir();
		of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
		of_log.close();
		return;
	}
	merge_stat_nonssd();
	print_stat();
	if(!gp.trim_fq1.empty()){
		merge_trim_data();
	}
	if(gp.output_clean<=0 && !(gp.l_total_reads_num>0 && gp.total_reads_num_random==false)){
		merge_clean_data();
		remove_tmpDir();
	}else{
		if(!gp.trim_fq1.empty()){
			remove_tmpDir();
		}
		if(limit_end==0 && (gp.output_clean>0 || gp.l_total_reads_num>0)){
			gzclose(gz_fq1);
			gzclose(gz_fq2);
		}
	}
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
}
void peProcess::run_extract_random(){
	if(gp.total_reads_num<=0 || gp.total_reads_num_random==false){
		cerr<<"Error:extract random clean reads error cuz parameters are wrong"<<endl;
		exit(1);
	}
	unsigned long long total_clean_reads(0);
	for(int i=0;i!=gp.threads_num;i++){
		total_clean_reads+=local_clean_stat1[i].gs.reads_number;
	}
	if(gp.f_total_reads_ratio>0){
		if(gp.f_total_reads_ratio>=1){
			cerr<<"Error:the ratio extract from clean fq file should not be more than 1"<<endl;
			exit(1);
		}
		if(gp.l_total_reads_num>0){
			cerr<<"Error:reads number and ratio should not be both assigned at the same time"<<endl;
			exit(1);
		}
		gp.l_total_reads_num=total_clean_reads*gp.f_total_reads_ratio;
	}
	if(total_clean_reads<gp.l_total_reads_num){
		cerr<<"Warning:the reads number in clean fastq file("<<total_clean_reads<<") is less than you assigned to output("<<gp.l_total_reads_num<<")"<<endl;
	}
	//cout<<gp.l_total_reads_num<<"\t"<<total_clean_reads<<endl;
	//vector<int> include_threads;
	int last_thread(0);
	int sticky_end(0);
	unsigned long long cur_total(0);
	for(int i=0;i!=gp.threads_num;i++){
		cur_total+=local_clean_stat1[i].gs.reads_number;
		if(cur_total>gp.l_total_reads_num){
			last_thread=i;
			sticky_end=local_clean_stat1[i].gs.reads_number-(cur_total-gp.l_total_reads_num);
			break;
		}
	}
	if(cur_total<=gp.l_total_reads_num)
		last_thread=gp.threads_num-1;
	//create the last patch clean fq file and stat
	if(sticky_end>0){
		process_some_reads(last_thread,sticky_end);
		merge_stat(last_thread);
		print_stat();
		merge_clean_data(last_thread);
	}else{
		merge_stat(last_thread);
		print_stat();
		merge_clean_data();
	}
	if(!gp.trim_fq1.empty()){
		merge_trim_data();
	}
}
void peProcess::process_some_reads(int index,int out_number){
	//open target fq file
	ostringstream target_file_fq1,target_file_fq2;
	target_file_fq1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r1.fq.gz";
	target_file_fq2<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r2.fq.gz";
	int file_reads_number=local_clean_stat1[index].gs.reads_number;
	local_clean_stat1[index].clear();
	local_clean_stat2[index].clear();
	//cout<<file_reads_number<<"\t"<<out_number<<endl;
	if(file_reads_number==out_number){
		return;
	}

	int times=(int)floor(file_reads_number/out_number);
	gzFile tmp_fq1=gzopen(target_file_fq1.str().c_str(),"rb");
	gzFile tmp_fq2=gzopen(target_file_fq2.str().c_str(),"rb");
	gzsetparams(tmp_fq1, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(tmp_fq1,2048*2048);
	gzsetparams(tmp_fq2, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(tmp_fq2,2048*2048);
	//set output file
	string last_file1=gp.output_dir+"/"+tmp_dir+"/last.r1.fq.gz";
	string last_file2=gp.output_dir+"/"+tmp_dir+"/last.r2.fq.gz";
	gzFile tmp_out_fq1=gzopen(last_file1.c_str(),"wb");
	gzFile tmp_out_fq2=gzopen(last_file2.c_str(),"wb");
	gzsetparams(tmp_out_fq1, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(tmp_out_fq1,2048*2048);
	gzsetparams(tmp_out_fq2, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(tmp_out_fq2,2048*2048);

	char buf1[READBUF],buf2[READBUF];
	C_fastq fastq1,fastq2;
	C_fastq_init(fastq1,fastq2);
	unsigned long long file1_line_num(0),file2_line_num(0);
	vector<C_fastq> fq1s,fq2s;
	int processed_number(0),rest_number(out_number);
	while(1){
		//cout<<"here"<<endl;
		if(rest_number<=0 || processed_number>=out_number){
			break;
		}
		if(gzgets(tmp_fq1,buf1,READBUF)!=NULL){
			//cout<<"here2\t"<<file1_line_num<<endl;
			if(file1_line_num%(times*4)<4){
				if(file1_line_num%4==0){
					fastq1.seq_id.assign(buf1);
					fastq1.seq_id.erase(fastq1.seq_id.size()-1);
				}
				if(file1_line_num%4==1){
					fastq1.sequence.assign(buf1);
					fastq1.sequence.erase(fastq1.sequence.size()-1);
				}
				if(file1_line_num%4==3){
					fastq1.qual_seq.assign(buf1);
					fastq1.qual_seq.erase(fastq1.qual_seq.size()-1);
				}
			}
			file1_line_num++;
		}
		if(gzgets(tmp_fq2,buf2,READBUF)!=NULL){
			//cout<<"here2.2\t"<<file2_line_num<<endl;
			if(file2_line_num%(times*4)<4){
				if(file2_line_num%4==0){
					fastq2.seq_id.assign(buf2);
					fastq2.seq_id.erase(fastq2.seq_id.size()-1);
				}else if(file2_line_num%4==1){
					fastq2.sequence.assign(buf2);
					fastq2.sequence.erase(fastq2.sequence.size()-1);
				}else if(file2_line_num%4==3){
					fastq2.qual_seq.assign(buf2);
					fastq2.qual_seq.erase(fastq2.qual_seq.size()-1);
					fq1s.push_back(fastq1);
					fq2s.push_back(fastq2);
					if(fq1s.size()==gp.patchSize || fq1s.size()==rest_number){
						of_log<<get_local_time()<<" last sticky thread processed reads:\t"<<file1_line_num/4<<endl;
						//thread_process_reads(index,fq1s,fq2s);
						processed_number+=fq1s.size();
						rest_number=out_number-processed_number;
						//cout<<"processed number\t"<<processed_number<<"\trest number\t"<<rest_number<<endl;
						limit_process_reads(index,fq1s,fq2s,tmp_out_fq1,tmp_out_fq2);
					}
				}
			}
			file2_line_num++;
		}else{
			if(fq1s.size()>0){
				limit_process_reads(index,fq1s,fq2s,tmp_out_fq1,tmp_out_fq2);
			}
			break;
		}
	}
	gzclose(tmp_fq1);
	gzclose(tmp_fq2);
	gzclose(tmp_out_fq1);
	gzclose(tmp_out_fq2);
}
void peProcess::limit_process_reads(int index,vector<C_fastq> &fq1s,vector<C_fastq> &fq2s,gzFile gzfq1,gzFile gzfq2){
	PEstatOption opt_clean;
	opt_clean.fq1s=&fq1s;
	opt_clean.stat1=&local_clean_stat1[index];
	opt_clean.fq2s=&fq2s;
	opt_clean.stat2=&local_clean_stat2[index];
	stat_pe_fqs(opt_clean);		//statistic raw fastqs
	peWrite(fq1s,fq2s,"clean",gzfq1,gzfq2);
	fq1s.clear();
	fq2s.clear();
}
void peProcess::remove_tmpDir(){
	string rm_dir=gp.output_dir+"/"+tmp_dir;
	int iter(0);
	while(1){
		DIR* tmpdir;
		struct dirent* ptr;
		tmpdir=opendir(rm_dir.c_str());
		bool empty_flag=true;
		if(tmpdir==NULL){
			cerr<<"Error:open directory error,"<<rm_dir<<endl;
			exit(1);
		}
		while((ptr=readdir(tmpdir))!=NULL){
			if(strcmp(ptr->d_name,"..")!=0 && strcmp(ptr->d_name,".")!=0){
				empty_flag=false;
				break;
			}
		}
		closedir(tmpdir);
		if(empty_flag){
			string remove="rm -r "+rm_dir;
			if(system(remove.c_str())==-1){
				cerr<<"Error:rmdir error,"<<remove<<endl;
				exit(1);
			}
			break;
		}else{
			sleep(2);
			iter++;
		}
		if(iter>30){
			break;
		}
	}
}
void peProcess::make_tmpDir(){
	srand(time(0));
	ostringstream tmp_str;
	for(int i=0;i!=6;i++){
		int tmp_rand=random(26)+'A';
		tmp_str<<(char)tmp_rand;
	}
	tmp_dir="TMP"+tmp_str.str();
	string mkdir_str="mkdir -p "+gp.output_dir+"/"+tmp_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir error,"<<mkdir_str<<endl;
		exit(1);
	}
}
void peProcess::process(){
	string mkdir_str="mkdir -p "+gp.output_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir fail"<<endl;
		exit(1);
	}
	of_log.open(gp.log.c_str());
	if(!of_log){
		cerr<<"Error:cannot open such file,"<<gp.log<<endl;
		exit(1);
	}
	of_log<<get_local_time()<<"\tAnalysis start!"<<endl;
	new_fq1_path=gp.fq1_path;
	new_fq2_path=gp.fq2_path;
	if(gp.fq1_path.rfind(".gz")==gp.fq1_path.size()-3){
		if(gp.threads_num==1){
			run_pigz(1);
			run_pigz(2);
		}else{
			thread t_pigz1(bind(&peProcess::run_pigz,this,1));
			thread t_pigz2(bind(&peProcess::run_pigz,this,2));
			t_pigz1.join();
			t_pigz2.join();
		}
		new_fq1_path=gp.output_dir+"/raw.r1.fq";
		new_fq2_path=gp.output_dir+"/raw.r2.fq";
		
	}
	if(gp.output_clean<=0){
		if(!(gp.l_total_reads_num>0 && gp.total_reads_num_random==false)){
			make_tmpDir();
		}else if(!gp.trim_fq1.empty()){
			make_tmpDir();
		}
	}else if(!gp.trim_fq1.empty()){
		make_tmpDir();
	}
	fq1fd=open(new_fq1_path.c_str(),O_RDONLY);
	fq2fd=open(new_fq2_path.c_str(),O_RDONLY);
	off_t file_size,file2_size;
	struct stat st,st2;
	fstat(fq1fd,&st);
	file_size=st.st_size;
	fstat(fq1fd,&st2);
	file2_size=st2.st_size;
	//cout<<file_size<<endl;
	int pieces=file_size/buffer;
	if(file_size%buffer!=0 || pieces==0)
		pieces+=1;
	int blocks=pieces/gp.threads_num;
	int real_block=blocks;
	if(blocks==0){
		real_block=1;
		used_threads_num=pieces;
	}else{
		used_threads_num=gp.threads_num;
	}
	for(int i=0;i!=used_threads_num;i++){
		t_start_pos[i]=i*real_block*buffer;
		if(i==used_threads_num-1){
			t_end_pos[i]=file_size;
		}else{
			t_end_pos[i]=(i+1)*real_block*buffer;
		}
	}	
	if(file_size!=file2_size){
		cerr<<"Error:input PE fastqs are different in char size"<<endl;
		exit(1);
	}
	close(fq1fd);
	close(fq2fd);
	thread t_array[used_threads_num];
	for(int i=0;i!=used_threads_num;i++){
		t_array[i]=thread(bind(&peProcess::sub_thread,this,i));
	}
	for(int i=0;i<used_threads_num;i++){
		t_array[i].join();
	}
	//cout<<limit_end<<endl;
	if(limit_end<=0){
		string fq1_unprocess,fq2_unprocess;
		for(int i=0;i<used_threads_num;i++){
			//cout<<sticky_head1[i]<<"\n"<<sticky_tail1<<endl;
			fq1_unprocess+=sticky_reads1[i];
			fq2_unprocess+=sticky_reads2[i];
		}
		//
		int line_num(0);
		C_fastq fastq1,fastq2;
		C_fastq_init(fastq1,fastq2);
		vector<C_fastq> fq1s,fq2s;
		for(int i=0;i!=fq1_unprocess.size();i++){
			if(fq1_unprocess[i]=='\n'){
				if(line_num%4==3){
					fq1s.push_back(fastq1);
					fastq1.seq_id.clear();
					fastq1.sequence.clear();
					fastq1.qual_seq.clear();
				}
				line_num++;
				continue;
			}
			if(line_num%4==0){
				fastq1.seq_id.insert(fastq1.seq_id.end(),fq1_unprocess[i]);
			}
			if(line_num%4==1){
				fastq1.sequence.insert(fastq1.sequence.end(),fq1_unprocess[i]);
			}
			if(line_num%4==3){
				fastq1.qual_seq.insert(fastq1.qual_seq.end(),fq1_unprocess[i]);
			}
		}
		line_num=0;
		for(int i=0;i!=fq2_unprocess.size();i++){
			if(fq2_unprocess[i]=='\n'){
				if(line_num%4==3){
					fq2s.push_back(fastq2);
					fastq2.seq_id.clear();
					fastq2.sequence.clear();
					fastq2.qual_seq.clear();
				}
				line_num++;
				continue;
			}
			if(line_num%4==0){
				fastq2.seq_id.insert(fastq2.seq_id.end(),fq2_unprocess[i]);
			}
			if(line_num%4==1){
				fastq2.sequence.insert(fastq2.sequence.end(),fq2_unprocess[i]);
			}
			if(line_num%4==3){
				fastq2.qual_seq.insert(fastq2.qual_seq.end(),fq2_unprocess[i]);
			}
		}
		long long clean_size(0);
		thread_process_reads(0,fq1s,fq2s);
	}
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1[0]);
		gzclose(gz_trim_out2[0]);
	}
	if(!gp.clean_fq1.empty()){
		if(gp.output_clean<=0 && !(gp.l_total_reads_num>0 && gp.total_reads_num_random==false)){
			gzclose(gz_clean_out1[0]);
			gzclose(gz_clean_out2[0]);
		}
	}
	if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
		string new_fq1_path=gp.output_dir+"/raw.r1.fq";
		string new_fq2_path=gp.output_dir+"/raw.r2.fq";
		string rm_cmd="rm "+new_fq1_path+";rm "+new_fq2_path;
		if(system(rm_cmd.c_str())){
			cerr<<"Error:rm error"<<endl;
			exit(1);
		} 
	}
	if(gp.total_reads_num_random==true && gp.total_reads_num>0){
		run_extract_random();
		remove_tmpDir();
		of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
		of_log.close();
		return;
	}
	if(!gp.trim_fq1.empty()){
		merge_trim_data();
	}
	//
	merge_stat();
	print_stat();
	if(gp.output_clean<=0 && !(gp.l_total_reads_num>0 && gp.total_reads_num_random==false)){
		merge_trim_data();
		merge_clean_data();
		remove_tmpDir();
	}else{
		if(!gp.trim_fq1.empty()){
			remove_tmpDir();
		}
		if(limit_end==0 && (gp.output_clean>0 || gp.l_total_reads_num>0)){
			gzclose(gz_fq1);
			gzclose(gz_fq2);
		}
	}
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
	
}
void peProcess::output_fastqs2(int type,vector<C_fastq> &fq1,ofstream& outfile){
	//m.lock();
	string out_content;
	
	for(int i=0;i!=fq1.size();i++){
	//for(vector<C_fastq>::iterator i=fq1->begin();i!=fq1->end();i++){
		if(gp.output_file_type=="fasta"){
			fq1[i].seq_id=fq1[i].seq_id.replace(fq1[i].seq_id.find("@"),1,">");
			out_content+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n";
		}else if(gp.output_file_type=="fastq"){
			if(gp.outputQualityPhred!=gp.qualityPhred){
				for(string::size_type ix=0;ix!=fq1[i].qual_seq.size();ix++){
					fq1[i].qual_seq[ix]=(char)(fq1[i].qual_seq[ix]-gp.outputQualityPhred);
				}
			}
			out_content+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
		}else{
			cerr<<"Error:output_file_type value error"<<endl;
			exit(1);
		}
	}
	if(gp.is_streaming){
		cout<<out_content;
	}else{
		outfile<<out_content;
	}
	//m.unlock();
}
void peProcess::output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile){
	//m.lock();
	string out_content,streaming_out;
	for(int i=0;i!=fq1.size();i++){
	//for(vector<C_fastq>::iterator i=fq1->begin();i!=fq1->end();i++){
		if(gp.output_file_type=="fasta"){
			fq1[i].seq_id=fq1[i].seq_id.replace(fq1[i].seq_id.find("@"),1,">");
			out_content+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n";
		}else if(gp.output_file_type=="fastq"){
			if(gp.outputQualityPhred!=gp.qualityPhred){
				for(string::size_type ix=0;ix!=fq1[i].qual_seq.size();ix++){
					int b_q=fq1[i].qual_seq[ix]-gp.qualityPhred;
					fq1[i].qual_seq[ix]=(char)(b_q+gp.outputQualityPhred);
				}
			}
			if(gp.is_streaming){
				string modify_id=fq1[i].seq_id;
				modify_id.erase(0,1);
				streaming_out+=">+\t"+modify_id+"\t"+type+"\t"+fq1[i].sequence+"\t"+fq1[i].qual_seq+"\n";
			}else{
				out_content+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
			}
		}else{
			cerr<<"Error:output_file_type value error"<<endl;
			exit(1);
		}
	}
	if(gp.is_streaming){
		cout<<streaming_out;
	}else{
		gzwrite(outfile,out_content.c_str(),out_content.size());
		//gzflush(outfile,1);
	}
	//m.unlock();
}
void peProcess::output_split_fastqs(string type,vector<C_fastq> &fq1){
	//m.lock();
	if(gp.clean_file_reads<=0){
		cerr<<"Error:output clean fastq file reads number should more than 0 when assigned -w or -L"<<endl;
		exit(1);
	}
	string streaming_out,out_content;
	int patch_idx,patch_mod;
	for(int i=0;i!=fq1.size();i++){
	//for(vector<C_fastq>::iterator i=fq1->begin();i!=fq1->end();i++){
		if(gp.output_file_type=="fastq"){
			if(gp.outputQualityPhred!=gp.qualityPhred){
				for(string::size_type ix=0;ix!=fq1[i].qual_seq.size();ix++){
					int b_q=fq1[i].qual_seq[ix]-gp.qualityPhred;
					fq1[i].qual_seq[ix]=(char)(b_q+gp.outputQualityPhred);
				}
			}
			if(gp.is_streaming){
				string modify_id=fq1[i].seq_id;
				modify_id.erase(0,1);
				streaming_out+=">+\t"+modify_id+"\t"+type+"\t"+fq1[i].sequence+"\t"+fq1[i].qual_seq+"\n";
			}else{
				out_content+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				
			}
		}
	}

	if(gp.is_streaming){
		cout<<streaming_out;
	}else{
		if(type=="1"){
			
			gp.have_output1+=fq1.size();
			//cout<<"here\t"<<gp.have_output1<<"\t"<<gp.output_clean<<endl;
			int idx=(gp.have_output1-1)/gp.clean_file_reads;
			pe1_out[idx]++;
			int mod=gp.have_output1%gp.clean_file_reads;
			if(pe1_out[idx]==1){
				if(limit_end>0){
					return;
				}
				int to_output=fq1.size()-mod;
				string sticky_tail,sticky_head;
				for(int i=0;i!=to_output;i++){
					sticky_tail+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				}
				for(int i=to_output;i!=fq1.size();i++){
					sticky_head+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				}
				if(!sticky_tail.empty()){
					gzwrite(gz_fq1,sticky_tail.c_str(),sticky_tail.size());
					gzclose(gz_fq1);
					if(gp.total_reads_num_random==false && gp.l_total_reads_num>0)
						goto LAB1;
				}else{
					if(idx>0){
						gzclose(gz_fq1);
						if(gp.total_reads_num_random==false && gp.l_total_reads_num>0){
							goto LAB1;
						}
					}
				}
				ostringstream out_fq1;
				out_fq1<<gp.output_dir<<"/split."<<idx<<"."<<gp.clean_fq1;
				if(gp.total_reads_num_random==false && gp.l_total_reads_num>0){
					string fq1_whole_path=gp.output_dir+"/"+gp.clean_fq1;
					gz_fq1=gzopen(fq1_whole_path.c_str(),"wb");
				}else{
					gz_fq1=gzopen(out_fq1.str().c_str(),"wb");
				}
				gzsetparams(gz_fq1, 2, Z_DEFAULT_STRATEGY);
				gzbuffer(gz_fq1,1024*1024*10);
				if(!sticky_head.empty()){
					gzwrite(gz_fq1,sticky_head.c_str(),sticky_head.size());
				}
			}else{
				//cout<<idx<<"\t"<<pe1_out[idx]<<"\t"<<gp.have_output1<<"\t"<<gp.clean_file_reads<<"\t"<<mod<<"\t"<<out_content.size()<<endl;
				gzwrite(gz_fq1,out_content.c_str(),out_content.size());
			}
		}
		LAB1:
		if(type=="2"){
			gp.have_output2+=fq1.size();
			int idx=(gp.have_output2-1)/gp.clean_file_reads;
			pe2_out[idx]++;
			int mod=gp.have_output2%gp.clean_file_reads;
			if(pe2_out[idx]==1){
				int to_output=fq1.size()-mod;
				string sticky_tail,sticky_head;
				for(int i=0;i!=to_output;i++){
					sticky_tail+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				}
				for(int i=to_output;i!=fq1.size();i++){
					sticky_head+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				}
				if(!sticky_tail.empty()){
					gzwrite(gz_fq2,sticky_tail.c_str(),sticky_tail.size());
					gzclose(gz_fq2);
					if(gp.total_reads_num_random==false && gp.l_total_reads_num>0){
						limit_end++;
						return;
					}
				}else{
					if(idx>0){
						if(gp.total_reads_num_random==false && gp.l_total_reads_num>0){
							gzclose(gz_fq2);
							limit_end++;
							return;
						}else{
							gzclose(gz_fq2);
						}
					}
				}
				ostringstream out_fq2;
				out_fq2<<gp.output_dir<<"/split."<<idx<<"."<<gp.clean_fq2;
				if(gp.total_reads_num_random==false && gp.l_total_reads_num>0){
					string fq2_whole_path=gp.output_dir+"/"+gp.clean_fq2;
					gz_fq2=gzopen(fq2_whole_path.c_str(),"wb");
				}else{
					gz_fq2=gzopen(out_fq2.str().c_str(),"wb");
				}
				gzsetparams(gz_fq2, 2, Z_DEFAULT_STRATEGY);
				gzbuffer(gz_fq2,1024*1024*10);
				if(!sticky_head.empty()){
					gzwrite(gz_fq2,sticky_head.c_str(),sticky_head.size());
				}
			}else{
				gzwrite(gz_fq2,out_content.c_str(),out_content.size());
			}
		}
	}
	//m.unlock();
}
void peProcess::peStreaming_stat(C_global_variable& local_gv){
	cout<<"#Total_statistical_information"<<"\n";
	/*int output_reads_num;
	int in_adapter_list_num;
	int include_adapter_seq_num;
	int include_contam_seq_num;
	int n_ratio_num;
	int highA_num,polyX_num;
	int tile_num,fov_num;
	int low_qual_base_ratio_num;
	int mean_quality_num;
	int short_len_num,long_len_num;
	int over_lapped_num;
	int no_3_adapter_num,int_insertNull_num;
	*/
	int total=local_gv.fs.include_adapter_seq_num+local_gv.fs.include_contam_seq_num+local_gv.fs.low_qual_base_ratio_num+local_gv.fs.mean_quality_num+local_gv.fs.n_ratio_num+local_gv.fs.over_lapped_num+local_gv.fs.highA_num+local_gv.fs.polyX_num;
	cout<<total<<" "<<local_gv.fs.include_adapter_seq_num<<" "<<local_gv.fs.include_contam_seq_num<<" "<<local_gv.fs.low_qual_base_ratio_num<<" "<<local_gv.fs.mean_quality_num<<" "<<local_gv.fs.n_ratio_num<<" "<<local_gv.fs.over_lapped_num<<" "<<local_gv.fs.highA_num<<" "<<local_gv.fs.polyX_num<<"\n";
	/*int read_max_length;
	int read_length;
	int reads_number;
	unsigned long long base_number;
	unsigned long long a_number,c_number,g_number,t_number,n_number;
	//unsigned long long a_ratio,c_ratio,g_ratio,t_ratio,n_ratio;
	unsigned long long q20_num,q30_num;
	*/
	cout<<"#Fq1_statistical_information"<<"\n";
	cout<<local_gv.raw1_stat.gs.read_length<<" "<<local_gv.clean1_stat.gs.read_length<<" "<<local_gv.raw1_stat.gs.reads_number<<" "<<local_gv.clean1_stat.gs.reads_number<<" "<<local_gv.raw1_stat.gs.base_number<<" "<<local_gv.clean1_stat.gs.base_number<<" "<<local_gv.raw1_stat.gs.a_number<<" "<<local_gv.clean1_stat.gs.a_number<<" "<<local_gv.raw1_stat.gs.c_number<<" "<<local_gv.clean1_stat.gs.c_number<<" "<<local_gv.raw1_stat.gs.g_number<<" "<<local_gv.clean1_stat.gs.g_number<<" "<<local_gv.raw1_stat.gs.t_number<<" "<<local_gv.clean1_stat.gs.t_number<<" "<<local_gv.raw1_stat.gs.n_number<<" "<<local_gv.clean1_stat.gs.n_number<<" "<<local_gv.raw1_stat.gs.q20_num<<" "<<local_gv.clean1_stat.gs.q20_num<<" "<<local_gv.raw1_stat.gs.q30_num<<" "<<local_gv.clean1_stat.gs.q30_num<<"\n";
	cout<<"#Base_distributions_by_read_position"<<"\n";
	//unsigned long long position_acgt_content[READ_MAX_LEN][5];
	for(int i=0;i!=local_gv.raw1_stat.gs.read_length;i++){
		for(int j=0;j!=4;j++){
			cout<<local_gv.raw1_stat.bs.position_acgt_content[i][j]<<" ";
		}
		cout<<local_gv.raw1_stat.bs.position_acgt_content[i][4]<<"\n";
	}
	for(int i=0;i!=local_gv.clean1_stat.gs.read_length;i++){
		for(int j=0;j!=4;j++){
			cout<<local_gv.clean1_stat.bs.position_acgt_content[i][j]<<" ";
		}
		cout<<local_gv.clean1_stat.bs.position_acgt_content[i][4]<<"\n";
	}
	cout<<"#Raw_Base_quality_value_distribution_by_read_position"<<"\n";
	//position_qual[READ_MAX_LEN][MAX_QUAL]
	for(int i=0;i!=local_gv.raw1_stat.gs.read_length;i++){
		for(int j=0;j!=40;j++){
			cout<<local_gv.clean1_stat.qs.position_qual[i][j]<<" ";
		}
		cout<<"0\n";
	}
	for(int i=0;i!=local_gv.clean1_stat.gs.read_length;i++){
		for(int j=0;j!=40;j++){
			cout<<local_gv.clean1_stat.qs.position_qual[i][j]<<" ";
		}
		cout<<"0\n";
	}
	cout<<"#Fq2_statistical_information"<<"\n";
	cout<<local_gv.raw2_stat.gs.read_length<<" "<<local_gv.clean2_stat.gs.read_length<<" "<<local_gv.raw2_stat.gs.reads_number<<" "<<local_gv.clean2_stat.gs.reads_number<<" "<<local_gv.raw2_stat.gs.base_number<<" "<<local_gv.clean2_stat.gs.base_number<<" "<<local_gv.raw2_stat.gs.a_number<<" "<<local_gv.clean2_stat.gs.a_number<<" "<<local_gv.raw2_stat.gs.c_number<<" "<<local_gv.clean2_stat.gs.c_number<<" "<<local_gv.raw2_stat.gs.g_number<<" "<<local_gv.clean2_stat.gs.g_number<<" "<<local_gv.raw2_stat.gs.t_number<<" "<<local_gv.clean2_stat.gs.t_number<<" "<<local_gv.raw2_stat.gs.n_number<<" "<<local_gv.clean2_stat.gs.n_number<<" "<<local_gv.raw2_stat.gs.q20_num<<" "<<local_gv.clean2_stat.gs.q20_num<<" "<<local_gv.raw2_stat.gs.q30_num<<" "<<local_gv.clean2_stat.gs.q30_num<<"\n";
	cout<<"#Base_distributions_by_read_position"<<"\n";
	//unsigned long long position_acgt_content[READ_MAX_LEN][5];
	for(int i=0;i!=local_gv.raw2_stat.gs.read_length;i++){
		for(int j=0;j!=4;j++){
			cout<<local_gv.raw2_stat.bs.position_acgt_content[i][j]<<" ";
		}
		cout<<local_gv.raw2_stat.bs.position_acgt_content[i][4]<<"\n";
	}
	for(int i=0;i!=local_gv.clean2_stat.gs.read_length;i++){
		for(int j=0;j!=4;j++){
			cout<<local_gv.clean2_stat.bs.position_acgt_content[i][j]<<" ";
		}
		cout<<local_gv.clean2_stat.bs.position_acgt_content[i][4]<<"\n";
	}
	cout<<"#Raw_Base_quality_value_distribution_by_read_position"<<"\n";
	//position_qual[READ_MAX_LEN][MAX_QUAL]
	for(int i=0;i!=local_gv.raw2_stat.gs.read_length;i++){
		for(int j=0;j!=41;j++){
			cout<<local_gv.raw2_stat.qs.position_qual[i][j]<<" ";
		}
		cout<<"0\n";
	}
	for(int i=0;i!=local_gv.clean2_stat.gs.read_length;i++){
		for(int j=0;j!=41;j++){
			cout<<local_gv.clean2_stat.qs.position_qual[i][j]<<" ";
		}
		cout<<"0\n";
	}
}
