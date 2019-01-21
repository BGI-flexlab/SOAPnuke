#include <iostream>
#include <string>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <thread>
#include <mutex>
#include <fstream>
#include <iomanip>
#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <math.h>
#include "seprocess.h"
#include "process_argv.h"
#include "zlib.h"
#include "gc.h"
#include "sequence.h"
using namespace::std;
#define READBUF 500
#define random(x) (rand()%x)
mutex se_read_m,se_stat_m,se_write_m;
C_filter_stat se_local_fs[max_thread];
C_fastq_file_stat se_local_raw_stat1[max_thread],se_local_trim_stat1[max_thread],se_local_clean_stat1[max_thread];
gzFile gz_fq_se;
string se_sticky_reads1[max_thread];
map<int,int> se_out;
int exceed_output_se(0);
string se_new_fq1_path;
unsigned long long se_buffer(1024*1024*2);
seProcess::seProcess(C_global_parameter m_gp){
	gp=m_gp;
	gv=C_global_variable();
	used_threads_num=0;
	processed_reads=0;
	srand((unsigned)time(NULL));
	ostringstream tmpstring;
	tmpstring<<rand()%100;
	random_num=tmpstring.str();
	limit_end=0;
}
void seProcess::print_stat(){
	string filter_out=gp.output_dir+"/Statistics_of_Filtered_Reads.txt";
	string general_out=gp.output_dir+"/Basic_Statistics_of_Sequencing_Quality.txt";
	string bs1_out=gp.output_dir+"/Base_distributions_by_read_position_1.txt";
	string qs1_out=gp.output_dir+"/Base_quality_value_distribution_by_read_position_1.txt";
	string q20_out1=gp.output_dir+"/Distribution_of_Q20_Q30_bases_by_read_position_1.txt";
	string trim_stat1=gp.output_dir+"/Statistics_of_Trimming_Position_of_Reads_1.txt";
	ofstream of_filter_stat(filter_out.c_str());
	ofstream of_general_stat(general_out.c_str());
	ofstream of_readPos_base_stat1(bs1_out.c_str());
	ofstream of_readPos_qual_stat1(qs1_out.c_str());
	ofstream of_q2030_stat1(q20_out1.c_str());
	ofstream of_trim_stat1(trim_stat1.c_str());
	if(!of_filter_stat){
		cerr<<"Error:cannot open such file,"<<filter_out<<endl;
		exit(1);
	}
	if(!of_general_stat){
		cerr<<"Error:cannot open such file,"<<general_out<<endl;
		exit(1);
	}
	if(!of_readPos_base_stat1){
		cerr<<"Error:cannot open such file,Base_distributions_by_read_position*.txt"<<endl;
		exit(1);
	}
	if(!of_readPos_qual_stat1){
		cerr<<"Error:cannot open such file,Base_quality_value_distribution_by_read_position*.txt"<<endl;
		exit(1);
	}
	if(!of_q2030_stat1){
		cerr<<"Error:cannot open such file,Distribution_of_Q20_Q30_bases_by_read_position*.txt"<<endl;
		exit(1);
	}
	of_filter_stat<<"Item\tTotal\tPercentage"<<endl;
	vector<string> filter_items;
	filter_items.push_back("Reads limited to output number");
	filter_items.push_back("Reads with filtered tile");
	filter_items.push_back("Reads with filtered fov");
	filter_items.push_back("Reads too short");
	filter_items.push_back("Reads too long");
	filter_items.push_back("Reads with contam sequence");
	filter_items.push_back("Reads with n rate exceed");
	filter_items.push_back("Reads with highA");
	filter_items.push_back("Reads with polyX");
	filter_items.push_back("Reads with low quality");
	filter_items.push_back("Reads with low mean quality");
	filter_items.push_back("Reads with adapter");
	filter_items.push_back("Reads with global contam sequence");
	map<string,unsigned long long> filter_number;
	filter_number["Reads with contam sequence"]=gv.fs.include_contam_seq_num;
	filter_number["Reads with global contam sequence"]=gv.fs.include_global_contam_seq_num;
	filter_number["Reads too short"]=gv.fs.short_len_num;
	filter_number["Reads with adapter"]=gv.fs.include_adapter_seq_num;
	filter_number["Reads with low quality"]=gv.fs.low_qual_base_ratio_num;
	filter_number["Reads with low mean quality"]=gv.fs.mean_quality_num;
	filter_number["Reads with n rate exceed"]=gv.fs.n_ratio_num;
	filter_number["Reads with highA"]=gv.fs.highA_num;
	filter_number["Reads with polyX"]=gv.fs.polyX_num;
	filter_number["Reads with filtered tile"]=gv.fs.tile_num;
	filter_number["Reads with filtered fov"]=gv.fs.fov_num;
	filter_number["Reads too long"]=gv.fs.long_len_num;
	//filter_number["Reads limited to output number"]=gv.fs.output_reads_num;
	unsigned long long total_filter_fq1_num=0;
	for(map<string,unsigned long long>::iterator ix=filter_number.begin();ix!=filter_number.end();ix++){
		total_filter_fq1_num+=ix->second;
	}
	//int total_filter_fq1_num=gv.fs.output_reads_num+gv.fs.include_contam_seq_num+gv.fs.include_adapter_seq_num+gv.fs.n_ratio_num+gv.fs.highA_num+gv.fs.tile_num+gv.fs.low_qual_base_ratio_num+gv.fs.mean_quality_num+gv.fs.short_len_num;
	of_filter_stat<<setiosflags(ios::fixed);
	of_filter_stat<<"Total filtered read pair number\t"<<total_filter_fq1_num<<"\t100.00%"<<endl;
	for(vector<string>::iterator ix=filter_items.begin();ix!=filter_items.end();ix++){
		if(filter_number[*ix]>0){
			of_filter_stat<<*ix<<"\t"<<filter_number[*ix]<<"\t";
			of_filter_stat<<setprecision(2)<<100*(float)filter_number[*ix]/total_filter_fq1_num<<"%"<<endl;
		}
	}
	/*of_filter_stat<<"Reads too short\t\t\t"<<gv.fs.short_len_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.short_len_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat<<"Reads with contam sequence\t"<<gv.fs.include_contam_seq_num<<"\t";
	if(gv.fs.include_contam_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.include_contam_seq_num/total_filter_fq1_num<<"%"<<endl;
	}

	of_filter_stat<<"Reads with adapter\t"<<gv.fs.include_adapter_seq_num<<"\t";
	if(gv.fs.include_adapter_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.include_adapter_seq_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat<<"Reads with low quality\t"<<gv.fs.low_qual_base_ratio_num<<"\t";
	if(gv.fs.include_adapter_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.low_qual_base_ratio_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat<<"Reads with low mean quality\t"<<gv.fs.mean_quality_num<<"\t";
	if(gv.fs.include_adapter_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.mean_quality_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat<<"Read with n rate exceed\t"<<gv.fs.n_ratio_num<<"\t";
	if(gv.fs.include_adapter_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.n_ratio_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat<<"Read with small insert size\t"<<gv.fs.over_lapped_num<<"\t";
	if(gv.fs.include_adapter_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.over_lapped_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat<<"Reads with highA\t"<<gv.fs.highA_num<<"\t";
	if(gv.fs.include_adapter_seq_num==0){
		of_filter_stat<<"0%"<<endl;
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.highA_num/total_filter_fq1_num<<"%"<<endl;
	}
	of_filter_stat.close();
	*/
	of_general_stat<<"Item\traw reads(fq1)\tclean reads(fq1)"<<endl;

	float raw1_rl(0),clean1_rl(0);
	char filter_r1_ratio[100];
	char raw_r1[7][100];
	char clean_r1[7][100];
	if(gv.raw1_stat.gs.reads_number!=0){
		raw1_rl=(float)gv.raw1_stat.gs.base_number/gv.raw1_stat.gs.reads_number;
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
		clean1_rl=(float)gv.clean1_stat.gs.base_number/gv.clean1_stat.gs.reads_number;
		sprintf(clean_r1[0],"%.2f",100*(float)gv.clean1_stat.gs.a_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[1],"%.2f",100*(float)gv.clean1_stat.gs.c_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[2],"%.2f",100*(float)gv.clean1_stat.gs.g_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[3],"%.2f",100*(float)gv.clean1_stat.gs.t_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[4],"%.2f",100*(float)gv.clean1_stat.gs.n_number/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[5],"%.2f",100*(float)gv.clean1_stat.gs.q20_num/gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[6],"%.2f",100*(float)gv.clean1_stat.gs.q30_num/gv.clean1_stat.gs.base_number);
	}
	of_general_stat<<setiosflags(ios::fixed)<<setprecision(1)<<"Read length\t"<<raw1_rl<<"\t"<<clean1_rl<<endl;
	of_general_stat<<"Total number of reads\t"<<setprecision(15)<<gv.raw1_stat.gs.reads_number<<" (100.00%)\t"<<gv.clean1_stat.gs.reads_number<<" (100.00%)"<<endl;
	of_general_stat<<"Number of filtered reads\t"<<total_filter_fq1_num<<" ("<<filter_r1_ratio<<"%)\t-"<<endl;
	unsigned long long filter_base1=total_filter_fq1_num*gv.raw1_stat.gs.read_length;
	of_general_stat<<"Total number of bases\t"<<setprecision(15)<<gv.raw1_stat.gs.base_number<<" (100.00%)\t"<<gv.clean1_stat.gs.base_number<<" (100.00%)"<<endl;
	of_general_stat<<"Number of filtered bases\t"<<setprecision(15)<<filter_base1<<" ("<<filter_r1_ratio<<"%)\t-"<<endl;
	of_general_stat<<"Number of base A\t"<<setprecision(15)<<gv.raw1_stat.gs.a_number<<" ("<<raw_r1[0]<<"%)\t"<<gv.clean1_stat.gs.a_number<<" ("<<clean_r1[0]<<"%)\t"<<endl;
	of_general_stat<<"Number of base C\t"<<setprecision(15)<<gv.raw1_stat.gs.c_number<<" ("<<raw_r1[1]<<"%)\t"<<gv.clean1_stat.gs.c_number<<" ("<<clean_r1[1]<<"%)\t"<<endl;
	of_general_stat<<"Number of base G\t"<<setprecision(15)<<gv.raw1_stat.gs.g_number<<" ("<<raw_r1[2]<<"%)\t"<<gv.clean1_stat.gs.g_number<<" ("<<clean_r1[2]<<"%)\t"<<endl;
	of_general_stat<<"Number of base T\t"<<setprecision(15)<<gv.raw1_stat.gs.t_number<<" ("<<raw_r1[3]<<"%)\t"<<gv.clean1_stat.gs.t_number<<" ("<<clean_r1[3]<<"%)\t"<<endl;
	of_general_stat<<"Number of base N\t"<<setprecision(15)<<gv.raw1_stat.gs.n_number<<" ("<<raw_r1[4]<<"%)\t"<<gv.clean1_stat.gs.n_number<<" ("<<clean_r1[4]<<"%)\t"<<endl;
	of_general_stat<<"Q20 number\t"<<setprecision(15)<<gv.raw1_stat.gs.q20_num<<" ("<<raw_r1[5]<<"%)\t"<<gv.clean1_stat.gs.q20_num<<" ("<<clean_r1[5]<<"%)"<<endl;
	//of_general_stat<<"Q20 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q20_num/gv.raw1_stat.gs.base_number<<"\t"<<(float)gv.clean1_stat.gs.q20_num/gv.clean1_stat.gs.base_number<<"\t"<<endl;
	of_general_stat<<"Q30 number\t"<<setprecision(15)<<gv.raw1_stat.gs.q30_num<<" ("<<raw_r1[6]<<"%)\t"<<gv.clean1_stat.gs.q30_num<<" ("<<clean_r1[6]<<"%)"<<endl;
	//of_general_stat<<"Q30 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q30_num/gv.raw1_stat.gs.base_number<<"\t"<<(float)gv.clean1_stat.gs.q30_num/gv.clean1_stat.gs.base_number<<"\t"<<endl;
	of_general_stat.close();

	of_readPos_base_stat1<<"Pos\tA\tC\tG\tT\tN\tclean A\tclean C\tclean G\tclean T\tclean N"<<endl;
	for(int i=0;i<gv.raw1_stat.gs.read_length;i++){
		of_readPos_base_stat1<<i+1<<"\t";
		string base_set="ACGTN";
		float raw1_cur_pos_total_base=0;
		float clean1_cur_pos_total_base=0;
		for(int j=0;j!=base_set.size();j++){
			raw1_cur_pos_total_base+=gv.raw1_stat.bs.position_acgt_content[i][j];
			clean1_cur_pos_total_base+=gv.clean1_stat.bs.position_acgt_content[i][j];
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
	}
	of_readPos_base_stat1.close();

	of_readPos_qual_stat1<<"#raw fastq1 quality distribution"<<endl;
	of_readPos_qual_stat1<<"Pos\t";
	int max_qual=0;
	for(int i=0;i<gv.raw1_stat.gs.read_length;i++){
		for(int j=1;j<=MAX_QUAL;j++){
			if(gv.raw1_stat.qs.position_qual[i][j]>0)
				max_qual=j;
		}
	}
	for(int i=0;i<=max_qual;i++){
		of_readPos_qual_stat1<<"Q"<<i<<"\t";
	}
	of_readPos_qual_stat1<<"Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile"<<endl;
	
	float raw1_q20[gv.raw1_stat.gs.read_length],raw1_q30[gv.raw1_stat.gs.read_length];
	float clean1_q20[gv.raw1_stat.gs.read_max_length],clean1_q30[gv.raw1_stat.gs.read_max_length];
	for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
		of_readPos_qual_stat1<<i+1<<"\t";
		unsigned long long raw1_q20_num(0),raw1_q30_num(0),raw1_total(0);
		for(int j=0;j<=max_qual;j++){
			if(j>=20){
				raw1_q20_num+=gv.raw1_stat.qs.position_qual[i][j];
			}
			if(j>=30){
				raw1_q30_num+=gv.raw1_stat.qs.position_qual[i][j];
			}
			raw1_total+=gv.raw1_stat.qs.position_qual[i][j];
			of_readPos_qual_stat1<<setiosflags(ios::fixed);
			of_readPos_qual_stat1<<setprecision(0)<<gv.raw1_stat.qs.position_qual[i][j]<<"\t";
		}
		raw1_q20[i]=(float)raw1_q20_num/raw1_total;
		raw1_q30[i]=(float)raw1_q30_num/raw1_total;
		quartile_result raw1_quar=cal_quar_from_array(gv.raw1_stat.qs.position_qual[i],max_qual+1);
		of_readPos_qual_stat1<<setiosflags(ios::fixed)<<setprecision(2)<<raw1_quar.mean<<"\t";
		of_readPos_qual_stat1<<setprecision(0)<<raw1_quar.median<<"\t"<<raw1_quar.lower_quar<<"\t"<<raw1_quar.upper_quar<<"\t"<<raw1_quar.first10_quar<<"\t"<<raw1_quar.last10_quar<<endl;
	}
	of_readPos_qual_stat1<<"#clean fastq1 quality distribution"<<endl;
	of_readPos_qual_stat1<<"Pos\t";
	for(int i=0;i<=max_qual;i++){
		of_readPos_qual_stat1<<"Q"<<i<<"\t";
	}
	of_readPos_qual_stat1<<"Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile"<<endl;


	of_q2030_stat1<<"Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\tPercentage of Clean Q20+\tPercentage of Clean Q30+"<<endl;
	for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
		of_readPos_qual_stat1<<i+1<<"\t";
		unsigned long long clean1_q20_num(0),clean1_q30_num(0),clean1_total(0);
		for(int j=0;j<=max_qual;j++){
			if(j>=20){
				clean1_q20_num+=gv.clean1_stat.qs.position_qual[i][j];
			}
			if(j>=30){
				clean1_q30_num+=gv.clean1_stat.qs.position_qual[i][j];
			}
			of_readPos_qual_stat1<<setiosflags(ios::fixed);
			clean1_total+=gv.clean1_stat.qs.position_qual[i][j];
			of_readPos_qual_stat1<<setprecision(0)<<gv.clean1_stat.qs.position_qual[i][j]<<"\t";
		}
		clean1_q20[i]=(float)clean1_q20_num/clean1_total;
		clean1_q30[i]=(float)clean1_q30_num/clean1_total;
		quartile_result clean1_quar=cal_quar_from_array(gv.clean1_stat.qs.position_qual[i],max_qual+1);
		of_readPos_qual_stat1<<setiosflags(ios::fixed)<<setprecision(2)<<clean1_quar.mean<<"\t";
		of_readPos_qual_stat1<<setprecision(0)<<clean1_quar.median<<"\t"<<clean1_quar.lower_quar<<"\t"<<clean1_quar.upper_quar<<"\t"<<clean1_quar.first10_quar<<"\t"<<clean1_quar.last10_quar<<endl;
		of_q2030_stat1<<i+1<<setiosflags(ios::fixed)<<setprecision(4)<<"\t"<<raw1_q20[i]<<"\t"<<raw1_q30[i]<<"\t"<<clean1_q20[i]<<"\t"<<clean1_q30[i]<<endl;
	}
	of_readPos_qual_stat1.close();
	of_q2030_stat1.close();
	of_trim_stat1<<"Pos\tHeadLowQual\tHeadFixLen\tTailAdapter\tTailLowQual\tTailFixLen\tCleanHeadLowQual\tCleanHeadFixLen\tCleanTailAdapter\tCleanTailLowQual\tCleanTailFixLen"<<endl;
	long long head_total1(0),tail_total1(0);
	long long head_total_clean1(0),tail_total_clean1(0);
	for(int i=0;i<gv.raw1_stat.gs.read_length;i++){
		head_total1+=gv.raw1_stat.ts.ht[i]+gv.raw1_stat.ts.hlq[i];
		tail_total1+=gv.raw1_stat.ts.ta[i]+gv.raw1_stat.ts.tlq[i]+gv.raw1_stat.ts.tt[i];
		head_total_clean1+=gv.clean1_stat.ts.ht[i]+gv.clean1_stat.ts.hlq[i];
		tail_total_clean1+=gv.clean1_stat.ts.ta[i]+gv.clean1_stat.ts.tlq[i]+gv.clean1_stat.ts.tt[i];
		//of_trim_stat1<<i+1<<"\t"<<gv.trim1_stat.ts.hlq[i]<<"\t"<<gv.trim1_stat.ts.ht[i]<<"\t"<<gv.trim1_stat.ts.ta[i]
	}
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
	}
	of_trim_stat1.close();
}
void seProcess::update_stat(C_fastq_file_stat& fq1s_stat,C_filter_stat& fs_stat,string type){
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
		string base_set="ACGTN";	//base content and quality along read position stat
		int max_qual=0;
		for(int i=0;i!=gv.raw1_stat.gs.read_length;i++){
			for(int j=0;j!=base_set.size();j++){
				gv.raw1_stat.bs.position_acgt_content[i][j]+=fq1s_stat.bs.position_acgt_content[i][j];
			}
		}
		for(int i=1;i<=gv.raw1_stat.gs.read_length;i++){
			gv.raw1_stat.ts.ht[i]+=fq1s_stat.ts.ht[i];
			gv.raw1_stat.ts.hlq[i]+=fq1s_stat.ts.hlq[i];
			gv.raw1_stat.ts.tt[i]+=fq1s_stat.ts.tt[i];
			gv.raw1_stat.ts.tlq[i]+=fq1s_stat.ts.tlq[i];
			gv.raw1_stat.ts.ta[i]+=fq1s_stat.ts.ta[i];
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
			}
		}
		
		//gv.fs.output_reads_num+=fs_stat.output_reads_num;	//filter stat
		gv.fs.in_adapter_list_num+=fs_stat.in_adapter_list_num;
		gv.fs.include_adapter_seq_num+=fs_stat.include_adapter_seq_num;
		gv.fs.n_ratio_num+=fs_stat.n_ratio_num;
		gv.fs.highA_num+=fs_stat.highA_num;
		gv.fs.polyX_num+=fs_stat.polyX_num;
		gv.fs.tile_num+=fs_stat.tile_num;
		gv.fs.low_qual_base_ratio_num+=fs_stat.low_qual_base_ratio_num;
		gv.fs.mean_quality_num+=fs_stat.mean_quality_num;
		gv.fs.short_len_num+=fs_stat.short_len_num;
		gv.fs.long_len_num+=fs_stat.long_len_num;
		gv.fs.include_contam_seq_num+=fs_stat.include_contam_seq_num;
		gv.fs.include_global_contam_seq_num+=fs_stat.include_global_contam_seq_num;
	}else if(type=="trim"){
		//gv.trim1_stat.gs.read_length=fq1s_stat.gs.read_length;	//generate stat
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

		int max_qual(0);
		string base_set="ACGTN";	//base content and quality along read position stat
		for(int i=0;i!=gv.trim1_stat.gs.read_max_length;i++){
			for(int j=0;j!=base_set.size();j++){
				gv.trim1_stat.bs.position_acgt_content[i][j]+=fq1s_stat.bs.position_acgt_content[i][j];
			}
		}
		for(int i=0;i!=gv.trim1_stat.gs.read_max_length;i++){
			gv.trim1_stat.ts.ht[i]+=fq1s_stat.ts.ht[i];
			gv.trim1_stat.ts.hlq[i]+=fq1s_stat.ts.hlq[i];
			gv.trim1_stat.ts.tt[i]+=fq1s_stat.ts.tt[i];
			gv.trim1_stat.ts.tlq[i]+=fq1s_stat.ts.tlq[i];
			gv.trim1_stat.ts.ta[i]+=fq1s_stat.ts.ta[i];
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
			}
		}
	}else if(type=="clean"){
		//gp.output_reads_num+=gv.clean1_stat.gs.reads_number;
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
		
		string base_set="ACGTN";	//base content and quality along read position stat
		int max_qual(0);
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
			for(int j=0;j!=base_set.size();j++){
				gv.clean1_stat.bs.position_acgt_content[i][j]+=fq1s_stat.bs.position_acgt_content[i][j];
			}
		}
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
			gv.clean1_stat.ts.ht[i]+=fq1s_stat.ts.ht[i];
			gv.clean1_stat.ts.hlq[i]+=fq1s_stat.ts.hlq[i];
			gv.clean1_stat.ts.tt[i]+=fq1s_stat.ts.tt[i];
			gv.clean1_stat.ts.tlq[i]+=fq1s_stat.ts.tlq[i];
			gv.clean1_stat.ts.ta[i]+=fq1s_stat.ts.ta[i];
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
			}
		}
	}else{
		cerr<<"Error:code error"<<endl;
		exit(1);
	}
	

}
void* seProcess::stat_se_fqs(SEstatOption opt){
	opt.stat1->gs.reads_number+=opt.fq1s->size();
	for(vector<C_fastq>::iterator ix=opt.fq1s->begin();ix!=opt.fq1s->end();ix++){
		if((*ix).head_hdcut>0 || (*ix).head_lqcut>0){
			if((*ix).head_hdcut>=(*ix).head_lqcut){
				opt.stat1->ts.ht[(*ix).head_hdcut]++;
			}else{
				opt.stat1->ts.hlq[(*ix).head_lqcut]++;
			}
		}
		if((*ix).tail_hdcut>0 || (*ix).tail_lqcut>0 || (*ix).adacut_pos>=0){
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
			if(base_quality>MAX_QUAL){
				cerr<<"Error:quality is too high,please check the quality system parameter or fastq file"<<endl;
				exit(1);
			}
			if(base_quality<MIN_QUAL){
				cerr<<"Error:quality is too low,please check the quality system parameter or fastq file"<<endl;
				exit(1);
			}
			opt.stat1->qs.position_qual[i][base_quality]++;
			if(base_quality>=20)
				opt.stat1->gs.q20_num++;
			if(base_quality>=30)
				opt.stat1->gs.q30_num++;
		}
		opt.stat1->gs.read_length=(*ix).sequence.size();
		opt.stat1->gs.base_number+=opt.stat1->gs.read_length;
	}
}

void seProcess::filter_se_fqs(SEcalOption opt){
	//C_reads_trim_stat cut_pos;
	for(vector<C_fastq>::iterator i=opt.fq1s->begin();i!=opt.fq1s->end();i++){
		C_single_fastq_filter se_fastq_filter=C_single_fastq_filter(*i,gp);
		se_fastq_filter.se_trim(gp);
		if(gp.adapter_discard_or_trim=="trim" || gp.contam_discard_or_trim=="trim" || !gp.trim.empty() || !gp.trimBadHead.empty() || !gp.trimBadTail.empty()){
			(*i).head_hdcut=se_fastq_filter.read.head_hdcut;
			(*i).head_lqcut=se_fastq_filter.read.head_lqcut;
			(*i).tail_hdcut=se_fastq_filter.read.tail_hdcut;
			(*i).tail_lqcut=se_fastq_filter.read.tail_lqcut;
			(*i).adacut_pos=se_fastq_filter.read.adacut_pos;
			//(*i).contam_pos=se_fastq_filter.read.contam_pos;
			//(*i).global_contam_pos=se_fastq_filter.read.global_contam_pos;
			//(*i).raw_length=se_fastq_filter.read.raw_length;	
		}
		//*i=se_fastq_filter.read;
		if(!gp.trim_fq1.empty()){
			preOutput(1,se_fastq_filter.read);
			opt.trim_result1->push_back(se_fastq_filter.read);
		}
		int whether_discard(0);
		if(gp.module_name=="filtersRNA"){
			whether_discard=se_fastq_filter.sRNA_discard(opt.se_local_fs,gp);
		}else{
			whether_discard=se_fastq_filter.se_discard(opt.se_local_fs,gp);
		}
		if(whether_discard!=1){
			if(!gp.clean_fq1.empty()){
				preOutput(1,se_fastq_filter.read);
				opt.clean_result1->push_back(se_fastq_filter.read);
			}
		}
	}
	//return cut_pos;
}

void  seProcess::preOutput(int type,C_fastq& a){
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
void seProcess::seWrite_split(vector<C_fastq>& se1){
	output_split_fastqs("1",se1);
}
void seProcess::seWrite(vector<C_fastq>& se1,string type,gzFile out1){
	output_fastqs("1",se1,out1);
}
void seProcess::C_fastq_init(C_fastq& a){
	a.seq_id="";
	a.sequence="";
	a.qual_seq="";
	a.adapter_seq=gp.adapter1_seq;
	a.contam_seq=gp.contam1_seq;
	//a.global_contams=gp.global_contams;
	a.head_hdcut=-1;
	a.head_lqcut=-1;
	a.tail_hdcut=-1;
	a.tail_lqcut=-1;
	a.adacut_pos=-1;
	a.contam_pos=-1;
	a.global_contam_pos=-1;
	a.raw_length=0;
	if(gp.module_name=="filtersRNA"){
		a.adapter_seq2=gp.adapter2_seq;
	}
	if(!gp.trim.empty()){
		vector<string> tmp_eles=get_pe_hard_trim(gp.trim);
		a.head_trim_len=tmp_eles[0];
		a.tail_trim_len=tmp_eles[1];
	}
}
int seProcess::read(vector<C_fastq>& pe1,ifstream& infile1){
	
	string buf1;
	int file1_line_num(0);
	C_fastq fastq1;
	C_fastq_init(fastq1);
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
				pe1.push_back(fastq1);
				//fq1s.push_back(fastq1);
			}
		}else{
			return -1;
		}
		
	}
}
int seProcess::read_gz(vector<C_fastq>& se1){
	char buf1[READBUF];
	C_fastq fastq1;
	C_fastq_init(fastq1);
	int file1_line_num(0);
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
				se1.push_back(fastq1);
			}
		}else{
			return -1;
		}
	}
}
void seProcess::create_thread_outputFile(int index){
	if(!gp.trim_fq1.empty()){	//create output trim files handle
		ostringstream trim_outfile1;
		trim_outfile1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".trim.r1.fq.gz";
		gz_trim_out1[index]=gzopen(trim_outfile1.str().c_str(),"wb");
		gzsetparams(gz_trim_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1[index],1024*1024*16);
	}
	if(!gp.clean_fq1.empty()){	//create output clean files handle
		if(gp.output_clean<=0){
			ostringstream outfile1;
			outfile1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r1.fq.gz";
			gz_clean_out1[index]=gzopen(outfile1.str().c_str(),"wb");
			gzsetparams(gz_clean_out1[index], 2, Z_DEFAULT_STRATEGY);
			gzbuffer(gz_clean_out1[index],1024*1024*16);
		}
	}
}
void seProcess::create_thread_trimoutputFile(int index){
	if(!gp.trim_fq1.empty()){	//create output trim files handle
		ostringstream trim_outfile1;
		trim_outfile1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".trim.r1.fq.gz";
		gz_trim_out1[index]=gzopen(trim_outfile1.str().c_str(),"wb");
		gzsetparams(gz_trim_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1[index],1024*1024*16);
	}
}
void seProcess::create_thread_cleanoutputFile(int index){
	if(!gp.clean_fq1.empty()){	//create output clean files handle
		if(gp.output_clean<=0){
			ostringstream outfile1;
			outfile1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r1.fq.gz";
			gz_clean_out1[index]=gzopen(outfile1.str().c_str(),"wb");
			gzsetparams(gz_clean_out1[index], 2, Z_DEFAULT_STRATEGY);
			gzbuffer(gz_clean_out1[index],1024*1024*16);
		}
	}
}
void* seProcess::sub_thread(int index){
	off_t start_pos=t_start_pos[index];
	off_t end_pos=t_end_pos[index];
	string head,tail;
	int flag(0);
	int self_fq1fd=open(se_new_fq1_path.c_str(),O_RDONLY);
	int min_len=1024*4*10;
	if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
		//create_thread_outputFile(index);
		create_thread_trimoutputFile(index);
		create_thread_cleanoutputFile(index);
	}else if(!gp.trim_fq1.empty()){
		create_thread_trimoutputFile(index);
	}
	//create_thread_outputFile(index);
	
	//char *fq1_buf,*fq2_buf;
	unsigned long long mmap_start(start_pos);
	unsigned long long copysz(se_buffer);
	int tmp_iter(0);
	char *buf1;
	C_fastq fastq1;
	C_fastq_init(fastq1);
	vector<C_fastq> fq1s;
	vector<C_fastq> trim_result1,clean_result1;
	while(1){
		if(mmap_start>=end_pos)
			break;
		if(end_pos - mmap_start<se_buffer){
			copysz=end_pos - mmap_start;
		}else{
			if(end_pos - mmap_start -se_buffer <se_buffer/4){
				copysz=end_pos - mmap_start;
			}else{
				copysz=se_buffer;
			}
		}
		if((buf1=(char*)mmap(0,copysz,PROT_READ,MAP_SHARED,self_fq1fd,mmap_start))==MAP_FAILED){
			cerr<<"Error:mmap error"<<endl;
			exit(1);
		}
		string tmp_head1,tmp_tail1;
		int head_idx=0;
		for(;head_idx<copysz;head_idx++){
			if(head_idx==0 || head_idx==1){
				tmp_head1.insert(tmp_head1.end(),buf1[head_idx]);
			}else{
				if(buf1[head_idx]=='@' && buf1[head_idx-1]=='\n' && buf1[head_idx-2]!='+'){
					break;
				}else{
					tmp_head1.insert(tmp_head1.end(),buf1[head_idx]);
				}
			}
		}
		int tail_idx=copysz-1;
		for(;tail_idx>1;tail_idx--){
			if(tail_idx==copysz-1){
				tmp_tail1.insert(tmp_tail1.begin(),buf1[tail_idx]);
			}else{
				if(buf1[tail_idx]=='\n' && buf1[tail_idx+1]=='@' && buf1[tail_idx-1]!='+'){
					break;
				}else{
					tmp_tail1.insert(tmp_tail1.begin(),buf1[tail_idx]);
				}
			}
		}
		se_sticky_reads1[index]+=tmp_head1;
		se_sticky_reads1[index]+=tmp_tail1;
		int line_num(0);
		
		
		for(int i=head_idx;i<=tail_idx;i++){
			if(buf1[i]=='\n'){
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
				fastq1.seq_id.insert(fastq1.seq_id.end(),buf1[i]);
			}
			if(line_num%4==1){
				fastq1.sequence.insert(fastq1.sequence.end(),buf1[i]);
			}
			if(line_num%4==3){
				fastq1.qual_seq.insert(fastq1.qual_seq.end(),buf1[i]);
			}
		}
		thread_process_reads(index,fq1s);
		if(limit_end>0){
			break;
		}
		mmap_start+=copysz;
		munmap(buf1,copysz);
		tmp_iter++;
	}
	close(self_fq1fd);
	if(!gp.trim_fq1.empty()){
		if(index!=0){
			gzclose(gz_trim_out1[index]);
		}
		
	}
	if(!gp.clean_fq1.empty()){
		if(index!=0){
			if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
				gzclose(gz_clean_out1[index]);
			}
		}
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}
void seProcess::run_pigz(){	//split raw files with pigz and "split" command
	ostringstream cmd1;
	int pigz_thread=gp.threads_num>2?gp.threads_num/2:1;
	cmd1<<"pigz -c -d -p "<<pigz_thread<<" "<<gp.fq1_path<<" > "<<gp.output_dir<<"/raw.r1.fq";
	if(system(cmd1.str().c_str())==-1){
		cerr<<"Error:run pigz error"<<endl;
		exit(1);
	}
}
void seProcess::run_pigz_split(int type){
	ostringstream cmd1;
	int pigz_thread=gp.threads_num>16?gp.threads_num:16;
	int split_line_num=gp.split_line*4;
	if(type==1){
		if(gp.fq1_path.rfind(".gz")==gp.fq1_path.size()-3){
			cmd1<<"pigz -c -d -p "<<pigz_thread<<" "<<gp.fq1_path<<" | split -l "<<split_line_num<<" -d - "<<gp.output_dir<<"/tmp"<<random_num<<".r1.fq.";
		}else{
			cmd1<<"split -l "<<split_line_num<<" -d "<<gp.fq1_path<<" "<<gp.output_dir<<"/tmp"<<random_num<<".r1.fq.";
		}
	}else if(type==2){
		cerr<<"Error:pigz&split error"<<endl;
		exit(1);
	}
	if(system(cmd1.str().c_str())==-1){
		cerr<<"Error:pigz&split error"<<endl;
		exit(1);
	}
}
void seProcess::merge_stat(){
	for(int i=0;i!=used_threads_num;i++){
		update_stat(se_local_raw_stat1[i],se_local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(se_local_trim_stat1[i],se_local_fs[i],"trim");
		}
		if(!gp.clean_fq1.empty()){
			update_stat(se_local_clean_stat1[i],se_local_fs[i],"clean");
		}
	}
}
void seProcess::merge_stat(int index){
	for(int i=0;i<=gp.threads_num;i++){
		update_stat(se_local_raw_stat1[i],se_local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(se_local_trim_stat1[i],se_local_fs[i],"trim");
		}
	}
	for(int i=0;i<=index;i++){
		if(!gp.clean_fq1.empty()){
			update_stat(se_local_clean_stat1[i],se_local_fs[i],"clean");
		}
	}
}
void seProcess::merge_stat_nonssd(){
	for(int i=0;i!=gp.threads_num;i++){
		update_stat(se_local_raw_stat1[i],se_local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(se_local_trim_stat1[i],se_local_fs[i],"trim");
		}
		if(!gp.clean_fq1.empty()){
			update_stat(se_local_clean_stat1[i],se_local_fs[i],"clean");
		}
	}
}
void seProcess::merge_data(){
	if(!gp.trim_fq1.empty()){
		string trim_file1;
		trim_file1=gp.output_dir+"/"+gp.trim_fq1;
		if(check_gz_empty(trim_file1)==1){
			string rm_cmd="rm "+trim_file1;
			system(rm_cmd.c_str());
		}
	}
	if(!gp.clean_fq1.empty()){
		string clean_file1;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream trim_out_fq1_tmp;
			trim_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r1.fq.gz";
			
			if(access(trim_out_fq1_tmp.str().c_str(),0)!=-1){

				ostringstream cat_cmd1;
				cat_cmd1<<"cat "<<trim_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq1<<";rm "<<trim_out_fq1_tmp.str();
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
			
		}
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp;
			clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1;
				cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
		}
	}
	
}
void seProcess::merge_trim_data(){
	if(!gp.trim_fq1.empty()){
		string trim_file1;
		trim_file1=gp.output_dir+"/"+gp.trim_fq1;
		if(check_gz_empty(trim_file1)==1){
			string rm_cmd="rm "+trim_file1;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream trim_out_fq1_tmp;
			trim_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r1.fq.gz";
			
			if(access(trim_out_fq1_tmp.str().c_str(),0)!=-1){

				ostringstream cat_cmd1;
				cat_cmd1<<"cat "<<trim_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq1<<";rm "<<trim_out_fq1_tmp.str();
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
			
		}
	}
	
}
void seProcess::merge_clean_data(){
	if(!gp.clean_fq1.empty()){
		string clean_file1;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp;
			clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1;
				cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
		}
	}
	
}
void seProcess::merge_data(int index){	//cat all output files to a single large file in limit output mode
	if(!gp.trim_fq1.empty()){
		string trim_file1;
		trim_file1=gp.output_dir+"/"+gp.trim_fq1;
		if(check_gz_empty(trim_file1)==1){
			string rm_cmd="rm "+trim_file1;
			system(rm_cmd.c_str());
		}
	}
	if(!gp.clean_fq1.empty()){
		string clean_file1;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream trim_out_fq1_tmp;
			trim_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".trim.r1.fq.gz";
			
			if(access(trim_out_fq1_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1;
				cat_cmd1<<"cat "<<trim_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.trim_fq1<<";rm "<<trim_out_fq1_tmp.str();
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
			
		}
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp;
			if(i==index){
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/last.r1.fq.gz";
			}else{
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			}
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1;
				if(i<=index){
					cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
					if(i==index){
						cat_cmd1<<";rm "<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
					}
				}else{
					cat_cmd1<<"rm "<<clean_out_fq1_tmp.str();
				}
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
		}
	}
}
void seProcess::merge_clean_data(int index){	//cat all output files to a single large file in limit output mode
	if(!gp.clean_fq1.empty()){
		string clean_file1;
		clean_file1=gp.output_dir+"/"+gp.clean_fq1;
		if(check_gz_empty(clean_file1)==1){
			string rm_cmd="rm "+clean_file1;
			system(rm_cmd.c_str());
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		if(!gp.clean_fq1.empty()){
			ostringstream clean_out_fq1_tmp;
			if(i==index){
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/last.r1.fq.gz";
			}else{
				clean_out_fq1_tmp<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			}
			if(access(clean_out_fq1_tmp.str().c_str(),0)!=-1){
				ostringstream cat_cmd1;
				if(i<=index){
					cat_cmd1<<"cat "<<clean_out_fq1_tmp.str()<<" >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<clean_out_fq1_tmp.str();
					if(i==index){
						cat_cmd1<<";rm "<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<i<<".clean.r1.fq.gz";
					}
				}else{
					cat_cmd1<<"rm "<<clean_out_fq1_tmp.str();
				}
				if(system(cat_cmd1.str().c_str())==-1){
					cerr<<"Error:cat file error,"<<cat_cmd1.str()<<endl;
				}
			}
		}
	}
}
void seProcess::create_thread_read(int index){
	multi_gzfq1[index]=gzopen((gp.fq1_path).c_str(),"rb");
	gzsetparams(gzfp1, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(gzfp1,2048*2048);
}
void* seProcess::sub_thread_nonssd_realMultiThreads(int index){
	of_log<<get_local_time()<<"\tthread "<<index<<" start"<<endl;
	if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false)){
		//create_thread_outputFile(index);
		create_thread_trimoutputFile(index);
		create_thread_cleanoutputFile(index);
	}else if(!gp.trim_fq1.empty()){
		create_thread_trimoutputFile(index);
	}
	//create_thread_outputFile(index);
	create_thread_read(index);
	
	char buf1[READBUF];
	C_fastq fastq1;
	C_fastq_init(fastq1);
	unsigned long long file1_line_num(0);
	unsigned long long block_line_num1(0);
	int thread_read_block=4*gp.patchSize;
	vector<C_fastq> fq1s;
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
					fq1s.push_back(fastq1);
					if(fq1s.size()==gp.patchSize){
						if(index==0)
							of_log<<get_local_time()<<" processed_reads:\t"<<file1_line_num/4<<endl;
						thread_process_reads(index,fq1s);
						if(limit_end>0){
							break;
						}
					}
				}
			}
			file1_line_num++;
		}else{
			gzclose(multi_gzfq1[index]);
			if(fq1s.size()>0){
				thread_process_reads(index,fq1s);
				if(limit_end>0){
					break;
				}
			}
			break;
		}
	}
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1[index]);
	}
	if(!gp.clean_fq1.empty()){
		if(gp.output_clean<=0 && !(gp.total_reads_num>0 && gp.total_reads_num_random==false))
			gzclose(gz_clean_out1[index]);
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}
void seProcess::process_nonssd(){
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
	//thread read_monitor(bind(&seProcess::monitor_read_thread,this));
	//sleep(10);
	for(int i=0;i<gp.threads_num;i++){
		//t_array[i]=thread(bind(&seProcess::sub_thread_nonssd_multiOut,this,i));
		t_array[i]=thread(bind(&seProcess::sub_thread_nonssd_realMultiThreads,this,i));
	}
	for(int i=0;i<gp.threads_num;i++){
		t_array[i].join();
	}

	if(gp.total_reads_num_random==true && gp.total_reads_num>0){
		run_extract_random();
		remove_tmpDir();
		of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
		of_log.close();
		return;
	}
	merge_stat_nonssd();
	if(!gp.trim_fq1.empty()){
		merge_trim_data();
	}
	print_stat();
	if(gp.output_clean<=0 && !(gp.l_total_reads_num>0 && gp.total_reads_num_random==false)){
		merge_clean_data();
		remove_tmpDir();
	}else{
		if(!gp.trim_fq1.empty()){
			remove_tmpDir();
		}
		if(gp.output_clean>0){
			gzclose(gz_fq_se);
		}
	}
	//read_monitor.join();
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
}
/*void seProcess::add_raw_trim(C_fastq_file_stat& a,C_reads_trim_stat& b){
	//unsigned long long hlq[READ_MAX_LEN],ht[READ_MAX_LEN];
	//unsigned long long ta[READ_MAX_LEN],tlq[READ_MAX_LEN],tt[READ_MAX_LEN];
	for(int i=1;i<=a.gs.read_length;i++){
		a.ts.hlq[i]+=b.hlq[i];
		a.ts.ht[i]+=b.ht[i];
		a.ts.ta[i]+=b.ta[i];
		a.ts.tlq[i]+=b.tlq[i];
		a.ts.tt[i]+=b.tt[i];
	}
}*/
void seProcess::thread_process_reads(int index,vector<C_fastq> &fq1s){
	vector<C_fastq> trim_result1,clean_result1;
	
	SEcalOption opt2;
	opt2.se_local_fs=&se_local_fs[index];
	opt2.fq1s=&fq1s;
	opt2.trim_result1=&trim_result1;
	opt2.clean_result1=&clean_result1;
	filter_se_fqs(opt2);		//filter raw fastqs by the given parameters
	SEstatOption opt_raw;
	opt_raw.fq1s=&fq1s;
	opt_raw.stat1=&se_local_raw_stat1[index];
	stat_se_fqs(opt_raw);		//statistic raw fastqs
	fq1s.clear();
	//add_raw_trim(se_local_raw_stat1[index],raw_cut);
	
	SEstatOption opt_trim,opt_clean;
	if(!gp.trim_fq1.empty()){	//trim means only trim but not discard.
		opt_trim.fq1s=&trim_result1;
		opt_trim.stat1=&se_local_trim_stat1[index];
		stat_se_fqs(opt_trim);	//statistic trim fastqs
	}
	/*
	if(!gp.clean_fq1.empty()){
		opt_clean.fq1s=&clean_result1;
		opt_clean.stat1=&se_local_clean_stat1[index];
		stat_se_fqs(opt_clean);	//statistic clean fastqs
	}
	*/
	//
	if(!gp.trim_fq1.empty()){
		seWrite(trim_result1,"trim",gz_trim_out1[index]);	//output trim files
		trim_result1.clear();
	}
	if(!gp.clean_fq1.empty()){
		opt_clean.stat1=&se_local_clean_stat1[index];
		if(gp.output_clean>0 || (gp.total_reads_num_random==false && gp.l_total_reads_num>0)){
			se_write_m.lock();
			if(limit_end>0){
				se_write_m.unlock();
				clean_result1.clear();
				return;
			}
			seWrite_split(clean_result1);
			if(limit_end==1){
				int to_remove=gp.have_output1-gp.clean_file_reads;
				//cout<<to_remove<<"\t"<<clean_result1.size()<<"\t"<<gp.have_output1<<endl;
				if(to_remove>=clean_result1.size()){
					se_write_m.unlock();
					clean_result1.clear();
					return;
				}else{
					clean_result1.erase(clean_result1.end()-to_remove,clean_result1.end());
				}
			}
			se_write_m.unlock();
		}else{
			seWrite(clean_result1,"clean",gz_clean_out1[index]);	//output clean files
		}
		opt_clean.fq1s=&clean_result1;
		stat_se_fqs(opt_clean);	//statistic clean fastqs
		clean_result1.clear();
		/*thread_write_m[index].lock();
		thread pewrite_t(bind(&seProcess::peWrite,this,clean_result1,clean_result2,"clean",gz_clean_out1[index],gz_clean_out2[index]));
		thread_write_m[index].unlock();*/
		if(gp.is_streaming){
			se_write_m.lock();
			C_global_variable tmp_gv;
			tmp_gv.fs=*(opt2.se_local_fs);
			tmp_gv.raw1_stat=*(opt_raw.stat1);
			//tmp_gv.trim1_stat=local_trim_stat1[index];
			//tmp_gv.trim2_stat=local_trim_stat2[index];
			tmp_gv.clean1_stat=*(opt_clean.stat1);
			seStreaming_stat(tmp_gv);
			se_write_m.unlock();
		}
	}
	//
}
void seProcess::run_extract_random(){
	if(gp.total_reads_num<=0 || gp.total_reads_num_random==false){
		cerr<<"Error:extract random clean reads error cuz parameters are wrong"<<endl;
		exit(1);
	}
	unsigned long long total_clean_reads(0);
	for(int i=0;i!=gp.threads_num;i++){
		total_clean_reads+=se_local_clean_stat1[i].gs.reads_number;
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
		cur_total+=se_local_clean_stat1[i].gs.reads_number;
		//cout<<se_local_clean_stat1[i].gs.reads_number<<"\t"<<cur_total<<"\t"<<gp.l_total_reads_num<<endl;
		if(cur_total>gp.l_total_reads_num){
			last_thread=i;
			sticky_end=se_local_clean_stat1[i].gs.reads_number-(cur_total-gp.l_total_reads_num);
			break;
		}
	}
	
	//create the last patch clean fq file and stat
	//cout<<"last thread\t"<<last_thread<<endl;
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
void seProcess::process_some_reads(int index,int out_number){
	//open target fq file
	ostringstream target_file_fq1;
	target_file_fq1<<gp.output_dir<<"/"<<tmp_dir<<"/thread."<<index<<".clean.r1.fq.gz";
	int file_reads_number=se_local_clean_stat1[index].gs.reads_number;
	se_local_clean_stat1[index].clear();
	if(file_reads_number==out_number){
		return;
	}
	int times=(int)floor(file_reads_number/out_number);
	gzFile tmp_fq1=gzopen(target_file_fq1.str().c_str(),"rb");
	gzsetparams(tmp_fq1, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(tmp_fq1,2048*2048);
	//set output file
	string last_file1=gp.output_dir+"/"+tmp_dir+"/last.r1.fq.gz";
	gzFile tmp_out_fq1=gzopen(last_file1.c_str(),"wb");
	gzsetparams(tmp_out_fq1, 2, Z_DEFAULT_STRATEGY);
	gzbuffer(tmp_out_fq1,2048*2048);

	char buf1[READBUF];
	C_fastq fastq1;
	C_fastq_init(fastq1);
	unsigned long long file1_line_num(0);
	vector<C_fastq> fq1s;
	int processed_number(0),rest_number(1000000);
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
					fq1s.push_back(fastq1);
					if(fq1s.size()==gp.patchSize || fq1s.size()==rest_number){
						of_log<<get_local_time()<<" last sticky thread processed reads:\t"<<file1_line_num/4<<endl;
						//thread_process_reads(index,fq1s,fq2s);
						processed_number+=fq1s.size();
						rest_number=out_number-processed_number;
						//cout<<"processed number\t"<<processed_number<<"\trest number\t"<<rest_number<<endl;
						limit_process_reads(index,fq1s,tmp_out_fq1);

					}
				}
			}
			file1_line_num++;
		}else{
			if(fq1s.size()>0){
				limit_process_reads(index,fq1s,tmp_out_fq1);
			}
			break;
		}
	}
	gzclose(tmp_fq1);
	gzclose(tmp_out_fq1);
}

void seProcess::limit_process_reads(int index,vector<C_fastq> &fq1s,gzFile gzfq1){
	SEstatOption opt_clean;
	opt_clean.fq1s=&fq1s;
	opt_clean.stat1=&se_local_clean_stat1[index];
	stat_se_fqs(opt_clean);		//statistic raw fastqs
	seWrite(fq1s,"clean",gzfq1);
	fq1s.clear();
}
void seProcess::process(){
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
	se_new_fq1_path=gp.fq1_path;
	if(gp.fq1_path.rfind(".gz")==gp.fq1_path.size()-3){
		run_pigz();
		se_new_fq1_path=gp.output_dir+"/raw.r1.fq";
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
	fq1fd=open(se_new_fq1_path.c_str(),O_RDONLY);
	off_t file_size;
	struct stat st;
	fstat(fq1fd,&st);
	file_size=st.st_size;
	int pieces=file_size/se_buffer;
	if(file_size%se_buffer!=0 || pieces==0)
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
		t_start_pos[i]=i*real_block*se_buffer;
		if(i==used_threads_num-1){
			t_end_pos[i]=file_size;
		}else{
			t_end_pos[i]=(i+1)*real_block*se_buffer;
		}
	}	
	close(fq1fd);
	thread t_array[used_threads_num];
	for(int i=0;i!=used_threads_num;i++){
		t_array[i]=thread(bind(&seProcess::sub_thread,this,i));
	}
	for(int i=0;i<used_threads_num;i++){
		t_array[i].join();
	}
	if(limit_end<=0){
		string fq1_unprocess;
		for(int i=0;i<used_threads_num;i++){
			fq1_unprocess+=se_sticky_reads1[i];
		}
		int line_num(0);
		C_fastq fastq1;
		C_fastq_init(fastq1);
		vector<C_fastq> fq1s;
		//cout<<fq1_unprocess<<endl;
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
		long long clean_size(0);
		thread_process_reads(0,fq1s);
		
	}
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1[0]);
	}
	if(!gp.clean_fq1.empty()){
		if(gp.output_clean<=0 && !(gp.l_total_reads_num>0 && gp.total_reads_num_random==false))
			gzclose(gz_clean_out1[0]);
	}

	//
	if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
		string se_new_fq1_path=gp.output_dir+"/raw.r1.fq";
		string rm_cmd="rm "+se_new_fq1_path;
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
		if(gp.output_clean>0){
			gzclose(gz_fq_se);
		}
	}

	
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
}
void seProcess::remove_tmpDir(){
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
void seProcess::make_tmpDir(){
	srand(time(0));
	ostringstream tmp_str;
	for(int i=0;i!=6;i++){
		int tmp_rand=random(26)+'A';
		tmp_str<<(char)tmp_rand;
	}
	tmp_dir=tmp_str.str();
	string mkdir_str="mkdir -p "+gp.output_dir+"/"+tmp_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir error,"<<mkdir_str<<endl;
		exit(1);
	}
}
void seProcess::output_fastqs2(int type,vector<C_fastq> &fq1,ofstream& outfile){
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
void seProcess::output_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile){
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
		gzflush(outfile,1);
	}
	//m.unlock();
}
void seProcess::output_split_fastqs(string type,vector<C_fastq> &fq1){
	//m.lock();
	
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
		//cout<<gp.output_clean<<endl;
		gp.have_output1+=fq1.size();
		int idx=(gp.have_output1-1)/gp.clean_file_reads;
		se_out[idx]++;
		int mod=gp.have_output1%gp.clean_file_reads;
		if(se_out[idx]==1){
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
				gzwrite(gz_fq_se,sticky_tail.c_str(),sticky_tail.size());
				gzclose(gz_fq_se);
				limit_end++;
				return;
			}
			ostringstream out_fq1;
			out_fq1<<gp.output_dir<<"/split."<<idx<<".clean.r1.fq.gz";
			if(gp.total_reads_num_random==false && gp.l_total_reads_num>0){
				string fq1_whole_path=gp.output_dir+"/"+gp.clean_fq1;
				gz_fq_se=gzopen(fq1_whole_path.c_str(),"wb");
			}else{
				gz_fq_se=gzopen(out_fq1.str().c_str(),"wb");
			}
			gzsetparams(gz_fq_se, 2, Z_DEFAULT_STRATEGY);
			gzbuffer(gz_fq_se,1024*1024*10);

			if(!sticky_head.empty()){
				gzwrite(gz_fq_se,sticky_head.c_str(),sticky_head.size());
			}
		}else{
			gzwrite(gz_fq_se,out_content.c_str(),out_content.size());
		}
	}
	//m.unlock();
}
void seProcess::seStreaming_stat(C_global_variable& local_gv){
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
}