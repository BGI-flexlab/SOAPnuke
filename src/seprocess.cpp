#include <iostream>
#include <string>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <thread>
#include <mutex>
#include <fstream>
#include <iomanip>
#include "seprocess.h"
#include "process_argv.h"
#include "zlib.h"
#include "gc.h"
#include "sequence.h"
using namespace::std;
#define READBUF 500
mutex se_read_m,se_stat_m,se_write_m;
C_filter_stat se_local_fs[max_thread];
C_fastq_file_stat se_local_raw_stat1[max_thread],se_local_raw_stat2[max_thread],se_local_trim_stat1[max_thread],se_local_trim_stat2[max_thread],se_local_clean_stat1[max_thread],se_local_clean_stat2[max_thread];
gzFile gz_fq_se;
int exceed_output_se(0);
seProcess::seProcess(C_global_parameter m_gp){
	gp=m_gp;
	gv=C_global_variable();
	used_threads_num=0;
	processed_reads=0;
	srand((unsigned)time(NULL));
	ostringstream tmpstring;
	tmpstring<<rand()%100;
	random_num=tmpstring.str();
	of_log.open(gp.log.c_str());
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

	map<string,unsigned long long> filter_number;
	filter_number["Reads with contam sequence"]=gv.fs.include_contam_seq_num;
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
	filter_number["Reads limited to output number"]=gv.fs.output_reads_num;
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
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
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
		
		gv.fs.output_reads_num+=fs_stat.output_reads_num;	//filter stat
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
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
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
		gp.output_reads_num+=gv.clean1_stat.gs.reads_number;
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
		if((*ix).tail_hdcut>0 || (*ix).tail_lqcut>0 || (*ix).adacut_pos>0){
			if((*ix).tail_hdcut>=(*ix).tail_lqcut){
				if((*ix).tail_hdcut>=(*ix).adacut_pos){
					opt.stat1->ts.tt[(*ix).sequence.size()-(*ix).tail_hdcut+1]++;
				}else{
					opt.stat1->ts.ta[(*ix).sequence.size()-(*ix).adacut_pos+1]++;
				}
			}else{
				if((*ix).tail_lqcut>=(*ix).adacut_pos){
					opt.stat1->ts.tlq[(*ix).sequence.size()-(*ix).tail_lqcut+1]++;
				}else{
					opt.stat1->ts.ta[(*ix).sequence.size()-(*ix).adacut_pos+1]++;
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
void* seProcess::filter_se_fqs(SEcalOption opt){
	for(vector<C_fastq>::iterator i=opt.fq1s->begin();i!=opt.fq1s->end();i++){
		C_single_fastq_filter se_fastq_filter=C_single_fastq_filter(*i,gp);
		se_fastq_filter.se_trim(gp);
		if(!gp.trim_fq1.empty()){
			preOutput(1,*i);
			opt.trim_result1->push_back(*i);
		}
		int whether_discard(0);
		if(gp.module_name=="filtersRNA"){
			whether_discard=se_fastq_filter.sRNA_discard(opt.se_local_fs,gp);
		}else{
			whether_discard=se_fastq_filter.se_discard(opt.se_local_fs,gp);
		}
		if(whether_discard!=1){
			if(!gp.clean_fq1.empty()){
				preOutput(1,*i);
				opt.clean_result1->push_back(*i);
			}
		}
	}
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
void seProcess::seWrite(vector<C_fastq>& pe1,string type,gzFile out1){
	if(type=="trim"){
		output_fastqs("1",pe1,out1);
	}else if(type=="clean"){
		if(gp.mode=="nonssd" && gp.output_clean>0){
			output_split_fastqs("1",pe1,out1);
		}else{
			output_fastqs("1",pe1,out1);
		}
	}
}
void seProcess::C_fastq_init(C_fastq& a){
	a.seq_id="";
	a.sequence="";
	a.qual_seq="";
	a.adapter_seq=gp.adapter1_seq;
	a.contam_seq=gp.contam1_seq;
	a.head_hdcut=-1;
	a.head_lqcut=-1;
	a.tail_hdcut=-1;
	a.tail_lqcut=-1;
	a.adacut_pos=-1;
	a.contam_pos=-1;
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
int seProcess::read_gz(vector<C_fastq>& pe1){
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
				pe1.push_back(fastq1);
			}
		}else{
			return -1;
		}
	}
}
void* seProcess::sub_thread_nonssd(int index){
	/*
	se_local_fs[index]=C_filter_stat();
	se_local_raw_stat1[index]=C_fastq_file_stat();
	if(!gp.trim_fq1.empty()){
		se_local_trim_stat1[index]=C_fastq_file_stat();
	}
	if(!gp.clean_fq1.empty()){
		se_local_clean_stat1[index]=C_fastq_file_stat();
	}
	*/
	vector<C_fastq> fq1s;
	vector<C_fastq> trim_result1,clean_result1;
	int done(0);

	while(1){
		se_read_m.lock();
		if(gp.fq1_path.rfind(".gz")==gp.fq1_path.size()-3){
			if(read_gz(fq1s)==-1){		//read fastqs from raw files
				done=1;
			}
		}else{
			if(read(fq1s,nongzfp1)==-1){		
				done=1;
			}
		}
		processed_reads+=fq1s.size();
		if(processed_reads%(gp.patchSize*gp.threads_num)==0){
			of_log<<get_local_time()<<"\tprocessed reads number:"<<processed_reads<<endl;
		}
		se_read_m.unlock();
		
		SEcalOption opt2;
		opt2.se_local_fs=&se_local_fs[index];
		opt2.fq1s=&fq1s;
		opt2.trim_result1=&trim_result1;
		opt2.clean_result1=&clean_result1;
		filter_se_fqs(opt2);
		SEstatOption opt_raw;
		opt_raw.fq1s=&fq1s;
		opt_raw.stat1=&se_local_raw_stat1[index];
		stat_se_fqs(opt_raw);
		fq1s.clear();
		SEstatOption opt_trim,opt_clean;
		if(!gp.trim_fq1.empty()){
			opt_trim.fq1s=&trim_result1;
			opt_trim.stat1=&se_local_trim_stat1[index];
			stat_se_fqs(opt_trim);
		}
		if(!gp.clean_fq1.empty()){
			opt_clean.fq1s=&clean_result1;
			opt_clean.stat1=&se_local_clean_stat1[index];
			stat_se_fqs(opt_clean);
		}
		se_write_m.lock();
		if(!gp.trim_fq1.empty()){
			seWrite(trim_result1,"trim",gz_trim_out1_nonssd);
		}
		if(!gp.clean_fq1.empty()){
			seWrite(clean_result1,"clean",gz_clean_out1_nonssd);
			if(gp.is_streaming){
				C_global_variable tmp_gv;
				tmp_gv.fs=*(opt2.se_local_fs);
				tmp_gv.raw1_stat=*(opt_raw.stat1);
				tmp_gv.clean1_stat=*(opt_clean.stat1);
				seStreaming_stat(tmp_gv);
			}
		}
		se_write_m.unlock();
		trim_result1.clear();
		clean_result1.clear();
		if(done==1){
			break;
		}
		if(gp.is_streaming){
			/*
			se_local_fs[index]=C_filter_stat();
			se_local_raw_stat1[index]=C_fastq_file_stat();
			se_local_clean_stat1[index]=C_fastq_file_stat();
			*/
		}
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}
void* seProcess::sub_thread(SEthreadOpt opt){
	
	int index=opt.index;
	if(!gp.trim_fq1.empty()){
		ostringstream trim_outfile1;
		trim_outfile1<<gp.output_dir<<"/thread."<<index<<".trim.r1.fq.gz";
		gz_trim_out1[index]=gzopen(trim_outfile1.str().c_str(),"wb");
		gzsetparams(gz_trim_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1[index],1024*1024*160);
	}
	if(!gp.clean_fq1.empty()){
		ostringstream outfile1;
		outfile1<<gp.output_dir<<"/thread."<<index<<".clean.r1.fq.gz";
		gz_clean_out1[index]=gzopen(outfile1.str().c_str(),"wb");
		gzsetparams(gz_clean_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out1[index],1024*1024*160);
	}
	/*
	se_local_fs[index]=C_filter_stat();
	se_local_raw_stat1[index]=C_fastq_file_stat();
	if(!gp.trim_fq1.empty()){
		se_local_trim_stat1[index]=C_fastq_file_stat();
	}
	if(!gp.clean_fq1.empty()){
		se_local_clean_stat1[index]=C_fastq_file_stat();
	}
	*/
	for(vector<string>::iterator ix=opt.files.begin();ix!=opt.files.end();ix++){
		ifstream fp1((*ix).c_str());
		int done=0;
		int iter=0;
		while(1){
			vector<C_fastq> fq1s;
			if(read(fq1s,fp1)==-1){
				done=1;
			}
			if(fq1s.size()==0)
				break;
			vector<C_fastq> trim_result1,clean_result1;
			
			
			SEcalOption opt2;
			
			opt2.se_local_fs=&se_local_fs[index];
			opt2.fq1s=&fq1s;
			opt2.trim_result1=&trim_result1;
			opt2.clean_result1=&clean_result1;
			filter_se_fqs(opt2);
			SEstatOption opt_raw;
			opt_raw.fq1s=&fq1s;
			opt_raw.stat1=&se_local_raw_stat1[index];
			stat_se_fqs(opt_raw);
			fq1s.clear();
			SEstatOption opt_trim,opt_clean;
			if(!gp.trim_fq1.empty()){
				opt_trim.fq1s=&trim_result1;
				opt_trim.stat1=&se_local_trim_stat1[index];
				stat_se_fqs(opt_trim);
			}
			if(!gp.clean_fq1.empty()){
				opt_clean.fq1s=&clean_result1;
				opt_clean.stat1=&se_local_clean_stat1[index];
				stat_se_fqs(opt_clean);
			}
			if(!gp.trim_fq1.empty()){
				seWrite(trim_result1,"trim",gz_trim_out1[index]);
				trim_result1.clear();
			}
			if(!gp.clean_fq1.empty()){
				seWrite(clean_result1,"clean",gz_clean_out1[index]);
				if(gp.is_streaming){
					C_global_variable tmp_gv;
					tmp_gv.fs=*(opt2.se_local_fs);
					tmp_gv.raw1_stat=*(opt_raw.stat1);
					tmp_gv.clean1_stat=*(opt_clean.stat1);
					seStreaming_stat(tmp_gv);
				}
				clean_result1.clear();
			}
			iter++;
			of_log<<get_local_time()<<"\tthread "<<index<<" : "<<"patch "<<iter<<" done"<<endl;
			if(done==1)
				break;
			if(gp.is_streaming){
				/*
				se_local_fs[index]=C_filter_stat();
				se_local_raw_stat1[index]=C_fastq_file_stat();
				se_local_clean_stat1[index]=C_fastq_file_stat();
				*/
			}
		}
		fp1.close();
		string rm_file="rm "+*ix;
		if(system(rm_file.c_str())==-1){
			cerr<<"Error:remove small fastq error"<<endl;
			exit(1);
		}
	}
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1[index]);
	}
	if(!gp.clean_fq1.empty()){
		gzclose(gz_clean_out1[index]);
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
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
	for(int i=0;i!=used_threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream cat_cmd1;
			cat_cmd1<<"cat "<<gp.output_dir<<"/thread."<<i<<".trim.r1.fq.gz >>"<<gp.trim_fq1<<";rm thread."<<i<<".trim.r1.fq.gz";
			if(system(cat_cmd1.str().c_str())==-1){
				cerr<<"Error:cat error"<<endl;
			}
		}
		if(!gp.clean_fq1.empty()){
			ostringstream cat_cmd1;
			cat_cmd1<<"cat "<<gp.output_dir<<"/thread."<<i<<".clean.r1.fq.gz >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<gp.output_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			if(system(cat_cmd1.str().c_str())==-1){
				cerr<<"Error:cat error"<<endl;
			}
		}
	}
}
void seProcess::process_nonssd(){
	of_log<<get_local_time()<<"\tAnalysis start!"<<endl;
	string mkdir_str="mkdir -p "+gp.output_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir fail"<<endl;
		exit(1);
	}
	if(!gp.trim_fq1.empty()){
		string out_file1=gp.output_dir+"/"+gp.trim_fq1;
		gz_trim_out1_nonssd=gzopen(out_file1.c_str(),"wb");
		gzsetparams(gz_trim_out1_nonssd, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1_nonssd,1024*1024*160);
	}
	if(!gp.clean_fq1.empty()){
		string out_file1=gp.output_dir+"/"+gp.clean_fq1;
		gz_clean_out1_nonssd=gzopen(out_file1.c_str(),"wb");
		gzsetparams(gz_clean_out1_nonssd, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out1_nonssd,1024*1024*160);

	}
	if(gp.fq1_path.rfind(".gz")==gp.fq1_path.size()-3){
		gzfp1=gzopen((gp.fq1_path).c_str(),"rb");
		gzsetparams(gzfp1, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gzfp1,2048*2048);
	}else{
		nongzfp1.open(gp.fq1_path.c_str());
	}
	thread t_array[gp.threads_num];
	for(int i=0;i<gp.threads_num;i++){
		t_array[i]=thread(bind(&seProcess::sub_thread_nonssd,this,i));
	}
	for(int i=0;i<gp.threads_num;i++){
		t_array[i].join();
	}
	merge_stat_nonssd();
	print_stat();
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1_nonssd);
	}
	if(!gp.clean_fq1.empty()){
		gzclose(gz_clean_out1_nonssd);
	}
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
	return;
}
void seProcess::process(){
	of_log<<get_local_time()<<"\tAnalysis start!"<<endl;
	string mkdir_str="mkdir -p "+gp.output_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir fail"<<endl;
		exit(1);
	}
	string clean_dir="rm "+gp.output_dir+"/tmp*fq*";
	if(system(clean_dir.c_str())==-1){
		cerr<<"Error:remove tmp fastq error"<<endl;
		exit(1);
	}
	int ii=0;
	run_pigz_split(1);
	vector<string> pe1_list;
	string find_cmd1="find "+gp.output_dir+" -name \"tmp"+random_num+".r1.fq.*\"";
	FILE* file_handle=popen(find_cmd1.c_str(),"r");
	char tmp_buf[200];
	while(fgets(tmp_buf,200,file_handle)!=NULL){
		string tmp_s(tmp_buf);
		tmp_s.erase(tmp_s.size()-1);
		pe1_list.push_back(tmp_s);
	}
	pclose(file_handle);
	vector<string> thread_files[gp.threads_num];
	int i=0;
	for(vector<string>::iterator ix=pe1_list.begin();ix!=pe1_list.end();ix++){
		thread_files[i].push_back(*ix);
		i++;
		if(i==gp.threads_num){
			i=0;
		}
	}
	for(int i=0;i!=gp.threads_num;i++){
		of_log<<"thread "<<i<<" is assigned "<<thread_files[i].size()<<" small files"<<endl;
	}
	thread t_array[gp.threads_num];
	for(int i=0;i<gp.threads_num;i++){
		SEthreadOpt curOpt;
		if(thread_files[i].size()==0){
			break;
		}
		curOpt.files=thread_files[i];
		curOpt.index=i;
		t_array[i]=thread(bind(&seProcess::sub_thread,this,curOpt));
		used_threads_num++;
	}
	if(used_threads_num>gp.threads_num){
		used_threads_num=gp.threads_num;
	}
	for(int i=0;i<used_threads_num;i++){
		t_array[i].join();
	}
	merge_stat();
	print_stat();
	merge_data();
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
	return;
	
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
void seProcess::output_split_fastqs(string type,vector<C_fastq> &fq1,gzFile outfile){
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
		if(exceed_output_se==0){
			gp.have_output1+=fq1.size();
			if(gp.have_output1>gp.output_clean){
				string out_fq1=gp.output_dir+"/rest."+gp.clean_fq1;
				gz_fq_se=gzopen(out_fq1.c_str(),"wb");
				gzsetparams(gz_fq_se, 2, Z_DEFAULT_STRATEGY);
				gzbuffer(gz_fq_se,1024*1024*160);
				exceed_output_se=1;
				int to_output=fq1.size()-(gp.have_output1-gp.output_clean);
				string sticky_tail,sticky_head;
				for(int i=0;i!=to_output;i++){
					sticky_tail+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				}
				for(int i=to_output;i!=fq1.size();i++){
					sticky_head+=fq1[i].seq_id+"\n"+fq1[i].sequence+"\n+\n"+fq1[i].qual_seq+"\n";
				}
				if(!sticky_tail.empty()){
					gzwrite(outfile,sticky_tail.c_str(),sticky_tail.size());
					gzclose(outfile);
				}
				if(!sticky_head.empty()){
					gzwrite(gz_fq_se,sticky_head.c_str(),sticky_head.size());
				}
			}else{
				gzwrite(outfile,out_content.c_str(),out_content.size());
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