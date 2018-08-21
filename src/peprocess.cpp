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
#include "peprocess.h"
#include "process_argv.h"
#include "zlib.h"
#include "gc.h"
#include "sequence.h"
using namespace::std;
#define READBUF 500
mutex pe_cal_m;
mutex read_m,stat_m,write_m;
C_filter_stat local_fs[max_thread];
C_fastq_file_stat local_raw_stat1[max_thread],local_raw_stat2[max_thread],local_trim_stat1[max_thread],local_trim_stat2[max_thread],local_clean_stat1[max_thread],local_clean_stat2[max_thread];
peProcess::peProcess(C_global_parameter m_gp){ //initialize
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
	int total_filter_fq1_num=gv.fs.output_reads_num+gv.fs.short_len_num+gv.fs.include_contam_seq_num+gv.fs.include_adapter_seq_num+gv.fs.n_ratio_num+gv.fs.polyA_num+gv.fs.tile_num+gv.fs.fov_num+gv.fs.low_qual_base_ratio_num+gv.fs.mean_quality_num+gv.fs.short_len_num+gv.fs.over_lapped_num;
	of_filter_stat<<setiosflags(ios::fixed);
	of_filter_stat<<"Total filtered read pair number\t"<<total_filter_fq1_num<<"\t100.00%\t\t"<<total_filter_fq1_num<<"\t"<<total_filter_fq1_num<<"\t"<<total_filter_fq1_num<<endl;
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
	of_filter_stat<<"Reads with PolyA\t\t"<<gv.fs.polyA_num<<"\t";
	if(total_filter_fq1_num==0){
		of_filter_stat<<"0%"<<"\t\t";
	}else{
		of_filter_stat<<setprecision(2)<<100*(float)gv.fs.polyA_num/total_filter_fq1_num<<"%"<<"\t\t";
	}
	of_filter_stat<<gv.fs.polyA_num1<<"\t"<<gv.fs.polyA_num2<<"\t"<<gv.fs.polyA_num_overlap<<endl;
	of_filter_stat.close();

	of_general_stat<<"Item\traw reads(fq1)\tclean reads(fq1)\traw reads(fq2)\tclean reads(fq2)"<<endl;
	float raw1_rl(0),raw2_rl(0),clean1_rl(0),clean2_rl(0);
	if(gv.raw1_stat.gs.reads_number!=0)
		raw1_rl=1.0*gv.raw1_stat.gs.base_number/gv.raw1_stat.gs.reads_number;
	if(gv.clean1_stat.gs.reads_number!=0)
		clean1_rl=1.0*gv.clean1_stat.gs.base_number/gv.clean1_stat.gs.reads_number;
	if(gv.raw2_stat.gs.reads_number!=0)
		raw2_rl=1.0*gv.raw2_stat.gs.base_number/gv.raw2_stat.gs.reads_number;
	if(gv.clean2_stat.gs.reads_number!=0)
		clean2_rl=1.0*gv.clean2_stat.gs.base_number/gv.clean2_stat.gs.reads_number;
	of_general_stat<<setiosflags(ios::fixed)<<setprecision(1)<<"Read length\t"<<raw1_rl<<"\t"<<clean1_rl<<"\t"<<raw2_rl<<"\t"<<clean2_rl<<endl;
	of_general_stat<<"Total number of reads\t"<<setprecision(15)<<gv.raw1_stat.gs.reads_number<<"\t"<<gv.clean1_stat.gs.reads_number<<"\t"<<gv.raw2_stat.gs.reads_number<<"\t"<<gv.clean2_stat.gs.reads_number<<endl;
	of_general_stat<<"Number of filtered reads\t"<<total_filter_fq1_num<<"\t-\t"<<total_filter_fq1_num<<"\t-"<<endl;
	of_general_stat<<"Total number of bases\t"<<setprecision(15)<<gv.raw1_stat.gs.base_number<<"\t"<<gv.clean1_stat.gs.base_number<<"\t"<<gv.raw2_stat.gs.base_number<<"\t"<<gv.clean2_stat.gs.base_number<<endl;
	of_general_stat<<"Number of filtered bases\t"<<setprecision(15)<<(float)total_filter_fq1_num*(float)gv.raw1_stat.gs.read_length<<"\t-\t"<<(float)total_filter_fq1_num*(float)gv.raw1_stat.gs.read_length<<"\t-"<<endl;
	of_general_stat<<"Number of base A\t"<<setprecision(15)<<gv.raw1_stat.gs.a_number<<"\t"<<gv.clean1_stat.gs.a_number<<"\t"<<gv.raw2_stat.gs.a_number<<"\t"<<gv.clean2_stat.gs.a_number<<endl;
	of_general_stat<<"Number of base C\t"<<setprecision(15)<<gv.raw1_stat.gs.c_number<<"\t"<<gv.clean1_stat.gs.c_number<<"\t"<<gv.raw2_stat.gs.c_number<<"\t"<<gv.clean2_stat.gs.c_number<<endl;
	of_general_stat<<"Number of base G\t"<<setprecision(15)<<gv.raw1_stat.gs.g_number<<"\t"<<gv.clean1_stat.gs.g_number<<"\t"<<gv.raw2_stat.gs.g_number<<"\t"<<gv.clean2_stat.gs.g_number<<endl;
	of_general_stat<<"Number of base T\t"<<setprecision(15)<<gv.raw1_stat.gs.t_number<<"\t"<<gv.clean1_stat.gs.t_number<<"\t"<<gv.raw2_stat.gs.t_number<<"\t"<<gv.clean2_stat.gs.t_number<<endl;
	of_general_stat<<"Number of base N\t"<<setprecision(15)<<gv.raw1_stat.gs.n_number<<"\t"<<gv.clean1_stat.gs.n_number<<"\t"<<gv.raw2_stat.gs.n_number<<"\t"<<gv.clean2_stat.gs.n_number<<endl;
	of_general_stat<<"Q20 number\t"<<setprecision(15)<<gv.raw1_stat.gs.q20_num<<"\t"<<gv.clean1_stat.gs.q20_num<<"\t"<<gv.raw2_stat.gs.q20_num<<"\t"<<gv.clean2_stat.gs.q20_num<<endl;
	of_general_stat<<"Q20 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q20_num/gv.raw1_stat.gs.base_number<<"\t"<<(float)gv.clean1_stat.gs.q20_num/gv.clean1_stat.gs.base_number<<"\t"<<(float)gv.raw2_stat.gs.q20_num/gv.raw2_stat.gs.base_number<<"\t"<<(float)gv.clean2_stat.gs.q20_num/gv.clean2_stat.gs.base_number<<endl;
	of_general_stat<<"Q30 number\t"<<setprecision(15)<<gv.raw1_stat.gs.q30_num<<"\t"<<gv.clean1_stat.gs.q30_num<<"\t"<<gv.raw2_stat.gs.q30_num<<"\t"<<gv.clean2_stat.gs.q30_num<<endl;
	of_general_stat<<"Q30 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q30_num/gv.raw1_stat.gs.base_number<<"\t"<<(float)gv.clean1_stat.gs.q30_num/gv.clean1_stat.gs.base_number<<"\t"<<(float)gv.raw2_stat.gs.q30_num/gv.raw2_stat.gs.base_number<<"\t"<<(float)gv.clean2_stat.gs.q30_num/gv.clean2_stat.gs.base_number<<endl;
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
		for(int i=0;i!=gv.clean1_stat.gs.read_max_length;i++){
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
		
		gv.fs.output_reads_num+=fs_stat.output_reads_num;	//filter stat
		gv.fs.in_adapter_list_num+=fs_stat.in_adapter_list_num;
		gv.fs.include_adapter_seq_num+=fs_stat.include_adapter_seq_num;
		gv.fs.include_adapter_seq_num1+=fs_stat.include_adapter_seq_num1;
		gv.fs.include_adapter_seq_num2+=fs_stat.include_adapter_seq_num2;
		gv.fs.include_adapter_seq_num_overlap+=fs_stat.include_adapter_seq_num_overlap;

		gv.fs.n_ratio_num+=fs_stat.n_ratio_num;
		gv.fs.n_ratio_num1+=fs_stat.n_ratio_num1;
		gv.fs.n_ratio_num2+=fs_stat.n_ratio_num2;
		gv.fs.n_ratio_num_overlap+=fs_stat.n_ratio_num_overlap;

		gv.fs.polyA_num+=fs_stat.polyA_num;
		gv.fs.polyA_num1+=fs_stat.polyA_num1;
		gv.fs.polyA_num2+=fs_stat.polyA_num2;
		gv.fs.polyA_num_overlap+=fs_stat.polyA_num_overlap;

		gv.fs.polyX_num+=fs_stat.polyX_num;
		gv.fs.polyX_num1+=fs_stat.polyX_num1;
		gv.fs.polyX_num2+=fs_stat.polyX_num2;
		gv.fs.polyX_num_overlap+=fs_stat.polyX_num_overlap;

		gv.fs.tile_num+=fs_stat.tile_num;

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

		gv.fs.over_lapped_num+=fs_stat.over_lapped_num;
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
void* peProcess::filter_pe_fqs(PEcalOption opt){
	vector<C_fastq>::iterator i2=opt.fq2s->begin();
	for(vector<C_fastq>::iterator i=opt.fq1s->begin();i!=opt.fq1s->end();i++){
		C_pe_fastq_filter pe_fastq_filter=C_pe_fastq_filter(*i,*i2,gp);
		/*if((*i).adacut_pos>0){
			cout<<"here\t"<<(*i).adacut_pos<<endl;
		}
		*/
		pe_fastq_filter.pe_trim(gp);
		if(!gp.trim_fq1.empty()){
			preOutput(1,*i);
			preOutput(2,*i2);
			opt.trim_result1->push_back(pe_fastq_filter.fq1);
			opt.trim_result2->push_back(pe_fastq_filter.fq2);
		}	
		if(pe_fastq_filter.pe_discard(opt.local_fs,gp)!=1){
			if(!gp.clean_fq1.empty()){
				preOutput(1,*i);
				preOutput(2,*i2);
				opt.clean_result1->push_back(pe_fastq_filter.fq1);
				opt.clean_result2->push_back(pe_fastq_filter.fq2);
			}
		}
		i2++;
	}
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
void peProcess::peWrite(vector<C_fastq>& pe1,vector<C_fastq>& pe2,string type,gzFile out1,gzFile out2){	//output the sequences to  files
	if(gp.threads_num==1){
		output_fastqs("1",pe1,out1);
		output_fastqs("2",pe2,out2);
	}else{
		thread write1(bind(&peProcess::output_fastqs,this,"1",pe1,out1));
		thread write2(bind(&peProcess::output_fastqs,this,"2",pe2,out2));
		write1.join();
		write2.join();
	}
}
void peProcess::C_fastq_init(C_fastq& a,C_fastq& b){
	a.adapter_seq=gp.adapter1_seq;
	b.adapter_seq=gp.adapter2_seq;
	a.contam_seq=gp.contam1_seq;
	b.contam_seq=gp.contam2_seq;
	a.head_hdcut=-1;
	a.head_lqcut=-1;
	a.tail_hdcut=-1;
	a.tail_lqcut=-1;
	a.adacut_pos=-1;
	a.contam_pos=-1;
	b.head_hdcut=-1;
	b.head_lqcut=-1;
	b.tail_hdcut=-1;
	b.tail_lqcut=-1;
	b.adacut_pos=-1;
	b.contam_pos=-1;
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
void* peProcess::sub_thread_nonssd(int index){	//sub thread process in non-ssd mode 
	of_log<<get_local_time()<<"\tthread "<<index<<" start"<<endl;
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
			local_fs[index]=C_filter_stat();
			local_raw_stat1[index]=C_fastq_file_stat();
			local_raw_stat2[index]=C_fastq_file_stat();
			local_clean_stat1[index]=C_fastq_file_stat();
			local_clean_stat2[index]=C_fastq_file_stat();
		}
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
}
void* peProcess::sub_thread(PEthreadOpt opt){	//sub thread in ssd mode
	
	int index=opt.index;
	if(opt.files.size()%2!=0){
		cerr<<"Error:pe files error"<<endl;
		exit(1);
	}
	if(!gp.trim_fq1.empty()){	//create output trim files handle
		ostringstream trim_outfile1,trim_outfile2;
		trim_outfile1<<gp.output_dir<<"/thread."<<index<<".trim.r1.fq.gz";
		trim_outfile2<<gp.output_dir<<"/thread."<<index<<".trim.r2.fq.gz";
		gz_trim_out1[index]=gzopen(trim_outfile1.str().c_str(),"wb");
		gzsetparams(gz_trim_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1[index],1024*1024*160);
		gz_trim_out2[index]=gzopen(trim_outfile2.str().c_str(),"wb");
		gzsetparams(gz_trim_out2[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out2[index],1024*1024*160);
	}
	if(!gp.clean_fq1.empty()){	//create output clean files handle
		ostringstream outfile1,outfile2;
		outfile1<<gp.output_dir<<"/thread."<<index<<".clean.r1.fq.gz";
		outfile2<<gp.output_dir<<"/thread."<<index<<".clean.r2.fq.gz";
		gz_clean_out1[index]=gzopen(outfile1.str().c_str(),"wb");
		gzsetparams(gz_clean_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out1[index],1024*1024*160);
		gz_clean_out2[index]=gzopen(outfile2.str().c_str(),"wb");
		gzsetparams(gz_clean_out2[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out2[index],1024*1024*160);
	}
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
	for(vector<string>::iterator ix=opt.files.begin();ix!=opt.files.end();ix+=2){	//process files assigned to the thread
		ifstream fp1((*ix).c_str());	//pe1
		ifstream fp2((*(ix+1)).c_str());	//pe2
		int done=0;
		int iter=0;
		while(1){
			vector<C_fastq> fq1s,fq2s;
			if(read(fq1s,fq2s,fp1,fp2)==-1){
				done=1;
			}
			if(fq1s.size()==0)
				break;
			vector<C_fastq> trim_result1,trim_result2,clean_result1,clean_result2;
			
			
			
			
			PEcalOption opt2;
			
			opt2.local_fs=&local_fs[index];
			opt2.fq1s=&fq1s;
			opt2.fq2s=&fq2s;
			opt2.trim_result1=&trim_result1;
			opt2.trim_result2=&trim_result2;
			opt2.clean_result1=&clean_result1;
			opt2.clean_result2=&clean_result2;
			filter_pe_fqs(opt2);

			PEstatOption opt_raw;
			opt_raw.fq1s=&fq1s;
			opt_raw.stat1=&local_raw_stat1[index];
			opt_raw.fq2s=&fq2s;
			opt_raw.stat2=&local_raw_stat2[index];
			stat_pe_fqs(opt_raw);
			fq1s.clear();
			fq2s.clear();
			PEstatOption opt_trim,opt_clean;
			if(!gp.trim_fq1.empty()){
				opt_trim.fq1s=&trim_result1;
				opt_trim.stat1=&local_trim_stat1[index];
				opt_trim.fq2s=&trim_result2;
				opt_trim.stat2=&local_trim_stat2[index];
				stat_pe_fqs(opt_trim);
			}
			if(!gp.clean_fq1.empty()){
				
				opt_clean.fq1s=&clean_result1;
				opt_clean.stat1=&local_clean_stat1[index];
				opt_clean.fq2s=&clean_result2;
				opt_clean.stat2=&local_clean_stat2[index];
				stat_pe_fqs(opt_clean);
			}
			if(!gp.trim_fq1.empty()){
				peWrite(trim_result1,trim_result2,"trim",gz_trim_out1[index],gz_trim_out2[index]);
				trim_result1.clear();
				trim_result2.clear();
			}
			if(!gp.clean_fq1.empty()){
				peWrite(clean_result1,clean_result2,"clean",gz_clean_out1[index],gz_clean_out2[index]);
				if(gp.is_streaming){
					C_global_variable tmp_gv;
					tmp_gv.fs=*(opt2.local_fs);
					tmp_gv.raw1_stat=*(opt_raw.stat1);
					tmp_gv.raw2_stat=*(opt_raw.stat2);
					tmp_gv.clean1_stat=*(opt_clean.stat1);
					tmp_gv.clean2_stat=*(opt_clean.stat2);
					peStreaming_stat(tmp_gv);
				}
				clean_result1.clear();
				clean_result2.clear();
			}
			iter++;
			of_log<<get_local_time()<<"\tthread "<<index<<" : "<<"patch "<<iter<<" done"<<endl;
			if(done==1)
				break;
			if(gp.is_streaming){
				local_fs[index]=C_filter_stat();
				local_raw_stat1[index]=C_fastq_file_stat();
				local_raw_stat2[index]=C_fastq_file_stat();
				local_clean_stat1[index]=C_fastq_file_stat();
				local_clean_stat2[index]=C_fastq_file_stat();
			}
		}
		fp1.close();
		fp2.close();
		string rm_file="rm "+*ix+";rm "+*(ix+1);
		if(system(rm_file.c_str())==-1){
			cerr<<"Error:remove small fastq error"<<endl;
			exit(1);
		}
	}
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1[index]);
		gzclose(gz_trim_out2[index]);
	}
	if(!gp.clean_fq1.empty()){
		gzclose(gz_clean_out1[index]);
		gzclose(gz_clean_out2[index]);
	}
	of_log<<get_local_time()<<"\tthread "<<index<<" done\t"<<endl;
	of_log.close();
}

void peProcess::run_pigz_split(int type){	//split raw files with pigz and "split" command
	ostringstream cmd1;
	int pigz_thread=gp.threads_num>16?gp.threads_num:16;
	int split_line_num=gp.split_line*4;
	if(type==1){
		if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
			cmd1<<"pigz -c -d -p "<<pigz_thread<<" "<<gp.fq1_path<<" | split -l "<<split_line_num<<" -d - "<<gp.output_dir<<"/tmp"<<random_num<<".r1.fq.";
		}else{
			cmd1<<"split -l "<<split_line_num<<" -d "<<gp.fq1_path<<" "<<gp.output_dir<<"/tmp"<<random_num<<".r1.fq.";
		}
	}else if(type==2){
		if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
			cmd1<<"pigz -c -d -p "<<pigz_thread<<" "<<gp.fq2_path<<" | split -l "<<split_line_num<<" -d - "<<gp.output_dir<<"/tmp"<<random_num<<".r2.fq.";
		}else{
			cmd1<<"split -l "<<split_line_num<<" -d "<<gp.fq2_path<<" "<<gp.output_dir<<"/tmp"<<random_num<<".r2.fq.";
		}
	}
	if(system(cmd1.str().c_str())==-1){
		cerr<<"Error:pigz&split error"<<endl;
		exit(1);
	}
}
void peProcess::merge_stat(){
	for(int i=0;i!=used_threads_num;i++){
		update_stat(local_raw_stat1[i],local_raw_stat2[i],local_fs[i],"raw");
		if(!gp.trim_fq1.empty()){
			update_stat(local_trim_stat1[i],local_trim_stat2[i],local_fs[i],"trim");
		}
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
void peProcess::merge_data(){	//cat all output files generated by multi-threads to a single large file
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
	for(int i=0;i!=used_threads_num;i++){
		if(!gp.trim_fq1.empty()){
			ostringstream cat_cmd1,cat_cmd2;
			cat_cmd1<<"cat "<<gp.output_dir<<"/thread."<<i<<".trim.r1.fq.gz >>"<<gp.trim_fq1<<";rm thread."<<i<<".trim.r1.fq.gz";
			cat_cmd2<<"cat "<<gp.output_dir<<"/thread."<<i<<".trim.r2.fq.gz >>"<<gp.trim_fq2<<";rm thread."<<i<<".trim.r2.fq.gz";
			if(system(cat_cmd1.str().c_str())==-1){
				cerr<<"Error:cat error"<<endl;
			}
			if(system(cat_cmd2.str().c_str())==-1){
				cerr<<"Error:cat error"<<endl;
			}
		}
		if(!gp.clean_fq1.empty()){
			ostringstream cat_cmd1,cat_cmd2;
			cat_cmd1<<"cat "<<gp.output_dir<<"/thread."<<i<<".clean.r1.fq.gz >>"<<gp.output_dir<<"/"<<gp.clean_fq1<<";rm "<<gp.output_dir<<"/thread."<<i<<".clean.r1.fq.gz";
			cat_cmd2<<"cat "<<gp.output_dir<<"/thread."<<i<<".clean.r2.fq.gz >>"<<gp.output_dir<<"/"<<gp.clean_fq2<<";rm "<<gp.output_dir<<"/thread."<<i<<".clean.r2.fq.gz";
			if(system(cat_cmd1.str().c_str())==-1){
				cerr<<"Error:cat error"<<endl;
			}
			if(system(cat_cmd2.str().c_str())==-1){
				cerr<<"Error:cat error"<<endl;
			}
		}
	}
}
void peProcess::process_nonssd(){
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
		string out_file2=gp.output_dir+"/"+gp.trim_fq2;
		gz_trim_out2_nonssd=gzopen(out_file2.c_str(),"wb");
		gzsetparams(gz_trim_out2_nonssd, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out2_nonssd,1024*1024*160);
	}
	if(!gp.clean_fq1.empty()){
		string out_file1=gp.output_dir+"/"+gp.clean_fq1;
		gz_clean_out1_nonssd=gzopen(out_file1.c_str(),"wb");
		gzsetparams(gz_clean_out1_nonssd, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out1_nonssd,1024*1024*160);
		string out_file2=gp.output_dir+"/"+gp.clean_fq2;

		gz_clean_out2_nonssd=gzopen(out_file2.c_str(),"wb");
		gzsetparams(gz_clean_out2_nonssd, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out2_nonssd,1024*1024*160);
	}
	if(gp.fq1_path.find(".gz")==gp.fq1_path.size()-3){
		gzfp1=gzopen((gp.fq1_path).c_str(),"rb");
		gzsetparams(gzfp1, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gzfp1,2048*2048);
		gzfp2=gzopen((gp.fq2_path).c_str(),"rb");
		gzsetparams(gzfp2, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gzfp2,2048*2048);
	}else{
		nongzfp1.open(gp.fq1_path.c_str());
		nongzfp2.open(gp.fq2_path.c_str());
	}
	thread t_array[gp.threads_num];
	for(int i=0;i<gp.threads_num;i++){
		t_array[i]=thread(bind(&peProcess::sub_thread_nonssd,this,i));
	}
	for(int i=0;i<gp.threads_num;i++){
		t_array[i].join();
	}
	merge_stat_nonssd();
	print_stat();
	if(!gp.trim_fq1.empty()){
		gzclose(gz_trim_out1_nonssd);
		gzclose(gz_trim_out2_nonssd);
	}
	if(!gp.clean_fq1.empty()){
		gzclose(gz_clean_out1_nonssd);
		gzclose(gz_clean_out2_nonssd);
	}
	of_log<<get_local_time()<<"\tAnalysis accomplished!"<<endl;
	of_log.close();
	return;
}
void peProcess::process(){

	of_log<<get_local_time()<<"\tAnalysis start!"<<endl;
	string mkdir_str="mkdir -p "+gp.output_dir;
	if(system(mkdir_str.c_str())==-1){
		cerr<<"Error:mkdir fail"<<endl;
		exit(1);
	}
	string clean_dir="rm -f "+gp.output_dir+"/tmp*fq*";
	if(system(clean_dir.c_str())==-1){
		cerr<<"Error:remove tmp fastq error"<<endl;
		exit(1);
	}
	if(gp.threads_num==1){
		run_pigz_split(1);
		run_pigz_split(2);
	}else if(gp.threads_num>1){
		thread pigz1(&peProcess::run_pigz_split,this,1);
		thread pigz2(&peProcess::run_pigz_split,this,2);
		pigz1.join();
		pigz2.join();
	}
	
	vector<string> pe1_list,pe2_list;
	string find_cmd1="find "+gp.output_dir+" -name \"tmp"+random_num+".r1.fq.*\"";
	FILE* file_handle=popen(find_cmd1.c_str(),"r");
	char tmp_buf[200];
	while(fgets(tmp_buf,200,file_handle)!=NULL){
		string tmp_s(tmp_buf);
		tmp_s.erase(tmp_s.size()-1);
		pe1_list.push_back(tmp_s);
	}
	pclose(file_handle);
	for(vector<string>::iterator ix=pe1_list.begin();ix!=pe1_list.end();ix++){
		string tmp_s2=*ix;
		tmp_s2.replace(tmp_s2.find("r1"),2,"r2");
		pe2_list.push_back(tmp_s2);
	}
	vector<string> thread_files[gp.threads_num];
	int i=0;
	vector<string>::iterator ix2=pe2_list.begin();
	for(vector<string>::iterator ix=pe1_list.begin();ix!=pe1_list.end();ix++){
		thread_files[i].push_back(*ix);
		thread_files[i].push_back(*ix2);
		ix2++;
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
		PEthreadOpt curOpt;
		if(thread_files[i].size()==0){
			break;
		}
		curOpt.files=thread_files[i];
		curOpt.index=i;
		t_array[i]=thread(bind(&peProcess::sub_thread,this,curOpt));
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
		gzflush(outfile,1);
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
	int polyA_num,polyX_num;
	int tile_num,fov_num;
	int low_qual_base_ratio_num;
	int mean_quality_num;
	int short_len_num,long_len_num;
	int over_lapped_num;
	int no_3_adapter_num,int_insertNull_num;
	*/
	int total=local_gv.fs.include_adapter_seq_num+local_gv.fs.include_contam_seq_num+local_gv.fs.low_qual_base_ratio_num+local_gv.fs.mean_quality_num+local_gv.fs.n_ratio_num+local_gv.fs.over_lapped_num+local_gv.fs.polyA_num+local_gv.fs.polyX_num;
	cout<<total<<" "<<local_gv.fs.include_adapter_seq_num<<" "<<local_gv.fs.include_contam_seq_num<<" "<<local_gv.fs.low_qual_base_ratio_num<<" "<<local_gv.fs.mean_quality_num<<" "<<local_gv.fs.n_ratio_num<<" "<<local_gv.fs.over_lapped_num<<" "<<local_gv.fs.polyA_num<<" "<<local_gv.fs.polyX_num<<"\n";
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
