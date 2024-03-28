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

using namespace ::std;
#define READBUF 1000
#define random(x) (rand() % x)

seProcess::seProcess(C_global_parameter m_gp)
{
	gp = m_gp;
	gv = C_global_variable(gp);
	used_threads_num = 0;
	srand((unsigned)time(NULL));
	ostringstream tmpstring;
	tmpstring << rand() % 100;
	random_num = tmpstring.str();
	gz_trim_out1 = new gzFile[gp.threads_num];
	gz_clean_out1 = new gzFile[gp.threads_num];
	nongz_clean_out1 = new FILE *[gp.threads_num];
	multi_gzfq1 = new gzFile[gp.threads_num];
	multi_Nongzfq1 = new FILE *[gp.threads_num];
	// gzFile* multi_gzfq1;

	se_local_fs = new C_filter_stat[gp.threads_num];
	se_local_raw_stat1 = new C_fastq_file_stat[gp.threads_num];
	se_local_trim_stat1 = new C_fastq_file_stat[gp.threads_num];
	se_local_clean_stat1 = new C_fastq_file_stat[gp.threads_num];
	for (int i = 0; i < gp.threads_num; i++)
	{
		se_local_raw_stat1[i] = C_fastq_file_stat(gp);
		se_local_trim_stat1[i] = C_fastq_file_stat(gp);
		se_local_clean_stat1[i] = C_fastq_file_stat(gp);
	}
	se_bq_check = 0;
	cur_cat_cycle = 0;
	nongz_trim_out1 = new FILE *[gp.threads_num];
	readyTrimFiles1 = new vector<string>[gp.threads_num];
	readyCleanFiles1 = new vector<string>[gp.threads_num];
	clean_file_readsNum = new vector<int>[gp.threads_num];
	nongz_trim_out1 = new FILE *[gp.threads_num];
	sub_thread_done = new int[gp.threads_num];
	for (int i = 0; i < gp.threads_num; i++)
	{
		sub_thread_done[i] = 0;
	}
	end_sub_thread = 0;
	patch = 160 / gp.threads_num;
	threadCurReadReadsNumIdx = new uint64_t[gp.threads_num];
	memset(
		threadCurReadReadsNumIdx, 0, sizeof(uint64_t) * gp.threads_num);
	if (gp.rmdup)
	{

		for (int i = 0; i < gp.threads_num; i++)
		{
			vector<uint64_t *> tmp;
			threadData.push_back(tmp);
			vector<size_t> tmp2;
			threadDataNum.push_back(tmp2);
		}
		threadReadsNum = new uint64_t[gp.threads_num];
		memset(
			threadReadsNum, 0, sizeof(uint64_t) * gp.threads_num);

		dupNum = 0;
		dupThreadOut1 = new gzFile[gp.threads_num];
		mkDir(gp.output_dir);
		for (int i = 0; i < gp.threads_num; i++)
		{
			dupThreadOut1[i] = gzopen((
										  gp.output_dir + "/dupReads." + to_string(i) + ".1.gz")
										  .c_str(),
									  "wb");
		}
	}
}

void seProcess::print_stat()
{
	string filter_out = gp.output_dir + "/Statistics_of_Filtered_Reads.txt";
	string general_out = gp.output_dir + "/Basic_Statistics_of_Sequencing_Quality.txt";
	string bs1_out = gp.output_dir + "/Base_distributions_by_read_position_1.txt";
	string qs1_out = gp.output_dir + "/Base_quality_value_distribution_by_read_position_1.txt";
	string q20_out1 = gp.output_dir + "/Distribution_of_Q20_Q30_bases_by_read_position_1.txt";
	string trim_stat1 = gp.output_dir + "/Statistics_of_Trimming_Position_of_Reads_1.txt";
	ofstream of_filter_stat(filter_out.c_str());
	ofstream of_general_stat(general_out.c_str());
	ofstream of_readPos_base_stat1(bs1_out.c_str());
	ofstream of_readPos_qual_stat1(qs1_out.c_str());
	ofstream of_q2030_stat1(q20_out1.c_str());
	ofstream of_trim_stat1(trim_stat1.c_str());
	if (!of_filter_stat)
	{
		cerr << "Error:cannot open such file," << filter_out << endl;
		exit(1);
	}
	if (!of_general_stat)
	{
		cerr << "Error:cannot open such file," << general_out << endl;
		exit(1);
	}
	if (!of_readPos_base_stat1)
	{
		cerr << "Error:cannot open such file,Base_distributions_by_read_position*.txt" << endl;
		exit(1);
	}
	if (!of_readPos_qual_stat1)
	{
		cerr << "Error:cannot open such file,Base_quality_value_distribution_by_read_position*.txt" << endl;
		exit(1);
	}
	if (!of_q2030_stat1)
	{
		cerr << "Error:cannot open such file,Distribution_of_Q20_Q30_bases_by_read_position*.txt" << endl;
		exit(1);
	}
	of_filter_stat << "Item\tTotal\tPercentage" << endl;
	vector<string> filter_items;
	filter_items.emplace_back("Reads are duplicate");
	filter_items.emplace_back("Reads limited to output number");
	filter_items.emplace_back("Reads with filtered tile");
	filter_items.emplace_back("Reads with filtered fov");
	filter_items.emplace_back("Reads too short");
	filter_items.emplace_back("Reads too long");
	filter_items.emplace_back("Reads with contam sequence");
	filter_items.emplace_back("Reads with n rate exceed");
	filter_items.emplace_back("Reads with highA");
	filter_items.emplace_back("Reads with polyX");
	filter_items.emplace_back("Reads with low quality");
	filter_items.emplace_back("Reads with low mean quality");
	filter_items.emplace_back("Reads with adapter");
	filter_items.emplace_back("Reads with global contam sequence");
	map<string, uint64_t> filter_number;
	filter_number["Reads are duplicate"] = gv.fs.dupReadsNum;
	filter_number["Reads with contam sequence"] = gv.fs.include_contam_seq_num;
	filter_number["Reads with global contam sequence"] = gv.fs.include_global_contam_seq_num;
	filter_number["Reads too short"] = gv.fs.short_len_num;
	filter_number["Reads with adapter"] = gv.fs.include_adapter_seq_num;
	filter_number["Reads with low quality"] = gv.fs.low_qual_base_ratio_num;
	filter_number["Reads with low mean quality"] = gv.fs.mean_quality_num;
	filter_number["Reads with n rate exceed"] = gv.fs.n_ratio_num;
	filter_number["Reads with highA"] = gv.fs.highA_num;
	filter_number["Reads with polyX"] = gv.fs.polyX_num;
	filter_number["Reads with filtered tile"] = gv.fs.tile_num;
	filter_number["Reads with filtered fov"] = gv.fs.fov_num;
	filter_number["Reads too long"] = gv.fs.long_len_num;
	// filter_number["Reads limited to output number"]=gv.fs.output_reads_num;
	uint64_t total_filter_fq1_num = 0;
	for (map<string, uint64_t>::iterator ix = filter_number.begin(); ix != filter_number.end(); ix++)
	{
		total_filter_fq1_num += ix->second;
	}
	// int total_filter_fq1_num=gv.fs.output_reads_num+gv.fs.include_contam_seq_num+gv.fs.include_adapter_seq_num+gv.fs.n_ratio_num+gv.fs.highA_num+gv.fs.tile_num+gv.fs.low_qual_base_ratio_num+gv.fs.mean_quality_num+gv.fs.short_len_num;
	of_filter_stat << setiosflags(ios::fixed);
	of_filter_stat << "Total filtered read pair number\t" << total_filter_fq1_num << "\t100.00%" << endl;
	for (vector<string>::iterator ix = filter_items.begin(); ix != filter_items.end(); ix++)
	{
		if (filter_number[*ix] > 0)
		{
			of_filter_stat << *ix << "\t" << filter_number[*ix] << "\t";
			of_filter_stat << setprecision(2) << 100 * (float)filter_number[*ix] / total_filter_fq1_num << "%" << endl;
		}
	}
	of_general_stat << "Item\traw reads(fq1)\tclean reads(fq1)" << endl;

	float raw1_rl(0), clean1_rl(0);
	char filter_r1_ratio[100];
	char raw_r1[7][100];
	char clean_r1[7][100];
	if (gv.raw1_stat.gs.reads_number != 0)
	{
		raw1_rl = (float)gv.raw1_stat.gs.base_number / gv.raw1_stat.gs.reads_number;
		sprintf(filter_r1_ratio, "%.2f", 100 * (float)total_filter_fq1_num / gv.raw1_stat.gs.reads_number);
		sprintf(raw_r1[0], "%.2f", 100 * (float)gv.raw1_stat.gs.a_number / gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[1], "%.2f", 100 * (float)gv.raw1_stat.gs.c_number / gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[2], "%.2f", 100 * (float)gv.raw1_stat.gs.g_number / gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[3], "%.2f", 100 * (float)gv.raw1_stat.gs.t_number / gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[4], "%.2f", 100 * (float)gv.raw1_stat.gs.n_number / gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[5], "%.2f", 100 * (float)gv.raw1_stat.gs.q20_num / gv.raw1_stat.gs.base_number);
		sprintf(raw_r1[6], "%.2f", 100 * (float)gv.raw1_stat.gs.q30_num / gv.raw1_stat.gs.base_number);
	}
	if (gv.clean1_stat.gs.reads_number != 0)
	{
		clean1_rl = (float)gv.clean1_stat.gs.base_number / gv.clean1_stat.gs.reads_number;
		sprintf(clean_r1[0], "%.2f", 100 * (float)gv.clean1_stat.gs.a_number / gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[1], "%.2f", 100 * (float)gv.clean1_stat.gs.c_number / gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[2], "%.2f", 100 * (float)gv.clean1_stat.gs.g_number / gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[3], "%.2f", 100 * (float)gv.clean1_stat.gs.t_number / gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[4], "%.2f", 100 * (float)gv.clean1_stat.gs.n_number / gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[5], "%.2f", 100 * (float)gv.clean1_stat.gs.q20_num / gv.clean1_stat.gs.base_number);
		sprintf(clean_r1[6], "%.2f", 100 * (float)gv.clean1_stat.gs.q30_num / gv.clean1_stat.gs.base_number);
	}
	of_general_stat << setiosflags(ios::fixed) << setprecision(1) << "Read length\t" << raw1_rl << "\t" << clean1_rl << endl;
	of_general_stat << "Total number of reads\t" << setprecision(15) << gv.raw1_stat.gs.reads_number << " (100.00%)\t"
					<< gv.clean1_stat.gs.reads_number << " (100.00%)" << endl;
	of_general_stat << "Number of filtered reads\t" << total_filter_fq1_num << " (" << filter_r1_ratio << "%)\t-" << endl;
	uint64_t filter_base1 = total_filter_fq1_num * gv.raw1_stat.gs.read_length;
	of_general_stat << "Total number of bases\t" << setprecision(15) << gv.raw1_stat.gs.base_number << " (100.00%)\t"
					<< gv.clean1_stat.gs.base_number << " (100.00%)" << endl;
	of_general_stat << "Number of filtered bases\t" << setprecision(15) << filter_base1 << " (" << filter_r1_ratio << "%)\t-" << endl;
	of_general_stat << "Number of base A\t" << setprecision(15) << gv.raw1_stat.gs.a_number << " (" << raw_r1[0] << "%)\t"
					<< gv.clean1_stat.gs.a_number << " (" << clean_r1[0] << "%)\t" << endl;
	of_general_stat << "Number of base C\t" << setprecision(15) << gv.raw1_stat.gs.c_number << " (" << raw_r1[1] << "%)\t"
					<< gv.clean1_stat.gs.c_number << " (" << clean_r1[1] << "%)\t" << endl;
	of_general_stat << "Number of base G\t" << setprecision(15) << gv.raw1_stat.gs.g_number << " (" << raw_r1[2] << "%)\t"
					<< gv.clean1_stat.gs.g_number << " (" << clean_r1[2] << "%)\t" << endl;
	of_general_stat << "Number of base T\t" << setprecision(15) << gv.raw1_stat.gs.t_number << " (" << raw_r1[3] << "%)\t"
					<< gv.clean1_stat.gs.t_number << " (" << clean_r1[3] << "%)\t" << endl;
	of_general_stat << "Number of base N\t" << setprecision(15) << gv.raw1_stat.gs.n_number << " (" << raw_r1[4] << "%)\t"
					<< gv.clean1_stat.gs.n_number << " (" << clean_r1[4] << "%)\t" << endl;
	of_general_stat << "Q20 number\t" << setprecision(15) << gv.raw1_stat.gs.q20_num << " (" << raw_r1[5] << "%)\t"
					<< gv.clean1_stat.gs.q20_num << " (" << clean_r1[5] << "%)" << endl;
	// of_general_stat<<"Q20 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q20_num/gv.raw1_stat.gs.base_number<<"\t"<<(float)gv.clean1_stat.gs.q20_num/gv.clean1_stat.gs.base_number<<"\t"<<endl;
	of_general_stat << "Q30 number\t" << setprecision(15) << gv.raw1_stat.gs.q30_num << " (" << raw_r1[6] << "%)\t"
					<< gv.clean1_stat.gs.q30_num << " (" << clean_r1[6] << "%)" << endl;
	// of_general_stat<<"Q30 ratio\t"<<setprecision(4)<<(float)gv.raw1_stat.gs.q30_num/gv.raw1_stat.gs.base_number<<"\t"<<(float)gv.clean1_stat.gs.q30_num/gv.clean1_stat.gs.base_number<<"\t"<<endl;
	of_general_stat.close();

	of_readPos_base_stat1 << "Pos\tA\tC\tG\tT\tN\tclean A\tclean C\tclean G\tclean T\tclean N" << endl;
	for (int i = 0; i < gv.raw1_stat.gs.read_length; i++)
	{
		of_readPos_base_stat1 << i + 1 << "\t";
		string base_set = "ACGTN";
		float raw1_cur_pos_total_base = 0;
		float clean1_cur_pos_total_base = 0;
		for (int j = 0; j != base_set.size(); j++)
		{
			raw1_cur_pos_total_base += gv.raw1_stat.bs.position_acgt_content[i][j];
			clean1_cur_pos_total_base += gv.clean1_stat.bs.position_acgt_content[i][j];
		}
		for (int j = 0; j != base_set.size(); j++)
		{
			of_readPos_base_stat1 << setiosflags(ios::fixed) << setprecision(2)
								  << 100 * (float)gv.raw1_stat.bs.position_acgt_content[i][j] / raw1_cur_pos_total_base << "%\t";
		}
		for (int j = 0; j != base_set.size(); j++)
		{
			of_readPos_base_stat1 << setiosflags(ios::fixed) << setprecision(2)
								  << 100 * (float)gv.clean1_stat.bs.position_acgt_content[i][j] / clean1_cur_pos_total_base << "%";
			if (base_set[j] != 'N')
			{
				of_readPos_base_stat1 << "\t";
			}
			else
			{
				of_readPos_base_stat1 << endl;
			}
		}
	}
	of_readPos_base_stat1.close();

	of_readPos_qual_stat1 << "#raw fastq1 quality distribution" << endl;
	of_readPos_qual_stat1 << "Pos\t";
	int max_qual = 0;
	for (int i = 0; i < gv.raw1_stat.gs.read_length; i++)
	{
		for (int j = 1; j <= gp.maxBaseQuality; j++)
		{
			if (gv.raw1_stat.qs.position_qual[i][j] > 0)
			{
				max_qual = max_qual > j ? max_qual : j;
			}
		}
	}
	for (int i = 0; i <= max_qual; i++)
	{
		of_readPos_qual_stat1 << "Q" << i << "\t";
	}
	of_readPos_qual_stat1 << "Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile" << endl;

	float *raw1_q20 = new float[gv.raw1_stat.gs.read_max_length];
	float *raw1_q30 = new float[gv.raw1_stat.gs.read_max_length];
	float *clean1_q20 = new float[gv.raw1_stat.gs.read_max_length];
	float *clean1_q30 = new float[gv.raw1_stat.gs.read_max_length];
	for (int i = 0; i != gv.raw1_stat.gs.read_length; i++)
	{
		of_readPos_qual_stat1 << i + 1 << "\t";
		uint64_t raw1_q20_num(0), raw1_q30_num(0), raw1_total(0);
		for (int j = 0; j <= max_qual; j++)
		{
			if (j >= 20)
			{
				raw1_q20_num += gv.raw1_stat.qs.position_qual[i][j];
			}
			if (j >= 30)
			{
				raw1_q30_num += gv.raw1_stat.qs.position_qual[i][j];
			}
			raw1_total += gv.raw1_stat.qs.position_qual[i][j];
			of_readPos_qual_stat1 << setiosflags(ios::fixed);
			of_readPos_qual_stat1 << setprecision(0) << gv.raw1_stat.qs.position_qual[i][j] << "\t";
		}
		raw1_q20[i] = (float)raw1_q20_num / raw1_total;
		raw1_q30[i] = (float)raw1_q30_num / raw1_total;
		quartile_result raw1_quar = cal_quar_from_array(gv.raw1_stat.qs.position_qual[i], max_qual + 1);
		of_readPos_qual_stat1 << setiosflags(ios::fixed) << setprecision(2) << raw1_quar.mean << "\t";
		of_readPos_qual_stat1 << setprecision(0) << raw1_quar.median << "\t" << raw1_quar.lower_quar << "\t" << raw1_quar.upper_quar
							  << "\t" << raw1_quar.first10_quar << "\t" << raw1_quar.last10_quar << endl;
	}
	of_readPos_qual_stat1 << "#clean fastq1 quality distribution" << endl;
	of_readPos_qual_stat1 << "Pos\t";
	for (int i = 0; i <= max_qual; i++)
	{
		of_readPos_qual_stat1 << "Q" << i << "\t";
	}
	of_readPos_qual_stat1 << "Mean\tMedian\tLower quartile\tUpper quartile\t10th percentile\t90th percentile" << endl;

	of_q2030_stat1
		<< "Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\tPercentage of Clean Q20+\tPercentage of Clean Q30+"
		<< endl;
	for (int i = 0; i != gv.clean1_stat.gs.read_max_length; i++)
	{
		of_readPos_qual_stat1 << i + 1 << "\t";
		uint64_t clean1_q20_num(0), clean1_q30_num(0), clean1_total(0);
		for (int j = 0; j <= max_qual; j++)
		{
			if (j >= 20)
			{
				clean1_q20_num += gv.clean1_stat.qs.position_qual[i][j];
			}
			if (j >= 30)
			{
				clean1_q30_num += gv.clean1_stat.qs.position_qual[i][j];
			}
			of_readPos_qual_stat1 << setiosflags(ios::fixed);
			clean1_total += gv.clean1_stat.qs.position_qual[i][j];
			of_readPos_qual_stat1 << setprecision(0) << gv.clean1_stat.qs.position_qual[i][j] << "\t";
		}
		clean1_q20[i] = (float)clean1_q20_num / clean1_total;
		clean1_q30[i] = (float)clean1_q30_num / clean1_total;
		quartile_result clean1_quar = cal_quar_from_array(gv.clean1_stat.qs.position_qual[i], max_qual + 1);
		of_readPos_qual_stat1 << setiosflags(ios::fixed) << setprecision(2) << clean1_quar.mean << "\t";
		of_readPos_qual_stat1 << setprecision(0) << clean1_quar.median << "\t" << clean1_quar.lower_quar << "\t"
							  << clean1_quar.upper_quar << "\t" << clean1_quar.first10_quar << "\t" << clean1_quar.last10_quar << endl;
		of_q2030_stat1 << i + 1 << setiosflags(ios::fixed) << setprecision(4) << "\t" << raw1_q20[i] << "\t" << raw1_q30[i] << "\t"
					   << clean1_q20[i] << "\t" << clean1_q30[i] << endl;
	}
	delete[] raw1_q20;
	delete[] raw1_q30;
	delete[] clean1_q20;
	delete[] clean1_q30;
	of_readPos_qual_stat1.close();
	of_q2030_stat1.close();
	of_trim_stat1
		<< "Pos\tHeadLowQual\tHeadFixLen\tTailAdapter\tTailLowQual\tTailFixLen\tCleanHeadLowQual\tCleanHeadFixLen\tCleanTailAdapter\tCleanTailLowQual\tCleanTailFixLen"
		<< endl;
	long long head_total1(0), tail_total1(0);
	long long head_total_clean1(0), tail_total_clean1(0);
	for (int i = 0; i < gv.raw1_stat.gs.read_length; i++)
	{
		head_total1 += gv.raw1_stat.ts.ht[i] + gv.raw1_stat.ts.hlq[i];
		tail_total1 += gv.raw1_stat.ts.ta[i] + gv.raw1_stat.ts.tlq[i] + gv.raw1_stat.ts.tt[i];
		head_total_clean1 += gv.clean1_stat.ts.ht[i] + gv.clean1_stat.ts.hlq[i];
		tail_total_clean1 += gv.clean1_stat.ts.ta[i] + gv.clean1_stat.ts.tlq[i] + gv.clean1_stat.ts.tt[i];
		// of_trim_stat1<<i+1<<"\t"<<gv.trim1_stat.ts.hlq[i]<<"\t"<<gv.trim1_stat.ts.ht[i]<<"\t"<<gv.trim1_stat.ts.ta[i]
	}
	for (int i = 1; i <= gv.raw1_stat.gs.read_length; i++)
	{
		of_trim_stat1 << i << "\t";
		if (head_total1 > 0)
		{
			of_trim_stat1 << gv.raw1_stat.ts.hlq[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.raw1_stat.ts.hlq[i] / head_total1 << "%\t";
			of_trim_stat1 << gv.raw1_stat.ts.ht[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.raw1_stat.ts.ht[i] / head_total1 << "%\t";
		}
		else
		{
			of_trim_stat1 << gv.raw1_stat.ts.hlq[i] << "\t0.00%\t";
			of_trim_stat1 << gv.raw1_stat.ts.ht[i] << "\t0.00%\t";
		}
		if (tail_total1 > 0)
		{
			of_trim_stat1 << gv.raw1_stat.ts.ta[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.raw1_stat.ts.ta[i] / tail_total1 << "%\t";
			of_trim_stat1 << gv.raw1_stat.ts.tlq[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.raw1_stat.ts.tlq[i] / tail_total1 << "%\t";
			of_trim_stat1 << gv.raw1_stat.ts.tt[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.raw1_stat.ts.tt[i] / tail_total1 << "%\t";
		}
		else
		{
			of_trim_stat1 << gv.raw1_stat.ts.ta[i] << "\t0.00%\t";
			of_trim_stat1 << gv.raw1_stat.ts.tlq[i] << "\t0.00%\t";
			of_trim_stat1 << gv.raw1_stat.ts.tlq[i] << "\t0.00%\t";
		}
		if (head_total_clean1 > 0)
		{
			of_trim_stat1 << gv.clean1_stat.ts.hlq[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.clean1_stat.ts.hlq[i] / head_total_clean1 << "%\t";
			of_trim_stat1 << gv.clean1_stat.ts.ht[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.clean1_stat.ts.ht[i] / head_total_clean1 << "%\t";
		}
		else
		{
			of_trim_stat1 << gv.clean1_stat.ts.hlq[i] << "\t0.00%\t";
			of_trim_stat1 << gv.clean1_stat.ts.ht[i] << "\t0.00%\t";
		}
		if (tail_total_clean1 > 0)
		{
			of_trim_stat1 << gv.clean1_stat.ts.ta[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.clean1_stat.ts.ta[i] / tail_total_clean1 << "%\t";
			of_trim_stat1 << gv.clean1_stat.ts.tlq[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.clean1_stat.ts.tlq[i] / tail_total_clean1 << "%\t";
			of_trim_stat1 << gv.clean1_stat.ts.tt[i] << "\t" << setiosflags(ios::fixed) << setprecision(2)
						  << 100 * (float)gv.clean1_stat.ts.tt[i] / tail_total_clean1 << "%" << endl;
		}
		else
		{
			of_trim_stat1 << gv.clean1_stat.ts.ta[i] << "\t0.00%\t";
			of_trim_stat1 << gv.clean1_stat.ts.tlq[i] << "\t0.00%\t";
			of_trim_stat1 << gv.clean1_stat.ts.tlq[i] << "\t0.00%" << endl;
		}
	}
	of_trim_stat1.close();
}

void seProcess::update_stat(C_fastq_file_stat &fq1s_stat, C_filter_stat &fs_stat, string type)
{
	if (type == "raw")
	{
		if (gv.raw1_stat.gs.read_length == 0)
			gv.raw1_stat.gs.read_length = fq1s_stat.gs.read_length; // generate stat
		if (gv.raw1_stat.gs.read_max_length < fq1s_stat.gs.read_length)
		{
			gv.raw1_stat.gs.read_max_length = fq1s_stat.gs.read_length;
		}
		gv.raw1_stat.gs.reads_number += fq1s_stat.gs.reads_number;
		gv.raw1_stat.gs.base_number += fq1s_stat.gs.base_number;
		gv.raw1_stat.gs.a_number += fq1s_stat.gs.a_number;
		gv.raw1_stat.gs.c_number += fq1s_stat.gs.c_number;
		gv.raw1_stat.gs.g_number += fq1s_stat.gs.g_number;
		gv.raw1_stat.gs.t_number += fq1s_stat.gs.t_number;
		gv.raw1_stat.gs.n_number += fq1s_stat.gs.n_number;
		gv.raw1_stat.gs.q20_num += fq1s_stat.gs.q20_num;
		gv.raw1_stat.gs.q30_num += fq1s_stat.gs.q30_num;
		string base_set = "ACGTN"; // base content and quality along read position stat
		int max_qual = 0;
		for (int i = 0; i != gv.raw1_stat.gs.read_max_length; i++)
		{
			for (int j = 0; j != base_set.size(); j++)
			{
				gv.raw1_stat.bs.position_acgt_content[i][j] += fq1s_stat.bs.position_acgt_content[i][j];
			}
		}
		for (int i = 1; i <= gv.raw1_stat.gs.read_max_length; i++)
		{
			gv.raw1_stat.ts.ht[i] += fq1s_stat.ts.ht[i];
			gv.raw1_stat.ts.hlq[i] += fq1s_stat.ts.hlq[i];
			gv.raw1_stat.ts.tt[i] += fq1s_stat.ts.tt[i];
			gv.raw1_stat.ts.tlq[i] += fq1s_stat.ts.tlq[i];
			gv.raw1_stat.ts.ta[i] += fq1s_stat.ts.ta[i];
		}
		for (int i = 0; i != gv.raw1_stat.gs.read_max_length; i++)
		{
			for (int j = 1; j <= gp.maxBaseQuality; j++)
			{
				if (fq1s_stat.qs.position_qual[i][j] > 0)
					max_qual = max_qual > j ? max_qual : j;
			}
		}

		for (int i = 0; i != gv.raw1_stat.gs.read_max_length; i++)
		{
			for (int j = 0; j <= max_qual; j++)
			{
				gv.raw1_stat.qs.position_qual[i][j] += fq1s_stat.qs.position_qual[i][j];
			}
		}

		// gv.fs.output_reads_num+=fs_stat.output_reads_num;	//filter stat
		gv.fs.in_adapter_list_num += fs_stat.in_adapter_list_num;
		gv.fs.include_adapter_seq_num += fs_stat.include_adapter_seq_num;
		gv.fs.n_ratio_num += fs_stat.n_ratio_num;
		gv.fs.highA_num += fs_stat.highA_num;
		gv.fs.polyX_num += fs_stat.polyX_num;
		gv.fs.tile_num += fs_stat.tile_num;
		gv.fs.low_qual_base_ratio_num += fs_stat.low_qual_base_ratio_num;
		gv.fs.mean_quality_num += fs_stat.mean_quality_num;
		gv.fs.short_len_num += fs_stat.short_len_num;
		gv.fs.long_len_num += fs_stat.long_len_num;
		gv.fs.include_contam_seq_num += fs_stat.include_contam_seq_num;
		gv.fs.include_global_contam_seq_num += fs_stat.include_global_contam_seq_num;
		if (gp.rmdup)
		{
			gv.fs.dupReadsNum += fs_stat.dupReadsNum;
		}
	}
	else if (type == "trim")
	{
		// gv.trim1_stat.gs.read_length=fq1s_stat.gs.read_length;	//generate stat
		gv.trim1_stat.gs.reads_number += fq1s_stat.gs.reads_number;
		gv.trim1_stat.gs.base_number += fq1s_stat.gs.base_number;
		if (gv.trim1_stat.gs.reads_number == 0)
		{
			gv.trim1_stat.gs.read_length = fq1s_stat.gs.read_length;
		}
		else
		{
			gv.trim1_stat.gs.read_length = gv.trim1_stat.gs.base_number / gv.trim1_stat.gs.reads_number;
		}
		if (gv.trim1_stat.gs.read_max_length < fq1s_stat.gs.read_length)
		{
			gv.trim1_stat.gs.read_max_length = fq1s_stat.gs.read_length;
		}
		gv.trim1_stat.gs.a_number += fq1s_stat.gs.a_number;
		gv.trim1_stat.gs.c_number += fq1s_stat.gs.c_number;
		gv.trim1_stat.gs.g_number += fq1s_stat.gs.g_number;
		gv.trim1_stat.gs.t_number += fq1s_stat.gs.t_number;
		gv.trim1_stat.gs.n_number += fq1s_stat.gs.n_number;
		gv.trim1_stat.gs.q20_num += fq1s_stat.gs.q20_num;
		gv.trim1_stat.gs.q30_num += fq1s_stat.gs.q30_num;

		int max_qual(0);
		string base_set = "ACGTN"; // base content and quality along read position stat
		for (int i = 0; i != gv.trim1_stat.gs.read_max_length; i++)
		{
			for (int j = 0; j != base_set.size(); j++)
			{
				gv.trim1_stat.bs.position_acgt_content[i][j] += fq1s_stat.bs.position_acgt_content[i][j];
			}
		}
		for (int i = 0; i != gv.trim1_stat.gs.read_max_length; i++)
		{
			gv.trim1_stat.ts.ht[i] += fq1s_stat.ts.ht[i];
			gv.trim1_stat.ts.hlq[i] += fq1s_stat.ts.hlq[i];
			gv.trim1_stat.ts.tt[i] += fq1s_stat.ts.tt[i];
			gv.trim1_stat.ts.tlq[i] += fq1s_stat.ts.tlq[i];
			gv.trim1_stat.ts.ta[i] += fq1s_stat.ts.ta[i];
		}
		for (int i = 0; i != gv.trim1_stat.gs.read_max_length; i++)
		{
			for (int j = 1; j <= gp.maxBaseQuality; j++)
			{
				if (fq1s_stat.qs.position_qual[i][j] > 0)
					max_qual = max_qual > j ? max_qual : j;
			}
		}
		for (int i = 0; i != gv.trim1_stat.gs.read_max_length; i++)
		{
			for (int j = 0; j <= max_qual; j++)
			{
				gv.trim1_stat.qs.position_qual[i][j] += fq1s_stat.qs.position_qual[i][j];
			}
		}
	}
	else if (type == "clean")
	{
		// gp.output_reads_num+=gv.clean1_stat.gs.reads_number;
		// if(gv.clean1_stat.gs.read_length==0)
		//	gv.clean1_stat.gs.read_length=fq1s_stat.gs.read_length;	//generate stat
		gv.clean1_stat.gs.reads_number += fq1s_stat.gs.reads_number;
		gv.clean1_stat.gs.base_number += fq1s_stat.gs.base_number;
		if (gv.clean1_stat.gs.reads_number == 0)
		{
			gv.clean1_stat.gs.read_length = fq1s_stat.gs.read_length;
		}
		else
		{
			gv.clean1_stat.gs.read_length = gv.clean1_stat.gs.base_number / gv.clean1_stat.gs.reads_number;
		}
		if (gv.clean1_stat.gs.read_max_length < fq1s_stat.gs.read_length)
		{
			gv.clean1_stat.gs.read_max_length = fq1s_stat.gs.read_length;
		}
		gv.clean1_stat.gs.a_number += fq1s_stat.gs.a_number;
		gv.clean1_stat.gs.c_number += fq1s_stat.gs.c_number;
		gv.clean1_stat.gs.g_number += fq1s_stat.gs.g_number;
		gv.clean1_stat.gs.t_number += fq1s_stat.gs.t_number;
		gv.clean1_stat.gs.n_number += fq1s_stat.gs.n_number;
		gv.clean1_stat.gs.q20_num += fq1s_stat.gs.q20_num;
		gv.clean1_stat.gs.q30_num += fq1s_stat.gs.q30_num;

		string base_set = "ACGTN"; // base content and quality along read position stat
		int max_qual(0);
		for (int i = 0; i != gv.clean1_stat.gs.read_max_length; i++)
		{
			for (int j = 0; j != base_set.size(); j++)
			{
				gv.clean1_stat.bs.position_acgt_content[i][j] += fq1s_stat.bs.position_acgt_content[i][j];
			}
		}
		for (int i = 0; i != gv.clean1_stat.gs.read_max_length; i++)
		{
			gv.clean1_stat.ts.ht[i] += fq1s_stat.ts.ht[i];
			gv.clean1_stat.ts.hlq[i] += fq1s_stat.ts.hlq[i];
			gv.clean1_stat.ts.tt[i] += fq1s_stat.ts.tt[i];
			gv.clean1_stat.ts.tlq[i] += fq1s_stat.ts.tlq[i];
			gv.clean1_stat.ts.ta[i] += fq1s_stat.ts.ta[i];
		}
		for (int i = 0; i != gv.clean1_stat.gs.read_max_length; i++)
		{
			for (int j = 1; j <= gp.maxBaseQuality; j++)
			{
				if (fq1s_stat.qs.position_qual[i][j] > 0)
					max_qual = max_qual > j ? max_qual : j;
			}
		}
		for (int i = 0; i != gv.clean1_stat.gs.read_max_length; i++)
		{
			for (int j = 0; j <= max_qual; j++)
			{
				gv.clean1_stat.qs.position_qual[i][j] += fq1s_stat.qs.position_qual[i][j];
			}
		}
	}
	else
	{
		cerr << "Error:code error" << endl;
		exit(1);
	}
}

void *seProcess::stat_se_fqs(SEstatOption opt, string dataType)
{
	opt.stat1->gs.reads_number += opt.fq1s->size();
	// int normal_bq(0);
	int qualityBase = 0;
	if (dataType == "clean")
	{
		qualityBase = gp.outputQualityPhred;
	}
	else
	{
		qualityBase = gp.qualityPhred;
	}
	for (vector<C_fastq>::iterator ix = opt.fq1s->begin(); ix != opt.fq1s->end(); ix++)
	{
		if ((*ix).head_hdcut > 0 || (*ix).head_lqcut > 0)
		{
			if ((*ix).head_hdcut >= (*ix).head_lqcut)
			{
				opt.stat1->ts.ht[(*ix).head_hdcut]++;
			}
			else
			{
				opt.stat1->ts.hlq[(*ix).head_lqcut]++;
			}
		}
		if ((*ix).tail_hdcut > 0 || (*ix).tail_lqcut > 0 || (*ix).adacut_pos >= 0)
		{
			if ((*ix).tail_hdcut >= (*ix).tail_lqcut)
			{
				if ((*ix).tail_hdcut >= (*ix).adacut_pos)
				{
					opt.stat1->ts.tt[(*ix).raw_length - (*ix).tail_hdcut + 1]++;
				}
				else
				{
					opt.stat1->ts.ta[(*ix).raw_length - (*ix).adacut_pos + 1]++;
				}
			}
			else
			{
				if ((*ix).tail_lqcut >= (*ix).adacut_pos)
				{
					opt.stat1->ts.tlq[(*ix).raw_length - (*ix).tail_lqcut + 1]++;
				}
				else
				{
					opt.stat1->ts.ta[(*ix).raw_length - (*ix).adacut_pos + 1]++;
				}
			}
		}
		for (string::size_type i = 0; i != (*ix).sequence.size(); i++)
		{
			switch (((*ix).sequence)[i])
			{
			case 'a':
			case 'A':
				opt.stat1->bs.position_acgt_content[i][0]++;
				opt.stat1->gs.a_number++;
				break;
			case 'c':
			case 'C':
				opt.stat1->bs.position_acgt_content[i][1]++;
				opt.stat1->gs.c_number++;
				break;
			case 'g':
			case 'G':
				opt.stat1->bs.position_acgt_content[i][2]++;
				opt.stat1->gs.g_number++;
				break;
			case 't':
			case 'T':
				opt.stat1->bs.position_acgt_content[i][3]++;
				opt.stat1->gs.t_number++;
				break;
			case 'n':
			case 'N':
				opt.stat1->bs.position_acgt_content[i][4]++;
				opt.stat1->gs.n_number++;
				break;
			default:
			{
				cerr << "Error:unrecognized sequence," << ((*ix).sequence) << endl;
				exit(1);
			}
			}
		}
		for (string::size_type i = 0; i != (*ix).qual_seq.size(); i++)
		{ // process quality sequence
			int base_quality = ((*ix).qual_seq)[i] - qualityBase;
			/*
			if(base_quality>gp.maxBaseQuality){
				cerr<<"Error:quality is too high,please check the quality system parameter or fastq file"<<endl;
				exit(1);
			}
			if(base_quality<MIN_QUAL){
				cerr<<"Error:quality is too low,please check the quality system parameter or fastq file"<<endl;
				exit(1);
			}
			*/
			opt.stat1->qs.position_qual[i][base_quality]++;
			if (base_quality >= 20)
				opt.stat1->gs.q20_num++;
			if (base_quality >= 30)
				opt.stat1->gs.q30_num++;
		}
		opt.stat1->gs.read_length = (*ix).sequence.size();
		opt.stat1->gs.base_number += opt.stat1->gs.read_length;
	}
	int q1_exceed(0), q1_normal_bq(0), q1_mean_bq(0);
	int q2_exceed(0), q2_normal_bq(0), q2_mean_bq(0);
	float q1_normal_ratio(0.0), q2_normal_ratio(0.0);
	float q1_mean(0.0), q2_mean(0.0);
	if (se_bq_check == 0)
	{
		se_stat_m.lock();
		for (vector<C_fastq>::iterator ix = opt.fq1s->begin(); ix != opt.fq1s->end(); ix++)
		{
			int qual_len = (*ix).qual_seq.size();
			for (string::size_type i = 0; i != qual_len; i++)
			{ // process quality sequence
				int base_quality = ((*ix).qual_seq)[i] - gp.qualityPhred;
				q1_mean_bq += base_quality;
				if (base_quality >= MIN_QUAL && base_quality <= gp.maxBaseQuality)
				{
					q1_normal_bq++;
				}
				else
				{
					if (base_quality < MIN_QUAL - 10 || base_quality > gp.maxBaseQuality + 10)
						q1_exceed++;
				}
				int another_q = gp.qualityPhred == 64 ? 33 : 64;
				int base_quality2 = ((*ix).qual_seq)[i] - another_q;
				q2_mean_bq += base_quality2;
				if (base_quality2 >= MIN_QUAL && base_quality2 <= gp.maxBaseQuality)
				{
					q2_normal_bq++;
				}
				else
				{
					if (base_quality2 < MIN_QUAL - 10 || base_quality2 > gp.maxBaseQuality + 10)
						q2_exceed++;
				}
			}
		}

		if (opt.stat1->gs.base_number == 0)
		{
			cerr << "Error:no data" << endl;
			exit(1);
		}
		q1_normal_ratio = (float)q1_normal_bq / opt.stat1->gs.base_number;
		q2_normal_ratio = (float)q2_normal_bq / opt.stat1->gs.base_number;
		q1_mean = (float)q1_mean_bq / opt.stat1->gs.base_number;
		q2_mean = (float)q2_mean_bq / opt.stat1->gs.base_number;
		// cout<<q1_exceed<<"\t"<<q1_normal_bq<<"\t"<<opt.stat1->gs.base_number<<"\t"<<q1_normal_ratio<<"\t"<<q1_mean<<endl;
		// cout<<q2_exceed<<"\t"<<q2_normal_bq<<"\t"<<opt.stat1->gs.base_number<<"\t"<<q2_normal_ratio<<"\t"<<q2_mean<<endl;
		int q1_score(0), q2_score(0);
		if (q1_exceed)
		{
			q1_score += 0;
		}
		else
		{
			q1_score += 1;
		}
		if (q2_exceed)
		{
			q2_score += 0;
		}
		else
		{
			q2_score += 1;
		}
		if (q1_normal_ratio > q2_normal_ratio)
		{
			q1_score += 3;
			q2_score += 0;
		}
		else if (q1_normal_ratio < q2_normal_ratio)
		{
			q1_score += 0;
			q2_score += 3;
		}
		else
		{
			q1_score += 3;
			q2_score += 3;
		}
		if (q1_mean < 10 || q1_mean > gp.maxBaseQuality)
		{
			q1_score += 0;
		}
		else
		{
			q1_score += 2;
		}
		if (q2_mean < 10 || q2_mean > gp.maxBaseQuality)
		{
			q2_score += 0;
		}
		else
		{
			q2_score += 2;
		}
		// cout<<q1_score<<"\t"<<q2_score<<endl;
		if (se_bq_check == 0)
		{
			if (q1_score - q2_score < -3)
			{
				cerr << "Error:base quality seems abnormal,please check the quality system parameter or fastq file" << endl;
				exit(1);
			}
			else if (q1_score - q2_score < 0)
			{
				cerr << "Warning:base quality seems abnormal,please check the quality system parameter or fastq file"
					 << endl;
			}
		}
		se_bq_check++;
		/*
		if(normal_ratio>1){
			cerr<<"Error:code error"<<endl;
			exit(1);
		}if(normal_ratio<0.5){
			cerr<<"Error:base quality seems abnormal,please check the quality system parameter or fastq file"<<endl;
			exit(1);
		}else if(normal_ratio<0.9){
			cerr<<"Warning:base quality seems abnormal,please check the quality system parameter or fastq file"<<endl;
		}else{
		}
		*/
		se_stat_m.unlock();
		// exit(1);
	}
	return &se_bq_check;
}

void seProcess::filter_se_fqs(SEcalOption opt)
{
	// C_reads_trim_stat cut_pos;
	// check dup
	int iter = 0;
	for (vector<C_fastq>::iterator i = opt.fq1s->begin(); i != opt.fq1s->end(); i++)
	{
		C_single_fastq_filter se_fastq_filter = C_single_fastq_filter(*i, gp);
		iter++;
		se_fastq_filter.se_trim(gp);
		if (gp.adapter_discard_or_trim == "trim" || gp.contam_discard_or_trim == "trim" || !gp.trim.empty() || !gp.trimBadHead.empty() || !gp.trimBadTail.empty())
		{
			(*i).head_hdcut = se_fastq_filter.read.head_hdcut;
			(*i).head_lqcut = se_fastq_filter.read.head_lqcut;
			(*i).tail_hdcut = se_fastq_filter.read.tail_hdcut;
			(*i).tail_lqcut = se_fastq_filter.read.tail_lqcut;
			(*i).adacut_pos = se_fastq_filter.read.adacut_pos;
			//(*i).contam_pos=se_fastq_filter.read.contam_pos;
			//(*i).global_contam_pos=se_fastq_filter.read.global_contam_pos;
			//(*i).raw_length=se_fastq_filter.read.raw_length;
		}
		//*i=se_fastq_filter.read;
		if (!gp.trim_fq1.empty())
		{
			preOutput(1, se_fastq_filter.read);
			opt.trim_result1->emplace_back(se_fastq_filter.read);
		}
		int whether_discard(0);
		if (gp.module_name == "filtersRNA")
		{
			whether_discard = se_fastq_filter.sRNA_discard(opt.se_local_fs, gp);
		}
		else
		{
			whether_discard = se_fastq_filter.se_discard(opt.se_local_fs, gp);
		}
		if (whether_discard != 1)
		{
			if (!gp.clean_fq1.empty())
			{
				preOutput(1, se_fastq_filter.read);
				opt.clean_result1->emplace_back(se_fastq_filter.read);
			}
		}
	}
	// return cut_pos;
}

void seProcess::preOutput(int type, C_fastq &a)
{
	if (!gp.base_convert.empty())
	{
		gp.base_convert = gp.base_convert.replace(gp.base_convert.find("TO"), 2, "");
		gp.base_convert = gp.base_convert.replace(gp.base_convert.find("2"), 1, "");
		if (gp.base_convert.size() != 2)
		{
			cerr << "Error:base_conver value format error" << endl;
			exit(1);
		}
		for (string::size_type ix = 0; ix != a.sequence.size(); ix++)
		{
			if (toupper(a.sequence[ix]) == toupper(gp.base_convert[0]))
				a.sequence[ix] = gp.base_convert[1];
		}
	}
}

void seProcess::seWrite(vector<C_fastq> &se1, gzFile out1)
{
	output_fastqs("1", se1, out1);
}

void seProcess::seWrite(vector<C_fastq> &se1, FILE *out1)
{
	output_fastqs("1", se1, out1);
}

void seProcess::C_fastq_init(C_fastq &a)
{
	a.seq_id = "";
	a.sequence = "";
	a.qual_seq = "";
	//	a.adapter_seq=gp.adapter1_seq;
	//	a.contam_seq=gp.contam1_seq;
	// a.global_contams=gp.global_contams;
	a.head_hdcut = -1;
	a.head_lqcut = -1;
	a.tail_hdcut = -1;
	a.tail_lqcut = -1;
	a.adacut_pos = -1;
	a.contam_pos = -1;
	a.global_contam_5pos = -1;
	a.global_contam_3pos = -1;
	a.raw_length = 0;
	a.pairPos = 1;
	//	if(gp.module_name=="filtersRNA"){
	//		a.adapter_seq2=gp.adapter2_seq;
	//	}
	if (!gp.trim.empty())
	{
		vector<string> tmp_eles = get_se_hard_trim(gp.trim);
		a.head_trim_len = tmp_eles[0];
		a.tail_trim_len = tmp_eles[1];
	}
}

int seProcess::read(vector<C_fastq> &pe1, ifstream &infile1)
{

	string buf1;
	int file1_line_num(0);
	C_fastq fastq1;
	C_fastq_init(fastq1);
	for (int i = 0; i < gp.patchSize * 4; i++)
	{
		if (getline(infile1, buf1))
		{
			file1_line_num++;
			// string s_line(buf1);
			// s_line.erase(s_line.size()-1);
			if (file1_line_num % 4 == 1)
			{
				fastq1.seq_id = buf1;
			}
			if (file1_line_num % 4 == 2)
			{
				fastq1.sequence = buf1;
			}
			if (file1_line_num % 4 == 0)
			{
				fastq1.qual_seq = buf1;
				pe1.emplace_back(fastq1);
				// fq1s.emplace_back(fastq1);
			}
		}
		else
		{
			return -1;
		}
	}
	return 0;
}

void *seProcess::sub_thread(int index)
{
	logLock.lock();
	of_log << get_local_time() << "\tthread " << index << " start" << endl;
	logLock.unlock();
	//	}
	create_thread_read(index);
	int thread_cycle = -1;
	char buf1[READBUF];
	C_fastq fastq1;
	C_fastq_init(fastq1);
	long long file1_line_num(0);
	long long block_line_num1(0);
	int thread_read_block = 4 * gp.patchSize * patch;
	vector<C_fastq> fq1s;
	gzFile tmpRead = gzopen((gp.fq1_path).c_str(), "rb");
	int spaceNum = 0;
	if (gzgets(tmpRead, buf1, READBUF) != NULL)
	{
		string tmpLine(buf1);
		while (isspace(tmpLine[tmpLine.size() - 1]))
		{
			spaceNum++;
			tmpLine.erase(tmpLine.size() - 1);
		}
	}
	gzclose(tmpRead);
	bool inputGzformat = true;
	if (gp.fq1_path.rfind(".gz") == gp.fq1_path.size() - 3)
	{
		inputGzformat = true;
	}
	else
	{
		inputGzformat = false;
	}
	if (inputGzformat)
	{
		while (1)
		{
			if (gzgets(multi_gzfq1[index], buf1, READBUF) != NULL)
			{
				if ((file1_line_num / thread_read_block) % gp.threads_num == index)
				{
					block_line_num1++;
					if (block_line_num1 % 4 == 1)
					{
						fastq1.seq_id.assign(buf1);
						fastq1.seq_id.erase(fastq1.seq_id.size() - spaceNum, spaceNum);
					}
					else if (block_line_num1 % 4 == 2)
					{
						fastq1.sequence.assign(buf1);
						fastq1.sequence.erase(fastq1.sequence.size() - spaceNum, spaceNum);
					}
					else if (block_line_num1 % 4 == 0)
					{
						fastq1.qual_seq.assign(buf1);
						fastq1.qual_seq.erase(fastq1.qual_seq.size() - spaceNum, spaceNum);
						fq1s.emplace_back(fastq1);
						if (fq1s.size() == gp.patchSize)
						{
							if (end_sub_thread == 1)
							{
								break;
							}
							int tmp_cycle = file1_line_num / (thread_read_block * gp.threads_num);
							if (tmp_cycle != thread_cycle && tmp_cycle > 0)
							{
								addCleanList(thread_cycle, index);
							}
							thread_cycle = tmp_cycle;
							threadCurReadReadsNumIdx[index] = file1_line_num / 4;
							thread_process_reads(
								index, thread_cycle, fq1s);
							if (index == 0)
							{
								of_log << get_local_time() << " processed_reads:\t" << file1_line_num / 4 << endl;
							}
						}
					}
				}
				file1_line_num++;
			}
			else
			{
				if (!fq1s.empty())
				{
					if (end_sub_thread == 1)
					{
						break;
					}
					int tmp_cycle = file1_line_num / (thread_read_block * gp.threads_num);
					if (tmp_cycle != thread_cycle && tmp_cycle > 0)
					{
						addCleanList(thread_cycle, index);
					}
					thread_cycle = tmp_cycle;
					threadCurReadReadsNumIdx[index] = file1_line_num / 4;
					thread_process_reads(
						index, thread_cycle, fq1s);
				}
				gzclose(multi_gzfq1[index]);
				break;
			}
		}
		if (thread_cycle >= 0)
			addCleanList(thread_cycle, index);
	}
	else
	{
		while (1)
		{
			if (fgets(buf1, READBUF, multi_Nongzfq1[index]) != NULL)
			{
				if ((file1_line_num / thread_read_block) % gp.threads_num == index)
				{
					block_line_num1++;
					if (block_line_num1 % 4 == 1)
					{
						fastq1.seq_id.assign(buf1);
						fastq1.seq_id.erase(fastq1.seq_id.size() - spaceNum, spaceNum);
					}
					else if (block_line_num1 % 4 == 2)
					{
						fastq1.sequence.assign(buf1);
						fastq1.sequence.erase(fastq1.sequence.size() - spaceNum, spaceNum);
					}
					else if (block_line_num1 % 4 == 0)
					{
						fastq1.qual_seq.assign(buf1);
						fastq1.qual_seq.erase(fastq1.qual_seq.size() - spaceNum, spaceNum);
						fq1s.emplace_back(fastq1);
						if (fq1s.size() == gp.patchSize)
						{
							if (end_sub_thread == 1)
							{
								break;
							}
							int tmp_cycle = file1_line_num / (thread_read_block * gp.threads_num);
							if (tmp_cycle != thread_cycle && tmp_cycle > 0)
							{
								addCleanList(thread_cycle, index);
							}
							thread_cycle = tmp_cycle;
							threadCurReadReadsNumIdx[index] = file1_line_num / 4;
							thread_process_reads(
								index, thread_cycle, fq1s);
							if (index == 0)
							{
								of_log << get_local_time() << " processed_reads:\t" << file1_line_num / 4 << endl;
							}
						}
					}
				}
				file1_line_num++;
			}
			else
			{
				if (!fq1s.empty())
				{
					if (end_sub_thread == 1)
					{
						break;
					}
					int tmp_cycle = file1_line_num / (thread_read_block * gp.threads_num);
					if (tmp_cycle != thread_cycle && tmp_cycle > 0)
					{
						addCleanList(thread_cycle, index);
					}
					thread_cycle = tmp_cycle;
					threadCurReadReadsNumIdx[index] = file1_line_num / 4;
					thread_process_reads(
						index, thread_cycle, fq1s);
				}
				fclose(multi_Nongzfq1[index]);
				break;
			}
		}
		if (thread_cycle >= 0)
			addCleanList(thread_cycle, index);
	}
	check_disk_available();
	sub_thread_done[index] = 1;
	logLock.lock();
	of_log << get_local_time() << "\tthread " << index << " done\t" << endl;
	logLock.unlock();
	return &se_bq_check;
}

void seProcess::addCleanList(int tmp_cycle, int index)
{
	ostringstream checkClean1;
	if (gp.cleanOutGzFormat)
	{
		checkClean1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << tmp_cycle << ".clean.r1.fq.gz";
	}
	else
	{
		checkClean1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << tmp_cycle << ".clean.r1.fq";
	}
	if (find(readyCleanFiles1[index].begin(), readyCleanFiles1[index].end(), checkClean1.str()) == readyCleanFiles1[index].end())
	{
		readyCleanFiles1[index].emplace_back(checkClean1.str());
	}
	if (!gp.trim_fq1.empty())
	{
		ostringstream readyTrimR1;
		if (gp.cleanOutGzFormat)
		{
			readyTrimR1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << tmp_cycle << ".trim.r1.fq.gz";
		}
		else
		{
			readyTrimR1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << tmp_cycle << ".trim.r1.fq";
		}
		if (find(readyTrimFiles1[index].begin(), readyTrimFiles1[index].end(), readyTrimR1.str()) == readyTrimFiles1[index].end())
			readyTrimFiles1[index].emplace_back(readyTrimR1.str());
	}
}

void seProcess::merge_stat()
{
	for (int i = 0; i != gp.threads_num; i++)
	{
		update_stat(se_local_raw_stat1[i], se_local_fs[i], "raw");
		if (!gp.trim_fq1.empty())
		{
			update_stat(se_local_trim_stat1[i], se_local_fs[i], "trim");
		}
		update_stat(se_local_clean_stat1[i], se_local_fs[i], "clean");
	}
}

void seProcess::create_thread_read(int index)
{
	if (gp.inputGzformat)
	{
		multi_gzfq1[index] = gzopen((gp.fq1_path).c_str(), "rb");
		if (!multi_gzfq1[index])
		{
			cerr << "Error:cannot open the file," << gp.fq1_path << endl;
			exit(1);
		}
		gzsetparams(multi_gzfq1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(multi_gzfq1[index], 2048 * 2048);
	}
	else
	{
		multi_Nongzfq1[index] = fopen((gp.fq1_path).c_str(), "r");
		if (!multi_Nongzfq1[index])
		{
			cerr << "Error:cannot open the file," << gp.fq1_path << endl;
			exit(1);
		}
	}
}

void seProcess::catRmFile(int index, int cycle, string type, bool gzFormat)
{
	ostringstream out_fq1_tmp, out_fq2_tmp;
	if (gzFormat)
	{
		out_fq1_tmp << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << cycle << "." << type << ".r1.fq.gz";
	}
	else
	{
		out_fq1_tmp << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << cycle << "." << type << ".r1.fq";
	}
	if (access(out_fq1_tmp.str().c_str(), 0) != -1)
	{
		ostringstream cat_cmd1;
		string outFile1;
		if (type == "trim")
		{
			outFile1 = gp.trim_fq1;
		}
		else if (type == "clean")
		{
			outFile1 = gp.clean_fq1;
		}
		else
		{
			cerr << "Error: code error!" << endl;
		}
		cat_cmd1 << "cat " << out_fq1_tmp.str() << " >>" << gp.output_dir << "/" << outFile1 << ";rm " << out_fq1_tmp.str();
		run_cmd(cat_cmd1.str());
	}
}

void seProcess::catRmFile(vector<int> indexes, int cycle, string type, bool gzFormat)
{
	ostringstream cmd1;
	cmd1 << "cat ";
	ostringstream rmCmd1;
	rmCmd1 << "rm ";
	for (vector<int>::iterator ix = indexes.begin(); ix != indexes.end(); ix++)
	{
		ostringstream out_fq1_tmp;
		if (gzFormat)
		{
			out_fq1_tmp << gp.output_dir << "/" << tmp_dir << "/thread." << *ix << "." << cycle << "." << type << ".r1.fq.gz";
		}
		else
		{
			out_fq1_tmp << gp.output_dir << "/" << tmp_dir << "/thread." << *ix << "." << cycle << "." << type << ".r1.fq";
		}
		if (access(out_fq1_tmp.str().c_str(), 0) != -1)
		{
			ostringstream cat_cmd1;
			cmd1 << out_fq1_tmp.str() << " ";
			rmCmd1 << out_fq1_tmp.str() << " ";
		}
	}
	string outFile1;
	if (type == "trim")
	{
		outFile1 = gp.trim_fq1;
	}
	else if (type == "clean")
	{
		outFile1 = gp.clean_fq1;
	}
	else
	{
		cerr << "Error: code error!" << endl;
	}
	cmd1 << " >>" << gp.output_dir << "/" << outFile1;
	if (indexes.size() > 0)
	{
		run_cmd(cmd1.str());
		run_cmd(rmCmd1.str());
	}
}

void seProcess::run_cmd(string cmd)
{
	if (system(cmd.c_str()) == -1)
	{
		cerr << "Error:when running " << cmd << endl;
		exit(1);
	}
}

void *seProcess::smallFilesProcess()
{
	if (!gp.trim_fq1.empty())
	{
		string trim_file1;
		trim_file1 = gp.output_dir + "/" + gp.trim_fq1;
		if (check_gz_empty(trim_file1) == 1)
		{
			string rm_cmd = "rm " + trim_file1;
			system(rm_cmd.c_str());
		}
	}
	if (!gp.clean_fq1.empty())
	{
		string clean_file1;
		clean_file1 = gp.output_dir + "/" + gp.clean_fq1;
		if (check_gz_empty(clean_file1) == 1)
		{
			string rm_cmd = "rm " + clean_file1;
			system(rm_cmd.c_str());
		}
	}
	int sleepTime = 5;
	if (gp.cleanOutSplit > 0)
	{
		int outputFileIndex = 0;
		int cur_avaliable_total_reads_number = 0;
		//        int sticky_tail_reads_number = 0;
		string tidyFile1 = gp.output_dir + "/split.0." + gp.clean_fq1;
		ostringstream rmCmd1;
		rmCmd1 << "rm " << gp.output_dir << "/split.*fq*";

		if (access(tidyFile1.c_str(), 0) != -1)
		{
			run_cmd(rmCmd1.str());
		}

		while (1)
		{
			bool subThreadAllDone = true;
			for (int i = 0; i < gp.threads_num; i++)
			{
				if (sub_thread_done[i] != 1)
				{
					subThreadAllDone = false;
					break;
				}
			}
			if (subThreadAllDone)
			{
				int ready_cycles = 0;
				for (int i = 0; i < gp.threads_num; i++)
				{
					if (readyCleanFiles1[i].size() > ready_cycles)
					{
						ready_cycles = readyCleanFiles1[i].size();
					}
				}
				for (int cycle = cur_cat_cycle; cycle < ready_cycles; cycle++)
				{
					for (int i = 0; i < gp.threads_num; i++)
					{
						if (cycle == ready_cycles - 1 && readyCleanFiles1[i].size() < ready_cycles)
						{
							break;
						}
						cur_avaliable_total_reads_number += clean_file_readsNum[i][cycle];
						if (cur_avaliable_total_reads_number >= gp.cleanOutSplit)
						{
							int toBeOutputReadsNumber = gp.cleanOutSplit - (cur_avaliable_total_reads_number - clean_file_readsNum[i][cycle]);
							extractReadsToFile(
								cycle, i, toBeOutputReadsNumber, "head", outputFileIndex, gp.cleanOutGzFormat);
							cur_avaliable_total_reads_number -= gp.cleanOutSplit;
						}
						else
						{
							ostringstream catFile1;
							if (gp.cleanOutGzFormat)
							{
								catFile1 << "cat " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle
										 << ".clean.r1.fq.gz >>" << gp.output_dir << "/split." << outputFileIndex << "."
										 << gp.clean_fq1 << ";rm " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle
										 << ".clean.r1.fq.gz";
							}
							else
							{
								catFile1 << "cat " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle
										 << ".clean.r1.fq >>" << gp.output_dir << "/split." << outputFileIndex << "."
										 << gp.clean_fq1 << ";rm " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle
										 << ".clean.r1.fq";
							}
							if (system(catFile1.str().c_str()) == -1)
							{
								cerr << "Error:cat file error" << endl;
								exit(1);
							}
						}
					}
				}
				break;
			}
			//            int stopIndex=0;

			int ready_cycles = readyCleanFiles1[0].size();
			for (int i = 1; i < gp.threads_num; i++)
			{
				if (readyCleanFiles1[i].size() < ready_cycles)
				{
					ready_cycles = readyCleanFiles1[i].size();
				}
			}

			for (int cycle = cur_cat_cycle; cycle < ready_cycles; cycle++)
			{
				for (int i = 0; i < gp.threads_num; i++)
				{

					if (clean_file_readsNum[i].size() < ready_cycles)
					{
						cerr << "Error:code error" << endl;
						exit(1);
					}
					cur_avaliable_total_reads_number += clean_file_readsNum[i][cycle];
					if (cur_avaliable_total_reads_number >= gp.cleanOutSplit)
					{
						// sticky_tail_reads_number=clean_file_readsNum[i][cycle]-cur_avaliable_total_reads_number;
						int toBeOutputReadsNumber =
							gp.cleanOutSplit - (cur_avaliable_total_reads_number - clean_file_readsNum[i][cycle]);
						extractReadsToFile(cycle, i, toBeOutputReadsNumber, "head", outputFileIndex, gp.cleanOutGzFormat);
						cur_avaliable_total_reads_number -= gp.cleanOutSplit;
					}
					else
					{
						ostringstream catFile1;
						if (gp.cleanOutGzFormat)
						{
							catFile1 << "cat " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle
									 << ".clean.r1.fq.gz >>" << gp.output_dir << "/split." << outputFileIndex << "." << gp.clean_fq1
									 << ";rm " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle << ".clean.r1.fq.gz";
						}
						else
						{
							catFile1 << "cat " << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle << ".clean.r1.fq >>"
									 << gp.output_dir << "/split." << outputFileIndex << "." << gp.clean_fq1 << ";rm "
									 << gp.output_dir << "/" << tmp_dir << "/thread." << i << "." << cycle << ".clean.r1.fq";
						}
						if (system(catFile1.str().c_str()) == -1)
						{
							cerr << "Error:cat file error" << endl;
							exit(1);
						}
					}
				}
			}
			if (ready_cycles > 0)
				cur_cat_cycle = ready_cycles;
			sleep(sleepTime);
		}
	}
	else
	{
		uint64_t total_merged_reads_number = 0;
		while (1)
		{ // merge small files by input order

			bool subThreadAllDone = true;
			for (int i = 0; i < gp.threads_num; i++)
			{
				if (sub_thread_done[i] != 1)
				{
					subThreadAllDone = false;
					break;
				}
			}
			if (subThreadAllDone)
			{
				int ready_cycles = 0;
				for (int i = 0; i < gp.threads_num; i++)
				{
					if (readyCleanFiles1[i].size() > ready_cycles)
					{
						ready_cycles = readyCleanFiles1[i].size();
					}
				}
				for (int cycle = cur_cat_cycle; cycle < ready_cycles; cycle++)
				{
					vector<int> readyCatFiles;
					for (int i = 0; i < gp.threads_num; i++)
					{
						if (cycle == ready_cycles - 1 && readyCleanFiles1[i].size() < ready_cycles)
						{
							break;
						}
						if (readyCleanFiles1[i].size() == 0)
						{
							break;
						}
						total_merged_reads_number += clean_file_readsNum[i][cycle];
						if (!gp.total_reads_num_random && gp.total_reads_num > 0)
						{
							if (total_merged_reads_number > gp.total_reads_num)
							{
								end_sub_thread = 1;
								catRmFile(readyCatFiles, cycle, "clean", gp.cleanOutGzFormat);
								// output some reads to the final file
								int toBeOutputReadsNumber =
									gp.total_reads_num - (total_merged_reads_number - clean_file_readsNum[i][cycle]);
								extractReadsToFile(cycle, i, toBeOutputReadsNumber, "head", gp.cleanOutGzFormat);
								rmTmpFiles();
								return &se_bq_check;
							}
						}
						readyCatFiles.emplace_back(i);
					}
					if (!gp.trim_fq1.empty())
					{
						catRmFile(readyCatFiles, cycle, "trim", gp.trimOutGzformat);
					}
					if (!gp.clean_fq1.empty())
					{
						catRmFile(readyCatFiles, cycle, "clean", gp.cleanOutGzFormat);
					}
				}
				break;
			}
			int ready_cycles = readyCleanFiles1[0].size();
			for (int i = 1; i < gp.threads_num; i++)
			{
				if (readyCleanFiles1[i].size() < ready_cycles)
				{
					ready_cycles = readyCleanFiles1[i].size();
				}
			}
			for (int cycle = cur_cat_cycle; cycle < ready_cycles; cycle++)
			{
				if (!gp.trim_fq1.empty())
				{
					for (int i = 0; i < gp.threads_num; i++)
					{
						catRmFile(i, cycle, "trim", gp.trimOutGzformat);
					}
				}
				if (!gp.clean_fq1.empty())
				{
					vector<int> readyCatFiles;
					for (int i = 0; i < gp.threads_num; i++)
					{
						if (readyCleanFiles1[i].size() == 0)
						{
							break;
						}
						total_merged_reads_number += clean_file_readsNum[i][cycle];
						if (!gp.total_reads_num_random && gp.total_reads_num > 0)
						{
							if (total_merged_reads_number > gp.total_reads_num)
							{
								end_sub_thread = 1;
								catRmFile(readyCatFiles, cycle, "clean", gp.cleanOutGzFormat);
								// output some reads to the final file
								int toBeOutputReadsNumber =
									gp.total_reads_num - (total_merged_reads_number - clean_file_readsNum[i][cycle]);
								extractReadsToFile(cycle, i, toBeOutputReadsNumber, "head", gp.cleanOutGzFormat);
								rmTmpFiles();
								return &se_bq_check;
							}
						}
						readyCatFiles.emplace_back(i);
					}
					catRmFile(readyCatFiles, cycle, "clean", gp.cleanOutGzFormat);
				}
			}
			if (ready_cycles > 0)
				cur_cat_cycle = ready_cycles;
			sleep(sleepTime);
			//        if(sleepTime<60)
			//            sleepTime+=10;
		}
	}
	return &se_bq_check;
}

void seProcess::create_thread_smalltrimoutputFile(int index, int cycle)
{
	ostringstream trim_outfile1;
	if (gp.trim_fq1.rfind(".gz") == gp.trim_fq1.size() - 3)
	{ // create output trim files handle
		trim_outfile1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << cycle << ".trim.r1.fq.gz";
		gz_trim_out1[index] = gzopen(trim_outfile1.str().c_str(), "ab");
		if (!gz_trim_out1[index])
		{
			cerr << "Error:cannot write to the file," << trim_outfile1.str() << endl;
			exit(1);
		}
		gzsetparams(gz_trim_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_trim_out1[index], 1024 * 1024 * 8);
	}
	else
	{
		trim_outfile1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << cycle << ".trim.r1.fq";
		if ((nongz_clean_out1[index] = fopen(trim_outfile1.str().c_str(), "a")) == NULL)
		{
			cerr << "Error:cannot write to the file," << trim_outfile1.str() << endl;
			exit(1);
		}
	}
}

void seProcess::create_thread_smallcleanoutputFile(int index, int cycle)
{
	ostringstream outfile1;
	if (gp.clean_fq1.rfind(".gz") == gp.clean_fq1.size() - 3)
	{
		outfile1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << cycle << ".clean.r1.fq.gz";
		gz_clean_out1[index] = gzopen(outfile1.str().c_str(), "ab");
		if (!gz_clean_out1[index])
		{
			cerr << "Error:cannot write to the file," << outfile1.str() << endl;
			exit(1);
		}
		gzsetparams(gz_clean_out1[index], 2, Z_DEFAULT_STRATEGY);
		gzbuffer(gz_clean_out1[index], 1024 * 1024 * 8);
	}
	else
	{
		outfile1 << gp.output_dir << "/" << tmp_dir << "/thread." << index << "." << cycle << ".clean.r1.fq";
		if ((nongz_clean_out1[index] = fopen(outfile1.str().c_str(), "a")) == NULL)
		{
			cerr << "Error:cannot write to the file," << outfile1.str() << endl;
			exit(1);
		}
	}
}

void seProcess::closeSmallTrimFileHandle(int index)
{
	if (gp.trim_fq1.rfind(".gz") == gp.trim_fq1.size() - 3)
	{
		gzclose(gz_trim_out1[index]);
	}
	else
	{
		fclose(nongz_trim_out1[index]);
	}
}

void seProcess::closeSmallCleanFileHandle(int index)
{
	if (gp.clean_fq1.rfind(".gz") == gp.clean_fq1.size() - 3)
	{
		gzclose(gz_clean_out1[index]);
	}
	else
	{
		fclose(nongz_clean_out1[index]);
	}
}

void *seProcess::sub_extract(string in, int mo, string out)
{
	int lineNum = 0;
	int readsNum = 0;
	if (gp.cleanOutGzFormat)
	{
		gzFile subFile = gzopen(out.c_str(), "wb");
		gzsetparams(subFile, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(subFile, 2048 * 2048);
		gzFile cleanFile = gzopen(in.c_str(), "rb");
		char buf[READBUF];
		while (gzgets(cleanFile, buf, READBUF))
		{
			if (lineNum % (4 * mo) >= 0 && lineNum % (4 * mo) <= 3)
			{
				string line(buf);
				gzwrite(subFile, line.c_str(), line.size());
				readsNum++;
				if (readsNum / 4 >= gp.l_total_reads_num && readsNum % 4 == 0)
				{
					break;
				}
			}
			lineNum++;
		}
		gzclose(subFile);
		gzclose(cleanFile);
	}
	else
	{
		FILE *subFile = fopen(out.c_str(), "w");
		FILE *cleanFile = fopen(in.c_str(), "r");
		char buf[READBUF];
		while (fgets(buf, READBUF, cleanFile))
		{
			if (lineNum % (4 * mo) >= 0 && lineNum % (4 * mo) <= 3)
			{
				string line(buf);
				fputs(line.c_str(), subFile);
				readsNum++;
				if (readsNum / 4 >= gp.l_total_reads_num && readsNum % 4 == 0)
				{
					break;
				}
			}
			lineNum++;
		}
		fclose(subFile);
		fclose(cleanFile);
	}
	return &se_bq_check;
}

void seProcess::rmTmpFiles()
{
	string cmd = "rm -f " + gp.output_dir + "/" + tmp_dir + "/*fq*";
	run_cmd(cmd);
}

void seProcess::extractReadsToFile(
	int cycle, int thread_index, int reads_number, string position, int &output_index, bool gzFormat)
{
	ostringstream outFile1;
	outFile1 << gp.output_dir << "/split." << output_index << "." << gp.clean_fq1;
	ostringstream cleanSmallFile1;

	if (gzFormat)
	{
		cleanSmallFile1 << gp.output_dir << "/" << tmp_dir << "/thread." << thread_index << "." << cycle << ".clean.r1.fq.gz";
		gzFile gzCleanSmall1;
		gzCleanSmall1 = gzopen(cleanSmallFile1.str().c_str(), "rb");

		gzFile splitGzFq1;
		splitGzFq1 = gzopen(outFile1.str().c_str(), "ab");
		gzsetparams(splitGzFq1, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(splitGzFq1, 1024 * 1024 * 10);
		char buf1[READBUF];

		if (position == "head")
		{
			for (int i = 0; i < reads_number * 4; i++)
			{
				if (gzgets(gzCleanSmall1, buf1, READBUF) != NULL)
				{
					string line(buf1);
					gzwrite(splitGzFq1, line.c_str(), line.size());
				}
			}
			gzclose(splitGzFq1);
			outFile1.str("");
			output_index++;
			outFile1 << gp.output_dir << "/split." << output_index << "." << gp.clean_fq1;
			splitGzFq1 = gzopen(outFile1.str().c_str(), "wb");
			gzsetparams(splitGzFq1, 2, Z_DEFAULT_STRATEGY);
			gzbuffer(splitGzFq1, 1024 * 1024 * 10);
			int readsNum = 0;
			while (gzgets(gzCleanSmall1, buf1, READBUF) != NULL)
			{
				readsNum++;
				string line(buf1);
				gzwrite(splitGzFq1, line.c_str(), line.size());
				if (readsNum == gp.cleanOutSplit * 4)
				{
					gzclose(splitGzFq1);
					output_index++;
					outFile1.str("");
					outFile1 << gp.output_dir << "/split." << output_index << "." << gp.clean_fq1;
					splitGzFq1 = gzopen(outFile1.str().c_str(), "wb");
					gzsetparams(splitGzFq1, 2, Z_DEFAULT_STRATEGY);
					gzbuffer(splitGzFq1, 1024 * 1024 * 10);
					readsNum = 0;
				}
			}
			gzclose(splitGzFq1);
			gzclose(gzCleanSmall1);
		}
	}
	else
	{
		cleanSmallFile1 << gp.output_dir << "/thread." << thread_index << "." << cycle << ".clean.r1.fq";
		FILE *nongzCleanSmall1;
		nongzCleanSmall1 = fopen(cleanSmallFile1.str().c_str(), "r");

		FILE *splitNonGzFq1;
		splitNonGzFq1 = fopen(outFile1.str().c_str(), "a");
		char buf1[READBUF];

		if (position == "head")
		{
			for (int i = 0; i < reads_number * 4; i++)
			{
				if (fgets(buf1, READBUF, nongzCleanSmall1) != NULL)
				{
					string line(buf1);
					fputs(line.c_str(), splitNonGzFq1);
				}
			}
			fclose(splitNonGzFq1);
			outFile1.str("");
			output_index++;
			outFile1 << gp.output_dir << "/split." << output_index << "." << gp.clean_fq1;
			splitNonGzFq1 = fopen(outFile1.str().c_str(), "w");
			int readsNum = 0;
			while (fgets(buf1, READBUF, nongzCleanSmall1) != NULL)
			{
				readsNum++;
				string line(buf1);
				fputs(line.c_str(), splitNonGzFq1);
				if (readsNum == gp.cleanOutSplit * 4)
				{
					fclose(splitNonGzFq1);
					output_index++;
					outFile1.str("");
					outFile1 << gp.output_dir << "/split." << output_index << "." << gp.clean_fq1;
					splitNonGzFq1 = fopen(outFile1.str().c_str(), "w");
					readsNum = 0;
				}
			}
			fclose(splitNonGzFq1);
			fclose(nongzCleanSmall1);
		}
	}
	string rmCmd1 = "rm " + cleanSmallFile1.str();
	run_cmd(rmCmd1);
}

void seProcess::extractReadsToFile(int cycle, int thread_index, int reads_number, string position, bool gzFormat)
{
	ostringstream cleanSmallFile1;

	if (gzFormat)
	{
		cleanSmallFile1 << gp.output_dir << "/" << tmp_dir << "/thread." << thread_index << "." << cycle << ".clean.r1.fq.gz";
		gzFile gzCleanSmall1;
		gzCleanSmall1 = gzopen(cleanSmallFile1.str().c_str(), "rb");

		gzFile splitGzFq1;
		string out1 = gp.output_dir + "/" + gp.clean_fq1;
		splitGzFq1 = gzopen(out1.c_str(), "ab");
		gzsetparams(splitGzFq1, 2, Z_DEFAULT_STRATEGY);
		gzbuffer(splitGzFq1, 1024 * 1024 * 10);
		char buf1[READBUF];

		if (position == "head")
		{
			for (int i = 0; i < reads_number * 4; i++)
			{
				if (gzgets(gzCleanSmall1, buf1, READBUF) != NULL)
				{
					string line(buf1);
					gzwrite(splitGzFq1, line.c_str(), line.size());
				}
			}
			gzclose(gzCleanSmall1);
			gzclose(splitGzFq1);
		}
	}
	else
	{
		cleanSmallFile1 << gp.output_dir << "/thread." << thread_index << "." << cycle << ".clean.r1.fq";
		FILE *nongzCleanSmall1;
		nongzCleanSmall1 = fopen(cleanSmallFile1.str().c_str(), "r");
		FILE *splitNonGzFq1;
		splitNonGzFq1 = fopen(gp.clean_fq1.c_str(), "a");
		char buf1[READBUF];

		if (position == "head")
		{
			for (int i = 0; i < reads_number * 4; i++)
			{
				if (fgets(buf1, READBUF, nongzCleanSmall1) != NULL)
				{
					string line(buf1);
					fputs(line.c_str(), splitNonGzFq1);
				}
			}
			fclose(splitNonGzFq1);
			fclose(nongzCleanSmall1);
		}
	}
	string rmCmd1 = "rm " + cleanSmallFile1.str();
	run_cmd(rmCmd1);
}

void seProcess::thread_process_reads(int index, int &cycle, vector<C_fastq> &fq1s)
{
	check_disk_available();
	create_thread_smallcleanoutputFile(index, cycle);
	if (!gp.trim_fq1.empty())
	{
		create_thread_smalltrimoutputFile(
			index, cycle);
	}
	vector<C_fastq> trim_result1, clean_result1;

	SEcalOption opt2;
	opt2.se_local_fs = &se_local_fs[index];
	opt2.fq1s = &fq1s;
	opt2.trim_result1 = &trim_result1;
	opt2.clean_result1 = &clean_result1;
	if (gp.rmdup && RMDUP == 2)
	{
		filter_se_fqs(
			opt2, index);
	}
	else
	{
		filter_se_fqs(opt2); // filter raw fastqs by the given parameters
	}
	SEstatOption opt_raw;
	opt_raw.fq1s = &fq1s;
	opt_raw.stat1 = &se_local_raw_stat1[index];
	stat_se_fqs(
		opt_raw, "raw"); // statistic raw fastqs
	fq1s.clear();
	// add_raw_trim(se_local_raw_stat1[index],raw_cut);

	SEstatOption opt_trim, opt_clean;
	if (!gp.trim_fq1.empty())
	{ // trim means only trim but not discard.
		opt_trim.fq1s = &trim_result1;
		opt_trim.stat1 = &se_local_trim_stat1[index];
		stat_se_fqs(opt_trim, "trim"); // statistic trim fastqs
	}
	/*
	if(!gp.clean_fq1.empty()){
		opt_clean.fq1s=&clean_result1;
		opt_clean.stat1=&se_local_clean_stat1[index];
		stat_se_fqs(opt_clean);	//statistic clean fastqs
	}
	*/
	//
	if (!gp.trim_fq1.empty())
	{
		seWrite(trim_result1, gz_trim_out1[index]); // output trim files
		trim_result1.clear();
		closeSmallTrimFileHandle(index);
	}
	opt_clean.fq1s = &clean_result1;
	if (!gp.clean_fq1.empty())
	{
		opt_clean.stat1 = &se_local_clean_stat1[index];
		if (gp.cleanOutGzFormat)
		{
			seWrite(clean_result1, gz_clean_out1[index]); // output clean files
		}
		else
		{
			seWrite(clean_result1, nongz_clean_out1[index]); // output clean files
		}
		closeSmallCleanFileHandle(index);

		stat_se_fqs(opt_clean, "clean"); // statistic clean fastqs
		if (gp.is_streaming)
		{
			se_write_m.lock();
			C_global_variable *tmp_gv = new C_global_variable();
			tmp_gv->fs = *(opt2.se_local_fs);
			tmp_gv->raw1_stat = *(opt_raw.stat1);
			// tmp_gv->trim1_stat=local_trim_stat1[index];
			// tmp_gv->trim2_stat=local_trim_stat2[index];
			tmp_gv->clean1_stat = *(opt_clean.stat1);
			seStreaming_stat(*tmp_gv);
			delete tmp_gv;
			se_write_m.unlock();
		}
	}

	if (clean_file_readsNum[index].size() < cycle + 1)
	{
		clean_file_readsNum[index].emplace_back(opt_clean.fq1s->size());
	}
	else
	{
		clean_file_readsNum[index][cycle] += opt_clean.fq1s->size();
	}

	clean_result1.clear();
	check_disk_available();
	//    cycle++;
}

void seProcess::run_extract_random()
{
	if (gp.total_reads_num <= 0 || gp.total_reads_num_random == false)
	{
		cerr << "Error:extract random clean reads error cuz parameters are wrong" << endl;
		exit(1);
	}
	uint64_t total_clean_reads(0);
	for (int i = 0; i != gp.threads_num; i++)
	{
		total_clean_reads += se_local_clean_stat1[i].gs.reads_number;
	}
	if (gp.f_total_reads_ratio > 0)
	{
		if (gp.f_total_reads_ratio >= 1)
		{
			cerr << "Error:the ratio extract from clean fq file should not be more than 1" << endl;
			exit(1);
		}
		if (gp.l_total_reads_num > 0)
		{
			cerr << "Error:reads number and ratio should not be both assigned at the same time" << endl;
			exit(1);
		}
		gp.l_total_reads_num = total_clean_reads * gp.f_total_reads_ratio;
	}
	if (total_clean_reads < gp.l_total_reads_num)
	{
		cerr << "Warning:the reads number in clean fastq file(" << total_clean_reads
			 << ") is less than you assigned to output(" << gp.l_total_reads_num << ")" << endl;
		return;
	}
	// cout<<gp.l_total_reads_num<<"\t"<<total_clean_reads<<endl;
	// vector<int> include_threads;
	if (gp.l_total_reads_num == 0)
	{
		cerr << "Error:assigned reads number should not be 0" << endl;
		return;
	}
	int interval = total_clean_reads / gp.l_total_reads_num;
	if (interval == 1)
		return;
	string in1 = gp.output_dir + "/" + gp.clean_fq1;
	string out1 = gp.cleanOutGzFormat ? gp.output_dir + "/cleanRandomExtractReads.r1.fq.gz" : gp.output_dir + "/cleanRandomExtractReads.r1.fq";
	sub_extract(in1, interval, out1);
	string cmd1 = "mv " + gp.output_dir + "/" + gp.clean_fq1 + " " + gp.output_dir + "/total." + gp.clean_fq1;
	cmd1 += "; mv " + out1 + " " + gp.output_dir + "/" + gp.clean_fq1;
	run_cmd(cmd1);
}

void seProcess::process()
{

	mkDir(gp.output_dir);
	//    string mkdir_str="mkdir -p "+gp.output_dir;
	//    if(system(mkdir_str.c_str())==-1){
	//        cerr<<"Error:mkdir fail"<<endl;
	//        exit(1);
	//    }
	of_log.open(
		gp.log.c_str());
	if (!of_log)
	{
		cerr << "Error:cannot open such file," << gp.log << endl;
		exit(1);
	}
	of_log << get_local_time() << "\tAnalysis start!" << endl;
	//    of_log<<"memSize used in rmdup:"<<(dupDB->realUseByteSize)/(1024*1024)<<"M"<<endl;
	make_tmpDir();
	if (gp.rmdup)
	{
		thread t_array[gp.threads_num];
		for (int i = 0; i < gp.threads_num; i++)
		{
			// t_array[i]=thread(bind(&peProcess::sub_thread_nonssd_multiOut,this,i));
			t_array[i] = thread(
				bind(
					&seProcess::sub_thread_rmdup_step1, this, i));
		}
		for (int i = 0; i < gp.threads_num; i++)
		{
			t_array[i].join();
		}
		int maxCycle = 0;
		uint64_t totalReadsNum = 0;
		for (int i = 0; i < gp.threads_num; i++)
		{
			totalReadsNum += threadReadsNum[i];
			if (threadData[i].size() > maxCycle)
			{
				maxCycle = threadData[i].size();
			}
		}
		delete[] threadReadsNum;
		if (totalReadsNum > (pow(
								 2, 32) -
							 1))
		{
			cerr << "Error:reads number is too large to do remove duplication," << totalReadsNum << endl;
			exit(1);
		}
		totalData = new uint64_t[totalReadsNum];
		memset(
			totalData, 0, sizeof(uint64_t) * totalReadsNum);
		//        int iter=0;
		uint64_t checkNum = 0;
		uint64_t *totalTmp = totalData;
		for (int i = 0; i < maxCycle; i += patch)
		{
			for (int j = 0; j < gp.threads_num; j++)
			{
				for (int k = 0; k < patch; k++)
				{
					if (threadData[j].size() > i + k)
					{
						checkNum += threadDataNum[j][i + k];
						if (checkNum > totalReadsNum)
						{
							cerr << "Error:code error," << __FILE__ << "," << __LINE__ << endl;
							exit(1);
						}
						memcpy(
							totalTmp, threadData[j][i + k], sizeof(uint64_t) * threadDataNum[j][i + k]);
						totalTmp += threadDataNum[j][i + k];
						delete[] threadData[j][i + k];
					}
				}
			}
		}
		for (int i = 0; i < gp.threads_num; i++)
		{
			vector<uint64_t *>().swap(threadData[i]);
			vector<size_t>().swap(threadDataNum[i]);
		}
		vector<vector<uint64_t *>>().swap(threadData);
		vector<vector<size_t>>().swap(threadDataNum);
		dupFlag = new bool[totalReadsNum];
		memset(
			dupFlag, 0, sizeof(bool) * totalReadsNum);
		rmdup *dormdup = new rmdup(
			totalData, totalReadsNum);
		dormdup->markDup(dupFlag);
		delete dormdup;

		for (int i = 0; i < totalReadsNum; i++)
		{
			if (dupFlag[i])
			{
				dupNum++;
			}
		}
		of_log << "duplicate reads number:\t" << dupNum << endl;
	}
	thread t_array[gp.threads_num];
	// thread read_monitor(bind(&peProcess::monitor_read_thread,this));
	// sleep(10);
	for (int i = 0; i < gp.threads_num; i++)
	{
		// t_array[i]=thread(bind(&peProcess::sub_thread_nonssd_multiOut,this,i));
		t_array[i] = thread(
			bind(
				&seProcess::sub_thread, this, i));
	}
	thread catFiles = thread(
		bind(
			&seProcess::smallFilesProcess, this));
	catFiles.join();

	for (int i = 0; i < gp.threads_num; i++)
	{
		t_array[i].join();
	}

	check_disk_available();
	if (gp.total_reads_num_random == true && gp.total_reads_num > 0)
	{
		run_extract_random();
	}
	merge_stat();
	print_stat();
	remove_tmpDir();
	check_disk_available();
	if (gp.rmdup)
	{
		for (int i = 0; i < gp.threads_num; i++)
		{
			if (dupThreadOut1[i] != NULL)
			{
				gzclose(dupThreadOut1[i]);
			}
		}
	}
	if (gp.rmdup)
	{
		of_log << "dup number:\t" << dupNum << endl;
	}
	of_log << get_local_time() << "\tAnalysis accomplished!" << endl;
	of_log.close();
}

void seProcess::remove_tmpDir()
{
	string rm_dir = gp.output_dir + "/" + tmp_dir;
	int iter(0);
	while (1)
	{
		DIR *tmpdir;
		struct dirent *ptr;
		tmpdir = opendir(rm_dir.c_str());
		bool empty_flag = true;
		if (tmpdir == NULL)
		{
			cerr << "Error:open directory error," << rm_dir << endl;
			exit(1);
		}
		while ((ptr = readdir(tmpdir)) != NULL)
		{
			if (strcmp(ptr->d_name, "..") != 0 && strcmp(ptr->d_name, ".") != 0)
			{
				empty_flag = false;
				break;
			}
		}
		closedir(tmpdir);
		if (empty_flag)
		{
			string remove = "rm -r " + rm_dir;
			if (system(remove.c_str()) == -1)
			{
				cerr << "Error:rmdir error," << remove << endl;
				exit(1);
			}
			break;
		}
		else
		{
			sleep(2);
			iter++;
		}
		if (iter > 30)
		{
			break;
		}
	}
}

void seProcess::make_tmpDir()
{
	srand(time(0));
	ostringstream tmp_str;
	for (int i = 0; i != 6; i++)
	{
		int tmp_rand = random(26) + 'A';
		tmp_str << (char)tmp_rand;
	}
	tmp_dir = "TMP" + tmp_str.str();
	string tmp_dir_abs = gp.output_dir + "/" + tmp_dir;
	mkDir(tmp_dir_abs);
	//	string mkdir_str="mkdir -p "+gp.output_dir+"/"+tmp_dir;
	//	if(system(mkdir_str.c_str())==-1){
	//		cerr<<"Error:mkdir error,"<<mkdir_str<<endl;
	//		exit(1);
	//	}
}

void seProcess::output_fastqs(string type, vector<C_fastq> &fq1, gzFile outfile)
{
	// m.lock();
	string out_content, streaming_out;
	for (int i = 0; i != fq1.size(); i++)
	{
		// for(vector<C_fastq>::iterator i=fq1->begin();i!=fq1->end();i++){
		if (gp.output_file_type == "fasta")
		{
			fq1[i].seq_id = fq1[i].seq_id.replace(fq1[i].seq_id.find("@"), 1, ">");
			out_content += fq1[i].seq_id + "\n" + fq1[i].sequence + "\n";
		}
		else if (gp.output_file_type == "fastq")
		{
			if (gp.outputQualityPhred != gp.qualityPhred)
			{
				for (string::size_type ix = 0; ix != fq1[i].qual_seq.size(); ix++)
				{
					int b_q = fq1[i].qual_seq[ix] - gp.qualityPhred;
					fq1[i].qual_seq[ix] = (char)(b_q + gp.outputQualityPhred);
				}
			}
			if (gp.is_streaming)
			{
				string modify_id = fq1[i].seq_id;
				modify_id.erase(0, 1);
				streaming_out += ">+\t" + modify_id + "\t" + type + "\t" + fq1[i].sequence + "\t" + fq1[i].qual_seq + "\n";
			}
			else
			{
				out_content += fq1[i].seq_id + "\n" + fq1[i].sequence + "\n+\n" + fq1[i].qual_seq + "\n";
			}
		}
		else
		{
			cerr << "Error:output_file_type value error" << endl;
			exit(1);
		}
	}
	if (gp.is_streaming)
	{
		cout << streaming_out;
	}
	else
	{
		gzwrite(outfile, out_content.c_str(), out_content.size());
		gzflush(outfile, 1);
	}
	// m.unlock();
}

void seProcess::output_fastqs(string type, vector<C_fastq> &fq1, FILE *outfile)
{
	// m.lock();
	string out_content, streaming_out;
	for (int i = 0; i != fq1.size(); i++)
	{
		// for(vector<C_fastq>::iterator i=fq1->begin();i!=fq1->end();i++){
		if (gp.output_file_type == "fasta")
		{
			fq1[i].seq_id = fq1[i].seq_id.replace(fq1[i].seq_id.find("@"), 1, ">");
			out_content += fq1[i].seq_id + "\n" + fq1[i].sequence + "\n";
		}
		else if (gp.output_file_type == "fastq")
		{
			if (gp.outputQualityPhred != gp.qualityPhred)
			{
				for (string::size_type ix = 0; ix != fq1[i].qual_seq.size(); ix++)
				{
					int b_q = fq1[i].qual_seq[ix] - gp.qualityPhred;
					fq1[i].qual_seq[ix] = (char)(b_q + gp.outputQualityPhred);
				}
			}
			if (gp.is_streaming)
			{
				string modify_id = fq1[i].seq_id;
				modify_id.erase(0, 1);
				streaming_out += ">+\t" + modify_id + "\t" + type + "\t" + fq1[i].sequence + "\t" + fq1[i].qual_seq + "\n";
			}
			else
			{
				out_content += fq1[i].seq_id + "\n" + fq1[i].sequence + "\n+\n" + fq1[i].qual_seq + "\n";
			}
		}
		else
		{
			cerr << "Error:output_file_type value error" << endl;
			exit(1);
		}
	}
	if (gp.is_streaming)
	{
		cout << streaming_out;
	}
	else
	{
		// gzwrite(outfile,out_content.c_str(),out_content.size());
		fputs(out_content.c_str(), outfile);
		// gzflush(outfile,1);
	}
	// m.unlock();
}

void seProcess::seStreaming_stat(C_global_variable &local_gv)
{
	cout << "#Total_statistical_information"
		 << "\n";
	int total = local_gv.fs.include_adapter_seq_num + local_gv.fs.include_contam_seq_num + local_gv.fs.low_qual_base_ratio_num + local_gv.fs.mean_quality_num + local_gv.fs.n_ratio_num + local_gv.fs.over_lapped_num + local_gv.fs.highA_num + local_gv.fs.polyX_num;
	cout << total << " " << local_gv.fs.include_adapter_seq_num << " " << local_gv.fs.include_contam_seq_num << " "
		 << local_gv.fs.low_qual_base_ratio_num << " " << local_gv.fs.mean_quality_num << " " << local_gv.fs.n_ratio_num << " "
		 << local_gv.fs.over_lapped_num << " " << local_gv.fs.highA_num << " " << local_gv.fs.polyX_num << "\n";
	cout << "#Fq1_statistical_information"
		 << "\n";
	cout << local_gv.raw1_stat.gs.read_length << " " << local_gv.clean1_stat.gs.read_length << " "
		 << local_gv.raw1_stat.gs.reads_number << " " << local_gv.clean1_stat.gs.reads_number << " "
		 << local_gv.raw1_stat.gs.base_number << " " << local_gv.clean1_stat.gs.base_number << " "
		 << local_gv.raw1_stat.gs.a_number << " " << local_gv.clean1_stat.gs.a_number << " " << local_gv.raw1_stat.gs.c_number
		 << " " << local_gv.clean1_stat.gs.c_number << " " << local_gv.raw1_stat.gs.g_number << " "
		 << local_gv.clean1_stat.gs.g_number << " " << local_gv.raw1_stat.gs.t_number << " "
		 << local_gv.clean1_stat.gs.t_number << " " << local_gv.raw1_stat.gs.n_number << " "
		 << local_gv.clean1_stat.gs.n_number << " " << local_gv.raw1_stat.gs.q20_num << " " << local_gv.clean1_stat.gs.q20_num
		 << " " << local_gv.raw1_stat.gs.q30_num << " " << local_gv.clean1_stat.gs.q30_num << "\n";
	cout << "#Base_distributions_by_read_position"
		 << "\n";
	// uint64_t position_acgt_content[READ_MAX_LEN][5];
	for (int i = 0; i != local_gv.raw1_stat.gs.read_length; i++)
	{
		for (int j = 0; j != 4; j++)
		{
			cout << local_gv.raw1_stat.bs.position_acgt_content[i][j] << " ";
		}
		cout << local_gv.raw1_stat.bs.position_acgt_content[i][4] << "\n";
	}
	for (int i = 0; i != local_gv.clean1_stat.gs.read_length; i++)
	{
		for (int j = 0; j != 4; j++)
		{
			cout << local_gv.clean1_stat.bs.position_acgt_content[i][j] << " ";
		}
		cout << local_gv.clean1_stat.bs.position_acgt_content[i][4] << "\n";
	}
	cout << "#Raw_Base_quality_value_distribution_by_read_position"
		 << "\n";
	// position_qual[READ_MAX_LEN][gp.maxBaseQuality]
	for (int i = 0; i != local_gv.raw1_stat.gs.read_length; i++)
	{
		for (int j = 0; j != 40; j++)
		{
			cout << local_gv.clean1_stat.qs.position_qual[i][j] << " ";
		}
		cout << "0\n";
	}
	for (int i = 0; i != local_gv.clean1_stat.gs.read_length; i++)
	{
		for (int j = 0; j != 40; j++)
		{
			cout << local_gv.clean1_stat.qs.position_qual[i][j] << " ";
		}
		cout << "0\n";
	}
}

void seProcess::check_disk_available()
{
	if (access(
			gp.fq1_path.c_str(), 0) == -1)
	{
		cerr << "Error:input raw fastq not exists suddenly, please check the disk" << endl;
		exit(1);
	}
	if (access(
			gp.output_dir.c_str(), 0) == -1)
	{
		cerr << "Error:output directory cannot open suddenly, please check the disk" << endl;
		exit(1);
	}
}

void *seProcess::sub_thread_rmdup_step1(int index)
{
	logLock.lock();
	of_log << get_local_time() << "\tthread " << index << " pre-rmdup start" << endl;
	logLock.unlock();
	create_thread_read(index);
	//    int thread_cycle=-1;
	char buf1[READBUF];
	C_fastq fastq1, fastq2;
	C_fastq_init(fastq1);
	long long file1_line_num(0);
	long long block_line_num1(0);
	int thread_read_block = 4 * gp.patchSize * patch;
	vector<C_fastq> fq1s;
	bool inputGzformat = true;
	gzFile tmpRead = gzopen((gp.fq1_path).c_str(), "rb");
	int spaceNum = 0;
	if (gzgets(
			tmpRead, buf1, READBUF) != NULL)
	{
		string tmpLine(buf1);
		while (isspace(
			tmpLine[tmpLine.size() - 1]))
		{
			spaceNum++;
			tmpLine.erase(
				tmpLine.size() - 1);
		}
	}
	gzclose(tmpRead);
	if (gp.fq1_path.rfind(".gz") == gp.fq1_path.size() - 3)
	{
		inputGzformat = true;
	}
	else
	{
		inputGzformat = false;
	}
	string fq1seq, fq2seq;
	vector<string> seqs;
	if (inputGzformat)
	{
		while (1)
		{
			if (gzgets(
					multi_gzfq1[index], buf1, READBUF) != NULL)
			{
				if ((
						file1_line_num / thread_read_block) %
						gp.threads_num ==
					index)
				{
					block_line_num1++;
					if (block_line_num1 % 4 == 2)
					{
						fq1seq.assign(buf1);
						fq1seq.erase(
							fq1seq.size() - spaceNum, spaceNum);
						seqs.emplace_back(fq1seq);
					}
					if (seqs.size() == gp.patchSize)
					{
						uint64_t *curData = new uint64_t[seqs.size()];
						for (int i = 0; i < seqs.size(); i++)
						{
							curData[i] = hash<string>()(seqs[i]);
							//                            MDString(seqs[i].c_str(),curData[i]);
						}
						threadData[index].emplace_back(curData);
						threadDataNum[index].emplace_back(seqs.size());
						threadReadsNum[index] += seqs.size();
						seqs.clear();
						if (index == 0)
						{
							of_log << get_local_time() << " pre-processed reads:\t" << file1_line_num / 4 << endl;
						}
					}
				}
				file1_line_num++;
			}
			else
			{
				if (!seqs.empty())
				{
					uint64_t *curData = new uint64_t[seqs.size()];
					//                    memset(curData,NULL,sizeof(uint64_t)*seqs.size());
					for (int i = 0; i < seqs.size(); i++)
					{
						curData[i] = hash<string>()(seqs[i]);
						//                        MDString(seqs[i].c_str(),curData[i]);
					}
					threadData[index].emplace_back(curData);
					threadReadsNum[index] += seqs.size();
					threadDataNum[index].emplace_back(seqs.size());
					seqs.clear();
				}
				if (multi_gzfq1[index] != NULL)
				{
					gzclose(multi_gzfq1[index]);
				}
				break;
			}
		}
	}
	else
	{
		while (1)
		{
			if (fgets(
					buf1, READBUF, multi_Nongzfq1[index]) != NULL)
			{
				if ((
						file1_line_num / thread_read_block) %
						gp.threads_num ==
					index)
				{
					block_line_num1++;
					if (block_line_num1 % 4 == 2)
					{
						fq1seq.assign(buf1);
						fq1seq.erase(
							fq1seq.size() - spaceNum, spaceNum);
						seqs.emplace_back(fq1seq);
					}
					if (seqs.size() == gp.patchSize)
					{
						uint64_t *curData = new uint64_t[seqs.size()];
						//                        memset(curData,NULL,sizeof(uint64_t)*seqs.size());
						for (int i = 0; i < seqs.size(); i++)
						{
							curData[i] = hash<string>()(seqs[i]);
							//                            MDString(seqs[i].c_str(),curData[i]);
						}
						threadData[index].emplace_back(curData);
						threadReadsNum[index] += seqs.size();
						threadDataNum[index].emplace_back(seqs.size());
						seqs.clear();
						if (index == 0)
						{
							of_log << get_local_time() << " pre-processed reads:\t" << file1_line_num / 4 << endl;
						}
					}
				}
				file1_line_num++;
			}
			else
			{
				if (!seqs.empty())
				{
					uint64_t *curData = new uint64_t[seqs.size()];
					//                    memset(curData,NULL,sizeof(uint64_t)*seqs.size());
					for (int i = 0; i < seqs.size(); i++)
					{
						curData[i] = hash<string>()(seqs[i]);
						//                        MDString(seqs[i].c_str(),curData[i]);
					}
					threadData[index].emplace_back(curData);
					threadReadsNum[index] += seqs.size();
					threadDataNum[index].emplace_back(seqs.size());
					seqs.clear();
				}
				if (multi_Nongzfq1[index] != NULL)
				{
					fclose(multi_Nongzfq1[index]);
				}
				break;
			}
		}
	}
	check_disk_available();
	//    sub_thread_done[index]=1;
	logLock.lock();
	of_log << get_local_time() << "\tthread " << index << " done\t" << endl;
	logLock.unlock();
	return &se_bq_check;
}

void seProcess::filter_se_fqs(
	SEcalOption opt, int index)
{
	// C_reads_trim_stat_2 cut_pos;
	vector<C_fastq>::iterator i_end = opt.fq1s->end();
	// check dup
	bool *dupFilter = new bool[opt.fq1s->size()];
	if (gp.rmdup)
	{
		checkDup.lock();
		memset(
			dupFilter, false, opt.fq1s->size());
		int iter = 0;
		for (vector<C_fastq>::iterator i = opt.fq1s->begin(); i != i_end; i++)
		{
			//            string checkSeq=(*i).sequence;
			//            if(checkDupMap.find(checkSeq)!=checkDupMap.end()){
			//                dupNum++;
			//                cout<<"real dup:\t"<<(*i).sequence<<endl;
			//            }else{
			//                checkDupMap.insert(checkSeq);
			//            }
			if (dupFlag[threadCurReadReadsNumIdx[index] - opt.fq1s->size() + iter])
			{
				//                    dupNum++;
				dupFilter[iter] = true;
				gzwrite(
					dupThreadOut1[index], (*i).toString().c_str(), (*i).toString().size());
			}
			iter++;
		}
		checkDup.unlock();
	}
	i_end = opt.fq1s->end();

	int iter = 0;
	for (vector<C_fastq>::iterator i = opt.fq1s->begin(); i != opt.fq1s->end(); i++)
	{
		C_single_fastq_filter se_fastq_filter = C_single_fastq_filter(
			*i, gp);
		if (dupFilter[iter])
		{
			se_fastq_filter.read_result.dup = true;
		}
		iter++;
		se_fastq_filter.se_trim(gp);
		if (gp.adapter_discard_or_trim == "trim" || gp.contam_discard_or_trim == "trim" || !gp.trim.empty() || !gp.trimBadHead.empty() || !gp.trimBadTail.empty())
		{
			(*i).head_hdcut = se_fastq_filter.read.head_hdcut;
			(*i).head_lqcut = se_fastq_filter.read.head_lqcut;
			(*i).tail_hdcut = se_fastq_filter.read.tail_hdcut;
			(*i).tail_lqcut = se_fastq_filter.read.tail_lqcut;
			(*i).adacut_pos = se_fastq_filter.read.adacut_pos;
			//(*i).contam_pos=se_fastq_filter.read.contam_pos;
			//(*i).global_contam_pos=se_fastq_filter.read.global_contam_pos;
			//(*i).raw_length=se_fastq_filter.read.raw_length;
		}
		//*i=se_fastq_filter.read;
		if (!gp.trim_fq1.empty())
		{
			preOutput(
				1, se_fastq_filter.read);
			opt.trim_result1->emplace_back(se_fastq_filter.read);
		}
		int whether_discard(0);
		if (gp.module_name == "filtersRNA")
		{
			whether_discard = se_fastq_filter.sRNA_discard(
				opt.se_local_fs, gp);
		}
		else
		{
			whether_discard = se_fastq_filter.se_discard(
				opt.se_local_fs, gp);
		}
		if (whether_discard != 1)
		{
			if (!gp.clean_fq1.empty())
			{
				preOutput(
					1, se_fastq_filter.read);
				opt.clean_result1->emplace_back(se_fastq_filter.read);
			}
		}
	}
	delete[] dupFilter;
	// return cut_pos;
}