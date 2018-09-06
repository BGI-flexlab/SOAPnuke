#include "global_variable.h"

C_general_stat::C_general_stat(){
	read_max_length=0;
	read_length=0;
	reads_number=0;
	base_number=0;
	a_number=0;
	c_number=0;
	g_number=0;
	t_number=0;
	n_number=0;
	//a_ratio,c_ratio,g_ratio,t_ratio,n_ratio=0;
	q20_num=0;
	q30_num=0;
}
C_reads_pos_base_stat::C_reads_pos_base_stat(){
	for(int i=0;i!=READ_MAX_LEN;i++){
		for(int j=0;j!=5;j++){
			position_acgt_content[i][j]=0;
		}
	}
}
C_reads_pos_qual_stat::C_reads_pos_qual_stat(){
	for(int i=0;i!=READ_MAX_LEN;i++){
		for(int j=0;j!=MAX_QUAL;j++){
			position_qual[i][j]=0;
		}
	}
}
C_reads_trim_stat::C_reads_trim_stat(){
	for(int i=0;i!=READ_MAX_LEN;i++){
		hlq[i]=0;
		ht[i]=0;
		ta[i]=0;
		tlq[i]=0;
		tt[i]=0;
	}
}
C_fastq_file_stat::C_fastq_file_stat(){
	gs=C_general_stat();
	//bs=C_reads_pos_base_stat();
	//qs=C_reads_pos_qual_stat();
	//ts=C_reads_trim_stat();
}
C_global_variable::C_global_variable(){
	fs=C_filter_stat();
	raw1_stat=C_fastq_file_stat();
	raw2_stat=C_fastq_file_stat();
	trim1_stat=C_fastq_file_stat();
	trim2_stat=C_fastq_file_stat();
	clean1_stat=C_fastq_file_stat();
	clean2_stat=C_fastq_file_stat();
}
