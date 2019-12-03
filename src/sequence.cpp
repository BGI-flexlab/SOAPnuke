#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <sstream>
#include <mutex>
#include "sequence.h"
#include "global_parameter.h"
#include "read_filter.h"
using namespace::std;

C_single_fastq_filter::C_single_fastq_filter(C_fastq& a,C_global_parameter& gp){
	read=a;
	read_result=stat_read(read,gp);
}
void C_single_fastq_filter::se_trim(C_global_parameter& gp){
	fastq_trim(read,gp);
}
int C_single_fastq_filter::sRNA_discard(C_filter_stat* fs,C_global_parameter& gp){
	int min_value=-1;
	if(gp.max_read_length!=-1){
		if(read.sequence.size()>gp.max_read_length){
			fs->long_len_num++;
			min_value=1;
			return min_value;
			//cout<<"max_read_length discard"<<endl;
		}
	}
	if(gp.lowQualityBaseRatio!=-1){
		if(read_result.low_qual_base_ratio>=gp.lowQualityBaseRatio){	//check low base quality ratio whether is too high
			fs->low_qual_base_ratio_num++;
			min_value=1;
			return min_value;
			//cout<<"low_qual_base_ratio discard"<<endl;
		}
	}
	if(read_result.include_3_adapter==-1){
		fs->no_3_adapter_num++;
		min_value=1;
		return min_value;
	}else if(read_result.include_3_adapter<=2){
		fs->int_insertNull_num++;
		min_value=1;
		return min_value;
	}
	if(read_result.include_5_adapter){
		fs->include_adapter_seq_num++;
		min_value=1;
		return min_value;
	}
	if(gp.highA_ratio!=-1){
		if(read_result.a_ratio>=gp.highA_ratio){	//check highA ratio whether is too high
			fs->highA_num++;
			min_value=1;
			return min_value;
			//cout<<"highA_ratio discard"<<endl;
		}
	}
	
	if(gp.polyX_num!=-1){
		if(read_result.contig_base>=gp.polyX_num){	//check polyX ratio whether is too high
			fs->polyX_num++;
			min_value=1;
			return min_value;
			//cout<<"polyX_num discard"<<endl;
		}
	}
	if(read.sequence.size()<gp.min_read_length){
		fs->short_len_num++;
		min_value=1;
		return min_value;
		//cout<<"min_read_length discard"<<endl;
	}
	return min_value;
}
int C_single_fastq_filter::se_discard(C_filter_stat* fs,C_global_parameter& gp){
	int min_value=-1;
	if(!gp.tile.empty()){	//check read tile whether in the given removal tile list
		if(check_tile_or_fov(read_result.read_tile,gp.tile)){
			fs->tile_num++;
			min_value=1;
			return min_value;
		}
	}
	
	if(!gp.fov.empty()){	//check read fov whether in the given removal fov list
		if(check_tile_or_fov(read_result.read_fov,gp.fov)){
			fs->fov_num++;
			min_value=1;
			return min_value;
		}
	}

	if(gp.min_read_length!=-1){
		if(read.sequence.size()<gp.min_read_length){
			fs->short_len_num++;
			min_value=1;
			return min_value;
			//cout<<"min_read_length discard"<<endl;
		}
	}
	if(gp.max_read_length!=-1){
		if(read.sequence.size()>gp.max_read_length){
			fs->long_len_num++;
			min_value=1;
			return min_value;
			//cout<<"max_read_length discard"<<endl;
		}
	}
	if(gp.contam_discard_or_trim=="discard"){
		if(read_result.include_contam==1){
			fs->include_contam_seq_num++;
			min_value=1;
			return min_value;
		}
		if(read_result.include_global_contam==1){
			fs->include_global_contam_seq_num++;
			min_value=1;
			return min_value;
		}
	}
	if(gp.n_ratio!=-1){
		if(read_result.n_ratio>=gp.n_ratio){	//check N base ratio whether is too high
			fs->n_ratio_num++;
			min_value=1;
			return min_value;
			//cout<<"n_ratio discard"<<endl;
		}
	}
	//cout<<"a_ratio\t"<<read_result.a_ratio<<"\t"<<read_result.fastq2_result.a_ratio<<endl;
	if(gp.highA_ratio!=-1){
		if(read_result.a_ratio>=gp.highA_ratio){	//check highA ratio whether is too high
			fs->highA_num++;
			min_value=1;
			return min_value;
			//cout<<"highA_ratio discard"<<endl;
		}
	}
	if(gp.polyX_num!=-1){
		if(read_result.contig_base>=gp.polyX_num){	//check polyX ratio whether is too high
			fs->polyX_num++;
			min_value=1;
			return min_value;
			//cout<<"polyX_num discard"<<endl;
		}
	}
	//cout<<"lbq\t"<<read_result.low_qual_base_ratio<<"\t"<<read_result.fastq2_result.low_qual_base_ratio<<"\t"<<gp.lowQualityBaseRatio<<endl;
	if(gp.lowQualityBaseRatio!=-1){
		if(read_result.low_qual_base_ratio>=gp.lowQualityBaseRatio){	//check low base quality ratio whether is too high
			fs->low_qual_base_ratio_num++;
			min_value=1;
			return min_value;
			//cout<<"low_qual_base_ratio discard"<<endl;
		}
	}
	//cout<<"mq\t"<<read_result.mean_quality<<"\t"<<read_result.fastq2_result.mean_quality<<endl;
	if(gp.meanQuality!=-1){
		if(read_result.mean_quality<gp.meanQuality){
			fs->mean_quality_num++;
			min_value=1;
			return min_value;
			//cout<<"mean_quality discard"<<endl;
		}
	}
	if(read_result.include_adapter_seq==1 && gp.adapter_discard_or_trim=="discard"){	//whether include adapter sequence
		fs->include_adapter_seq_num++;
		min_value=1;
		return min_value;
		//cout<<"include_adapter_seq discard"<<endl;
	}
	return min_value;
}

//C_pe_fastq_filter::C_pe_fastq_filter(C_single_fastq_filter a,C_single_fastq_filter b,C_global_parameter gp){

C_pe_fastq_filter::C_pe_fastq_filter(C_fastq& a,C_fastq& b,C_global_parameter& gp){
	C_global_parameter gp2=gp;
	gp2.adaMis=gp2.adaMis2;
	gp2.adaMR=gp2.adaMR2;
	gp2.adaEdge=gp2.adaEdge2;
	C_single_fastq_filter fastq1=C_single_fastq_filter(a,gp);
	C_single_fastq_filter fastq2=C_single_fastq_filter(b,gp2);
	fq1=fastq1.read;
	fq2=fastq2.read;
	//reads_result.fastq1_result=stat_read(fastq1.read,gp);
	//reads_result.fastq2_result=stat_read(fastq2.read,gp);
	reads_result.fastq1_result=fastq1.read_result;
	reads_result.fastq2_result=fastq2.read_result;
	reads_result.over_lapped=false;
	if(gp.overlap_length!=-1)
		reads_result.over_lapped=whether_over_overlapped(a,b,gp);
}

int C_pe_fastq_filter::pe_discard(C_filter_stat* fs,C_global_parameter& gp){
	int min_value=-1;
	if(!gp.tile.empty()){	//check read tile whether in the given removal tile list
		if(check_tile_or_fov(reads_result.fastq1_result.read_tile,gp.tile)){
			fs->tile_num++;
			min_value=1;
			return min_value;
			//cout<<"tile discard"<<endl;
		}
	}
	
	if(!gp.fov.empty()){	//check read fov whether in the given removal fov list
		//cout<<reads_result.fastq1_result.read_fov<<endl;
		//cout<<gp.fov<<endl;
		if(check_tile_or_fov(reads_result.fastq1_result.read_fov,gp.fov)){
			fs->fov_num++;
			min_value=1;
			return min_value;
		}
	}

	if(gp.min_read_length!=-1){
		int v=pe_dis(fq1.sequence.size()<gp.min_read_length,fq2.sequence.size()<gp.min_read_length);
		if(v>0){
			switch(v){
				case 1:fs->short_len_num1++;break;
				case 2:fs->short_len_num2++;break;
				case 3:fs->short_len_num1++;fs->short_len_num2++;fs->short_len_num_overlap++;break;
				default:break;
			}
			fs->short_len_num++;
			return 1;
		}
			//cout<<"min_read_length discard"<<endl;
	}else{
		if(fq1.sequence.empty() || fq2.sequence.empty()) {
			return 1;
		}
	}
	if(gp.max_read_length!=-1){
		int v=pe_dis(fq1.sequence.size()>gp.max_read_length,fq2.sequence.size()>gp.max_read_length);
		if(v>0){
			switch(v){
				case 1:fs->long_len_num1++;break;
				case 2:fs->long_len_num2++;break;
				case 3:fs->long_len_num1++;fs->long_len_num2++;fs->long_len_num_overlap++;break;
				default:break;
			}
			fs->long_len_num++;
			return 1;
		}
			//cout<<"max_read_length discard"<<endl;
	}
	if(gp.contam_discard_or_trim=="discard"){
		int v=pe_dis(reads_result.fastq1_result.include_global_contam==1,reads_result.fastq2_result.include_global_contam==1);
		if(v>0){
			switch(v){
				case 1:fs->include_global_contam_seq_num1++;break;
				case 2:fs->include_global_contam_seq_num2++;break;
				case 3:fs->include_global_contam_seq_num1++;fs->include_global_contam_seq_num2++;fs->include_global_contam_seq_num_overlap++;break;
				default:break;
			}
			fs->include_global_contam_seq_num++;
			return 1;
		}
	}
	if(gp.contam_discard_or_trim=="discard"){
		int v=pe_dis(reads_result.fastq1_result.include_contam==1,reads_result.fastq2_result.include_contam==1);
		if(v>0){

			switch(v){
				case 1:fs->include_contam_seq_num1++;break;
				case 2:fs->include_contam_seq_num2++;break;
				case 3:fs->include_contam_seq_num1++;fs->include_contam_seq_num2++;fs->include_contam_seq_num_overlap++;break;
				default:break;
			}
			fs->include_contam_seq_num++;
			return 1;
		}
	}
	if(gp.n_ratio!=-1){
		int v=pe_dis(reads_result.fastq1_result.n_ratio>=gp.n_ratio,reads_result.fastq2_result.n_ratio>=gp.n_ratio);	//check N base ratio whether is too high
		if(v>0){
			switch(v){
				case 1:fs->n_ratio_num1++;break;
				case 2:fs->n_ratio_num2++;break;
				case 3:fs->n_ratio_num1++;fs->n_ratio_num2++;fs->n_ratio_num_overlap++;break;
				default:break;
			}
			fs->n_ratio_num++;
			return 1;
		}
	}
	//cout<<"a_ratio\t"<<reads_result.fastq1_result.a_ratio<<"\t"<<reads_result.fastq2_result.a_ratio<<endl;
	if(gp.highA_ratio!=-1){
		int v=pe_dis(reads_result.fastq1_result.a_ratio>=gp.highA_ratio,reads_result.fastq2_result.a_ratio>=gp.highA_ratio);	//check highA ratio whether is too high
		if(v>0){
			switch(v){
				case 1:fs->highA_num1++;break;
				case 2:fs->highA_num2++;break;
				case 3:fs->highA_num1++;fs->highA_num2++;fs->highA_num_overlap++;break;
				default:break;
			}
			fs->highA_num++;
			return 1;
		}
	}
	if(gp.polyX_num!=-1){
		int v=pe_dis(reads_result.fastq1_result.contig_base>=gp.polyX_num,reads_result.fastq2_result.contig_base>=gp.polyX_num);	//check polyX ratio whether is too high
		if(v>0){
			switch(v){
				case 1:fs->polyX_num1++;break;
				case 2:fs->polyX_num2++;break;
				case 3:fs->polyX_num1++;fs->polyX_num2++;fs->polyX_num_overlap++;break;
				default:break;
			}
			fs->polyX_num++;
			return 1;
		}
	}
	//cout<<"lbq\t"<<reads_result.fastq1_result.low_qual_base_ratio<<"\t"<<reads_result.fastq2_result.low_qual_base_ratio<<"\t"<<gp.lowQualityBaseRatio<<endl;
	if(gp.lowQualityBaseRatio!=-1){
		int v=pe_dis(reads_result.fastq1_result.low_qual_base_ratio>=gp.lowQualityBaseRatio,reads_result.fastq2_result.low_qual_base_ratio>=gp.lowQualityBaseRatio);	//check low base quality ratio whether is too high
		if(v>0){
			if(reads_result.fastq1_result.low_qual_base_ratio>1 || reads_result.fastq2_result.low_qual_base_ratio>1){
				cerr<<"Error:low quality base ratio stat error,"<<fq1.seq_id<<endl;
				exit(1);
			}
			//cout<<fq1.seq_id<<"\t"<<fq2.seq_id<<"\t"<<reads_result.fastq1_result.low_qual_base_ratio<<"\t"<<reads_result.fastq2_result.low_qual_base_ratio<<"\t"<<gp.lowQualityBaseRatio<<endl;
			switch(v){
				case 1:fs->low_qual_base_ratio_num1++;break;
				case 2:fs->low_qual_base_ratio_num2++;break;
				case 3:fs->low_qual_base_ratio_num1++;fs->low_qual_base_ratio_num2++;fs->low_qual_base_ratio_num_overlap++;break;
				default:break;
			}
			fs->low_qual_base_ratio_num++;
			return 1;
		}
	}
	//cout<<"mq\t"<<reads_result.fastq1_result.mean_quality<<"\t"<<reads_result.fastq2_result.mean_quality<<endl;
	if(gp.meanQuality!=-1){
		int v=pe_dis(reads_result.fastq1_result.mean_quality<gp.meanQuality,reads_result.fastq2_result.mean_quality<gp.meanQuality);
		if(v>0){
			switch(v){
				case 1:fs->mean_quality_num1++;break;
				case 2:fs->mean_quality_num2++;break;
				case 3:fs->mean_quality_num1++;fs->mean_quality_num2++;fs->mean_quality_num_overlap++;break;
				default:break;
			}
			fs->mean_quality_num++;
			return 1;
		}
	}
	if(gp.overlap_length!=-1){
		if(reads_result.over_lapped){
			fs->over_lapped_num++;
			min_value=1;
			return min_value;
			//cout<<"over_lapped discard"<<endl;
		}
	}
	if(gp.adapter_discard_or_trim=="discard"){
		int v=pe_dis(reads_result.fastq1_result.include_adapter_seq==1,reads_result.fastq2_result.include_adapter_seq==1);	//whether include adapter sequence
		if(v>0){
			switch(v){
				case 1:fs->include_adapter_seq_num1++;break;
				case 2:fs->include_adapter_seq_num2++;break;
				case 3:fs->include_adapter_seq_num1++;fs->include_adapter_seq_num2++;fs->include_adapter_seq_num_overlap++;break;
				default:break;
			}
			fs->include_adapter_seq_num++;
			return 1;
		}
	}
	
	return min_value;
}
void C_pe_fastq_filter::pe_trim(C_global_parameter& gp){
	fastq_trim(fq1,gp);
	fastq_trim(fq2,gp);
}
int C_pe_fastq_filter::pe_dis(bool a,bool b){
	int return_value=0;
	if(a)
		return_value+=1;
	if(b)
		return_value+=2;
	return return_value;
}