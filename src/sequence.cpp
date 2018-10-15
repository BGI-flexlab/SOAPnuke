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

mutex m2;
void C_fastq::clear(){
	seq_id="";
	sequence="";
	qual_seq="";
}
void C_fastq::output(){
	cout<<seq_id<<endl;
	cout<<sequence<<endl;
	cout<<qual_seq<<endl;
}
void C_fastq::output2(int type,C_global_parameter gp,gzFile outfile){
	string out_content=seq_id+"\n"+sequence+"\n+\n"+qual_seq+"\n";
	gzwrite(outfile,out_content.c_str(),out_content.size());
}
C_single_fastq_filter::C_single_fastq_filter(C_fastq& a,C_global_parameter& gp){
	read=a;
	read_result=stat_read(a,gp);
}
void C_single_fastq_filter::se_trim(C_global_parameter& gp){
	fastq_trim(read,gp);
}
int C_single_fastq_filter::sRNA_discard(C_filter_stat* fs,C_global_parameter& gp){
	int min_value=-1;
	if(gp.total_reads_num!=-1){	//global reads number limit
		if(gp.output_reads_num>gp.total_reads_num){
			fs->output_reads_num++;
			min_value=1;
			return min_value;
			//cout<<"total_reads_num discard"<<endl;
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
		if(read_result.contig_base>gp.polyX_num){	//check polyX ratio whether is too high
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
	if(gp.total_reads_num!=-1){	//global reads number limit
		if(gp.output_reads_num>gp.total_reads_num){
			fs->output_reads_num++;
			min_value=1;
			return min_value;
			//cout<<"total_reads_num discard"<<endl;
		}
	}
	if(!gp.tile.empty()){	//check read tile whether in the given removal tile list
		if(check_tile_or_fov(read_result.read_tile,gp.tile)){
			fs->tile_num++;
			min_value=1;
			return min_value;
			//cout<<"tile discard"<<endl;
		}
	}
	
	if(!gp.fov.empty()){	//check read fov whether in the given removal fov list
		if(check_tile_or_fov(read_result.read_fov,gp.fov)){
			fs->fov_num++;
			min_value=1;
			return min_value;
		}
	}
	/*
	if(read_result.in_adapter_list){	//whether in adapter list
		fs->in_adapter_list_num++;
		min_value=1;
		return min_value;
	}
	*/

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
		if(read_result.contig_base>gp.polyX_num){	//check polyX ratio whether is too high
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
	C_single_fastq_filter fastq1=C_single_fastq_filter(a,gp);
	C_single_fastq_filter fastq2=C_single_fastq_filter(b,gp);
	fq1=a;
	fq2=b;
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
	//cout<<"n_ratio\t"<<reads_result.fastq1_result.n_ratio<<"\t"<<reads_result.fastq2_result.n_ratio<<endl;
	if(gp.total_reads_num!=-1){	//global reads number limit
		if(gp.output_reads_num>gp.total_reads_num){
			fs->output_reads_num++;
			min_value=1;
			return min_value;
			//cout<<"total_reads_num discard"<<endl;
		}
	}
	if(!gp.tile.empty()){	//check read tile whether in the given removal tile list
		if(check_tile_or_fov(reads_result.fastq1_result.read_tile,gp.tile)){
			fs->tile_num++;
			min_value=1;
			return min_value;
			//cout<<"tile discard"<<endl;
		}
	}
	
	if(!gp.fov.empty()){	//check read fov whether in the given removal fov list
		if(check_tile_or_fov(reads_result.fastq1_result.read_fov,gp.fov)){
			fs->fov_num++;
			min_value=1;
			return min_value;
		}
	}
	/*
	if(reads_result.fastq1_result.in_adapter_list){	//whether in adapter list
		fs->in_adapter_list_num++;
		min_value=1;
		return min_value;
	}
	*/

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
		int v=pe_dis(reads_result.fastq1_result.contig_base>gp.polyX_num,reads_result.fastq2_result.contig_base>gp.polyX_num);	//check polyX ratio whether is too high
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

C_pe_fastqs_filter::C_pe_fastqs_filter(vector<string> seq_id1,vector<string> sequence1,vector<string> qual_seq1,vector<string> seq_id2,vector<string> sequence2,vector<string> qual_seq2,C_global_parameter gp,C_global_variable& gv){
	/*vector<int> trimmed_index;
	vector<int> filtered_index;
	vector<int> output_index;
	int reads_num;
	*/
	v_seq_id1=seq_id1;
	v_sequence1=sequence1;
	v_qual_seq1=qual_seq1;
	v_seq_id2=seq_id2;
	v_sequence2=sequence2;
	v_qual_seq2=qual_seq2;
	for(int i=0;i!=seq_id1.size();i++){
		C_fastq a,b;
		a.seq_id=seq_id1[i];
		a.sequence=sequence1[i];
		a.qual_seq=qual_seq1[i];
		b.seq_id=seq_id2[i];
		b.sequence=sequence2[i];
		b.qual_seq=qual_seq2[i];
		C_pe_fastq_filter pe_fastq_filter=C_pe_fastq_filter(a,b,gp);
		pe_fastq_filter.pe_trim(gp);
		/*
		if(pe_fastq_filter.pe_discard(gv,gp)!=1){
			filtered_index.push_back(1);
		}else{
			filtered_index.push_back(0);
		}
		*/
	}
}
int C_pe_fastq_filter::pe_dis(bool a,bool b){
	int return_value=0;
	if(a)
		return_value+=1;
	if(b)
		return_value+=2;
	return return_value;
}
void C_pe_fastqs_filter::pe_output(C_global_parameter gp,gzFile outfile1,gzFile outfile2){
	string out_content1,out_content2;
	for(int i=0;i!=filtered_index.size();i++){
		out_content1+=v_seq_id1[i]+"\n"+v_sequence1[i]+"\n+\n"+v_qual_seq1[i]+"\n";
		out_content2+=v_seq_id2[i]+"\n"+v_sequence2[i]+"\n+\n"+v_qual_seq2[i]+"\n";
		
		/*
		if(filtered_index[i]==1){
			if(gp.whether_add_pe_info){
				v_seq_id1[i]=v_seq_id1[i]+"/1";
				v_seq_id2[i]=v_seq_id2[i]+"/2";
			}
			if(!gp.base_convert.empty()){
				gp.base_convert=gp.base_convert.replace(gp.base_convert.find("TO"),2,"");
				gp.base_convert=gp.base_convert.replace(gp.base_convert.find("2"),1,"");
				if(gp.base_convert.size()!=2){
					cerr<<"Error:base_conver value format error"<<endl;
					exit(1);
				}
				for(string::size_type ix=0;ix!=v_sequence1[i].size();ix++){
					if(toupper(v_sequence1[i][ix])==toupper(gp.base_convert[0]))
						v_sequence1[i][ix]=gp.base_convert[1];
				}
				for(string::size_type ix=0;ix!=v_sequence2[i].size();ix++){
					if(toupper(v_sequence2[i][ix])==toupper(gp.base_convert[0]))
						v_sequence2[i][ix]=gp.base_convert[1];
				}
			}
			if(gp.output_file_type=="fasta"){
				v_seq_id1[i]=v_seq_id1[i].replace(v_seq_id1[i].find("@"),1,"");
				out_content1=">"+v_seq_id1[i]+"\n"+v_sequence2[i]+"\n";
				v_seq_id1[i]=v_seq_id2[i].replace(v_seq_id2[i].find("@"),1,"");
				out_content2=">"+v_seq_id2[i]+"\n"+v_sequence2[i]+"\n";
			}else if(gp.output_file_type=="fastq"){
				if(gp.outputQualityPhred!=gp.qualityPhred){
					for(string::size_type ix=0;ix!=v_qual_seq1[i].size();ix++){
						v_qual_seq1[i][ix]=(char)(v_qual_seq1[i][ix]-atoi(gp.outputQualityPhred.c_str()));
					}
				}
				
				//out_content1=v_seq_id1[i]+"\n"+v_sequence1[i]+"\n+\n"+v_qual_seq1[i]+"\n";
				//out_content2=v_seq_id2[i]+"\n"+v_sequence2[i]+"\n+\n"+v_qual_seq2[i]+"\n";
				//gzwrite(outfile1,out_content1.c_str(),out_content1.size());
				//gzwrite(outfile2,out_content2.c_str(),out_content2.size());
			}else{
				cerr<<"Error:output_file_type value error"<<endl;
				exit(1);
			}
				
		}
		*/
	}
	if(gp.is_streaming){
		cout<<out_content1;
	}else{
		gzwrite(outfile1,out_content1.c_str(),out_content1.size());
		gzwrite(outfile2,out_content2.c_str(),out_content2.size());
	}

}