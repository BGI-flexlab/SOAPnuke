#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <map>
#include "read_filter.h"
#include "global_parameter.h"
#include "gc.h"
using namespace::std;

#define MAX_LENGTH 500
bool check_tile_or_fov(string tile,string& tile_parameter){
	if(tile_parameter.find("C")==string::npos){
		if(tile_parameter.find(",")==string::npos){
			if(tile_parameter.find("-")==string::npos){
				if(tile!=tile_parameter){
					return -1;
				}else{
					return 1;
				}
			}else{
				vector<string> line_eles;
				line_split(tile_parameter,'-',line_eles);
				if(line_eles.size()!=2){
					cerr<<"Error:input tile parameter format error,"<<tile_parameter<<endl;
					exit(1);
				}
				for(int target_tile=atoi(line_eles[0].c_str());target_tile<=atoi(line_eles[1].c_str());target_tile++){
					if(tile==tile_parameter){
						return 1;
					}
				}
				return -1;
			}
		}else{
			vector<string> line_eles;
			line_split(tile_parameter,',',line_eles);
			for(vector<string>::iterator ix=line_eles.begin();ix!=line_eles.end();ix++){
				if((*ix).find("-")!=string::npos){
					vector<string> line_eles2;
					line_split(*ix,'-',line_eles2);
					if(line_eles2.size()!=2){
						cerr<<"Error:input tile parameter format error,"<<*ix<<endl;
						exit(1);
					}
					for(int target_tile=atoi(line_eles2[0].c_str());target_tile<=atoi(line_eles2[1].c_str());target_tile++){
						if(tile==*ix){
							return 1;
						}
					}
				}else{
					if(*ix==tile){
						return 1;
					}
				}
			}
			return -1;
		}
	}else{
		if(tile_parameter.find(",")==string::npos){
			if(tile!=tile_parameter){
				return -1;
			}else{
				return 1;
			}
		}else{
			vector<string> line_eles;
			line_split(tile_parameter,',',line_eles);
			for(vector<string>::iterator ix=line_eles.begin();ix!=line_eles.end();ix++){
				if(*ix==tile){
					return 1;
				}
			}
			return -1;
		}
	}
}
C_fastq_stat_result stat_read(C_fastq& fq_read,C_global_parameter& gp){ //stat some feature of read,not include adapter-related information
	C_fastq_stat_result return_value;
	return_value.in_adapter_list=0;
	return_value.include_adapter_seq=-1;
	fq_read.raw_length=fq_read.sequence.size();
	if(gp.tile.empty()){
		return_value.read_tile="";
	}else{
		int i=0;
		int num=0;
		if(gp.seq_type=="0"){
			for(;i<fq_read.seq_id.size();i++){
				if(fq_read.seq_id[i]==':')
					num++;
				if(num>=2)
					break;
				if(num<2 || i+4>=fq_read.seq_id.size()){
					cerr<<"Warning:input fastq maybe not include tile information\n"<<endl;
				}
			}
		}else{
			for(;i<fq_read.seq_id.size();i++){
				if(fq_read.seq_id[i]==':')
					num++;
				if(num>=4)
					break;
				if(num<4 || i+4>=fq_read.seq_id.size()){
					cerr<<"Warning:input fastq maybe not include tile information\n"<<endl;
				}
			}
		}
		for(int j=0; j!=4; j++){
            if(fq_read.seq_id[i+j+1]>='0' && fq_read.seq_id[i+j+1] <= '9')
            {
                return_value.read_tile.insert(return_value.read_tile.end(),fq_read.seq_id[i+j+1]);
            }else{
                cerr<<"Warning:input fastq maybe not include tile information\n"<<endl;
            }
		}
		if(return_value.read_tile.size()!=4){
			cerr<<"Warning:input fastq maybe not include tile information\n"<<endl;
		}
	}
	if(gp.fov.empty()){
		return_value.read_fov="";
	}else{
		int i=0;
        if(gp.seq_type=="0"){
            for(i=0;i<fq_read.seq_id.size();i++){
                if(fq_read.seq_id[i]=='C')
                    if(i+8<fq_read.seq_id.size() && fq_read.seq_id[i+4]=='R' && fq_read.seq_id[i+8]=='_')
                        break;
            }
        }else{
            cerr<<"Warning:Zebra-500 data(--fov), --seqType is 0"<<endl;
            exit(1);
        }
        for(int j=0;j!=8;j++){
            return_value.read_fov.insert(return_value.read_fov.end(),fq_read.seq_id[i+j]);
        }
        if(return_value.read_fov.size()!=8){
        	cerr<<"Warning:input fastq maybe not include fov information\n"<<endl;
        }
	}
	
/*
	if(!gp.adapter_file.empty()){	//check whether in adapter list;
		set<string> adapter_list;	
		ifstream if_adapter_file(gp.adapter_file.c_str());
		if(!if_adapter_file){
			cerr<<"Error:cannot open file,"<<adapter_list<<endl;
			exit(1);
		}else{
			string tmp_line;
			while(getline(if_adapter_file,tmp_line)){
				adapter_list.insert(tmp_line);
			}
		}
		if(adapter_list.find(fq_read.seq_id)!=adapter_list.end()){
			return_value.in_adapter_list=1;
		}
	}
*/	
    if(gp.module_name=="filtersRNA"){
    	int ada_pos=sRNA_findAdapter(fq_read.sequence,fq_read.adapter_seq2,gp);
    	return_value.include_3_adapter=ada_pos;
    	return_value.include_5_adapter=sRNA_hasAdapter(fq_read.sequence,fq_read.adapter_seq,gp);
    }else{
    	int ada_pos=adapter_pos(fq_read.sequence,fq_read.adapter_seq,gp);
	    if(ada_pos>=0){
	    	return_value.include_adapter_seq=1;
	    	fq_read.adacut_pos=fq_read.sequence.size()-ada_pos;
	    }
	    if(fq_read.contam_seq.find(",")==string::npos){
	    	fq_read.contam_pos=hasContam(fq_read.sequence,fq_read.contam_seq,gp);
	    	if(fq_read.contam_pos>=0)
	    		return_value.include_contam=1;
	    }else{
	    	vector<int> contam_poses=hasContams(fq_read.sequence,fq_read.contam_seq,gp);
	    	for(vector<int>::iterator ix=contam_poses.begin();ix!=contam_poses.end();ix++){
	    		if(*ix>=0){
	    			return_value.include_contam=1;
	    			if(fq_read.contam_pos!=-1){
	    				if(*ix<=fq_read.contam_pos){
			    			fq_read.contam_pos=*ix;
			    		}
	    			}else{
	    				fq_read.contam_pos=*ix;
	    			}
	    		}
	    	}
	    }
	    if(!gp.global_contams.empty()){
	    	vector<int> global_contam_poses=hasGlobalContams(fq_read.sequence,gp);
	    	for(vector<int>::iterator ix=global_contam_poses.begin();ix!=global_contam_poses.end();ix++){
	    		if(*ix>=0){
	    			return_value.include_global_contam=1;
	    			if(fq_read.global_contam_pos!=-1){
	    				if(*ix<=fq_read.global_contam_pos){
			    			fq_read.global_contam_pos=*ix;
			    		}
	    			}else{
	    				fq_read.global_contam_pos=*ix;
	    			}
	    		}
	    	}
	    }
    }
	if((fq_read.sequence).size()==0){
		cerr<<"Error:empty sequence"<<endl;
		exit(1);
	}
	return_value.seq_len=(fq_read.sequence).size();
	char last_char('Q');
	int contig_base(0);
	int max_contig(1);
	for(string::size_type ix=0;ix!=return_value.seq_len;ix++){	//process base sequence
		//switch(tolower((fq_read.sequence)[ix])){
		if(gp.polyX_num!=-1){
			if(fq_read.sequence[ix]==last_char){
				contig_base++;
				if(max_contig<contig_base)
					max_contig=contig_base;
			}else{
				contig_base=1;
			}
		}
		switch((fq_read.sequence)[ix]){
		//switch((fq_read.sequence)[ix]){
			case 'a':
			case 'A':return_value.a_num++;break;
			case 'c':
			case 'C':return_value.c_num++;break;
			case 'g':
			case 'G':return_value.g_num++;break;
			case 't':
			case 'T':return_value.t_num++;break;
			case 'n':
			case 'N':return_value.n_num++;break;
			default:{
				cerr<<"Error:unrecognized sequence,"<<(fq_read.sequence)<<endl;
				exit(1);
			}
		}
	}

	return_value.contig_base=max_contig;
	return_value.a_ratio=float(return_value.a_num)/(fq_read.sequence).size();
	return_value.c_ratio=float(return_value.c_num)/(fq_read.sequence).size();
	return_value.g_ratio=float(return_value.g_num)/(fq_read.sequence).size();
	return_value.t_ratio=float(return_value.t_num)/(fq_read.sequence).size();
	return_value.n_ratio=float(return_value.n_num)/(fq_read.sequence).size();
	return_value.gc_ratio=float((return_value.g_num+return_value.c_num))/(fq_read.sequence).size();
	int qual_basis=gp.qualityPhred;
	int i_low_qual_threshold=gp.lowQual;
	int total_base_quality=0;
	for(string::size_type ix=0;ix!=return_value.seq_len;ix++){	//process quality sequence
		int base_quality=(fq_read.qual_seq)[ix]-qual_basis;
		total_base_quality+=base_quality;
		if(base_quality<=i_low_qual_threshold)
			return_value.low_qual_base_num++;
		if(base_quality>=20)
			return_value.q20_num++;
		if(base_quality>=30)
			return_value.q30_num++;
	}
	return_value.qual_len=(fq_read.qual_seq).size();
	return_value.low_qual_base_ratio=float(return_value.low_qual_base_num)/(fq_read.qual_seq).size();
	return_value.mean_quality=(float)total_base_quality/(fq_read.qual_seq).size();
	return return_value;
}

bool whether_over_overlapped(C_fastq fastq1,C_fastq fastq2,C_global_parameter& gp){
        	//code is similar to SOAPnuke
    string seq2 = reversecomplementary(fastq2.sequence);
    int mismatch, max_mismatch;
    int max_match_length = (fastq1.sequence.size() > fastq2.sequence.size()) ? fastq2.sequence.size() : fastq1.sequence.size();
    for (int i = gp.overlap_length; i <= max_match_length; i++){
        max_mismatch = (int)(gp.peMismatchRatio * (float)i);
        mismatch = 0;
        for (int j = 0; j < i; j++){
            if (fastq1.sequence[fastq1.sequence.size() - i + j] == 'N' || seq2[j] == 'N')
            {
                mismatch++;
            }
            else if (fastq1.sequence[fastq1.sequence.size() - i + j] != seq2[j])
            {
                mismatch++;
            }
        }
        if (mismatch <= max_mismatch)
            return true;
    }
    return false;
}
void fastq_trim(C_fastq& read,C_global_parameter& gp){	//	1.index_remove	2.adapter_trim	3.hard_trim	4.qual_trim

	string hard_trim=gp.trim;
	string low_qual_head_trim=gp.trimBadHead;
	string low_qual_tail_trim=gp.trimBadTail;
	bool ht_flag(0),lqt_flag(0),index_flag(0),ada_trim_flag(0),contam_trim_flag(0);
	if(!hard_trim.empty())
		ht_flag=1;
	if(!low_qual_head_trim.empty() && !low_qual_tail_trim.empty())
		lqt_flag=1;
	if(gp.index_remove)
		index_flag=1;
	if(gp.adapter_discard_or_trim=="trim")
		ada_trim_flag=1;
	if(gp.contam_discard_or_trim=="trim")
		contam_trim_flag=1;
	if(!ht_flag && !lqt_flag && !index_flag && !ada_trim_flag && !contam_trim_flag && gp.polyG_tail==-1){
		return;
	}else{
		if(index_flag){	//remove index. Code is similar to SOAP nuke
			if(gp.seq_type=="0"){
			/* old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2\n";
        cout << "\t                                  new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC\n";
        */	
                read.seq_id=read.seq_id.substr(0,read.seq_id.find("#")-0);
            }else{
            	read.seq_id=read.seq_id.substr(0,read.seq_id.find_last_of("#")-0);
            }
		}
		int head_cut(0),tail_cut(0);
		if(ht_flag){	//hard trim
			read.head_hdcut=atoi(read.head_trim_len.c_str());
			read.tail_hdcut=atoi(read.tail_trim_len.c_str());
			head_cut=read.head_hdcut;
			tail_cut=read.tail_hdcut;
		}
		if(lqt_flag){	//low quality end trim
			vector<string> head_eles,tail_eles;
			int head_low_qual_threshold,head_low_qual_trim_length_limit,tail_low_qual_threshold,tail_low_qual_trim_length_limit;
			line_split(low_qual_head_trim,',',head_eles);
			line_split(low_qual_tail_trim,',',tail_eles);
			if(head_eles.size()!=2 || tail_eles.size()!=2){
				cerr<<"Error:low quality base at end format error,"<<low_qual_head_trim<<" "<<low_qual_head_trim<<endl;
				exit(1);
			}
			head_low_qual_threshold=atoi(head_eles[0].c_str());
			head_low_qual_trim_length_limit=atoi(head_eles[1].c_str());
			tail_low_qual_threshold=atoi(tail_eles[0].c_str());
			tail_low_qual_trim_length_limit=atoi(tail_eles[1].c_str());
			int qual_basis=gp.qualityPhred;
			int head_ix(0),tail_ix(0);
			for(int ix=0;ix!=head_low_qual_trim_length_limit;ix++){
				int base_quality=read.qual_seq[ix]-qual_basis;
				if(base_quality<head_low_qual_threshold){
					head_ix++;
				}else{
					break;
				}
			}
			for(int ix=0;ix!=tail_low_qual_trim_length_limit;ix++){
				int base_quality=read.qual_seq[read.qual_seq.size()-ix-1]-qual_basis;
				if(base_quality<head_low_qual_threshold){
					tail_ix++;
				}else{
					break;
				}
			}
			read.head_lqcut=head_ix;
			read.tail_lqcut=tail_ix;
			head_cut=head_cut>=head_ix?head_cut:head_ix;
			tail_cut=tail_cut>=tail_ix?tail_cut:tail_ix;
		}
		if(ada_trim_flag){
			//cout<<read.sequence<<endl;
			if(gp.module_name=="filtersRNA"){
				int ada_pos=sRNA_findAdapter(read.sequence,read.adapter_seq2,gp);
				if(ada_pos>2 && ada_pos<read.sequence.size()){
					read.sequence=read.sequence.substr(0,ada_pos);
					read.qual_seq=read.qual_seq.substr(0,ada_pos);
				}
			}
			if(read.adacut_pos>0){
				tail_cut=tail_cut>=read.adacut_pos?tail_cut:read.adacut_pos;
			}
		}
		if(contam_trim_flag){
			if(read.global_contam_pos>=0 && read.sequence.size()-read.global_contam_pos>tail_cut){
				tail_cut=read.sequence.size()-read.global_contam_pos;
			}
			if(read.contam_pos>=0 && read.sequence.size()-read.contam_pos>tail_cut){
				tail_cut=read.sequence.size()-read.contam_pos;
			}
		}
		if(gp.polyG_tail!=-1){
			int polyG_n=polyG_number(read.sequence);
			if(polyG_n>=gp.polyG_tail){
				if(polyG_n>tail_cut){
					tail_cut=polyG_n;
				}
			}
		}
		read.sequence=read.sequence.substr(head_cut,read.sequence.size()-head_cut-tail_cut);
		read.qual_seq=read.qual_seq.substr(head_cut,read.qual_seq.size()-head_cut-tail_cut);
	}
}
int polyG_number(string& ref_sequence){
	int polyG_tail_number(0);
	for(int i=ref_sequence.size()-1;i>=0;i--){
		if(ref_sequence[i]=='G' || ref_sequence[i]=='g'){
			polyG_tail_number++;
		}else{
			break;
		}
	}
	return polyG_tail_number;
}
vector<int> hasContams(string& ref_sequence,string& contam,C_global_parameter& gp){
	if(gp.ctMatchR.find(",")==string::npos){
		cerr<<"Error:the number of ctMatchR value should equal to that of contam sequences"<<endl;
		exit(1);
	}else{
		vector<string> contams,mrs;
		line_split(contam,',',contams);
		line_split(gp.ctMatchR,',',mrs);
		if(contams.size()!=mrs.size()){
			cerr<<"Error:the number of ctMatchR value should equal to that of contam sequences,"<<contams.size()<<".vs."<<mrs.size()<<endl;
			exit(1);
		}
		vector<int> return_poses;
		for(int i=0;i!=contams.size();i++){
			float tmp_mr=atof(mrs[i].c_str());
			int tmp_pos=hasContam(ref_sequence,contams[i],gp,tmp_mr);
			return_poses.push_back(tmp_pos);
			if(tmp_pos>=0 && tmp_pos<gp.min_read_length){
				break;
			}
		}
		return return_poses;
	}
}
int hasContam(string& ref_sequence,string& contam,C_global_parameter& gp,float mr){
	int readLen=ref_sequence.size();
	int contamLen=contam.size();
	if(contamLen==0){
		return -1;
	}
	float misGrad = (contamLen-gp.adaEdge)/(gp.adaMis+1);
    int r1, mis, maxSegMatch;
    int segMatchThr = (int)ceil(contamLen * mr);
    float segGrad;
    if(segMatchThr-7+1==0){
    	segGrad=0;
    }else{
    	segGrad = (contamLen-gp.adaEdge)/(segMatchThr-7+1);
    }
    int misMatchTemp, segMatchTemp;
    for (r1 = 0; r1 < contamLen - gp.adaEdge; ++r1)
    {
        mis = 0;
        maxSegMatch = 0;
        misMatchTemp = r1/misGrad;
        if(segGrad!=0){
        	segMatchTemp = 7 + r1/segGrad;
        }else{
        	segMatchTemp=7;
        }
        for (int c = 0; c < r1+gp.adaEdge; ++c)
        {
            if (contam[contamLen-r1-gp.adaEdge+c] == ref_sequence[c]){
                maxSegMatch++;
                if (maxSegMatch >= segMatchTemp)
                    return 0;
            }
            else{
                if (ref_sequence[c] != 'N'){
                    mis++;
                    maxSegMatch = 0;
                    if (mis > misMatchTemp)
                        break;
                }
            }
        }
        if (mis <= misMatchTemp)
            return 0;
    }

    for (r1 = 0; r1 <= readLen - contamLen; ++r1)
    {
        maxSegMatch = 0;
        mis = 0;

        for (int c = 0; c < contamLen; ++c)
        {
            if (contam[c] == ref_sequence[r1 + c]){
                maxSegMatch ++;
                if (maxSegMatch >= segMatchThr)
                    return r1;
            }
            else{
                if (ref_sequence[r1 + c] != 'N'){
                    mis++;
                    maxSegMatch = 0;
                    if (mis > gp.adaMis)
                        break;
                }
            }
        }
        if (mis <= gp.adaMis)
            return r1;
    }

    for (r1 = 0; r1 < contamLen - gp.adaEdge; ++r1)
    {
        mis = 0;
        maxSegMatch = 0;
        misMatchTemp = r1/misGrad;
        segMatchTemp = 7 + r1/segGrad;

        for (int c = 0; c < r1+gp.adaEdge; ++c)
        {
            if (contam[c] == ref_sequence[readLen-r1-gp.adaEdge+c]){
                maxSegMatch++;
                if (maxSegMatch >= segMatchTemp)
                    return readLen-r1-gp.adaEdge;
            }
            else{
                if (ref_sequence[readLen-r1-gp.adaEdge+c] != 'N'){
                    mis++;
                    maxSegMatch = 0;
                    if (mis > misMatchTemp)
                        break;
                }
            }
        }
        if (mis <= misMatchTemp)
            return readLen-r1-gp.adaEdge;
    }
    return -1;
}
int hasContam(string& ref_sequence,string& contam,C_global_parameter& gp){
	int readLen=ref_sequence.size();
	int contamLen=contam.size();
	if(contamLen==0){
		return -1;
	}
	float misGrad = (contamLen-gp.adaEdge)/(gp.adaMis+1);
    int r1, mis, maxSegMatch;
    int segMatchThr = (int)ceil(contamLen * atof(gp.ctMatchR.c_str()));
    float segGrad(0);
    if(segMatchThr-7+1==0){
    	segGrad=0;
    }else{
    	segGrad = (contamLen-gp.adaEdge)/(segMatchThr-7+1);
    }
    int misMatchTemp, segMatchTemp;
    for (r1 = 0; r1 < contamLen - gp.adaEdge; ++r1)
    {
        mis = 0;
        maxSegMatch = 0;
        misMatchTemp = r1/misGrad;
        if(segGrad!=0){
        	segMatchTemp = 7 + r1/segGrad;
        }else{
        	segMatchTemp=7;
        }
        

        for (int c = 0; c < r1+gp.adaEdge; ++c)
        {
            if (contam[contamLen-r1-gp.adaEdge+c] == ref_sequence[c]){
                maxSegMatch++;
                if (maxSegMatch >= segMatchTemp)
                    return 0;
            }
            else{
                if (ref_sequence[c] != 'N'){
                    mis++;
                    maxSegMatch = 0;
                    if (mis > misMatchTemp)
                        break;
                }
            }
        }
        if (mis <= misMatchTemp)
            return 0;
    }

    for (r1 = 0; r1 <= readLen - contamLen; ++r1)
    {
        maxSegMatch = 0;
        mis = 0;

        for (int c = 0; c < contamLen; ++c)
        {
            if (contam[c] == ref_sequence[r1 + c]){
                maxSegMatch ++;
                if (maxSegMatch >= segMatchThr)
                    return r1;
            }
            else{
                if (ref_sequence[r1 + c] != 'N'){
                    mis++;
                    maxSegMatch = 0;
                    if (mis > gp.adaMis)
                        break;
                }
            }
        }
        if (mis <= gp.adaMis)
            return r1;
    }

    for (r1 = 0; r1 < contamLen - gp.adaEdge; ++r1)
    {
        mis = 0;
        maxSegMatch = 0;
        misMatchTemp = r1/misGrad;
        segMatchTemp = 7 + r1/segGrad;

        for (int c = 0; c < r1+gp.adaEdge; ++c)
        {
            if (contam[c] == ref_sequence[readLen-r1-gp.adaEdge+c]){
                maxSegMatch++;
                if (maxSegMatch >= segMatchTemp)
                    return readLen-r1-gp.adaEdge;
            }
            else{
                if (ref_sequence[readLen-r1-gp.adaEdge+c] != 'N'){
                    mis++;
                    maxSegMatch = 0;
                    if (mis > misMatchTemp)
                        break;
                }
            }
        }
        if (mis <= misMatchTemp)
            return readLen-r1-gp.adaEdge;
    }
    return -1;
}
int adapter_pos(string& ref_sequence,string& adapter,C_global_parameter& gp){
	if(gp.module_name=="filter" || gp.module_name=="filterMeta"){
		int adptLen=adapter.size();
		if(adptLen==0){
			return -1;
		}
		int readLen=ref_sequence.size();
		float misGrad = (adptLen-gp.adaEdge)/(gp.adaMis+1);
        int r1, mis, maxSegMatch;
        int segMatchThr = (int)ceil(adptLen * gp.adaMR);
        int misMatchTemp, segMatchTemp;
        int minEdge5 = adptLen - 5;
        for (r1 = 0; r1 < adptLen - minEdge5; ++r1)
        {
            mis = 0;
            maxSegMatch = 0;
            misMatchTemp = r1/misGrad;

            for (int c = 0; c < r1+minEdge5; ++c)
            {
                if (adapter[adptLen-r1-minEdge5+c] == ref_sequence[c]){
                    maxSegMatch++;
                    if (maxSegMatch >= segMatchThr)
                        return 0;
                }
                else{
                    mis++;
                    maxSegMatch = 0;
                    if (mis > misMatchTemp)
                        break;
                }
            }
            if (mis <= misMatchTemp)
                return 0;
        }
        for (r1 = 0; r1 <= readLen - adptLen; ++r1)
        {
            maxSegMatch = 0;
            mis = 0;

            for (int c = 0; c < adptLen; ++c)
            {
                if (adapter[c] == ref_sequence[r1 + c]){
                    maxSegMatch ++;
                    if (maxSegMatch >= segMatchThr)
                        return r1;
                }
                else{
                    mis++;
                    maxSegMatch = 0;
                    if (mis > gp.adaMis)
                        break;
                }
            }
            if (mis <= gp.adaMis)
                return r1;
        }
        for (r1 = 0; r1 < adptLen - gp.adaEdge; ++r1)
        {
            mis = 0;
            maxSegMatch = 0;
            misMatchTemp = r1/misGrad;

            for (int c = 0; c < r1+gp.adaEdge; ++c)
            {
                if (adapter[c] == ref_sequence[readLen-r1-gp.adaEdge+c]){
                    maxSegMatch++;
                    if (maxSegMatch >= segMatchThr)
                        return readLen-r1-gp.adaEdge;
                }
                else{
                    mis++;
                    maxSegMatch = 0;
                    if (mis > misMatchTemp)
                        break;
                }
            }

            if (mis <= misMatchTemp)
                return readLen-r1-gp.adaEdge;
        }
        return -1;
	}
}
int sRNA_findAdapter(string sequence,string adapter,C_global_parameter& gp)
{
	int startPos = -1;
	int readLen = sequence.size();
	int adptLen = adapter.size();
	if(adptLen==0){
		return -1;
	}
	int a1 = 2;
	int r1 = 0;
	int len;
	int mis;
	bool flagType = false;
	int misTmp = 0;
	int totalMapTmp = 0;
	for (r1 = 0; r1 <= readLen - gp.adaRMa;)
	{
		int len1 = adptLen - a1;
		int len2 = readLen - r1;
		len = (len1 < len2) ? len1 : len2;
		mis = 0;
		int map[MAX_LENGTH];
		map[0] = 0;
		int totalMap = 0;
		for (int c = 0; c < len; c++)
		{
			if (sequence[r1 + c] =='N'){
				continue;
			}
			if (adapter[a1 + c] == sequence[r1 + c])
			{
				map[mis]++;
				totalMap++;
			}
			else
			{
				mis++;
				map[mis] = 0;
			}
		}
		int misAndMap = mis + totalMap;
		float rate =1.0 * mis / totalMap;
		if (mis <= gp.adaRMm && misAndMap >= gp.adaRMa && rate <= gp.adaREr)
		{
			if(flagType)
			{
				if (mis <= misTmp && totalMap >= totalMapTmp)
				{
					startPos = r1;
					misTmp = mis;
					totalMapTmp = totalMap;
				}
			}
			else
			{
				startPos = r1;
				flagType = true;
				misTmp = mis;
				totalMapTmp = totalMap;
			}
		}
		if (a1 > 0)
		{
			a1--;
		}
		else
		{
			r1++;
		}
	}
	return startPos;
}
bool sRNA_hasAdapter(string sequence, string adapter,C_global_parameter& gp)
{
	bool find = false;
	int readLen = sequence.size();
	int adptLen = adapter.size();
	if(adptLen==0){
		return find;
	}
	int a1 = adptLen - gp.adaRCtg;
	int r1 = 0;
	int len;
	int mis;
	int readLenSmall = (readLen - gp.adaRCtg < 0) ? 0 : (readLen - gp.adaRCtg);
	for (r1 = 0; r1 <= readLenSmall;)
	{
		int len1 = adptLen - a1;
		int len2 = readLen - r1;
		len = (len1 < len2) ? len1 : len2;
		mis = 0;
		int map[MAX_LENGTH];
		map[0] = 0;
		int totalMap=0;
		for (int c = 0; c < len; c++)
		{
			if (adapter[a1 + c] == sequence[r1 + c])
			{
				map[mis]++;
				totalMap++;
			}
			else
			{
				mis++;
				map[mis] = 0;
			}
		}
		int max_map = 0;
		for (int c = 0; c <= mis; c++)
		{
			if (map[c] > max_map)
			{
				max_map = map[c];
			}
		}
		if (mis <= 4 && (max_map >= gp.adaRCtg || readLen < 12) && (1.0 * totalMap / readLen >= gp.adaRAr || 1.0 * totalMap / adptLen >= gp.adaRAr))
		{
			find = true;
			break;
		}
		/*		if (max_map >= 15){
				find = true;
				break;
				}*/
		if (a1 > 0)
		{
			a1--;
		}
		else
		{
			r1++;
		}
	}

	return find;
}
vector<int> hasGlobalContams(string& ref_sequence,C_global_parameter& gp){
	vector<int> return_poses;
	vector<string> g_contams,g_mr,g_mm;
	line_split(gp.global_contams,',',g_contams);
	line_split(gp.g_mrs,',',g_mr);
	line_split(gp.g_mms,',',g_mm);
	if(g_contams.size()!=g_mr.size() || g_contams.size()!=g_mm.size()){
		cerr<<"Error:the number of global contamination sequences should equal to that of related parameters"<<endl;
		exit(1);
	}
	for(int i=0;i!=g_contams.size();i++){
		float mr=atof(g_mr[i].c_str());
		int mm=atoi(g_mm[i].c_str());
		int pos=global_contam_pos(ref_sequence,g_contams[i],mr,mm);
		string reverse_contam=reversecomplementary(g_contams[i]);
		int reverse_pos=global_contam_pos(ref_sequence,reverse_contam,mr,mm);
		int push_pos;
		if(pos>=0){
			if(reverse_pos>=0){
				push_pos=pos<reverse_pos?pos:reverse_pos;
			}else{
				push_pos=pos;
			}
		}else{
			push_pos=reverse_pos;
		}
		return_poses.push_back(push_pos);
		if(push_pos>=0 && push_pos<gp.min_read_length){
			break;
		}
	}
	return return_poses;
}
int global_contam_pos(string& ref_sequence,string& global_contam,float min_matchRatio,int mismatch_number){
	int mismatch_score=-200;
	int match_score=1;
	//float min_matchRatio=0.3;
	//int mismatch_number=1;
	int rl=ref_sequence.size();
	int cl=global_contam.size();
	int min_match_len=int(cl*min_matchRatio);
	int lower_score=(min_match_len-mismatch_number)-mismatch_number*mismatch_score;
	int total_score(-1000),overlap(0);
	//contam sequence is at front of read sequence
	for(int i=cl-min_match_len;i>=0;i--){
		int j_max=cl-i>rl?rl:cl-i;
		for(int j=0;j!=j_max;j++){
			if(ref_sequence[j]==global_contam[i+j]){
				if(total_score>mismatch_number*mismatch_score){
					total_score+=match_score;
					overlap++;
				}else{
					total_score=match_score;
					overlap=1;
				}
			}else{
				if(total_score>mismatch_number*mismatch_score){
					total_score+=mismatch_score;
					overlap++;
				}
			}
			if(total_score>=lower_score && overlap>=min_match_len){
				return 0;
			}
		}
	}
	//contam sequence is at middle of read sequence
	total_score=-1000;
	overlap=0;
	for(int i=0;i<=rl-cl;i++){
		for(int j=0;j!=cl;j++){
			if(ref_sequence[i+j]==global_contam[j]){
				if(total_score>mismatch_number*mismatch_score){
					total_score+=match_score;
					overlap++;
				}else{
					total_score=match_score;
					overlap=1;
				}
			}else{
				if(total_score>mismatch_number*mismatch_score){
					total_score+=mismatch_score;
					overlap++;
				}
			}
			if(total_score>=lower_score && overlap>=min_match_len){
				return i+j-overlap+1;
			}
		}
	}
	//contam sequence is at tail of read sequence
	total_score=-1000;
	overlap=0;
	int i_min=cl>rl?cl-rl:0;
	for(int i=i_min;i<=cl-min_match_len;i++){
		for(int j=0;j!=cl-i;j++){
			if(ref_sequence[rl-(cl-i)+j]==global_contam[j]){
				if(total_score>mismatch_number*mismatch_score){
					total_score+=match_score;
					overlap++;
				}else{
					total_score=match_score;
					overlap=1;
				}
			}else{
				if(total_score>mismatch_number*mismatch_score){
					total_score+=mismatch_score;
					overlap++;
				}
			}
			if(total_score>=lower_score && overlap>=min_match_len){
				return rl-cl+i+j-overlap+1;
			}
		}
	}
	return -1;
}
int smithWatermanAlign(string query, string target) { // Smith-Waterman algorithm
		// the match in query must be started on base 0
	float match_[5][5];
	string chr_("ACGTN");
	int queryLen = query.length();
	int targetLen = target.length();
	int i, j;
	int open = 100;

	float **scoreMatrix, **directMatrix; // 0: up, 1: left, 2: northwest, 3: itself
	int maxcol, maxrow;
	float maxScore;
	scoreMatrix = (float **) malloc(sizeof(float *) * (targetLen + 1));
	for (i = 0; i <= targetLen; i++) {
		scoreMatrix[i] = (float *) malloc(sizeof(float) * (targetLen + 1));
	}
	directMatrix = (float **) malloc(sizeof(float *) * (targetLen + 1));
	for (i = 0; i <= targetLen; i++) {
		directMatrix[i] = (float *) malloc(sizeof(float) * (targetLen + 1));
	}

	for (i = 0; i <= targetLen; i++) {
		scoreMatrix[0][i] = 0;
		directMatrix[0][i] = 1;
	}
	for (i = 0; i <= queryLen; i++) {
		scoreMatrix[i][0] = i * open * (-1);
		directMatrix[i][0] = 0;
	}
	directMatrix[0][0] = 3;
	for (i = 1; i <= queryLen; i++) {
		for (j = 1; j <= targetLen; j++) {
			float a = scoreMatrix[i - 1][j - 1] + match_[chr_.find(toupper(query[i - 1]))][chr_.find(toupper(target[j - 1]))];
			float b = scoreMatrix[i - 1][j] - open;
			float c = scoreMatrix[i][j - 1] - open;
			float temp1;
			if (b >= c) {
				directMatrix[i][j] = 0;
				temp1 = b;
			} else {
				directMatrix[i][j] = 1;
				temp1 = c;
			}
			if (a >= temp1) {
				directMatrix[i][j] = 2;
				scoreMatrix[i][j] = a;
			} else {
				scoreMatrix[i][j] = temp1;
			}
			if (maxScore < scoreMatrix[i][j]) {
				maxScore = scoreMatrix[i][j];
				maxcol = j;
				maxrow = i;
			}
		}
	}
	// trace back start point
	int row = queryLen;
	int col = targetLen;

	// output
	float colmaxS = scoreMatrix[0][targetLen];
	row = 0;
	for (j = 1; j <= queryLen; j++) {
		if (scoreMatrix[j][targetLen] >= colmaxS) {
			row = j;
			colmaxS = scoreMatrix[j][targetLen];
		}
	}

	int queryS = -1, queryE, targetS = -1, targetE;
	targetE = col;
	queryE = row;
	while (directMatrix[row][col] != 3) {
		if (directMatrix[row][col] == 0) {
			row--;
		} else if (directMatrix[row][col] == 1) {
			col--;
		} else if (directMatrix[row][col] == 2) {
			targetS = col;
			queryS = row;
			col--;
			row--;
		}
	}

	// free memory
	for (i = 0; i <= targetLen; i++) {
		free(scoreMatrix[i]);
	}
	for (i = 0; i <= targetLen; i++) {
		free(directMatrix[i]);
	}
	free(scoreMatrix);
	free(directMatrix);

	//cout << queryS << "," << queryE << "," << targetS << "," << targetE;
	return targetS;
	
}
string reversecomplementary(string& a){	//get reverse complementary sequence
	string b;
	map<char,char> dna_base_pair;
	dna_base_pair['A']='T';
	dna_base_pair['T']='A';
	dna_base_pair['G']='C';
	dna_base_pair['C']='G';
	for(string::size_type ix=0;ix!=a.size();ix++){
		char tmp_base=toupper(a[ix]);
		if(tmp_base=='N'){
			b.insert(b.end(),tmp_base);
		}else if(dna_base_pair.find(tmp_base)==dna_base_pair.end()){
			cerr<<"Error:unrecognized base,"<<a<<endl;
			exit(1);
		}else{
			b.insert(b.end(),dna_base_pair[tmp_base]);
		}
	}
	return b;
}
