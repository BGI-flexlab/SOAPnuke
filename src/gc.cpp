//#include<iostream>
//#include<algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include <set>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include "zlib.h"
#include "gc.h"
//#include<stdio.h>

using namespace::std;
int MAXLEN=1000;
//char separator='\t';
quartile_result cal_quar_from_map(map<int,int> data){
	quartile_result return_value;
	float total_num(0);
	int data_num(0);
	for(map<int,int>::iterator ix=data.begin();ix!=data.end();ix++){
		total_num+=(ix->first)*(ix->second);
		data_num+=ix->second;
	}
	return_value.mean=total_num/data_num;
	int lower_pos=data_num/4;
	int upper_pos=data_num*3/4;
	int first10_pos=data_num/10;
	int last10_pos=data_num*9/10;
	int median_pos=data_num/2;
	int last_cur_pos(0),cur_pos(0);
	for(map<int,int>::iterator ix=data.begin();ix!=data.end();ix++){
		cur_pos+=ix->second;
		if(lower_pos>=last_cur_pos && lower_pos<=cur_pos){
			return_value.lower_quar=ix->first;
		}
		if(upper_pos>=last_cur_pos && upper_pos<=cur_pos){
			return_value.upper_quar=ix->first;
		}
		if(first10_pos>=last_cur_pos && first10_pos<=cur_pos){
			return_value.first10_quar=ix->first;
		}
		if(last10_pos>=last_cur_pos && last10_pos<=cur_pos){
			return_value.last10_quar=ix->first;
		}
		if(median_pos>=last_cur_pos && median_pos<=cur_pos){
			return_value.median=ix->first;
		}
		last_cur_pos=cur_pos;
	}
	return return_value;
}
quartile_result cal_quar_from_array(long long data[],int len){
	//cout<<len<<"\there"<<endl;
	quartile_result return_value;
	unsigned long long total_num(0);
	int data_num(0);
	for(int i=0;i<=len;i++){
		total_num+=i*data[i];
		data_num+=data[i];
	}
	if(data_num==0){
		return_value.mean=0;
	}else{
		return_value.mean=total_num/data_num;
	}
	int lower_pos=data_num/4;
	int upper_pos=data_num*3/4;
	int first10_pos=data_num/10;
	int last10_pos=data_num*9/10;
	int median_pos=data_num/2;
	int last_cur_pos(0),cur_pos(0);
	for(int i=0;i<=len;i++){
		cur_pos+=data[i];
		if(lower_pos>=last_cur_pos && lower_pos<=cur_pos){
			return_value.lower_quar=i;
		}
		if(upper_pos>=last_cur_pos && upper_pos<=cur_pos){
			return_value.upper_quar=i;
		}
		if(first10_pos>=last_cur_pos && first10_pos<=cur_pos){
			return_value.first10_quar=i;
		}
		if(last10_pos>=last_cur_pos && last10_pos<=cur_pos){
			return_value.last10_quar=i;
		}
		if(median_pos>=last_cur_pos && median_pos<=cur_pos){
			return_value.median=i;
		}
		last_cur_pos=cur_pos;
	}
	return return_value;
}
void check_gz_file(string a){
	gzFile gzfp=gzopen((a).c_str(),"rb");
	if(!gzfp){
		cerr<<"Error:cannot open file,"<<a<<endl;
		exit(1);
	}
	char line_info[MAXLEN];
	if(gzgets(gzfp,line_info,MAXLEN)==NULL){
		cerr<<"Error:empty file"<<endl;
		exit(1);
	}
	gzclose(gzfp);
}
int check_gz_empty(string a){
	gzFile gzfp=gzopen((a).c_str(),"rb");
	if(!gzfp){
		return -1;
	}
	char line_info[MAXLEN];
	if(gzgets(gzfp,line_info,MAXLEN)!=NULL){
		return 1;
	}
    gzclose(gzfp);
	return -1;
	
}
void remove_space(string &a){
	string b;
	for(string::size_type ix=0;ix!=a.size();ix++){
		if(!isspace(a[ix])){
			b+=a[ix];
		}
	}
	a=b;
}
int file_exist_and_not_empty(string file_name){
	int return_value=1;
	ifstream if_input;
	if_input.open(file_name.c_str());
	if(!if_input){
		return_value=0;
		return return_value;
	}
	if(if_input.peek()==EOF){
		return_value=0;
		return return_value;
	}
	return return_value;
}
void uniq_vector(vector<string> &a){
	sort(a.begin(),a.end());
	vector<string>::iterator end_iter=unique(a.begin(),a.end());
	a.erase(end_iter,a.end());
}
string get_local_time(){
	time_t a;
	time(&a);
	struct tm* b=localtime(&a);
	int cur_year=b->tm_year+1900;
	int cur_mon=b->tm_mon+1;
	int cur_day=b->tm_mday;
	ostringstream out_time;
	out_time<<cur_year<<"-"<<cur_mon<<"-"<<cur_day<<"  "<<b->tm_hour<<":"<<b->tm_min<<":"<<b->tm_sec<<endl;
	string out=out_time.str();
	out.erase(out.end()-1);
	return out;
}
void line_split(string line_info,char sep,vector<string> &elements){
	elements.clear();
	string element;
	for(string::size_type ix=0;ix!=line_info.size();ix++){
		if(line_info[ix]!=sep){
			element+=line_info[ix];
		}else{
			elements.emplace_back(element);
			element="";
		}
	}
	elements.emplace_back(element);
}
void line_split(string line_info,vector<string> &elements){
	elements.clear();
	string element;
	for(string::size_type ix=0;ix!=line_info.size();ix++){
		if(!isspace(line_info[ix])){
			element+=line_info[ix];
		}else{
			elements.emplace_back(element);
			element="";
		}
	}
	elements.emplace_back(element);
}
void line_split(string line_info,char sep,set<string> &elements){
	elements.clear();
	string element;
	for(string::size_type ix=0;ix!=line_info.size();ix++){
		if(line_info[ix]!=sep){
			element+=line_info[ix];
		}else{
			elements.insert(element);
			element="";
		}
	}
	elements.insert(element);
}
int count_gc(string a){
	string::size_type t_index=0,t_index2;
	int gc_counts=0;
	while((t_index2=a.find_first_of("GC",t_index))!=string::npos){
		t_index=t_index2+1;
		gc_counts++;
	}
	return gc_counts;
}
void int2string(int &a,string &b){
        stringstream ss; 
        ss<<a;
        b=ss.str();
}
void float2string(float &a,string &b){
		stringstream ss;
		ss<<a;
		b=ss.str();
}
void double2string(double &a,string &b){
	stringstream ss;
	ss<<a;
	b=ss.str();
}
void string2double(string &a,double &b){
	stringstream ss;
	ss<<a;
	ss>>b;
}
void string2upper(string &a){
	for(string::iterator ix=a.begin();ix!=a.end();ix++){
		*ix=toupper(*ix);
	}
}
string join_vector(vector<string> a,char sep){
	string result;
	for(vector<string>::iterator ix=a.begin();ix!=a.end();ix++){
		if(ix!=a.end()-1){
			result+=*ix+sep;
		}else{
			result+=*ix;
		}
	}
	return result;
}
string join_vector(vector<int> a,char sep){
	string result;
	string tmp;
	for(vector<int>::iterator ix=a.begin();ix!=a.end();ix++){
		int2string(*ix,tmp);
		if(ix!=a.end()-1){
			result+=tmp+sep;
		}else{
			result+=tmp;
		}
	}
	return result;
}
string join_vector(vector<double> a,char sep){
	string result;
	for(vector<double>::iterator ix=a.begin();ix!=a.end();ix++){
		string tmp_a;
		double2string(*ix,tmp_a);
		if(ix!=a.end()-1){
			result+=tmp_a+sep;
		}else{
			result+=tmp_a;
		}
	}
	return result;
}
string join_vector(set<string> a,char sep){
	string result;
	for(set<string>::iterator ix=a.begin();ix!=a.end();ix++){
		result+=*ix+sep;
	}
	result.erase(result.size()-1,1);
	return result;
}
string join_vector(vector<string> a,string sep){
	string result;
	for(vector<string>::iterator ix=a.begin();ix!=a.end();ix++){
		if(ix!=a.end()-1){
			result+=*ix+sep;
		}else{
			result+=*ix;
		}
	}
	return result;
}
string link_dir_file(string dir,string file){
	//string::size_type index;
	// remove null char in the start or end;
	dir.erase(0,dir.find_first_not_of(" "));  
	dir.erase(dir.find_last_not_of(" ")+1);
/*	while((index=dir.find(" ",index))!=string::npos){       //remove all null char
		dir.erase(index,1);
	}
*/
	if(dir[dir.size()-1]!='/'){
		dir=dir+"/";
	}
	if(file.find("/")!=string::npos){
		string::size_type ix=file.find("/");
		file.erase(0,ix+1);
	}
	string out=dir+file;
	return out;
}
string get_absolute_path(string path){
	char* _path=getenv("PATH");
	string s_path(_path);
	vector<string> path_eles;
	line_split(s_path,':',path_eles);
	struct stat buf;
	for(vector<string>::iterator ix=path_eles.begin();ix!=path_eles.end();ix++){
		string process_path=link_dir_file(*ix,path);
		if(lstat(process_path.c_str(),&buf)<0){
			continue;
		}
		if(S_ISREG(buf.st_mode) && S_IXUSR&(buf.st_mode)){
			return process_path;
		}
	}
	cerr<<"no such file: "<<path<<endl;
	exit(1);
}
int check_bam_ok(string bin_path,string bam_file){
	ifstream if_pre_check(bam_file.c_str());
	if(!if_pre_check)
		return 0;
	if_pre_check.close();
	string check_cmd=bin_path+"/samtools quickcheck "+bam_file+" && echo 'all ok' || echo 'fail!'";
	FILE* s_check=popen(check_cmd.c_str(),"r");
	char line[100];
	while(fgets(line,100,s_check)!=NULL){
		string s_line(line);
		if(s_line.find("ok")!=string::npos){
			return 1;
		}
		if(s_line.find("fail")!=string::npos){
			return -1;
		}
	}
	cerr<<"ERROR:bam check error"<<endl;
	exit(1);
}
int check_bai_ok(string bin_path,string bam_file){
	string check_cmd=bin_path+"/samtools tview -d T "+bam_file+" 2>&1";
	FILE* s_check=popen(check_cmd.c_str(),"r");
	char line[1000];
	while(fgets(line,1000,s_check)!=NULL){
		string s_line(line);
		if(s_line.find("Cannot read index for")!=string::npos){
			return 0;
		}
	}
	return 1;
}
string get_exe_path(string path){
	if(path.find("/")!=string::npos){
		return path.substr(0,path.find_last_of("/")+1);
	}else{
		string abs_path=get_absolute_path(path);
		return path.substr(0,abs_path.find_last_of("/")+1);
	}
}
void chomp_space(string& a,string type){
	if(type=="head"){
		int s_num=0;
		for(string::size_type ix=0;ix!=a.size();ix++){
			if(!isspace(a[ix])){
				break;
			}else{
				s_num++;
			}
		}
		a.erase(0,s_num);
	}else if(type=="tail"){
		string::size_type s_iter=0;;
		int s_num=0;
		for(string::size_type ix=a.size()-1;ix>=0;ix++){
			if(!isspace(a[ix])){
				s_iter=ix;
				break;
			}else{
				s_num++;
			}
		}
		a.erase(s_iter,s_num);
	}else{
		int s_num=0;
		for(string::size_type ix=0;ix!=a.size();ix++){
			if(!isspace(a[ix])){
				break;
			}else{
				s_num++;
			}
		}
		a.erase(0,s_num);
		string::size_type s_iter=0;;
		int s_num2=0;
		for(string::size_type ix=a.size()-1;ix>=0;ix++){
			if(!isspace(a[ix])){
				s_iter=ix;
				break;
			}else{
				s_num++;
			}
		}
		a.erase(s_iter,s_num2);
	}
}
vector<string> get_pe_hard_trim(string a){
	vector<string> tmp_eles;
	line_split(a,',',tmp_eles);
	if(tmp_eles.size()!=4){
		cerr<<"Error:trim value format error"<<endl;
		exit(1);
	}
	return tmp_eles;
}
vector<string> get_se_hard_trim(string a){
	vector<string> tmp_eles;
	line_split(a,',',tmp_eles);
	if(tmp_eles.size()!=2){
		cerr<<"Error:trim value format error"<<endl;
		exit(1);
	}
	return tmp_eles;
}

void mkDir(string dir) {
    struct stat s;
    if(stat(dir.c_str(),&s)==0){
        if(S_ISDIR(s.st_mode)){
            return;
        }
    }
    if(mkdir(dir.c_str(),0755)!=0){
        cerr<<"Error:mkdir fail,"<<dir<<endl;
        cerr<<strerror(errno)<<endl;
        exit(errno);
    }
}

void mkDir(string dir, mode_t mode) {
    struct stat s;
    if(stat(dir.c_str(),&s)==0){
        if(S_ISDIR(s.st_mode)){
            return;
        }
    }
    if(mkdir(dir.c_str(),mode)!=0){
        cerr<<"Error:mkdir fail,"<<dir<<endl;
        cerr<<strerror(errno)<<endl;
        exit(errno);
    }
}

long long guessReadsNum(string fq1) {
    struct stat s;
    long long fileSize=0;
    float conserveFlow=1.2;
    if(stat(fq1.c_str(),&s)==0){
        fileSize=s.st_size;
    }
    gzFile inFq1=gzopen(fq1.c_str(),"r");
    int bufSize=1024*1024*10;
    if(bufSize>fileSize){
        bufSize=fileSize/2;
    }
    char* buf=new char[bufSize];
    int readSize=0;
    int lfNum=0;
    int cSize=0;
    if((readSize=gzread(inFq1,buf,bufSize))>0){
        for(int i=0;i<readSize;i++){
            if(buf[i]=='\n'){
                lfNum++;
            }
        }
        string tmpFq=fq1+".tmp";
        gzFile outTmpFq=gzopen(tmpFq.c_str(),"wb");
        gzwrite(outTmpFq,buf,readSize);
        gzclose(outTmpFq);
        struct stat s;
        if(stat(tmpFq.c_str(),&s)==0){
            cSize=s.st_size;
        }
        if(remove(tmpFq.c_str())!=0){
            cerr<<"Error:cannot remove file,"<<tmpFq<<endl;
            exit(1);
        }
    }
    if(lfNum==0){
        cerr<<"Error:no reads found in input file,"<<fq1<<endl;
        exit(1);
    }
    long long readsNum=round(lfNum/4);
    if(cSize==0){
        cerr<<"Error:empty file,"<<fq1<<endl;
        exit(1);
    }
    long long guessedReadsNum=(readsNum*fileSize/cSize)*conserveFlow;
    gzclose(inFq1);


    delete[] buf;
    return guessedReadsNum;
}
