#ifndef GC_H
#define GC_H

#include <string>
#include <string.h>
#include <vector>
#include <sstream>
#include <cctype>
#include <set>
#include <map>
#include <math.h>
#include <cstdint>
using namespace ::std;
class quartile_result
{
public:
	float mean;
	float median;
	float lower_quar, upper_quar;
	float first10_quar, last10_quar;
};

// quartile_result cal_quar_from_map(map<int,int> data);
quartile_result cal_quar_from_array(uint64_t data[], int len);
vector<string> get_pe_hard_trim(string a);
vector<string> get_se_hard_trim(string a);
void check_gz_file(string a);
int check_gz_empty(string a);
// void remove_space(string &a);
int file_exist_and_not_empty(string file_name);
// void uniq_vector(vector<string> &a);
string get_local_time();
void line_split(string line_info, char sep, vector<string> &elements);
void line_split(string line_info, char sep, set<string> &elements);
void line_split(string line_info, vector<string> &elements);
// int count_gc(string a);
void int2string(int &a, string &b);
void double2string(double &a, string &b);
void float2string(float &a, string &b);
void string2double(string &a, double &b);

void string2upper(string &a);
string join_vector(vector<string> a, char sep);
string join_vector(vector<int> a, char sep);
string join_vector(set<string> a, char sep);
string join_vector(vector<double> a, char sep);
string join_vector(vector<string> a, string seq);
string get_absolute_path(string path);
int check_bam_ok(string bin_path, string bam_file);
int check_bai_ok(string bin_path, string bam_file);
string get_exe_path(string path);
string link_dir_file(string dir, string file);
void chomp_space(string &a, string type);
void mkDir(string dir);
void mkDir(string dir, mode_t mode);
long long guessReadsNum(string fq1);
long long guessReadsNum(string fq1, string fq2);
#endif
