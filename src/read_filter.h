#ifndef _READ_FILTER_H
#define _READ_FILTER_H

#include <string>
#include <vector>
#include <math.h>
#include "process_argv.h"
#include "sequence.h"
#include "global_parameter.h"
using namespace::std;


bool check_tile_or_fov(string tile,string& tile_parameter);
string reversecomplementary(string a);
C_fastq_stat_result stat_read(C_fastq& fq_read,C_global_parameter& gp);
bool whether_over_overlapped(C_fastq fastq1,C_fastq fastq2,C_global_parameter& gp);
void fastq_trim(C_fastq& read,C_global_parameter& gp);
int adapter_pos(string& ref_sequence,string& adapter,C_global_parameter& gp);
int hasContam(string& ref_sequence,string& contam,C_global_parameter& gp);
int sRNA_findAdapter(string sequence,string adapter,C_global_parameter& gp);
bool sRNA_hasAdapter(string sequence, string adapter,C_global_parameter& gp);
#endif