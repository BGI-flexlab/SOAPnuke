# -*- coding: UTF-8 -*-
# !/usr/bin/env python
# Filename = SOAPnuke_streaming.py

import cmd
import commands
import argparse
import os
from fnmatch import fnmatch
import re
import sys


class SOAPnuke(cmd.Cmd):
    intro = commands.getoutput('./SOAPnuke')            # The path should be ensured correct.
    intro += '\nType "help <ModuleName>" (such as "help filter") to list corresponding commands.\nType "<ModuleName> <parameters>" to start the program.\n'
    prompt = '(SOAPnuke)'
    file = None

    # ----- Four modules to choose -----
    def help_filter(self):
        print commands.getoutput('./SOAPnuke filter -h')
        print '\n'
        SOAPnuke.args_streaming().print_help()

    def do_filter(self, arg):
        (args1, args) = SOAPnuke.args2cmd(arg, 'filter')
        os.system('%s fs -rm -r -skipTrash %s' % (args1['hadoop'], args1['outdir_hdfs']))
        os.system(args)
        SOAPnuke.merge(args1)

    def help_filtersrna(self):
        print commands.getoutput('./SOAPnuke filtersrna -h')
        print '\n'
        SOAPnuke.args_streaming().print_help()

    def do_filtersrna(self, arg):
        (args1, args) = SOAPnuke.args2cmd(arg, 'filtersrna')
        os.system('%s fs -rm -r -skipTrash %s' % (args1['hadoop'], args1['outdir_hdfs']))
        os.system(args)
        SOAPnuke.merge(args1)

    def help_filtermeta(self):
        print commands.getoutput('./SOAPnuke filtermeta -h')
        print '\n'
        SOAPnuke.args_streaming().print_help()

    def do_filtermeta(self, arg):
        (args1, args) = SOAPnuke.args2cmd(arg, 'filtermeta')
        os.system('%s fs -rm -r -skipTrash %s' % (args1['hadoop'], args1['outdir_hdfs']))
        os.system(args)
        SOAPnuke.merge(args1)

    def help_filterdge(self):
        print commands.getoutput('./SOAPnuke filterdge -h')
        print '\n'
        SOAPnuke.args_streaming().print_help()

    def do_filterdge(self, arg):
        (args1, args) = SOAPnuke.args2cmd(arg, 'filterdge')
        os.system('%s fs -rm -r -skipTrash %s' % (args1['hadoop'], args1['outdir_hdfs']))
        os.system(args)
        SOAPnuke.merge(args1)

    # ----- Args Parser -----
    @staticmethod
    def args_streaming(*parser):  # The default paras below should be ensured correct.
        if not parser:
            parser = argparse.ArgumentParser()
        parser.add_argument("--hadoop", help="hadoop Path, default as /usr/bin/hadoop",                   type=str, default='/usr/bin/hadoop')
        parser.add_argument("--java",   help="java Path, default as /usr/bin/java",                       type=str, default='/usr/bin/java')
        parser.add_argument("--jar",    help="Streaming_fq.jar, default to be in the streaming folder",   type=str, default='./streaming/Streaming_fq.jar')
        parser.add_argument("--lj", dest='libjars',      help='SuffixMultipleTextOutputFormat.jar, default to be in the streaming folder', type=str,
                            default='./streaming/SuffixMultipleTextOutputFormat.jar')
        parser.add_argument("--jn", dest='jobname',      help='jobname, default as SOAPnuke_streaming',   type=str, default='SOAPnuke_streaming')
        parser.add_argument("--mt", dest='maptasks',     help='number of map_tasks',                      type=int)
        parser.add_argument("--rt", dest='reducetasks',  help='number of reduce_tasks',                   type=int)
        parser.add_argument("--input",                   help='input_directory',                          type=str)
        parser.add_argument("--oh", dest='outdir_hdfs',  help='output_directory_on_hdfs',                 type=str)
        parser.add_argument("--ol", dest='outdir_local', help='output_directory_on_local',                type=str)
        parser.add_argument("--merge", action='store_true', help='whether merge fastq pieces or not')
        return parser

    @staticmethod
    def args2cmd(arg, module):
        args1 = vars(SOAPnuke.args_streaming().parse_known_args(arg.split())[0])         # args_dict for streaming
        for i in args1:
            if args1[i] is None:
                print 'Error: %s required.' % i
                sys.exit()
            if i in ['java', 'hadoop', 'jar', 'libjars'] and not os.path.exists(args1[i]):
                print 'Error: %s not exists. Please check it again.' % i
                sys.exit()

        args2 = ' '.join(SOAPnuke.args_streaming().parse_known_args(arg.split())[1][0])  # args_str for SOAPnuke
        args = '''%s jar %s -D mapred.job.name="%s"
                    -D mapred.map.tasks=%s -D mapred.reduce.tasks=%s -D
                    -D stream.num.map.output.key.fields=3 -D num.key.fields.for.partition=2 -libjars %s
                    -partitioner org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner
                    -inputformat org.bgi.flexlab.hadoop.format.fastq.FqText
                    -outputformat org.bgi.flexlab.gaea.data.mapreduce.output.text.SuffixMultipleTextOutputFormat
                    -input %s -output %s -mapper "/bin/cat" -reducer "./SOAPnuke %s %s" ''' \
               % (args1['hadoop'],       args1['jar'],     args1['jobname'], args1['maptasks'], args1['reducetasks'],
                  args1['libjars'], args1['input'],   args1['outdir_hdfs'], module, args2)
        
        return args1, args.replace('\n', ' ')

    # ----- Merge streaming result -----
    @staticmethod
    def merge(args):
        # MapReduce
        (indir, outdir_h, outdir_l, javapath, Rpath, hadooppath, merge_mark) = \
            (args[i] for i in ('input', 'outdir_hdfs', 'outdir_local', 'java', 'R', 'hadoop', 'merge'))

        # Move statistics files from HDFS to local
        os.system('%s fs -copyToLocal "%s/part-*-S" %s' % (hadooppath, outdir_h, outdir_l))
        os.system('%s fs -rm "%s/part-*-S"' % (hadooppath, outdir_h))

        # Load stats from files
        def adds(num1, num2):
            try:
                return num1 + num2
            except:
                if num1:
                    return num1
                else:
                    return num2

        def matrix_add(matrix1, matrix2):
            if len(matrix1) < len(matrix2):
                matrix1.extend(['' for x in xrange(len(matrix2) - len(matrix1))])
            elif len(matrix1) > len(matrix2):
                matrix2.extend(['' for x in xrange(len(matrix1) - len(matrix2))])
            return map(lambda x, y: map(lambda m, n: adds(m, n), x, y), matrix1, matrix2)

        pe_mark = False
        for File in os.listdir(outdir_l):
            if fnmatch(File, '*-S'):
                with open(outdir_l + '/' + File, 'r') as F:
                    lines = F.readlines()
                    header_list = []
                    for i in xrange(len(lines)):
                        if re.match('#', lines[i].strip()):
                            header_list.append([lines[i].strip(), i])
                            if re.match('#Fq2', lines[i].strip()):
                                pe_mark = True
                    header_list.append(['null', len(lines)])   # The last one corresponds to the bottom of file

                    qc_dict = {}
                    for i in xrange(len(header_list)-1):
                        if 'statistical_information' in header_list[i][0]:
                            qc_dict[header_list[i][0]] = \
                                map(lambda x, y: adds(x, y), map(int, lines[header_list[i][1] + 1].strip().split()),
                                    qc_dict[header_list[i][0]])
                        if header_list[i][0] not in qc_dict:
                            qc_dict[header_list[i][0]] = lines[header_list[i][1]+1: header_list[i+1][1]]
                        else:
                            qc_dict[header_list[i][0][:1]+'Fq2_'+header_list[i][0][1:]] = \
                                lines[header_list[i][1]+1: header_list[i+1][1]]
                os.remove(File)

        # Stats Computation
        # -- Basic stats
        def div_perc(num1, num2, round_num=2):
            return '(' + str(round(float(num1)/float(num2)*100, round_num)) + ')'

        temp1 = qc_dict['#Fq1_statistical_information']
        bs = [['Item', 'rawreads(fq1)', 'cleanreads(fq1)'],
              ['Read length', temp1[0], temp1[1]],
              ['Total number of reads', (temp1[2], '(100.00%)'), (temp1[3], '(100.00%)')],
              ['Number of filtered reads (%)',
               (int(temp1[2]) - int(temp1[3]), div_perc(int(temp1[2]) - int(temp1[3]), int(temp1[2]))), '-'],
              '',
              ['Total number of bases', (temp1[5], '(100.00%)'), (temp1[6], '(100.00%)')],
              ['Number of filtered bases (%)',
               (int(temp1[5]) - int(temp1[6]), div_perc(int(temp1[5]) - int(temp1[6]), int(temp1[5]))), '-'],
              '',
              ['Number of base A (%)',
               (temp1[7], div_perc(temp1[7], temp1[5])), (temp1[8], div_perc(temp1[8], temp1[5]))],
              ['Number of base C (%)',
               (temp1[9], div_perc(temp1[9], temp1[5])), (temp1[10], div_perc(temp1[10], temp1[5]))],
              ['Number of base G (%)',
               (temp1[11], div_perc(temp1[11], temp1[5])), (temp1[12], div_perc(temp1[12], temp1[5]))],
              ['Number of base T (%)',
               (temp1[13], div_perc(temp1[13], temp1[5])), (temp1[14], div_perc(temp1[14], temp1[5]))],
              ['Number of base N (%)',
               (temp1[15], div_perc(temp1[15], temp1[5])), (temp1[16], div_perc(temp1[16], temp1[5]))],
              '',
              ['Number of base calls with quality value of 20 or higher (Q20+) (%)',
               (temp1[17], div_perc(temp1[17], temp1[5])), (temp1[18], div_perc(temp1[18], temp1[5]))],
              ['Number of base calls with quality value of 30 or higher (Q30+) (%)',
               (temp1[19], div_perc(temp1[19], temp1[5])), (temp1[20], div_perc(temp1[20], temp1[5]))]]

        if pe_mark:
            temp2 = qc_dict['#Fq2_statistical_information']
            bs[0].extend(['rawreads(fq2)', 'cleanreads(fq2)'])
            bs[1].extend([temp2[0], temp2[1]])
            bs[2].extend([(temp2[2], '(100.00%)'), (temp2[3], '(100.00%)')])
            bs[3].extend([(int(temp2[2])-int(temp2[3]), div_perc(int(temp2[2])-int(temp2[3]), int(temp2[2]))), '-'])
            bs[5].extend([(temp2[5], '(100.00%)'), (temp2[6], '(100.00%)')])
            bs[6].extend([(int(temp2[5])-int(temp2[6]), div_perc(int(temp2[5])-int(temp2[6]), int(temp2[5]))), '-'])
            bs[8].extend([(temp2[7], div_perc(temp2[7], temp2[5])), (temp2[8], div_perc(temp2[8], temp2[5]))])
            bs[9].extend([(temp2[9], div_perc(temp2[9], temp2[5])), (temp2[10], div_perc(temp2[10], temp2[5]))])
            bs[10].extend([(temp2[11], div_perc(temp2[11], temp2[5])), (temp2[12], div_perc(temp2[12], temp2[5]))])
            bs[11].extend([(temp2[13], div_perc(temp2[13], temp2[5])), (temp2[14], div_perc(temp2[14], temp2[5]))])
            bs[12].extend([(temp2[15], div_perc(temp2[15], temp2[5])), (temp2[16], div_perc(temp2[16], temp2[5]))])
            bs[14].extend([(temp2[17], div_perc(temp2[17], temp2[5])), (temp2[18], div_perc(temp2[18], temp2[5]))])
            bs[15].extend([(temp2[19], div_perc(temp2[19], temp2[5])), (temp2[20], div_perc(temp2[20], temp2[5]))])

        # -- Filter info          Caution that blank lines exist in those content.
        temp1 = qc_dict['#Fq1_statistical_information']
        flr_total_sum = sum([temp1[28], temp1[29], temp1[30], temp1[31], temp1[32], temp1[33], temp1[34]])
        flr = [['Item', 'Total', 'Percentage'],
               ['Total filtered reads (%)', flr_total_sum, '100.00%'],
               ['Reads with adapter (%)', temp1[29], div_perc(temp1[29], flr_total_sum)],
               ['Reads with low quality (%)', temp1[31], div_perc(temp1[31], flr_total_sum)],
               ['Reads with low mean quality (%)', temp1[32], div_perc(temp1[32], flr_total_sum)],
               ['Reads with duplications (%)', temp1[28], div_perc(temp1[28], flr_total_sum)],
               ['Read with n rate exceed: (%)', temp1[30], div_perc(temp1[30], flr_total_sum)],
               ['Read with small insert size: (%)', temp1[33], div_perc(temp1[33], flr_total_sum)],
               ['Reads with PolyA (%)', temp1[34], div_perc(temp1[34], flr_total_sum)]]

        if pe_mark:
            temp2 = qc_dict['#Fq2_statistical_information']
            flr1_sum = sum([temp1[21], temp1[22], temp1[23], temp1[24], temp1[25], temp1[26], temp1[27]])
            flr2_sum = sum([temp2[21], temp2[22], temp2[23], temp2[24], temp2[25], temp2[26], temp2[27]])
            flr[0].extend(['Counts(fq1)', 'Percentage', 'Counts(fq2)', 'Percentage'])
            flr[1].extend([flr1_sum, '100.00%', flr2_sum, '100.00%'])
            flr[2].extend([temp1[29], div_perc(temp1[29], flr1_sum), temp2[29], div_perc(temp2[29], flr2_sum)])
            flr[3].extend([temp1[31], div_perc(temp1[31], flr1_sum), temp2[31], div_perc(temp2[31], flr2_sum)])
            flr[4].extend([temp1[32], div_perc(temp1[32], flr1_sum), temp2[32], div_perc(temp2[32], flr2_sum)])
            flr[5].extend([temp1[28], div_perc(temp1[28], flr1_sum), temp2[28], div_perc(temp2[28], flr2_sum)])
            flr[6].extend([temp1[30], div_perc(temp1[30], flr1_sum), temp2[30], div_perc(temp2[30], flr2_sum)])
            flr[7].extend([temp1[33], div_perc(temp1[33], flr1_sum), temp2[33], div_perc(temp2[33], flr2_sum)])
            flr[8].extend([temp1[34], div_perc(temp1[34], flr1_sum), temp2[34], div_perc(temp2[34], flr2_sum)])

        # -- Per base- Base content distribution
        temp1 = qc_dict['#Base_distributions_by_read_position']
        pbbd1 = [['Pos', 'A', 'C', 'G',	'T', 'N', 'Clean A', 'Clean C', 'Clean G', 'Clean T', 'Clean N']]
        for i in temp1:
            (raw_sum, clean_sum) = (sum(map(int, i[:5])), sum(map(int, i[5:])))
            pbbl = [str(i+1)]    # A list for temporary storage of each line
            for ii in i[:5]:
                pbbl.append(div_perc(ii, raw_sum))
            for ii in i[5:]:
                pbbl.append(div_perc(ii, clean_sum))
            pbbd1.append(pbbl)

        if pe_mark:
            temp2 = qc_dict['#Fq2_Base_distributions_by_read_position']
            pbbd2 = [['Pos', 'A', 'C', 'G', 'T', 'N', 'Clean A', 'Clean C', 'Clean G', 'Clean T', 'Clean N']]
            for i in temp2:
                (raw_sum, clean_sum) = (sum(map(int, i[:5])), sum(map(int, i[5:])))
                pbbl = [str(i+1)]  # A list for temporary storage of each line
                for ii in i[:5]:
                    pbbl.append(div_perc(ii, raw_sum))
                for ii in i[5:]:
                    pbbl.append(div_perc(ii, clean_sum))
                pbbd2.append(pbbl)

        # -- Per base- Quality scores--
        temp1 = qc_dict['#Raw_Base_quality_value_distribution_by_read_position']
        raw_pbq1 = [['Pos']]
        raw_pbq1[0].extend(['Q%s' % i for i in xrange(len(temp1[0]))])
        raw_pbq1.extend([temp1[i].insert(0, str(i+1)) for i in xrange(len(temp1))])

        temp1 = qc_dict['#Clean_Base_quality_value_distribution_by_read_position']
        clean_pbq1 = [['Pos']]
        clean_pbq1[0].extend(['Q%s' % i for i in xrange(len(temp1[0]))])
        clean_pbq1.extend([temp1[i].insert(0, str(i+1)) for i in xrange(len(temp1))])

        if pe_mark:
            temp2 = qc_dict['#Raw_Base_quality_value_distribution_by_read_position']
            raw_pbq2 = [['Pos']]
            raw_pbq2[0].extend(['Q%s' % i for i in xrange(len(temp2[0]))])
            raw_pbq2.extend([temp1[i].insert(0, str(i + 1)) for i in xrange(len(temp2))])

            temp2 = qc_dict['#Clean_Base_quality_value_distribution_by_read_position']
            clean_pbq2 = [['Pos']]
            clean_pbq2[0].extend(['Q%s' % i for i in xrange(len(temp2[0]))])
            clean_pbq2.extend([temp1[i].insert(0, str(i + 1)) for i in xrange(len(temp2))])

        # -- Per base- Quality scores Distribution--
        def statis(data_list, *locs):   # Locs Must be sorted before.
            import math

            data_sum = sum(map(int, data_list))
            idx = []
            for i in xrange(len(locs)):
                idx.append(float(locs[i]) * (int(data_sum) + 1))

            def verify_int(num):
                return (False, True)[round(float(num)) == float(num)]

            res = []
            num_count = 0
            for i in xrange(len(data_list)):
                cut_count = 0
                for ii in xrange(len(idx)):
                    if num_count + int(data_list[i]) > int(idx[ii]):
                        res.append(i)
                        cut_count += 1
                    elif num_count + int(data_list[i]) == int(idx[ii]):
                        if verify_int(idx[ii]):
                            res.append(i)
                        else:
                            dec_num = math.modf(idx[ii])[0]
                            for iii in xrange(len(idx)-ii-1):
                                if idx[ii+iii] != 0:
                                    break
                            res.append(int(i+dec_num*(iii-i)))
                        cut_count += 1
                del idx[:cut_count]
                num_count += int(data_list[i])

            return res

        def mean(data_list):
            rsum = qsum = 0
            for i in xrange(len(data_list)):
                qsum += int(data_list[i]) * i
                rsum += data_list[i]
            return round(float(qsum)/rsum, 2)

        temp1 = qc_dict['#Raw_Base_quality_value_distribution_by_read_position']
        temp2 = qc_dict['#Clean_Base_quality_value_distribution_by_read_position']

        fq1_pbqd = [['Position in reads', 'Clean_10%_lower Q', 'Clean_25%_lower Q', 'Clean_Median Q',
                     'Clean_Mean Q', 'Clean_25%_higherQ', 'Clean_10%_higher Q',
                     'Clean_10%_lower Q', 'Clean_25%_lower Q', 'Clean_Median Q',
                     'Clean_Mean Q', 'Clean_25%_higherQ', 'Clean_10%_higher Q']]
        fq1_pbq2030 = [['Position in reads', 'Percentage of Q20+ bases', 'Percentage of Q30+ bases',
                        'Percentage of Clean Q20+', 'Percentage of Clean Q30+']]

        for i in xrange(max(len(temp1), len(temp2))):
            temp_q2030 = [str(i+1)]
            temp_qd = [str(i+1)]

            if i < len(temp1):
                temp_q2030.append(div_perc(sum(map(int, temp1[i][20:])), sum(map(int, temp1[i]))))
                temp_q2030.append(div_perc(sum(map(int, temp1[i][30:])), sum(map(int, temp1[i]))))
                temp_qd.extend(statis(temp1[i], 0.1, 0.25, 0.5, 0.75, 0.9).insert(3, mean(temp1[i])))
            else:
                temp_q2030.extend(['/', '/'])
                temp_qd.extend(['/', '/'])

            if i < len(temp2):
                temp_q2030.append(div_perc(sum(map(int, temp2[i][20:])), sum(map(int, temp2[i]))))
                temp_q2030.append(div_perc(sum(map(int, temp2[i][30:])), sum(map(int, temp2[i]))))
                temp_qd.extend(statis(temp2[i], 0.1, 0.25, 0.5, 0.75, 0.9).insert(3, mean(temp2[i])))
            else:
                temp_q2030.extend(['/', '/'])
                temp_qd.extend(['/', '/'])

            fq1_pbq2030.append(temp_q2030)
            fq1_pbqd.append(temp_qd)

        if pe_mark:
            temp1 = qc_dict['#Fq2_Raw_Base_quality_value_distribution_by_read_position']
            temp2 = qc_dict['#Fq2_Clean_Base_quality_value_distribution_by_read_position']

            fq2_pbqd = [['Position in reads', 'Clean_10%_lower Q', 'Clean_25%_lower Q', 'Clean_Median Q',
                         'Clean_Mean Q', 'Clean_25%_higherQ', 'Clean_10%_higher Q',
                         'Clean_10%_lower Q', 'Clean_25%_lower Q', 'Clean_Median Q',
                         'Clean_Mean Q', 'Clean_25%_higherQ', 'Clean_10%_higher Q']]

            fq2_pbq2030 = [['Position in reads', 'Percentage of Q20+ bases', 'Percentage of Q30+ bases',
                            'Percentage of Clean Q20+', 'Percentage of Clean Q30+']]

            for i in xrange(max(len(temp1), len(temp2))):
                temp_q2030 = [str(i + 1)]
                temp_qd = [str(i + 1)]

                if i < len(temp1):
                    temp_q2030.append(div_perc(sum(map(int, temp1[i][20:])), sum(map(int, temp1[i]))))
                    temp_q2030.append(div_perc(sum(map(int, temp1[i][30:])), sum(map(int, temp1[i]))))
                    temp_qd.extend(statis(temp1[i], 0.1, 0.25, 0.5, 0.75, 0.9).insert(3, mean(temp1[i])))
                else:
                    temp_q2030.extend(['/', '/'])
                    temp_qd.extend(['/', '/'])

                if i < len(temp2):
                    temp_q2030.append(div_perc(sum(map(int, temp2[i][20:])), sum(map(int, temp2[i]))))
                    temp_q2030.append(div_perc(sum(map(int, temp2[i][30:])), sum(map(int, temp2[i]))))
                    temp_qd.extend(statis(temp2[i], 0.1, 0.25, 0.5, 0.75, 0.9).insert(3, mean(temp2[i])))
                else:
                    temp_q2030.extend(['/', '/'])
                    temp_qd.extend(['/', '/'])

                fq2_pbq2030.append(temp_q2030)
                fq2_pbqd.append(temp_qd)

        # Stats Output
        def qc_format(data_list):
            col_list = [0 for i in xrange(max([len(i) for i in data_list]))]
            for line in data_list:
                for i in xrange(len(line)):
                    if isinstance(line[i], list):
                        i_length = len(line[i][0])+len(line[i][1])+1 # min_insert_size of 2 items in one col
                    else:
                        i_length = len(line[i])

                    if i_length > col_list[i]:
                        col_list[i] = i_length

            res = ''
            for line in data_list:
                res_line = ''
                for i in xrange(len(line)):
                    if isinstance(line[i], list):
                        res_line += line[i][0] + (col_list[i]-len(line[i][0])-len(line[i][1]))*' ' + line[i][1] + '\t'
                    else:
                        res_line += eval("'{:<%s}'.format('%s')" % (col_list[i], line[i])) + '\t'
                res += res_line[:-1] + '\n'
            return res

        # 然后Write到文件即可
        with open('%s/Basic_Statistics_of_Sequencing_Quality.txt' % outdir_l, 'w') as F:
            F.write(qc_format(bs))

        with open('%s/Statistics_of_filtered_reads.txt' % outdir_l, 'w') as F:
            F.write(qc_format(flr))

        with open('%s/Base_distribution_by_read_position_1.txt' % outdir_l, 'w') as F:
            F.write(qc_format(pbbd1))

        with open('%s/Base_quality_value_by_read_position_raw_1.txt' % outdir_l, 'w') as F:
            F.write(qc_format(raw_pbq1))

        with open('%s/Base_quality_value_by_read_position_clean_1.txt' % outdir_l, 'w') as F:
            F.write(qc_format(clean_pbq1))

        with open('%s/Base_quality_value_distribution_by_read_position_1.txt' % outdir_l, 'w') as F:
            F.write(qc_format(fq1_pbqd))

        with open('%s/Distribution_of_Q20_Q30_bases_by_read_position_1.txt' % outdir_l, 'w') as F:
            F.write(qc_format(fq1_pbq2030))

        if pe_mark:
            with open('%s/Base_distribution_by_read_position_2.txt' % outdir_l, 'w') as F:
                F.write(qc_format(pbbd2))

            with open('%s/Base_quality_value_by_read_position_raw_2.txt' % outdir_l, 'w') as F:
                F.write(qc_format(raw_pbq2))

            with open('%s/Base_quality_value_by_read_position_clean_2.txt' % outdir_l, 'w') as F:
                F.write(qc_format(clean_pbq2))

            with open('%s/Base_quality_value_distribution_by_read_position_2.txt' % outdir_l, 'w') as F:
                F.write(qc_format(fq2_pbqd))

            with open('%s/Distribution_of_Q20_Q30_bases_by_read_position_2.txt' % outdir_l, 'w') as F:
                F.write(qc_format(fq2_pbq2030))

        if merge_mark:
            if len(os.listdir(indir)) == 1:
                os.system('%s -jar ./FqLoad.jar -i %s -o %s -C %s' % (javapath, outdir_h, outdir_h, 'clean.fastq'))

            elif len(os.listdir(indir)) == 2:
                os.system('%s -jar ./FqLoad.jar -i %s -o %s -C %s -D %s' % (javapath, outdir_h, outdir_h,
                                                                            'clean.fastq1', 'clean.fastq2'))

if __name__ == '__main__':
    a = SOAPnuke()
    a.cmdloop()
