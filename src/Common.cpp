/*
 * Common.cpp
 *
 *  Created on: 2012-6-14
 *      Author: Haosen Chen
 * 		Mail  : chenhaosen@genomics.cn
 */
#include "Common.h"
#include "Logger.h"
#include <sstream>

#define GZ_BUF_SIZE 1048576


namespace PreProcessTool {

	bool isSequence(const string &seq)
	{
		if (seq.empty())
			return false;
		for (int i=seq.size()-1; i>=0; --i)
		{
			switch (seq[i])
			{
				case 'A':
				case 'a':
				case 'C':
				case 'c':
				case 'G':
				case 'g':
				case 'T':
				case 't':
					break;
				default:
					return false;
			}
		}
		return true;
	}

	int adapterType(bool isPE, const string &adapter1, const string &adapter2)
	{
		bool isSeq1 = isSequence(adapter1);
		bool isSeq2 = isSequence(adapter2);

		if (adapter2.empty())
		{
			if (isPE)
			{
				return 3; //error
			}
			if (isSeq1)
			{
				return 1; //adapter sequence
			}
			else
			{
				return 2; //adapter list
			}
		}

		if (isSeq1 && isSeq2)
			return 1;   //sequence
		else if (!isSeq1 && !isSeq2)
			return 2;  //adapter list
		else
			return 3;
	}

	void printFqInfo(const FqInfo *fqInfo)
	{
		cout << fqInfo->rawReadLength << "\t" << fqInfo->cleanReadLength << "\t" << fqInfo->rawTotalReadNum << "\t" << fqInfo->cleanTotalReadNum << "\t" << fqInfo->greyTotalReadNum;
		cout << "\t" << fqInfo->rawTotalBaseNum << "\t" << fqInfo->cleanTotalBaseNum << "\t" << fqInfo->rawBaseA << "\t" << fqInfo->cleanBaseA << "\t" << fqInfo->rawBaseC;
		cout << "\t" << fqInfo->cleanBaseC << "\t" << fqInfo->rawBaseG << "\t" << fqInfo->cleanBaseG << "\t" << fqInfo->rawBaseT << "\t" << fqInfo->cleanBaseT;
		cout << "\t" << fqInfo->rawBaseN << "\t" << fqInfo->cleanBaseN;
		cout << "\t" << fqInfo->rawQ20 << "\t" << fqInfo->cleanQ20 << "\t" << fqInfo->rawQ30 << "\t" << fqInfo->cleanQ30;
		cout << "\t" << fqInfo->duplicationNum << "\t" << fqInfo->adapterNum << "\t" << fqInfo->nExceedNum << "\t" << fqInfo->lowQualNum;
		cout << "\t" << fqInfo->lowMeanNum << "\t" << fqInfo->smallInsertNum << "\t" << fqInfo->polyANum;
		cout << "\t" << fqInfo->totalDuplicationNum << "\t" << fqInfo->totalAdapterNum << "\t" << fqInfo->totalNExceedNum << "\t" << fqInfo->totalLowQualNum;
		cout << "\t" << fqInfo->totalSmallInsertNum << "\t" << fqInfo->totalPolyANum << "\t" << fqInfo->totalCutAdaptorNum;
		cout << "\t" << fqInfo->maxQualityValue << "\t" << fqInfo->shortNum << "\t" << fqInfo->totalShortNum;
		cout << "\t" << fqInfo->polyXNum << "\t" << fqInfo->totalPolyXNum << "#S" << endl;  //这里的shortNum插入位置有点不规范，留待修改

		//base distributions by read position
		cout << "#Base_distributions_by_read_position\t#S" << endl;
		for (unsigned int i=0; i<fqInfo->rawReadLength; ++i)
		{
			for (int k=0; k<5; k++)
			{
				cout << fqInfo->base[i][k] << "\t" ;
			}

			cout << fqInfo->clean_base[i][0];
			for (int k=1; k<5; k++)
			{
				cout << "\t" << fqInfo->clean_base[i][k];
			}
			cout << "#S\n";
		}

		//Raw Base_quality_value_distribution_by_read_position && Distribution_of_Q20_Q30_bases_by_read_position
		cout << "#Raw_Base_quality_value_distribution_by_read_position\t#S" << endl;
		for (unsigned int i=0; i<fqInfo->rawReadLength; ++i)
		{
			cout << fqInfo->q20q30[i][0] << "\t" << fqInfo->q20q30[i][1];
			for (int k=0; k<=fqInfo->maxQualityValue; k++)
			{
				cout << "\t" << fqInfo->qual[i][k];
			}
			cout << "#S\n";
		}

		//Clean Base_quality_value_distribution_by_read_position && Distribution_of_Q20_Q30_bases_by_read_position
		cout << "#Clean_Base_quality_value_distribution_by_read_position\t#S" << endl;
		for (unsigned int i=0; i<fqInfo->cleanReadLength; ++i)
		{
			cout << fqInfo->clean_q20q30[i][0] << "\t" << fqInfo->clean_q20q30[i][1];
			for (int k=0; k<=fqInfo->maxQualityValue; k++)
			{
				cout << "\t" << fqInfo->clean_qual[i][k];
			}
			cout << "#S\n";
		}


	}

	void printFqInfo(const FqInfo *fqInfo1, const FqInfo *fqInfo2){
		cout << "#Fq1_statistical_information\t#S" << endl;
		printFqInfo(fqInfo1);
		if (fqInfo2 != NULL){
			cout << "#Fq2_statistical_information\t#S" << endl;
			printFqInfo(fqInfo2);
		}
	}

	void printFqInfo(const string outDir, const string prefix, const FqInfo *fqInfo1, const FqInfo *fqInfo2)
	{
		string outFile = outDir + "/" + prefix  + SEQUENCING_QUALITY + ".txt";
		ofstream out(outFile.c_str());
		if (!out)
		{
			LOG(ERROR, "open output file: " << outFile);
			exit(1);
		}

		// set style;
		int setwTmp = 1;
		long baseTmp = fqInfo1->rawTotalBaseNum;
		do
		{
			setwTmp++;
			baseTmp /= 10;
		}while(baseTmp > 0);
		if (setwTmp < 3)
			setwTmp = 3;
		int setwTmp2 = setwTmp + 8;
		int setwTmp3 = 5;
		int precisionTmp = 2;
		//

		out << fixed;
		out << "Item                                                              \t"
			<< std::left << setw(setwTmp2) << "raw reads(fq1)" << "\t" << std::left << setw(setwTmp2) << "clean reads(fq1)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp2) << "raw reads(fq2)" << "\t" << std::left << setw(setwTmp2) << "clean reads(fq2)";
		}
		out << "\n";

		out << "Read length                                                       \t"
			<< std::left << setw(setwTmp2) << fqInfo1->rawReadLength << "\t" <<std::left << setw(setwTmp2) << fqInfo1->cleanReadLength;
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp2) << fqInfo2->rawReadLength << "\t" <<std::left << setw(setwTmp2) << fqInfo2->cleanReadLength;
		}
		out << "\n";

		out << "Total number of reads                                             \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawTotalReadNum << "(100.00%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanTotalReadNum << "(100.00%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawTotalReadNum << "(100.00%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanTotalReadNum << "(100.00%)";
		}
		out << "\n";

		out << "Number of filtered reads (%)                                     \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawTotalReadNum - fqInfo1->cleanTotalReadNum
			<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * (fqInfo1->rawTotalReadNum - fqInfo1->cleanTotalReadNum) / fqInfo1->rawTotalReadNum << "%)"
			<< "\t" << std::left << setw(setwTmp2) << "-";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawTotalReadNum - fqInfo2->cleanTotalReadNum
				<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * (fqInfo2->rawTotalReadNum - fqInfo2->cleanTotalReadNum) / fqInfo2->rawTotalReadNum << "%)"
				<< "\t" << std::left << setw(setwTmp2) << "-";
		}
		out << "\n";

		out << "\n";

		out << "Total number of bases                                            \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawTotalBaseNum << "(100.00%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanTotalBaseNum << "(100.00%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawTotalBaseNum << "(100.00%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanTotalBaseNum << "(100.00%)";
		}
		out << "\n";

		out << "Number of filtered bases (%)                                     \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawTotalBaseNum - fqInfo1->cleanTotalBaseNum
			<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * (fqInfo1->rawTotalBaseNum - fqInfo1->cleanTotalBaseNum) / fqInfo1->rawTotalBaseNum << "%)"
			<< "\t" << std::left << setw(setwTmp2) << "-";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawTotalBaseNum - fqInfo2->cleanTotalBaseNum
				<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * (fqInfo2->rawTotalBaseNum - fqInfo2->cleanTotalBaseNum) / fqInfo2->rawTotalBaseNum << "%)"
				<< "\t" << std::left << setw(setwTmp2) << "-";
		}
		out << "\n";

		out << "\n";

		out << "Reads related to Adapter and Trimmed (%)                           \t"
			<< std::left << setw(setwTmp) << fqInfo1->totalCutAdaptorNum << "(" <<std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->totalCutAdaptorNum / fqInfo1->rawTotalReadNum << "%)" << "\t" << std::left << setw(setwTmp2) << "-";
		if(fqInfo2 != NULL){
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->totalCutAdaptorNum << "(" <<std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo2->totalCutAdaptorNum / fqInfo2->rawTotalReadNum << "%)" << "\t" << std::left << setw(setwTmp2) << "-";
		}
		out << "\n";

		out<< "\n";

		out << "Number of base A (%)                                             \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawBaseA << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawBaseA / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanBaseA << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanBaseA / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawBaseA << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawBaseA / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanBaseA << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanBaseA / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out << "Number of base C (%)                                             \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawBaseC << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawBaseC / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanBaseC << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanBaseC / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawBaseC << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawBaseC / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanBaseC << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanBaseC / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out << "Number of base G (%)                                             \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawBaseG<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawBaseG / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanBaseG << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanBaseG / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawBaseG<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawBaseG / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanBaseG << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanBaseG / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out << "Number of base T (%)                                             \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawBaseT << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawBaseT / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanBaseT << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanBaseT / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawBaseT << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawBaseT / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanBaseT << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanBaseT / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out << "Number of base N (%)                                             \t"
			<< std::left << setw(setwTmp) << fqInfo1->rawBaseN << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawBaseN / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanBaseN << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanBaseN / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawBaseN << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawBaseN / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanBaseN << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanBaseN / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out << "\n";

		out << "Number of base calls with quality value of 20 or higher (Q20+) (%)\t"
			<< std::left << setw(setwTmp) << fqInfo1->rawQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawQ20 / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanQ20 / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawQ20 / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanQ20 / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out << "Number of base calls with quality value of 30 or higher (Q30+) (%)\t"
			<< std::left << setw(setwTmp) << fqInfo1->rawQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->rawQ30 / fqInfo1->rawTotalBaseNum << "%)" << "\t"
			<< std::left << setw(setwTmp) << fqInfo1->cleanQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo1->cleanQ30 / fqInfo1->cleanTotalBaseNum << "%)";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->rawQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->rawQ30 / fqInfo2->rawTotalBaseNum << "%)" << "\t"
				<< std::left << setw(setwTmp) << fqInfo2->cleanQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
				<< 100.0 * fqInfo2->cleanQ30 / fqInfo2->cleanTotalBaseNum << "%)";
		}
		out << "\n";

		out.close();
		out.clear();


		outFile = outDir + "/" + prefix + FILTERED_READS + ".txt";

		out.open(outFile.c_str());
		if (!out)
		{
			LOG(ERROR, "open output file: " << outFile);
			exit(1);
		}


		int long filteredReads = fqInfo1->rawTotalReadNum - fqInfo1->cleanTotalReadNum;
		baseTmp = filteredReads;
		setwTmp = 1;
		do
		{
			setwTmp++;
			baseTmp /= 10;
		}while(baseTmp > 0);

		if (setwTmp < 6)
			setwTmp = 6;
		out << fixed;
		out << "Item                                    \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << "Total" << "\t" << std::right << setw(setwTmp) << "Percentage" << "\t";
		}
		out << std::left << setw(setwTmp) << "Counts(fq1)" << "\t" << std::right << setw(setwTmp) << "Percentage";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << "Counts(fq2)" << "\t" << std::right << setw(setwTmp) << "Percentage";
		}
		out << "\n";

		out << "Total filtered reads (%)                \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << filteredReads * 2 << "\t" << std::right << setw(setwTmp) << "100.00" << "%\t";
		}
		out << std::left << setw(setwTmp) << filteredReads << "\t" << std::right << setw(setwTmp) << "100.00" << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << filteredReads << "\t" << std::right << setw(setwTmp) << "100.00" << "%";
		}
		out << "\n";

		out << "Reads too short (%)                     \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalShortNum*2  << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalShortNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->shortNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->shortNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->shortNum  << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->shortNum / filteredReads << "%";
		}
		out << "\n";

		out << "Reads with adapter (%)                  \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalAdapterNum*2  << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalAdapterNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->adapterNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->adapterNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->adapterNum  << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->adapterNum / filteredReads << "%";
		}
		out << "\n";

		out << "Reads with low quality (%)               \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalLowQualNum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalLowQualNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->lowQualNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->lowQualNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->lowQualNum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->lowQualNum / filteredReads << "%";
		}
		out << "\n";

		out << "Reads with low mean quality (%)         \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalLowMeanNum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalLowMeanNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->lowMeanNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->lowMeanNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->lowMeanNum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->lowMeanNum / filteredReads << "%";
		}
		out << "\n";

		out << "Reads with duplications (%)             \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalDuplicationNum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalDuplicationNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->duplicationNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->duplicationNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->duplicationNum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->duplicationNum / filteredReads << "%";
		}
		out << "\n";

		out << "Read with n rate exceed: (%)            \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalNExceedNum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalNExceedNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->nExceedNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->nExceedNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->nExceedNum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->nExceedNum / filteredReads << "%";
		}
		out << "\n";

		out << "Read with small insert size: (%)        \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalSmallInsertNum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalSmallInsertNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->smallInsertNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->smallInsertNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->smallInsertNum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->smallInsertNum / filteredReads << "%";
		}
		out << "\n";

		out << "Reads with PolyA (%)                    \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalPolyANum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalPolyANum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->polyANum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->polyANum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->polyANum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->polyANum / filteredReads << "%";
		}
		out << "\n";

		out << "Reads with PolyX (%)                    \t";
		if (fqInfo2 != NULL)
		{
			out << std::left << setw(setwTmp) << fqInfo1->totalPolyXNum * 2 << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo1->totalPolyXNum / filteredReads << "%\t";
		}
		out << std::left << setw(setwTmp) << fqInfo1->polyXNum << "\t" << std::right << setw(setwTmp)
			<< setprecision(precisionTmp) << 100.0 * fqInfo1->polyXNum / filteredReads << "%";
		if (fqInfo2 != NULL)
		{
			out << "\t" << std::left << setw(setwTmp) << fqInfo2->polyXNum << "\t" << std::right << setw(setwTmp)
				<< setprecision(precisionTmp) << 100.0 * fqInfo2->polyXNum / filteredReads << "%";
		}
		out << "\n";

		out.close();
		out.clear();

		outFile = outDir + "/" + prefix + BASE_DISTRIBUTIONS + "_1.txt";
		out.open(outFile.c_str());
		if (!out)
		{
			LOG(ERROR, "open output file: " << outFile);
			exit(1);
		}

		setwTmp = 6;
		out << std::left << setw(setwTmp) << "Pos" << "\t" << std::right << setw(setwTmp) << "A" << "\t"
			<< setw(setwTmp) << "C" << "\t" <<  setw(setwTmp) << "G" << "\t"
			<< setw(setwTmp) << "T" <<  setw(setwTmp) << "N" << "\t"
			<< setw(setwTmp) << "Clean A" << "\t"	<< setw(setwTmp) << "Clean C" << "\t"
			<< setw(setwTmp) << "Clean G" << "\t" << setw(setwTmp) << "Clean T" << "\t"
			<< setw(setwTmp) << "Clean N" << "\n";

		setwTmp = 5;

		for (unsigned int i=0; i<fqInfo1->rawReadLength; ++i)
		{
			out << std::left << setw(setwTmp) << i+1 << "\t" << std::right
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->base[i][0] / fqInfo1->rawTotalReadNum << "%\t"
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->base[i][1] / fqInfo1->rawTotalReadNum << "%\t"
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->base[i][2] / fqInfo1->rawTotalReadNum << "%\t"
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->base[i][3] / fqInfo1->rawTotalReadNum << "%\t"
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->base[i][4] / fqInfo1->rawTotalReadNum << "%\t";
			if (i < fqInfo1->cleanReadLength)
			{
				out << std::right
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_base[i][0] / fqInfo1->cleanTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_base[i][1] / fqInfo1->cleanTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_base[i][2] / fqInfo1->cleanTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_base[i][3] / fqInfo1->cleanTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_base[i][4] / fqInfo1->cleanTotalReadNum << "%\n";
			}
			else
			{
				out << "\n";
			}
		}
		out.close();
		out.clear();

		outFile = outDir + "/" + prefix + DISTRIBUTION_OF_Q20_Q30 + "_1.txt";
		out.open(outFile.c_str());
		if (!out)
		{
			LOG(ERROR, "open output file: " << outFile);
			exit(1);
		}

		out <<"Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\tPercentage of Clean Q20+\tPercentage of Clean Q30+\n";
		for (unsigned int i=0; i<fqInfo1->rawReadLength; i++)
		{
			out << std::left << setw(setwTmp) << i+1 << "\t" << std::right
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->q20q30[i][0] / fqInfo1->rawTotalReadNum << "%\t"
				<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->q20q30[i][1] / fqInfo1->rawTotalReadNum << "%\t";

			if (i < fqInfo1->cleanReadLength)
			{
				out << std::right
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_q20q30[i][0] / fqInfo1->cleanTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo1->clean_q20q30[i][1] / fqInfo1->cleanTotalReadNum << "%\n";
			}
			else
			{
				out << "\n";
			}
		}
		out.close();
		out.clear();

		outFile = outDir + "/" + prefix + BASE_QUALITY_VALUE_DISTRIBUTION + "_1.txt";
		out.open(outFile.c_str());
		if (!out)
		{
			LOG(ERROR, "open output file: " << outFile);
			exit(1);
		}
		baseTmp = fqInfo1->rawTotalReadNum;
		setwTmp = 1;
		do
		{
			setwTmp++;
			baseTmp /= 10;
		}while(baseTmp > 0);
		if (setwTmp < 3)
			setwTmp = 3;


		float qBox[MAX_LENGTH][6];
		memset(qBox, 0, sizeof(float) * 6 * MAX_LENGTH);
		unsigned long qBoxPos[2][5];//median lower upper 10 90;
		memset(qBoxPos, 0, sizeof(long) * 5 * 2);
		if (fqInfo1->rawTotalReadNum != 0)
		{
			if (fqInfo1->rawTotalReadNum % 2 == 0)
			{
				qBoxPos[0][0] = fqInfo1->rawTotalReadNum / 2;
				if (qBoxPos[0][0] != 0)
					qBoxPos[1][0] = qBoxPos[0][0] + 1;
			}
			else
			{
				qBoxPos[0][0] = (fqInfo1->rawTotalReadNum + 1) / 2;
				qBoxPos[1][0] = qBoxPos[0][0];
			}

			if (qBoxPos[0][0] % 2 == 0)
			{
				qBoxPos[0][1] = qBoxPos[0][0] / 2;
				qBoxPos[1][1] = qBoxPos[0][1] + 1;
			}
			else
			{
				qBoxPos[0][1] = (qBoxPos[0][0] + 1) / 2;
				qBoxPos[1][1] = qBoxPos[0][1];
			}

			if ((fqInfo1->rawTotalReadNum - qBoxPos[1][0] + 1) % 2 == 0)
			{
				qBoxPos[0][2] = (fqInfo1->rawTotalReadNum - qBoxPos[1][0] + 1) / 2 + qBoxPos[1][0] - 1;
				qBoxPos[1][2] = qBoxPos[0][2] + 1;
			}else{
				qBoxPos[0][2] = (fqInfo1->rawTotalReadNum - qBoxPos[1][0]) / 2 + qBoxPos[1][0];
				qBoxPos[1][2] = qBoxPos[0][2];
			}

			qBoxPos[0][3] = (int)(fqInfo1->rawTotalReadNum * 0.1 + 0.5);
			qBoxPos[1][3] = qBoxPos[0][3];

			qBoxPos[0][4] = (int)(fqInfo1->rawTotalReadNum * 0.9 + 0.5);
			qBoxPos[1][4] = qBoxPos[0][4];
		}

		unsigned long qBoxTmp[2][6];
		for (unsigned int i=0; i<fqInfo1->rawReadLength; ++i)
		{
			memset(qBoxTmp, 0, sizeof(long) * 2 * 6);
			for (int j=0; j<=fqInfo1->maxQualityValue; j++)
			{
				qBoxTmp[1][5] += fqInfo1->qual[i][j] * j;
				qBoxTmp[0][5] += fqInfo1->qual[i][j];
				for (int k=0; k<5; k++)
				{
					if(qBoxTmp[0][5] < qBoxPos[0][k])
						qBoxTmp[0][k] = j;
					if(qBoxTmp[0][5] < qBoxPos[1][k])
						qBoxTmp[1][k] = j;
				}
			}

			qBox[i][0] = 1.0 * qBoxTmp[1][5] / fqInfo1->rawTotalReadNum;
			for (int k=0; k<5; k++)
			{
				qBox[i][k+1] = (qBoxTmp[0][k] + qBoxTmp[1][k]) / 2 + 1;
			}
		}

		out << std::left << setw(setwTmp) << "Pos";
		for (int i=0; i<=fqInfo1->maxQualityValue; i++)
		{
			out << "\tQ" << setw(setwTmp-1) << i;
		}
		out << "\tMean\tMedian\tLower quartile\tUpper quartile\t10thpercentile\t90thpercentile";
		out << "\n";

		for (unsigned int i=0; i<fqInfo1->rawReadLength; ++i)
		{
			out << setw(setwTmp) << (i+1);
			for (int j=0; j<=fqInfo1->maxQualityValue; j++)
			{
				out << "\t" << setw(setwTmp) << fqInfo1->qual[i][j];
			}
			for (int j=0; j<6; j++)
			{
				out << "\t" << setw(setwTmp) << setprecision(precisionTmp) << qBox[i][j];
			}
			out << "\n";
		}

		/////////////////////////////////////
		//clean base quality distribute
		memset(qBox, 0, sizeof(float) * 6 * MAX_LENGTH);
		memset(qBoxPos, 0, sizeof(long) * 5 * 2);
		if (fqInfo1->cleanTotalReadNum != 0)
		{
			if (fqInfo1->cleanTotalReadNum % 2 == 0)
			{
				qBoxPos[0][0] = fqInfo1->cleanTotalReadNum / 2;
				if (qBoxPos[0][0] != 0)
					qBoxPos[1][0] = qBoxPos[0][0] + 1;
			}
			else
			{
				qBoxPos[0][0] = (fqInfo1->cleanTotalReadNum + 1) / 2;
				qBoxPos[1][0] = qBoxPos[0][0];
			}

			if (qBoxPos[0][0] % 2 == 0)
			{
				qBoxPos[0][1] = qBoxPos[0][0] / 2;
				qBoxPos[1][1] = qBoxPos[0][1] + 1;
			}
			else
			{
				qBoxPos[0][1] = (qBoxPos[0][0] + 1) / 2;
				qBoxPos[1][1] = qBoxPos[0][1];
			}

			if ((fqInfo1->cleanTotalReadNum - qBoxPos[1][0] + 1) % 2 == 0)
			{
				qBoxPos[0][2] = (fqInfo1->cleanTotalReadNum - qBoxPos[1][0] + 1) / 2 + qBoxPos[1][0] - 1;
				qBoxPos[1][2] = qBoxPos[0][2] + 1;
			}else{
				qBoxPos[0][2] = (fqInfo1->cleanTotalReadNum - qBoxPos[1][0]) / 2 + qBoxPos[1][0];
				qBoxPos[1][2] = qBoxPos[0][2];
			}

			qBoxPos[0][3] = (int)(fqInfo1->cleanTotalReadNum * 0.1 + 0.5);
			qBoxPos[1][3] = qBoxPos[0][3];

			qBoxPos[0][4] = (int)(fqInfo1->cleanTotalReadNum * 0.9 + 0.5);
			qBoxPos[1][4] = qBoxPos[0][4];
		}

		for (unsigned int i=0; i<fqInfo1->cleanReadLength; ++i)
		{
			memset(qBoxTmp, 0, sizeof(long) * 2 * 6);
			for (int j=0; j<=fqInfo1->maxQualityValue; j++)
			{
				qBoxTmp[1][5] += fqInfo1->clean_qual[i][j] * j;
				qBoxTmp[0][5] += fqInfo1->clean_qual[i][j];
				for (int k=0; k<5; k++)
				{
					if(qBoxTmp[0][5] < qBoxPos[0][k])
						qBoxTmp[0][k] = j;
					if(qBoxTmp[0][5] < qBoxPos[1][k])
						qBoxTmp[1][k] = j;
				}
			}

			qBox[i][0] = 1.0 * qBoxTmp[1][5] / fqInfo1->cleanTotalReadNum;
			for (int k=0; k<5; k++)
			{
				qBox[i][k+1] = (qBoxTmp[0][k] + qBoxTmp[1][k]) / 2 + 1;
			}
		}

		out << "Clean Quality Value Distribute\n";
		out << std::left << setw(setwTmp) << "Pos";
		for (int i=0; i<=fqInfo1->maxQualityValue; i++)
		{
			out << "\tQ" << setw(setwTmp-1) << i;
		}

		out << "\tMean\tMedian\tLower quartile\tUpper quartile\t10thpercentile\t90thpercentile";
		out << "\n";

		for (unsigned int i=0; i<fqInfo1->cleanReadLength; ++i)
		{
			out << setw(setwTmp) << (i+1);
			for (int j=0; j<=fqInfo1->maxQualityValue; j++)
			{
				out << "\t" << setw(setwTmp) << fqInfo1->clean_qual[i][j];
			}
			for (int j=0; j<6; j++)
			{
				out << "\t" << setw(setwTmp) << setprecision(precisionTmp) << qBox[i][j];
			}
			out << "\n";
		}
		out.close();
		out.clear();

		//output fq2 info
		if (fqInfo2 != NULL)
		{
			outFile = outDir + "/" + prefix + BASE_DISTRIBUTIONS + "_2.txt";
			out.open(outFile.c_str());
			if (!out)
			{
				LOG(ERROR, "open output file: " << outDir + "/" + BASE_DISTRIBUTIONS + ".txt");
				exit(1);
			}

			setwTmp = 6;
			out << std::left << setw(setwTmp) << "Pos" << "\t" << std::right << setw(setwTmp) << "A" << "\t"
				<< setw(setwTmp) << "C" << "\t" <<  setw(setwTmp) << "G" << "\t"
				<< setw(setwTmp) << "T" <<  "\t" << setw(setwTmp) << "N" << "\t"
				<< setw(setwTmp) << "Clean A" << "\t"
				<< setw(setwTmp) << "Clean C" << "\t" << setw(setwTmp) << "Clean G" << "\t"
				<< setw(setwTmp) << "Clean T" << "\t" << setw(setwTmp) << "Clean N" << "\n";


			setwTmp = 5;

			for (unsigned int i=0; i<fqInfo2->rawReadLength; ++i)
			{
				out << std::left << setw(setwTmp) << i+1 << "\t" << std::right
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->base[i][0] / fqInfo2->rawTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->base[i][1] / fqInfo2->rawTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->base[i][2] / fqInfo2->rawTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->base[i][3] / fqInfo2->rawTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->base[i][4] / fqInfo2->rawTotalReadNum << "%\t";
				if (i < fqInfo2->cleanReadLength)
				{
					out << std::right
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_base[i][0] / fqInfo2->cleanTotalReadNum << "%\t"
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_base[i][1] / fqInfo2->cleanTotalReadNum << "%\t"
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_base[i][2] / fqInfo2->cleanTotalReadNum << "%\t"
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_base[i][3] / fqInfo2->cleanTotalReadNum << "%\t"
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_base[i][4] / fqInfo2->cleanTotalReadNum << "%\n";
				}
				else
				{
					out << "\n";
				}
			}
			out.close();
			out.clear();

			outFile = outDir + "/" + prefix + DISTRIBUTION_OF_Q20_Q30 + "_2.txt";
			out.open(outFile.c_str());
			if (!out)
			{
				LOG(ERROR, "open output file: " << outFile);
				exit(1);
			}

			out <<"Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\tPercentage of Clean Q20\tPercentage of Clean Q30\n";
			for (unsigned int i=0; i<fqInfo2->rawReadLength; i++)
			{
				out << std::left << setw(setwTmp) << i+1 << "\t" << std::right
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->q20q30[i][0] / fqInfo2->rawTotalReadNum << "%\t"
					<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->q20q30[i][1] / fqInfo2->rawTotalReadNum << "%\t";
				if (i < fqInfo2->cleanReadLength)
				{
					out << std::right
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_q20q30[i][0] / fqInfo2->cleanTotalReadNum << "%\t"
						<< setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo2->clean_q20q30[i][1] / fqInfo2->cleanTotalReadNum << "%\n";
				}
				else
				{
					out << "\n";
				}
			}
			out.close();
			out.clear();

			outFile = outDir + "/" + prefix + BASE_QUALITY_VALUE_DISTRIBUTION + "_2.txt";
			out.open(outFile.c_str());
			if (!out)
			{
				LOG(ERROR, "open output file: " << outFile);
				exit(1);
			}
			baseTmp = fqInfo2->rawTotalReadNum;
			setwTmp = 1;
			do
			{
				setwTmp++;
				baseTmp /= 10;
			}while(baseTmp > 0);
			if (setwTmp < 3)
				setwTmp = 3;


			float qBox[MAX_LENGTH][6];
			memset(qBox, 0, sizeof(float) * 6 * MAX_LENGTH);
			unsigned long qBoxPos[2][5];//median lower upper 10 90;
			memset(qBoxPos, 0, sizeof(long) * 5 * 2);
			if (fqInfo2->rawTotalReadNum != 0)
			{
				if (fqInfo2->rawTotalReadNum % 2 == 0)
				{
					qBoxPos[0][0] = fqInfo2->rawTotalReadNum / 2;
					if (qBoxPos[0][0] != 0)
						qBoxPos[1][0] = qBoxPos[0][0] + 1;
				}
				else
				{
					qBoxPos[0][0] = (fqInfo2->rawTotalReadNum + 1) / 2;
					qBoxPos[1][0] = qBoxPos[0][0];
				}

				if (qBoxPos[0][0] % 2 == 0)
				{
					qBoxPos[0][1] = qBoxPos[0][0] / 2;
					qBoxPos[1][1] = qBoxPos[0][1] + 1;
				}
				else
				{
					qBoxPos[0][1] = (qBoxPos[0][0] + 1) / 2;
					qBoxPos[1][1] = qBoxPos[0][1];
				}

				if ((fqInfo2->rawTotalReadNum - qBoxPos[1][0] + 1) % 2 == 0)
				{
					qBoxPos[0][2] = (fqInfo2->rawTotalReadNum - qBoxPos[1][0] + 1) / 2 + qBoxPos[1][0] - 1;
					qBoxPos[1][2] = qBoxPos[0][2] + 1;
				}else{
					qBoxPos[0][2] = (fqInfo2->rawTotalReadNum - qBoxPos[1][0]) / 2 + qBoxPos[1][0];
					qBoxPos[1][2] = qBoxPos[0][2];
				}

				qBoxPos[0][3] = (int)(fqInfo2->rawTotalReadNum * 0.1 + 0.5);
				qBoxPos[1][3] = qBoxPos[0][3];

				qBoxPos[0][4] = (int)(fqInfo2->rawTotalReadNum * 0.9 + 0.5);
				qBoxPos[1][4] = qBoxPos[0][4];
			}

			unsigned long qBoxTmp[2][6];
			for (unsigned int i=0; i<fqInfo2->rawReadLength; ++i)
			{
				memset(qBoxTmp, 0, sizeof(long) * 2 * 6);
				for (int j=0; j<=fqInfo2->maxQualityValue; j++)
				{
					qBoxTmp[1][5] += fqInfo2->qual[i][j] * j;
					qBoxTmp[0][5] += fqInfo2->qual[i][j];
					for (int k=0; k<5; k++)
					{
						if(qBoxTmp[0][5] < qBoxPos[0][k])
							qBoxTmp[0][k] = j;
						if(qBoxTmp[0][5] < qBoxPos[1][k])
							qBoxTmp[1][k] = j;
					}
				}

				qBox[i][0] = 1.0 * qBoxTmp[1][5] / fqInfo2->rawTotalReadNum;
				for (int k=0; k<5; k++)
				{
					qBox[i][k+1] = (qBoxTmp[0][k] + qBoxTmp[1][k]) / 2 + 1;
				}
			}

			out << std::left << setw(setwTmp) << "Pos";
			for (int i=0; i<=fqInfo2->maxQualityValue; i++)
			{
				out << "\tQ" << setw(setwTmp-1) << i;
			}
			out << "\tMean\tMedian\tLower quartile\tUpper quartile\t10thpercentile\t90thpercentile";
			out << "\n";

			for (unsigned int i=0; i<fqInfo2->rawReadLength; ++i)
			{
				out << setw(setwTmp) << (i+1);
				for (int j=0; j<=fqInfo2->maxQualityValue; j++)
				{
					out << "\t" << setw(setwTmp) << fqInfo2->qual[i][j];
				}
				for (int j=0; j<6; j++)
				{
					out << "\t" << setw(setwTmp) << setprecision(precisionTmp) << qBox[i][j];
				}
				out << "\n";
			}

			/////////////////////////////////
			//clean base quality value distribute
			memset(qBox, 0, sizeof(float) * 6 * MAX_LENGTH);
			memset(qBoxPos, 0, sizeof(long) * 5 * 2);
			if (fqInfo2->cleanTotalReadNum != 0)
			{
				if (fqInfo2->cleanTotalReadNum % 2 == 0)
				{
					qBoxPos[0][0] = fqInfo2->cleanTotalReadNum / 2;
					if (qBoxPos[0][0] != 0)
						qBoxPos[1][0] = qBoxPos[0][0] + 1;
				}
				else
				{
					qBoxPos[0][0] = (fqInfo2->cleanTotalReadNum + 1) / 2;
					qBoxPos[1][0] = qBoxPos[0][0];
				}

				if (qBoxPos[0][0] % 2 == 0)
				{
					qBoxPos[0][1] = qBoxPos[0][0] / 2;
					qBoxPos[1][1] = qBoxPos[0][1] + 1;
				}
				else
				{
					qBoxPos[0][1] = (qBoxPos[0][0] + 1) / 2;
					qBoxPos[1][1] = qBoxPos[0][1];
				}

				if ((fqInfo2->cleanTotalReadNum - qBoxPos[1][0] + 1) % 2 == 0)
				{
					qBoxPos[0][2] = (fqInfo2->cleanTotalReadNum - qBoxPos[1][0] + 1) / 2 + qBoxPos[1][0] - 1;
					qBoxPos[1][2] = qBoxPos[0][2] + 1;
				}else{
					qBoxPos[0][2] = (fqInfo2->cleanTotalReadNum - qBoxPos[1][0]) / 2 + qBoxPos[1][0];
					qBoxPos[1][2] = qBoxPos[0][2];
				}

				qBoxPos[0][3] = (int)(fqInfo2->cleanTotalReadNum * 0.1 + 0.5);
				qBoxPos[1][3] = qBoxPos[0][3];

				qBoxPos[0][4] = (int)(fqInfo2->cleanTotalReadNum * 0.9 + 0.5);
				qBoxPos[1][4] = qBoxPos[0][4];
			}

			for (unsigned int i=0; i<fqInfo2->cleanReadLength; ++i)
			{
				memset(qBoxTmp, 0, sizeof(long) * 2 * 6);
				for (int j=0; j<=fqInfo2->maxQualityValue; j++)
				{
					qBoxTmp[1][5] += fqInfo2->clean_qual[i][j] * j;
					qBoxTmp[0][5] += fqInfo2->clean_qual[i][j];
					for (int k=0; k<5; k++)
					{
						if(qBoxTmp[0][5] < qBoxPos[0][k])
							qBoxTmp[0][k] = j;
						if(qBoxTmp[0][5] < qBoxPos[1][k])
							qBoxTmp[1][k] = j;
					}
				}

				qBox[i][0] = 1.0 * qBoxTmp[1][5] / fqInfo2->cleanTotalReadNum;
				for (int k=0; k<5; k++)
				{
					qBox[i][k+1] = (qBoxTmp[0][k] + qBoxTmp[1][k]) / 2 + 1;
				}
			}

			out << "Clean Quality Value Distribute\n";
			out << std::left << setw(setwTmp) << "Pos";
			for (int i=0; i<=fqInfo2->maxQualityValue; i++)
			{
				out << "\tQ" << setw(setwTmp-1) << i;
			}
			out << "\tMean\tMedian\tLower quartile\tUpper quartile\t10thpercentile\t90thpercentile";
			out << "\n";

			for (unsigned int i=0; i<fqInfo2->cleanReadLength; ++i)
			{
				out << setw(setwTmp) << (i+1);
				for (int j=0; j<=fqInfo2->maxQualityValue; j++)
				{
					out << "\t" << setw(setwTmp) << fqInfo2->clean_qual[i][j];
				}
				for (int j=0; j<6; j++)
				{
					out << "\t" << setw(setwTmp) << setprecision(precisionTmp) << qBox[i][j];
				}
				out << "\n";
			}
			out.close();
			out.clear();
		}

	}

	string getOutputFileName(string filename, string prefix, string path)
	{
		string outFileName;
		if (filename.substr(filename.size() - 2, 2) != "gz")
		{
			int pos = filename.rfind('/');
			if (pos == -1)
			{
				outFileName = path + "/" + prefix + filename + ".gz";
			}
			else
			{
				outFileName = path + "/" + prefix
					+ filename.substr(pos + 1) + ".gz";
			}
		}
		else
		{
			int pos = filename.rfind('/');
			if (pos == -1)
			{
				outFileName = path + "/" + prefix + filename;
			}
			else
			{
				outFileName = path + "/" + prefix
					+ filename.substr(pos + 1);
			}
		}

		return outFileName;
	}

	//对字符串进行压缩
	void encode(string &seq, ReadSeq &base)
	{
		int x = 0;
		for (int m = 0; m < base.size; ++m)
		{
			x = m + 1;
			for (int i = (m<<5); i < (x<<5); ++i)
			{
				int k = ((seq[i]&0x06)>>1);
				base.readName[m] = base.readName[m] << 2;
				base.readName[m] += k;
			}
		}
		for (unsigned int i= (base.size<<5); i<seq.size(); ++i)
		{
			int k = ((seq[i] & 0x06)>>1);
			base.readName[base.size] = base.readName[base.size]<<2;
			base.readName[base.size] += k;
		}
	}

	string reverseComplementary(string &seq)
	{
		string reverse(seq);
		int len = reverse.size();
		char temp;
		for (int i=0; i<len/2; ++i)
		{
			temp = reverse[i];
			reverse[i] = reverse[len - 1 - i];
			reverse[len - 1 - i] = temp;
		}

		for (int i=reverse.size()-1; i>=0; --i)
		{
			switch (reverse[i])
			{
				case 'A':
				case 'a':
					reverse[i] = 'T';
					break;
				case 'C':
				case 'c':
					reverse[i] = 'G';
					break;
				case 'T':
				case 't':
					reverse[i] = 'A';
					break;
				case 'G':
				case 'g':
					reverse[i] = 'C';
					break;
			}
		}

		return reverse;
	}

	void turnBase(char* str, char from, char to)
	{
		while (*str)
		{
			if (*str == from)
			{
				*str = to;
			}
			++str;
		}
	}

	void maskLowQualBase(char* qual, char* seq, char qualityThreshold)
	{
		while (*qual)
		{
			if (*qual <= qualityThreshold)
			{
				*seq = 'N';
			}
			++qual;
			++seq;
		}
	}

	void upper(char* str)
	{
		while (*str)
		{
			if (*str >= 'a' && *str <= 'z')
			{
				*str -= 32; //'a' - 'A'
			}
			++str;
		}
	}

	void getTiles(string tiles, set<string> &tileSet)
	{
        stringstream strStream;
		int len = tiles.size();
		if (len < 0)
			return;

		char *s = new char[len+1];
		vector<int> pos;
		pos.push_back(0);
		for(int i=0; i<len; i++)
		{
			s[i] = tiles[i];
			if(s[i]==',')
			{
				pos.push_back(i+1);
				s[i] = '\0';
			}
		}
		s[len]='\0';

		for(size_t i=0; i<pos.size();i++)
		{
			if(strlen(s+pos[i])==4)
			{
				tileSet.insert(s+pos[i]);
			}
			else if(strlen(s+pos[i])==9)
			{
				s[pos[i]+4]='\0';
				int begin = atoi(s+pos[i]);
				int end = atoi(s+pos[i]+5);
				for(int j=begin; j<=end; j++)
				{
					tileSet.insert(i2s(j));
				}
			}
		}
		delete []s;
	}
    
    void getFovs(string fovs, set<string> &fovSet)
    {
        long len = fovs.size();
        if (len < 0)
            return;
        
        char *s = new char[len+1];
        vector<int> pos;
        pos.push_back(0);
        for(int i=0; i<len; i++)
        {
            s[i] = fovs[i];
            if(s[i]==',')
            {
                pos.push_back(i+1);
                s[i] = '\0';
            }
        }
        s[len]='\0';
        
        for(size_t i=0; i<pos.size();i++)
        {
            if(strlen(s+pos[i])==8)
            {
                fovSet.insert(s+pos[i]);
            }
            else
            {
                LOG(ERROR, "--fov parameter format: " + fovs + " error");
				exit(1);
            }
        }
        delete []s;
    }
    
    string i2s(int i)
    {
        stringstream ss;
        ss << i;
        return ss.str();
    }
	
	int CopyFile(const char *in, const char *out)  
    {  
       FILE *fp1;  
       fp1 = fopen(in, "r");  
       FILE *fp2;  
       fp2 = fopen(out, "w");  
       char buff[1024] = {'\0'};  
       int count = 0;  
       while((count = fread(buff, 1, 1024, fp1)) != 0)  
       {  
          fwrite(buff, 1, count, fp2);  
       }  
       fclose(fp1);  
       fclose(fp2);
	   //delete []buff;
       return 0;  
    }
	
	bool gzLoad(const char *gzfn, const char *out)
	{
	    //open .gz file
	    gzFile gzfp = gzopen(gzfn,"rb");
	    if(!gzfp)
	    {
	        return false;
	    }
		
		ofstream outfile;
		outfile.open(out);
 
	    //read and add it to out
	    unsigned char buf[GZ_BUF_SIZE];
	    int have;
	    while( (have = gzread(gzfp,buf,GZ_BUF_SIZE)) > 0)
	    {
			outfile << buf;
	       // outfile.append((const char*)buf,have);
	    }
 
	    //close .gz file
	    gzclose(gzfp);
		//delete []buf;
	    return true;
	}

	vector<string> split(const std::string &s, char delim) {
		    stringstream ss(s);
		    string item;
		    vector<string> elems;
		    while(getline(ss, item, delim)) {
		        elems.push_back(item);
		    }
		    return elems;
	}

	vector<string> split(char *s, char delim) {
	    stringstream ss(s);
	    string item;
	    vector<string> elems;
	    while(getline(ss, item, delim)) {
	        elems.push_back(item);
	    }
	    return elems;
	}

    
}  // namespace PreProcessTool


