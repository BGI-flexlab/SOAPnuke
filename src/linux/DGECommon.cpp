/*
 * Common.cpp
 *
 *  Created on: 2012-6-3
 *      Author: Haosen Chen
 *      Editor: Shuai JIANG
 * 		Mail  : chenhaosen@genomics.cn
 *              jaingshuai@genomics.cn
 */

#include "DGECommon.h"

namespace DGEProcessTool {

bool isSequence(const string &seq)
{
	if (seq.empty())
		return false;

	for (int i=seq.size() - 1; i>=0; --i)
	{
		switch (seq[i])
		{
		case 'A':
			break;
		case 'a':
			break;
		case 'C':
			break;
		case 'c':
			break;
		case 'G':
			break;
		case 'g':
			break;
		case 'T':
			break;
		case 't':
			break;
		default:
			return false;
		}
	}
	return true;
}

int adapterType(const string &adapter1, const string &adapter2)
{
	if (adapter1.empty() && adapter2.empty())
	{
		return 3;
	}

	bool isSeq1 = isSequence(adapter1);
	bool isSeq2 = isSequence(adapter2);

	if (adapter2.empty())
	{
		if (adapter1.empty())
		{
			return 3; //error
		}
		else
		{
			if (isSeq1)
			{
				return 1; //adapter sequence
			}
			else
			{
				return 2; //adapter list
			}
		}
	}

	if (isSeq1 && isSeq2)
		return 1;
	else if (!isSeq1 && !isSeq2)
		return 2;
	else
		return 3;
}

string getOutputFileName(string filename, string path)
{
	string outFileName;
	if (filename.substr(filename.size() - 2, 2) != "gz")
	{
		int pos = filename.rfind('/');
		if (pos == -1)
		{
			outFileName = path + "/" + CLEAN_FQ_PREFIX + filename + ".gz";
		}
		else
		{
			outFileName = path + "/" + CLEAN_FQ_PREFIX
					+ filename.substr(pos + 1);
		}
	}
	else
	{
		int pos = filename.rfind('/');
		if (pos == -1)
		{
			outFileName = path + "/" + CLEAN_FQ_PREFIX + filename;
		}
		else
		{
			outFileName = path + "/" + CLEAN_FQ_PREFIX
					+ filename.substr(pos + 1);
		}
	}

	return outFileName;
}


void printFqInfo(const string outDir, FqInfo &fqInfo)
{
	
	string outFile = outDir + "/" + SEQUENCING_QUALITY;
	ofstream out(outFile.c_str());
	if (!out)
	{
		LOG(ERROR, "open output file: " << outDir + "/" + SEQUENCING_QUALITY);
		exit(1);
	}

	// set style;
	int setwTmp = 1; long baseTmp = fqInfo.rawTotalBases;
	do
	{
		setwTmp++;
		baseTmp /= 10;
	}while(baseTmp > 0);
	if (setwTmp < 3)
		setwTmp = 3;
	int setwTmp2 = setwTmp + 8;
	int setwTmp3 = 5;
	int precisionTmp = 1;
	//


	out << fixed;
	out << "Item                                                              \t"
		<< std::left << setw(setwTmp2) << "raw reads" << "\t" << std::left << setw(setwTmp2) << "clean reads" << "\n";
	
	out << "Read length                                                       \t"
		<< std::left << setw(setwTmp2) << fqInfo.readslength << "\t" <<std::left << setw(setwTmp2) << "-" << "\n";
	
	out << "Total number of reads                                             \t"
		<< std::left << setw(setwTmp) << fqInfo.rawTotalReads << "(100.0%)" << "\t"
		<< std::left << setw(setwTmp) << fqInfo.cleanTotalReads << "(100.0%)" << "\n";
	
	out << "Number of filtered reads (%)                                     \t"
		<< std::left << setw(setwTmp) << fqInfo.rawTotalReads - fqInfo.cleanTotalReads
		<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp) << 100.0 * (fqInfo.rawTotalReads - fqInfo.cleanTotalReads) / fqInfo.rawTotalReads << "%)"
		<< "\t" << std::left << setw(setwTmp2) << "-" << "\n";
	
	out << "\n";
	
	out << "Total number of bases                                            \t"
		<< std::left << setw(setwTmp) << fqInfo.rawTotalBases << "(100.0%)" << "\t"
		<< std::left << "-\n";
//		<< std::left << setw(setwTmp) << fqInfo.cleanTotalBases << "(100.0%)" << "\n";
	
//	out << "Number of filtered bases (%)                                     \t"
//		<< std::left << setw(setwTmp) << fqInfo.rawTotalBases - fqInfo.cleanTotalBases
//		<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp) << 100.0 * (fqInfo.rawTotalBases - fqInfo.cleanTotalBases) / fqInfo.rawTotalBases << "%)"
//		<< "\t" << std::left << setw(setwTmp2) << "-" << "\n";

	out << "\n";

	out << "Number of base A (%)                                             \t"
		<< std::left << setw(setwTmp) << fqInfo.rawBaseA << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
		<< 100.0 * fqInfo.rawBaseA / fqInfo.rawTotalBases << "%)" << "\t"
		<< std::left << "-\n";
//		<< std::left << setw(setwTmp) << fqInfo.cleanBaseA << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//		<< 100.0 * fqInfo.cleanBaseA / fqInfo.cleanTotalBases << "%)" << "\n";
	
	out << "Number of base C (%)                                             \t"
		<< std::left << setw(setwTmp) << fqInfo.rawBaseC << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
		<< 100.0 * fqInfo.rawBaseC / fqInfo.rawTotalBases << "%)" << "\t"
		<< std::left << "-\n";
//		<< std::left << setw(setwTmp) << fqInfo.cleanBaseC << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//		<< 100.0 * fqInfo.cleanBaseC / fqInfo.cleanTotalBases << "%)" << "\n";
	
	out << "Number of base G (%)                                             \t"
		<< std::left << setw(setwTmp) << fqInfo.rawBaseG<< "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
		<< 100.0 * fqInfo.rawBaseG / fqInfo.rawTotalBases << "%)" << "\t"
		<< std::left << "-\n";
//		<< std::left << setw(setwTmp) << fqInfo.cleanBaseG << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//		<< 100.0 * fqInfo.cleanBaseG / fqInfo.cleanTotalBases << "%)" << "\n";
	
	out << "Number of base T (%)                                             \t"
		<< std::left << setw(setwTmp) << fqInfo.rawBaseT << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
		<< 100.0 * fqInfo.rawBaseT / fqInfo.rawTotalBases << "%)" << "\t"
		<< std::left << "-\n";
//		<< std::left << setw(setwTmp) << fqInfo.cleanBaseT << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//		<< 100.0 * fqInfo.cleanBaseT / fqInfo.cleanTotalBases << "%)" << "\n";
	
	out << "Number of base N (%)                                             \t"
		<< std::left << setw(setwTmp) << fqInfo.rawBaseN << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
		<< 100.0 * fqInfo.rawBaseN / fqInfo.rawTotalBases << "%)" << "\t"
		<< std::left << "-\n";
//		<< std::left << setw(setwTmp) << fqInfo.cleanBaseN << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//		<< 100.0 * fqInfo.cleanBaseN / fqInfo.cleanTotalBases << "%)" << "\n";
	
	out << "\n";

//	if (fqInfo.cleanQ20 == 0 && fqInfo.cleanQ30 == 0)
//	{
		out << "Number of base calls with quality value of 20 or higher (Q20+) (%)\t"
			<< std::left << setw(setwTmp) << fqInfo.rawQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo.rawQ20 / fqInfo.rawTotalBases << "%)" << "\t"
			<< std::left << setw(setwTmp2) << "-" << "\n";
		
		out << "Number of base calls with quality value of 30 or higher (Q30+) (%)\t"
			<< std::left << setw(setwTmp) << fqInfo.rawQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
			<< 100.0 * fqInfo.rawQ30 / fqInfo.rawTotalBases << "%)" << "\t"
			<< std::left << setw(setwTmp2) << "-" << "\n";
//	}
//	else
//	{
//		out << "Number of base calls with quality value of 20 or higher (Q20+) (%)\t"
//			<< std::left << setw(setwTmp) << fqInfo.rawQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//			<< 100.0 * fqInfo.rawQ20 / fqInfo.rawTotalBases << "%)" << "\t"
//			<< std::left << setw(setwTmp) << fqInfo.cleanQ20 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//			<< 100.0 * fqInfo.cleanQ20 / fqInfo.cleanTotalBases << "%)" << "\n";
//	
//		out << "Number of base calls with quality value of 30 or higher (Q30+) (%)\t"
//			<< std::left << setw(setwTmp) << fqInfo.rawQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//			<< 100.0 * fqInfo.rawQ30 / fqInfo.rawTotalBases << "%)" << "\t"
//			<< std::left << setw(setwTmp) << fqInfo.cleanQ30 << "(" << std::right << setw(setwTmp3) << setprecision(precisionTmp)
//			<< 100.0 * fqInfo.cleanQ30 / fqInfo.cleanTotalBases << "%)" << "\n";
//	}
	
	out.close();

/*
	outFile = outDir + "/" + FILTERED_READS;
	
	out.open(outFile.c_str());
	if (!out)
	{
		LOG(ERROR, "open output file: " << outDir + "/" + FILTERED_READS);
		exit(1);
	}


	int long filteredReads = fqInfo.rawTotalReads - fqInfo.cleanTotalReads;
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
	out << "Item                                    \t" << std::left << setw(setwTmp) << "Counts" << "\t" << "Percentage\n";
	out << "Total filtered reads (%)                \t" << std::left << setw(setwTmp) << filteredReads << "\t100.0%\n";
	out << "Reads with low quallity bases (%)       \t" << std::left << setw(setwTmp) << fqInfo.readsWithLowQual << "\t"
		<< setprecision(precisionTmp) << 100.0 * fqInfo.readsWithLowQual / filteredReads << "%\n";
	out << "Reads with invalid adapter sequences (%)\t" << std::left << setw(setwTmp) << fqInfo.readsWithAdapter << "\t"
		<< setprecision(precisionTmp) << 100.0 * fqInfo.readsWithAdapter / filteredReads << "%\n";
	out << "Reads with PolyA (%)                    \t" << std::left << setw(setwTmp) << fqInfo.readWithPolyA << "\t"
		<< setprecision(precisionTmp) << 100.0 * fqInfo.readWithPolyA / filteredReads << "%\n";
	out << "Reads with short valid length (%)       \t" << std::left << setw(setwTmp) << fqInfo.readWithShortValidLength << "\t"
		<< setprecision(precisionTmp) << 100.0 * fqInfo.readWithShortValidLength / filteredReads << "%\n";


	out.close();
*/
	outFile = outDir + "/" + BASE_DISTRIBUTIONS + ".txt";
	out.open(outFile.c_str());
	if (!out)
	{
		LOG(ERROR, "open output file: " << outDir + "/" + BASE_DISTRIBUTIONS + ".txt");
		exit(1);
	}

	setwTmp = 6;
	out << std::left << setw(setwTmp) << "Pos" << "\t" << std::right << setw(setwTmp) << "A" << "\t"
		<<  setw(setwTmp) << "C" << "\t" <<  setw(setwTmp) << "G" << "\t"
		<<  setw(setwTmp) << "T" << "\t" <<  setw(setwTmp) << "N" << "\n";
	setwTmp = 5;

	 for (unsigned int i=0; i<fqInfo.readslength; ++i)
	 {
		 out << std::left << setw(setwTmp) << i+1 << "\t" << std::right
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.base[i][0] / fqInfo.rawTotalReads << "%\t"
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.base[i][1] / fqInfo.rawTotalReads << "%\t"
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.base[i][2] / fqInfo.rawTotalReads << "%\t"
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.base[i][3] / fqInfo.rawTotalReads << "%\t"
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.base[i][4] / fqInfo.rawTotalReads << "%\n";
	 }
	 out.close();

	 outFile = outDir + "/" + DISTRIBUTION_OF_Q20_Q30 + ".txt";
	 out.open(outFile.c_str());
	 if (!out)
	 {
		 LOG(ERROR, "open output file: " << outDir + "/" + DISTRIBUTION_OF_Q20_Q30 + ".txt");
		 exit(1);
	 }

	 out <<"Position in reads\tPercentage of Q20+ bases\tPercentage of Q30+ bases\n";
	 for (unsigned int i=0; i<fqInfo.readslength; i++)
	 {
		 out << std::left << setw(setwTmp) << i+1 << "\t" << std::right
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.q20q30[i][0] / fqInfo.rawTotalReads << "%\t"
			 << setw(setwTmp) << setprecision(precisionTmp) << 100.0 * fqInfo.q20q30[i][1] / fqInfo.rawTotalReads << "%\n";
	 }
	 out.close();


	 outFile = outDir + "/" + BASE_QUALITY_VALUE_DISTRIBUTION + ".txt";
	 out.open(outFile.c_str());
	 if (!out)
	 {
		 LOG(ERROR, "open output file: " << outDir + "/" + BASE_QUALITY_VALUE_DISTRIBUTION + ".txt");
		 exit(1);
	 }
	 baseTmp = fqInfo.rawTotalReads;
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
	 if (fqInfo.rawTotalReads != 0)
	 {
		 if (fqInfo.rawTotalReads % 2 == 0)
		 {
			 qBoxPos[0][0] = fqInfo.rawTotalReads / 2;
			 if (qBoxPos[0][0] != 0)
				 qBoxPos[1][0] = qBoxPos[0][0] + 1;
		 }
		 else
		 {
			 qBoxPos[0][0] = (fqInfo.rawTotalReads + 1) / 2;
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

		 if ((fqInfo.rawTotalReads - qBoxPos[1][0] + 1) % 2 == 0)
		 {
			 qBoxPos[0][2] = (fqInfo.rawTotalReads - qBoxPos[1][0] + 1) / 2 + qBoxPos[1][0] - 1;
			 qBoxPos[1][2] = qBoxPos[0][2] + 1;
		 }else{
			 qBoxPos[0][2] = (fqInfo.rawTotalReads - qBoxPos[1][0]) / 2 + qBoxPos[1][0];
			 qBoxPos[1][2] = qBoxPos[0][2];
		 }

		 qBoxPos[0][3] = (int)(fqInfo.rawTotalReads * 0.1 + 0.5);
		 qBoxPos[1][3] = qBoxPos[0][3];

		 qBoxPos[0][4] = (int)(fqInfo.rawTotalReads * 0.9 + 0.5);
		 qBoxPos[1][4] = qBoxPos[0][4];
	 }

	unsigned long qBoxTmp[2][6];
	for (unsigned int i=0; i<fqInfo.readslength; ++i)
	{
		memset(qBoxTmp, 0, sizeof(long) * 2 * 6);
		for (int j=0; j<MAX_QUALITY; j++)
		{
			qBoxTmp[1][5] += fqInfo.qual[i][j] * j;
			qBoxTmp[0][5] += fqInfo.qual[i][j];
			for (int k=0; k<5; k++)
			{
				if(qBoxTmp[0][5] < qBoxPos[0][k])
					qBoxTmp[0][k] = j;
				if(qBoxTmp[0][5] < qBoxPos[1][k])
					qBoxTmp[1][k] = j;
			}
		}
		
		qBox[i][0] = 1.0 * qBoxTmp[1][5] / fqInfo.rawTotalReads;
		for (int k=0; k<5; k++)
		{
			qBox[i][k+1] = (qBoxTmp[0][k] + qBoxTmp[1][k]) / 2 + 1;
		}
	}
	 
	 out << std::left << setw(setwTmp) << "Pos";
	 for (int i=2; i<MAX_QUALITY; i++)
	 {
		 out << "\tQ" << setw(setwTmp-1) << i;
	 }
	 out << "\tMean\tMedian\tLower quartile\tUpper quartile\t10thpercentile\t90thpercentile";
	 out << "\n";

	 for (unsigned int i=0; i<fqInfo.readslength; ++i)
	 {
		 out << setw(setwTmp) << i+1;
		 for (int j=2; j<MAX_QUALITY; j++)
		 {
		     out << "\t" << setw(setwTmp) << fqInfo.qual[i][j];
		 }
		 for (int j=0; j<6; j++)
		 {
			 out << "\t" << setw(setwTmp) << setprecision(precisionTmp) << qBox[i][j];
		 }
		 out << "\n";
	 }
	 out.close();
/*
	 outFile = outDir + "/" + LENGTH_DISTRIBUTION;
	 out.open(outFile.c_str());if (!out)
	 {
		 LOG(ERROR, "open output file: " << outDir + "/" + LENGTH_DISTRIBUTION);
		 exit(1);
	 }

	 out<<"Length\tCounts  \tPercentage\n";
	 for (int i = fqInfo.lengthStart; i <= fqInfo.lengthEnd; i++)
	 {
		 out << setw(6) << i << "\t" << setw(8) << fqInfo.lengthDis[i] << "\t" 
			 << setprecision(precisionTmp) << 100.0 * fqInfo.lengthDis[i] /fqInfo.cleanTotalReads
			 << "%\n";
	 }
	 out.close();
	 */
}

void clearFqInfo(FqInfo &fqInfo)
{
	fqInfo.readslength = 0;
	fqInfo.rawTotalReads = 0;
	fqInfo.cleanTotalReads = 0;
	fqInfo.rawTotalBases = 0;
	fqInfo.cleanTotalBases = 0;
	fqInfo.rawBaseA = 0;
	fqInfo.cleanBaseA = 0;
	fqInfo.rawBaseC = 0;
	fqInfo.cleanBaseC = 0;
	fqInfo.rawBaseG = 0;
	fqInfo.cleanBaseG = 0;
	fqInfo.rawBaseT = 0;
	fqInfo.cleanBaseT = 0;
	fqInfo.rawBaseN = 0;
	fqInfo.cleanBaseN = 0;
	fqInfo.rawQ20 = 0;
	fqInfo.cleanQ20 = 0;
	fqInfo.rawQ30 = 0;
	fqInfo.cleanQ30 = 0;

	fqInfo.readsWithAdapter = 0;
	fqInfo.readsWithNrate = 0;
	fqInfo.readsWithLowQual = 0;
	fqInfo.readWithShortValidLength = 0;
	fqInfo.readWithPolyA = 0;
	fqInfo.adapter3Null = 0;
	fqInfo.insertNull = 0;
	fqInfo.adapter5Pollute = 0;

	memset(fqInfo.base, 0, sizeof(long) * 5 * MAX_LENGTH);
	memset(fqInfo.q20q30, 0, sizeof(long) * 2 * MAX_LENGTH);
	memset(fqInfo.qual, 0, sizeof(long) * MAX_QUALITY * MAX_LENGTH);
	memset(fqInfo.lengthDis, 0, sizeof(long) * MAX_LENGTH);
}

}  // namespace PreProcessTool




