#!/usr/bin/python
import sys
import pandas as pd

def main():
	f = pd.read_csv(MetaAAFile, sep = '\t')
	f3 = pd.read_csv(RAPIDAnnoFile, sep = '\t')
	f3 = f3[f3['chain_type'].str.contains('Heavy_Chain')].reset_index(drop = True)

	cdr3aa = list(set(f['CDR3aa']))
	cdr3aa1 = list(set(f[f['Donor count'] > 1]['CDR3aa']))
	cdr3aa2 = list(set(cdr3aa1) & set(f3['cdr3_aa_seq']))
	cdr3aa3 = list(set(f3['cdr3_aa_seq']))

	LenData = [5 if len(i) <= 5 else 25 if len(i) >= 25 else len(i) for i in cdr3aa]
	LenData1 = [5 if len(i) <= 5 else 25 if len(i) >= 25 else len(i) for i in cdr3aa1]
	LenData2 = [5 if len(i) <= 5 else 25 if len(i) >= 25 else len(i) for i in cdr3aa2]
	LenData3 = [5 if len(i) <= 5 else 25 if len(i) >= 25 else len(i) for i in cdr3aa3]

	lenList = sorted(list(set(LenData)))
	lenList1 = sorted(list(set(LenData1)))
	lenList2 = sorted(list(set(LenData2)))
	lenList3 = sorted(list(set(LenData3)))

	lenRatio = [100*(float(LenData.count(i))/len(LenData)) for i in lenList]
	lenRatio1 = [100*(float(LenData1.count(i))/len(LenData1)) for i in lenList1]
	lenRatio2 = [100*(float(LenData2.count(i))/len(LenData2)) for i in lenList2]
	lenRatio3 = [100*(float(LenData3.count(i))/len(LenData3)) for i in lenList3]

	df = pd.DataFrame()
	df['Length'] = lenList + lenList1 + lenList2 + lenList3
	df['Ratio'] = lenRatio + lenRatio1 + lenRatio2 + lenRatio3
	df['Group'] = ['Total' for i in lenList] + ['Inter-individual' for i in lenList1] +\
					['Inter-individual-RAPID' for i in lenList2] + ['RAPID' for i in lenList3]
	df.to_csv(outputFile, sep = '\t', index = False)

if __name__ == '__main__':
	MetaAAFile = sys.argv[1]
	RAPIDAnnoFile = sys.argv[2]
	outputFile = sys.argv[3]
	main()
