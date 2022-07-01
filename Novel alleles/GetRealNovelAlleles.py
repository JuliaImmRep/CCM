#!/usr/bin/python
import pandas as pd
from Bio import SeqIO
import sys,glob

def get_fasta_df(infile, cut):
	df1 = pd.DataFrame()
	df1['Alleles'] = [i.id for i in SeqIO.parse(infile,"fasta")]
	if cut == 1:
		df1['Seq'] = [str(i.seq)[24:] for i in SeqIO.parse(infile,"fasta")]
	else:
		df1['Seq'] = [str(i.seq) for i in SeqIO.parse(infile,"fasta")]#
	df1['Gene'] = df1['Alleles'].str.split('*', expand = True)[0]
	return df1

def if_unique(seq, gene, subdf):
	uniquedf = subdf[subdf['Seq'].str.contains(seq)]
	seq_status = 0
	if uniquedf.empty == True:
		seq_status += 1
	return seq_status

def compare_with_615_imgt(f, df_merge):
	result = pd.DataFrame()
	for index,line in f.iterrows():
		s = if_unique(line['Seq'], line['Gene'], df_merge)
		if s == 1:
			result = result.append(pd.Series({'Alleles':line['Alleles'], 'Gene':line['Gene'], 'Seq':line['Seq']}), ignore_index = True)
	return result

def compare_with_self(selfdf):
	selfdf = selfdf.drop_duplicates('Seq')
	selfdf = selfdf.reset_index(drop = True)
	subselfdf = pd.DataFrame()
	for n in range(len(selfdf)):
		sdf = selfdf.drop([n])
		sudf = sdf[sdf['Seq'].str.contains(selfdf['Seq'][n])]
		if sudf.empty == True:
			subselfdf = subselfdf.append(selfdf.iloc[n], ignore_index = True)
	subselfdf = subselfdf.reset_index(drop = True)
	return subselfdf

def main():
	df_615 = get_fasta_df(KIMDB, 1)
	df_imgt = get_fasta_df(IMGTDB, 1)
	DF_merge = pd.concat([df_615, df_imgt])
	DF_merge.to_csv('data/IMGT_615_DF.txt', sep = '\t', index = False)

	data = []
	for donor in [1,2,3]:
		cm1 = get_fasta_df('data/IgDiscover-CCM%s_V.fasta' % donor, 0)
		gn1 = pd.read_csv('data/Database-CCM%s-V.txt' %donor, sep = '\t')
		goalcm1 = cm1[cm1['Alleles'].isin(list(gn1['Allele']))].reset_index(drop = True)
		data.append(compare_with_615_imgt(goalcm1, DF_merge))
	data = pd.concat(data)
	uniqueGene = list(set(data['Gene']))
	df = compare_with_self(data)
	df = df.sort_values('Alleles')
	df.to_csv(outFile, sep = '\t', index = False)

if __name__ == '__main__':
	KIMDB = sys.argv[1]
	IMGTDB = sys.argv[2]
	outFile = sys.argv[3]
	main()
