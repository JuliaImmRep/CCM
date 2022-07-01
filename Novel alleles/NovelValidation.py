#!/usr/bin/python
import pandas as pd
import sys,re,glob
from Bio import SeqIO
import multiprocessing as mp
from Bio.Seq import Seq

def get_fasta_df(infile):
	df1 = pd.DataFrame()
	df1['Name'] = [i.id for i in SeqIO.parse(infile,"fasta")]
	df1['Seq'] = [str(i.seq) for i in SeqIO.parse(infile,"fasta")]
	return df1 

def validation(ss):
	sangerDF = sangerdf.copy()
	goaldf = sangerDF[sangerDF['Seq'].str.contains(ss) | sangerDF['Reverse'].str.contains(ss)].reset_index(drop = True)
	if len(goaldf) > 0:
		posData,nameData,seqData = [],[],[]
		for n in range(len(goaldf)):
			san = goaldf['Seq'][n]
			subdata = [i.start() for i in re.finditer(ss, san)]
			if subdata == []:
				san = goaldf['Reverse'][n]
				subdata = [i.start() for i in re.finditer(ss, san)]
			pos_s = subdata[0]
			pos_e = subdata[0] + len(ss)
			posData.append('%d:%d' % (pos_s,pos_e))
			seqData.append(san)
			nameData.append(goaldf['Name'][n])
		posData = '|'.join(posData)
		seqData = '|'.join(seqData)
		nameData = ':'.join(nameData)
	else:
		posData,nameData,seqData = '0','0','0'
	return [posData,nameData,seqData]

def compare_with_self(selfdf):
    selfdf = selfdf.drop_duplicates('Seq_val').reset_index(drop = True)
    subselfdf = pd.DataFrame()
    for n in range(len(selfdf)):
        sdf = selfdf.drop([n])
        sudf = [i for i in sdf['Seq_val'] if selfdf['Seq_val'][n] in i or i in selfdf['Seq_val'][n]]
        if len(sudf) == 0:
        #sudf = sdf[sdf['Seq_val'].str.contains(selfdf['Seq_val'][n])]
        #if sudf.empty == True:
            subselfdf = subselfdf.append(selfdf.iloc[n], ignore_index = True)
    subselfdf = subselfdf.reset_index(drop = True)
    return subselfdf

def get_donor_message(alleles):
	donorFiles = sorted(glob.glob('data/Database-CCM*-V.txt'))
	donorAl = [pd.read_csv(dfa, sep = '\t') for dfa in donorFiles]
	all_data = []
	for al in alleles:
		ALss = []
		if al in list(donorAl[0]['Allele']):
			ALss.append('CM1')
		if al in list(donorAl[1]['Allele']):
			ALss.append('CM2')
		if al in list(donorAl[2]['Allele']):
			ALss.append('CM3')
		all_data.append('|'.join(ALss))
	return all_data

def main():
	fpre = pd.read_csv('data/Real-novel-V-BCR.txt', sep = '\t')
	fpre['Donor_message'] = get_donor_message(list(fpre['Alleles']))
	fpre['Seq_val'] = fpre['Seq'].str[1:-25]
	
	f_imgt615 = pd.read_csv('data/IMGT_615_DF.txt', sep = '\t')
	fpre['IfIGMT615'] = ['Y' if len(f_imgt615[f_imgt615['Seq'].str.contains(i)]) > 0 else 'N' for i in fpre['Seq_val']]
	fpre = fpre[fpre['IfIGMT615'] == 'N'].reset_index(drop = True)
	
	fpre = fpre.drop_duplicates('Seq_val', keep = False)
	fpre = compare_with_self(fpre)
	fpre.index = list(fpre['Seq'])

	pool = mp.Pool()
	results = pool.map(validation, list(fpre['Seq_val']))
	df2 = pd.DataFrame(results, columns = ['Pos', 'Name', 'Seq_sanger'], index = list(fpre['Seq']))

	df = pd.concat([fpre, df2], axis = 1)
	df['Validation_in_donor'] = ['|'.join([i for i in df['Donor_message'][k].split('|') if i in df['Name'][k]]) \
								if df['Name'][k] != 'na' else '0' for k in range(len(df))]
	df = df.drop_duplicates('Seq_val', keep = False)##########
	df.to_csv(outfile, sep = '\t', index = False)

if __name__ == '__main__':
	outfile = sys.argv[1]
	sangerList = glob.glob('data/CM*-Sanger-*.fasta')
	sangerdf = pd.concat([get_fasta_df(i) for i in sangerList])
	sangerdf['Reverse'] = [str((Seq(i).reverse_complement())) for i in sangerdf['Seq']]
	main()
