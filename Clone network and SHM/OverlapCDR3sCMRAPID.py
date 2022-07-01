#!/usr/bin/python
import sys,glob
import pandas as pd
import multiprocessing as mp

def getdf(infile):
	subdf = pd.read_csv(infile, sep = '\t', usecols = ['CDR3aa','clone_id','sample_list'])
	subdf['Donor count'] = [len(list(set([i.split('-')[0] for i in j.split('|')]))) for j in subdf['sample_list']]
	subdf = subdf[['clone_id','Donor count','CDR3aa']]
	return subdf


def main():
	annoFile = pd.read_csv(RAPIDAnnoFile, sep = '\t')
	annodf = annoFile[annoFile['chain_type'].str.contains('Heavy_Chain')].reset_index(drop = True)
	# Get the CDR3aa sequences in Known Abs dataset
	annoData = annodf[~(annodf['source_id'].str.contains('TheraSAbDab'))].reset_index(drop = True)
	annoCDR3aa = list(annoData['cdr3_aa_seq'])
	# Get the CDR3aa sequences in Therapeutic Abs dataset
	therasData = annodf[annodf['source_id'].str.contains('TheraSAbDab')].reset_index(drop = True)
	therasCDR3aa = list(therasData['cdr3_aa_seq'])

	infiles = glob.glob('BCRExampleClone/BCR-CCM*-BM-IgM-ChangeO-GenoNovel.txt')
	pool = mp.Pool()
	results = pool.map(getdf, infiles)
	f = pd.concat(results)
	f = f.drop_duplicates(['clone_id','Donor count','CDR3aa']).reset_index(drop = True)
	
	f = f.sort_values('clone_id')
	f.to_csv('example-meta-aa.txt', sep = '\t', index = False)
	
	if pub_total == 'public':
		subdf = f[f['Donor count'] > 1].reset_index(drop = True)
		cloneList = list(set(subdf['CDR3aa']))
	else:
		cloneList = list(set(f['CDR3aa']))
	
	goalannoCDR3aa = list(set(cloneList) & set(annoCDR3aa))
	goalAnnoData = annoData[annoData['cdr3_aa_seq'].isin(goalannoCDR3aa)]
	goalAnnoData.to_csv('example_Known_CDR3_aa_in_' + pub_total + '.txt', sep = '\t', index = False)

	goaltheraCDR3aa = list(set(cloneList) & set(therasCDR3aa))
	goalTheraData = therasData[therasData['cdr3_aa_seq'].isin(goaltheraCDR3aa)]
	goalTheraData.to_csv('example_Thera_CDR3_aa_in_' + pub_total + '.txt', sep = '\t', index = False)

	annotheraCDR3aa = list(set(annoCDR3aa) & set(therasCDR3aa))
	goalannothera = list(set(annoCDR3aa) & set(therasCDR3aa) & set(cloneList))

	Clone = len(list(set(cloneList) - set(annoCDR3aa) - set(therasCDR3aa)))
	Anno = len(list(set(annoCDR3aa) - set(cloneList) - set(therasCDR3aa)))
	Thera = len(list(set(therasCDR3aa) - set(cloneList) - set(annoCDR3aa)))
	CAT = len(goalannothera)
	CA = len(goalannoCDR3aa) - CAT
	CT = len(goaltheraCDR3aa) - CAT
	AT = len(annotheraCDR3aa) - CAT
	data = [Clone,Anno,Thera,CA,CT,AT,CAT]
	title = 'Clone\tKnown\tThera\tCK\tCT\tKT\tCKT\n'
	message = title + '\t'.join([str(i) for i in data])
	with open('example_Anno_venn_data_in_' + pub_total + '.txt','w') as fw:
		fw.write(message)

if __name__ == '__main__':
	RAPIDAnnoFile = sys.argv[1]
	pub_total = sys.argv[2]
	main()
