#!/usr/bin/python
import pandas as pd
import sys,glob
import multiprocessing as mp

def parse(cdr3aaLen, isotype):
    LenData = [5 if int(i) <= 5 else 25 if int(i) >= 25 else int(i) for i in cdr3aaLen]
    lenList = sorted(list(set(LenData)))
    lenRatio = [100*(float(LenData.count(i))/len(LenData)) for i in lenList]
    subdf = pd.DataFrame()
    subdf['Length'] = lenList
    subdf['Ratio'] = lenRatio
    subdf['Group'] = [isotype for i in lenRatio]
    return subdf

def getdf(infile):
	subdf = pd.read_csv(infile, sep = '\t', usecols = ['CDR3aa','clone_id','sample_list'])
	subdf['Donor count'] = [len(list(set([i.split('-')[0] for i in j.split('|')]))) for j in subdf['sample_list']]
	subdf['IsotypeList'] = ['|'.join(list(set([i.split('-')[2] for i in j.split('|')]))) for j in subdf['sample_list']]
	subdf['CDR3aa len'] = subdf['CDR3aa'].apply(len)
	subdf = subdf[['clone_id','IsotypeList','Donor count','CDR3aa len']]
	return subdf

def main():
	infiles = glob.glob('BCRExampleClone/BCR-CCM*-BM-IgM-ChangeO-GenoNovel.txt')
	pool = mp.Pool()
	results = pool.map(getdf, infiles)
	df = pd.concat(results)
	df = df.drop_duplicates(['clone_id','IsotypeList','Donor count','CDR3aa len']).reset_index(drop = True)
	df = df.sort_values('clone_id').reset_index(drop = True)
	df.to_csv('example-Meta-clonedata-isotype.txt', sep = '\t', index = False)
	
	for donor in [0,1]:
		isotypeList = ['IgM','IgD','IgG','IgA','IgE']
		data = df[df['Donor count'] > donor].reset_index(drop = True)
		lenData = []
		for iso in isotypeList:
			ssdf = data[data['IsotypeList'].str.contains(iso)].reset_index(drop = True)
			lenData.append(parse(list(ssdf['CDR3aa len']), iso))
		df1 = pd.concat(lenData)
		df1.to_csv('example-Len-isotype-%s.txt' % group[donor], sep = '\t', index = False)

if __name__ == '__main__':
	group = ['All', 'Public']
	main()
