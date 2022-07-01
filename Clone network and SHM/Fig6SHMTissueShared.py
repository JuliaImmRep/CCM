#!/usr/bin/python
import pandas as pd
import numpy as np
import multiprocessing as mp
import sys,glob

def main():
	f = pd.read_csv('example-BCR-SHM-meta-data-6H.txt', sep = '\t')
	allclone = max(f['clone_id'])
	example_x = []
	for i in range(int(allclone/10000)+1):
		subdf = f[(f['clone_id'] >= i*10000) & (f['clone_id'] < (i+1)*10000)]
		example_x.append(subdf)
		#subdf.to_csv('example_x_%s' % str(i).zfill(4), sep = '\t', index = False)
	return example_x

def main_1(f):
        #f = pd.read_csv(infile, sep = '\t', usecols = ['clone_id','SHM','Donor','Tissue','Isotype'])
        f['SHM'] = f['SHM']*100

        pool = mp.Pool()
        results = pool.map(sub_main, list(f.groupby('clone_id')))
        df = pd.concat(results)
        df = df.sort_values('clone_id')
        return df
        #df.to_csv('example_shm_%s' % infile, sep = '\t', index = False)

def sub_main(groupdf):
        index,group = groupdf[0],groupdf[1]
        subdf = group.groupby(['Donor','Tissue'], as_index = False).mean().reset_index(drop = True)
        subdf = subdf.groupby(['Tissue'], as_index = False).mean().reset_index(drop = True)
        tissueList = list(set(group['Tissue']))
        if len(tissueList) == 1:
                subdf['Group'] = 'Private'
        else:
                subdf['Group'] = 'Shared'
        subdf['TissueList'] = '|'.join(sorted(tissueList))
        subdf['Tissue Count'] = len(tissueList)
        subdf = subdf[["clone_id","SHM","Tissue","Group","TissueList","Tissue Count"]]
        return subdf

def main_2():
	df = pd.concat(example_shm)
	df = df.sort_values('clone_id')
	df.to_csv('example-SHM-tissue-shared.txt', sep = '\t', index = False)

if __name__ == '__main__':
	example_x_data = main()
	example_shm = [main_1(sf) for sf in example_x_data]
	main_2()
