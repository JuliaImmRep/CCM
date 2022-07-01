#!/usr/bin/python
import pandas as pd
import sys,glob
import re
import multiprocessing as mp
import numpy as np

def get_donor_count(ss):
	subData = len(list(set([i.split('-')[0] for i in ss.split('|')])))
	return subData

def merge_SHM(fileName):
	donor = fileName.split('/')[1].split('-')[1]
	tissue = fileName.split('/')[1].split('-')[2]
	isotype = fileName.split('/')[1].split('-')[3]
	f1 = pd.read_csv(fileName, sep = '\t', usecols = ['CDR3aa','SHMinfo','Vinfo','clone_id','sample_list'])
	f1 = f1[~f1['CDR3aa'].str.contains('\*')]
	f1['Alignment length'] = [int(i.split(';')[1]) for i in f1['Vinfo']]
	f1['Top24bp'] = [25-int(i.split(';')[-2]) if int(i.split(';')[-2]) <= 24 else 0 for i in f1['Vinfo']]
	f1['V Length'] = f1['Alignment length'] - f1['Top24bp']
	f1 = f1[f1['Alignment length'] >= 200]
	f1['SHM count'] = [len([i for i in list(re.findall('\d+', str(s))) if int(i) > 24]) for s in f1['SHMinfo']]
	f1['SHM'] = f1['SHM count']/f1['V Length']
	f1['Donor count'] = f1['sample_list'].apply(lambda x: get_donor_count(x))
	f = f1[['clone_id','SHM','SHM count','Donor count']]#############
	f = f.groupby('clone_id', as_index = False).mean().reset_index(drop = True)
	f = f.drop_duplicates().reset_index(drop = True)
	subdf = f.copy()
	subdf['Donor'] = donor
	subdf['Tissue'] = tissue
	subdf['Isotype'] = isotype
	return subdf
	
def main_1():
	fileNameList = glob.glob('BCRExampleClone/BCR-CCM*-BM-IgM-ChangeO-GenoNovel.txt')
	pool = mp.Pool()
	data = pool.map(merge_SHM, fileNameList)
	df = pd.concat(data)
	df.to_csv('Partial-BCR-SHM-meta-data.txt', sep = '\t', index = False)

def main_2():
	f = pd.read_csv('example-BCR-SHM-meta-data.txt', sep = '\t', usecols = ['clone_id','SHM','Tissue','Isotype','Donor count'])
	if clonetype == 'public':
		df = f[f['Donor count'] > 1].reset_index(drop = True)
	else:
		df = f.copy()
	df['SHM'] = df['SHM']*100
	df = df.groupby(['clone_id','Tissue','Isotype'], as_index = False).mean().reset_index()
	df = df[['clone_id','SHM','Tissue','Isotype']]
	df.to_csv('example-BCR-SHM-tissue-isotype-%sclone.txt' % clonetype, sep = '\t', index = False)

def main_3():
	f = pd.read_csv('example-BCR-SHM-meta-data.txt', sep = '\t', usecols = ['clone_id','SHM','Donor','Tissue','Isotype','Donor count'])
	if clonetype == 'public':
		df = f[f['Donor count'] > 1].reset_index(drop = True)
	else:
		df = f.copy()
	df['SHM'] = df['SHM']*100
	df = df.groupby(['clone_id','Donor','Isotype'], as_index = False).mean().reset_index()
	df = df.groupby(['clone_id','Isotype'], as_index = False).mean().reset_index()
	df = df[['clone_id','SHM','Isotype']]
	df.to_csv('example-BCR-SHM-only-isotype-%sclone.txt' % clonetype, sep = '\t', index = False)

def main_4():
	f = pd.read_csv('example-BCR-SHM-meta-data.txt', sep = '\t', usecols = ['clone_id','SHM','Donor','Tissue','Isotype','Donor count'])
	if clonetype == 'public':
		df = f[f['Donor count'] > 1].reset_index(drop = True)
	else:
		df = f.copy()
	df['SHM'] = df['SHM']*100

	TissueList = ['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN',
		'Liver', 'Gallbladder', 'Parotid_gland',
		'Gastric_cardia', 'Gastric_body', 'Gastric_pylorus',
		'Duodenum', 'Jejunum', 'Ileum', 'Cecum',
		'Ascending_colon', 'Transverse_colon', 'Descending_colon', 'Sigmoid_colon', 'Rectum']
	group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],\
		'Others':['Liver', 'Gallbladder', 'Parotid_gland'],\
		'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],\
		'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],\
		'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
	new_dict = {}
	for k,v in group_dict.items():
		for n in v:
			new_dict[n] = k

	df['Tissue group'] = df['Tissue'].map(new_dict)
	df = df.groupby(['clone_id','Donor','Isotype','Tissue group'], as_index = False).mean().reset_index()
	df = df.groupby(['clone_id','Isotype','Tissue group'], as_index = False).mean().reset_index()
	df = df[['clone_id','SHM','Isotype','Tissue group']]
	df.to_csv('example-BCR-SHM-tissuegroup-isotype-%sclone.txt' % clonetype, sep = '\t', index = False)

if __name__ == '__main__':
	step,clonetype = str(sys.argv[1]),sys.argv[2]
	if step == '1':
		main_1()
	elif step == '2':
		main_2()
	elif step == '3':
		main_3()
	else:
		main_4()
