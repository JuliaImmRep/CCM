#!/usr/bin/python
import pandas as pd
import multiprocessing as mp
import sys

def get_tissue(line):
	subTissueList = sorted(list(set([i.split('-')[1] for i in line.split('|')])))
	return [len(subTissueList), '|'.join(subTissueList)]

def get_shared_data(subdf):
	sharedData = []
	for tis in TissueList:
		ssdf = subdf[subdf['TissueList'].str.contains(tis)]
		sharedNum = len(ssdf[ssdf['TissueNum'] > 1])
		unsharedNum = len(ssdf[ssdf['TissueNum'] == 1])
		sharedData.append([tis, sharedNum, unsharedNum, len(ssdf)])
	shareddf = pd.DataFrame(sharedData, columns = ['Tissue', '# Shared clone','# Unshared clone', '# All clone'])
	shareddf.to_csv('example-%sCR-Fig5A-sharedData.txt' % BT, sep = '\t', index = False)

def get_percentage(subdf):
	Tissuemessage = list(subdf['TissueList'])
	if BT == 'B':
		group_order = ['Immune related','Gastric','Small intestine','Large intestine','Others']
	else:
		group_order = ['Immune related','Gastric','Small intestine','Large intestine']
	data = [[0 for i in group_order] for j in TissueList]
	for mes in Tissuemessage:
		tissueList = mes.split('|')
		for tis in tissueList:
			restTis = list(set(tissueList) - set([tis]))
			goalGroup = list(set([new_dict[i] for i in restTis]))
			p1 = TissueList.index(tis)
			for group in goalGroup:
				p2 = group_order.index(group)
				data[p1][p2] += 1/len(goalGroup)
	df = pd.DataFrame(data, columns = group_order, index = TissueList)
	df.to_csv('example-%sCR-Fig5A-percentage.txt' % BT, sep = '\t')

def main():
	f1 = f[['clone_id','sample_list']]
	pool = mp.Pool()
	TissueListData = pool.map(get_tissue, list(f1['sample_list']))
	f1['TissueNum'] = [i[0] for i in TissueListData]
	f1['TissueList'] = [i[1] for i in TissueListData]
	
	get_shared_data(f1)

	f2 = f1[f1['TissueNum'] > 1].reset_index(drop = True)
	get_percentage(f2)

if __name__ == '__main__':
	TissueList = ['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN',
		'Gastric_cardia', 'Gastric_body', 'Gastric_pylorus',
		'Duodenum', 'Jejunum', 'Ileum', 'Cecum',
		'Ascending_colon', 'Transverse_colon', 'Descending_colon', 'Sigmoid_colon', 'Rectum']
	group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],\
		'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],\
		'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],\
		'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
	BT = sys.argv[1]
	if BT == 'B':
		TissueList += ['Liver', 'Gallbladder', 'Parotid_gland']
		group_dict['Others'] = ['Liver', 'Gallbladder', 'Parotid_gland']
	
	new_dict = {}
	for k,v in group_dict.items():
		for n in v:
			new_dict[n] = k
	if BT == 'B':
		f = pd.read_csv('example-BCR.AA.0.06.Overlapped.cloneNum.txt', sep = '\t')
	else:
		f = pd.read_csv('example-TCR.AA.0.Overlapped.cloneNum.txt', sep = '\t')
	main()
