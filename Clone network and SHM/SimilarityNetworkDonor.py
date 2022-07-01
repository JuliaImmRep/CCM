#!/usr/bin/python
import pandas as pd
import multiprocessing as mp
import sys

def get_clone(tissue):
	subCloneList = list(set(f[f['sample_list'].str.contains('%s-%s' % (donor, tissue))]['clone_id']))
	return subCloneList

def main():
	pool = mp.Pool()
	cloneList = pool.map(get_clone, TissueList)
	data = []
	for n1 in range(len(TissueList)-1):
		for n2 in range(n1+1,len(TissueList)):
			nclone1 = len(cloneList[n1])
			nclone2 = len(cloneList[n2])
			nshareclone = len(list(set(cloneList[n1]) & set(cloneList[n2])))
			score = (2 * nshareclone)/(nclone1 + nclone2)
			if TissueList[n1] == 'PB':
				data.append(['PBMC', nclone1, TissueList[n2], nclone2, score])
			elif TissueList[n2] == 'PB':
				data.append([TissueList[n1], nclone1, 'PBMC', nclone2, score])
			else:
				data.append([TissueList[n1], nclone1, TissueList[n2], nclone2, score])
	df = pd.DataFrame(data, columns = ['Tissue_1','Tissue_1_nCDR3aaVJ','Tissue_2','Tissue_2_nCDR3aaVJ','Similarity score'])
	df.to_csv('example-Network-figure-%sCR-CCM%s.txt' % (BT, donor), index = False)

if __name__ == '__main__':
	TissueList = ['PB', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN',
		 'Gastric_cardia', 'Gastric_body', 'Gastric_pylorus',
		 'Duodenum', 'Jejunum', 'Ileum', 'Cecum',
		 'Ascending_colon', 'Transverse_colon', 'Descending_colon', 'Sigmoid_colon', 'Rectum']
	BT, donor = sys.argv[1], str(sys.argv[2])
	if BT == 'B':
		TissueList += ['Liver', 'Gallbladder', 'Parotid_gland']
		f = pd.read_csv('example-BCR.AA.0.06.Overlapped.cloneNum.txt', sep = '\t')
	else:
		f = pd.read_csv('example-TCR.AA.0.Overlapped.cloneNum.txt', sep = '\t')
	f['sample_list'] = f['sample_list'].apply(lambda x: x.replace('PBMC', 'PB'))
	main()
