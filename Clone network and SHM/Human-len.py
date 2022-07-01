#!/usr/bin/python
import pandas as pd
from Bio import Seq
import multiprocessing as mp
import sys

def get_len(ss):
	if len(ss) <= 5:
		return 5
	elif len(ss) >= 25:
		return 25
	else:
		return len(ss)

def parse(cdr3aaList, group):
    pool_1 =  mp.Pool()
    LenData = pool_1.map(get_len, cdr3aaList)
    lenList = sorted(list(set(LenData)))
    lenRatio = [100*(float(LenData.count(i))/len(LenData)) for i in lenList]
    subdf = pd.DataFrame()
    subdf['Length'] = ['<= 5'] + lenList[1:-1] + ['>= 25']
    subdf['Ratio'] = lenRatio
    subdf['Group'] = group
    return subdf

def tran(seq):
	return Seq.Seq(seq).transcribe().translate()

def main():
	fcm = pd.read_csv(MetaAAFile, sep = '\t')
	cmall = list(set(fcm['CDR3aa']))
	cmpub = list(set(fcm[fcm['Donor count'] > 1]['CDR3aa']))
	
	fhumanall = pd.read_csv(HuamanAllClone,sep = '\t', header = None)
	fhumanall.columns = ['Clone','CDR3aa','a']
	humanall = list(set(list(fhumanall['CDR3aa'])))

	fhumanpub = pd.read_csv(HuamanPublicClone,sep = '\t', header = None)
	fhumanpub.columns = ['Clone','a','b']
	fhumanpub['nt'] = fhumanpub['Clone'].str.split('_', expand = True)[2]
	pool = mp.Pool()
	humanpub = pool.map(tran, list(set(list(fhumanpub['nt']))))

	lendf = pd.concat([parse(humanall,'All'), parse(humanpub, 'Public')])
	lendf.to_csv('example-len-distribution-human.txt', sep = '\t', index = False)

	data = [[len(set(cmall) & set(humanall)), len(set(cmall) & set(humanpub))],\
		[len(set(cmpub) & set(humanall)), len(set(cmpub) & set(humanpub))]]
	df = pd.DataFrame(data, columns = ['Human all', 'Human public'], index = ['CM all', 'CM public'])
	df.to_csv('example-CDR3aa-overlap-CM-human.txt', sep = '\t')

if __name__ == '__main__':
	MetaAAFile = sys.argv[1]
	HuamanAllClone = sys.argv[2]
	HuamanPublicClone = sys.argv[3]
	main()
