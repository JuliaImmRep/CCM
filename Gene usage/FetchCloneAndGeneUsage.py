import re,csv,glob,sys,os
import pandas as pd


def CalGeneUsage(df):
	mydf = df['Vgene'].value_counts(normalize=True).to_frame().reset_index()
	mydf.columns = ['Gene', 'Usage']
	return mydf

def FetchRandomNumber(infile, clone_num, isotype):
	sample = infile.split('/')[-3]
	df = pd.read_csv(infile, sep='\t')
	df = df[['VHits', 'clone_id']]
	df['Vgene'] = df['VHits'].str.split('*', expand=True)[0]
	df = df[['Vgene', 'clone_id']]
	df.drop_duplicates(keep='first', inplace=True, ignore_index=True)
	sample_list, gene_list = [], []
	if df.shape[0] >= int(clone_num):
		for t in range(100):
			df2 = df.sample(n=int(clone_num), random_state=t)
			df2['Isotype'] = isotype
			df2['Sample'] = sample+"-"+isotype
			df2['Time'] = t
			df2['Min_num'] = clone_num
			
			gene_df = CalGeneUsage(df2)
			gene_df['Isotype'] = isotype
			gene_df['Sample'] = sample+"-"+isotype
			gene_df['Time'] = t
			gene_df['Min_num'] = clone_num
			
			sample_list.append(df2)
			gene_list.append(gene_df)
			
	return sample_list, gene_list

def MergeSamplingDataAndGeneUsage(infiles, clone_num, isotype, outdir):
	sample_list, gene_list = [], []
	for infile in infiles:
		sample_list1, gene_list1 = FetchRandomNumber(infile, clone_num, isotype)
		sample_list = sample_list + sample_list1
		gene_list = gene_list + gene_list1
	final_sample = pd.concat(sample_list)
	final_gene = pd.concat(gene_list)
	final_sample.to_csv('%s/%s.Random.%s.origin-clone.txt'%(outdir, isotype, str(clone_num)), sep='\t', index=False)
	final_gene.to_csv('%s/%s.Random.%s.origin-geneusage.txt'%(outdir, isotype, str(clone_num)), sep='\t', index=False)
			
	return 0


def GetAllGeneUsage(infile, isotype):
	sample = infile.split('/')[-3]
	df = pd.read_csv(infile, sep='\t')
	df = df[['VHits', 'clone_id']]
	df['Vgene'] = df['VHits'].str.split('*', expand=True)[0]
	df = df[['Vgene', 'clone_id']]
	df.drop_duplicates(keep='first', inplace=True, ignore_index=True)
	
	gene_df = CalGeneUsage(df)
	gene_df['Isotype'] = isotype
	gene_df['Sample'] = sample+"-"+isotype
	gene_df['Time'] = -1
	gene_df['Min_num'] = -1
	
	return gene_df


def MergeAllCloneGU(infiles, outdir, isotype):
	all_list = []
	for file in infiles:
		gene_df = GetAllGeneUsage(file, isotype)
		all_list.append(gene_df)
	all_df = pd.concat(all_list)
	all_df.to_csv('%s/%s.All.geneusage.txt'%(outdir, isotype), sep='\t', index=False)

	return 0


def main():
	isotype = sys.argv[1] #IgM
	clone_num = sys.argv[2] #100
	outdir = sys.argv[3]
	Clone_files = sys.argv[4] #Clone_files
	print (outdir)
	infiles = glob.glob('%s'%(Clone_files))
	MergeAllCloneGU(infiles, outdir, isotype)
	MergeSamplingDataAndGeneUsage(infiles, clone_num, isotype, outdir)


if __name__ == '__main__':
	main()
