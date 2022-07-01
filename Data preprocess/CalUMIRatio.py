import re,csv,glob,sys,os
import pandas as pd

infile = sys.argv[1]
outdir = sys.argv[2]
mydf = pd.read_csv(infile, sep='\t')

err_fw, err_rev, err_dual = 0, 0, 0
sum_fw, sum_rev, sum_dual = 0, 0, 0

for fw_umi, fw_df in mydf.groupby(['UMI5']):
	sum_fw += 1
	if len(fw_df['Prim3ID'].unique()) > 1:
		err_fw += 1
	
		
for rev_umi, rev_df in mydf.groupby(['UMI3']):
	sum_rev += 1
	if len(rev_df['Prim3ID'].unique()) > 1:
		err_rev += 1


for dual_umi, dual_df in mydf.groupby(['UMI3', 'UMI5']):
	sum_dual += 1
	if len(dual_df['Prim3ID'].unique()) > 1:
		err_dual += 1

out = open('%s/UMI.isotype.txt'%(outdir), 'w')
out.write('Group\tTotal_group\tError_group\n')
out.write('Fw\t%s\t%s\n'%(sum_fw, err_fw))
out.write('Rev\t%s\t%s\n'%(sum_rev, err_rev))
out.write('Dual\t%s\t%s\n'%(sum_dual, err_dual))
