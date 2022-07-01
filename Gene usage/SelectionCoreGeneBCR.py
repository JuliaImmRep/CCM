#!/usr/bin/env python
# coding: utf-8

# In[1]:

import re,csv,glob,sys,os
from Bio import SeqIO
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn import linear_model
import pandas as pd
import numpy as np
from scipy import stats
from scipy . spatial . distance import pdist


# In[2]:


def CountData(infile, outdir):
    x, y, medium_usage = [], [], []
    num = 0
    for rec in csv.reader(open(infile, 'r'), delimiter="\t"):
        if "Occur" in rec:
            continue
        else:
            num += 1
            x.append(num)
            arrlist = [float(i) for i in rec[1:-2]]
            y.append(arrlist)
            medium_usage.append(np.median(np.array(arrlist)))

    return medium_usage


# In[13]:


def CountCoef(numlist, outdir):
    mynum = np.cumsum(np.array(numlist))
    reg = linear_model.LinearRegression()
    mycoef = []
    for i in range(5, len(mynum)-6):
        data = mynum[i-5:i+6]
        reg.fit(np.array(range(len(data))).reshape(-1, 1), np.array(data).reshape(-1, 1))
        mycoef.append(reg.coef_[0])

    coregene_num = 0
    for ind, k in enumerate(mycoef):
        if k < 0.001:
            coregene_num = ind + 6
            break

    mycoefnum = []
    coregene_number = 0
    for i in range(5, len(mynum)-6):
        min_data = mynum[i-5]
        max_data = mynum[i+5]
        coef_num = (max_data - min_data)/11
        mycoefnum.append(coef_num)

    for ind2, k2 in enumerate(mycoefnum):
        if k2 < 0.0001:
            coregene_number = ind2 + 6
            break

    return coregene_number


# In[4]:


def WriteCore(number, core_num):
    if number <= core_num:
        return 'Core'
    else:
        return 'Non-core'


# In[18]:


outdir = sys.argv[1]
gene_file = sys.argv[2]
gene_name = gene_file.split('/')[-1].split('.')[1]


# In[19]:


df = pd.read_csv(gene_file, sep='\t',index_col=0)
df['Mean'] = df.mean(axis=1)
df['Occur'] = df.count(axis=1)-1
#df.sort_values(by=['Mean', 'Occur'], inplace=True, ascending=False)
df.sort_values(by=['Occur', 'Mean'], inplace=True, ascending=False)
df.fillna(0, inplace=True)
df.to_csv('%s/%s.Vgene.input.core.selection.txt'%(outdir, 'CM'), sep='\t')


# In[15]:


infile = '%s/%s.Vgene.input.core.selection.txt'%(outdir, 'CM')
medium_usage = CountData(infile, outdir)
Coregene_directNum = CountCoef(medium_usage, outdir)

mydf = pd.read_csv('%s/%s.Vgene.input.core.selection.txt'%(outdir, 'CM'), sep='\t')
rank_df = mydf[['Gene', 'Mean']]
rank_df['Rank'] = range(1, rank_df.shape[0]+1)
rank_df['Type'] = rank_df.apply(lambda x:WriteCore(x['Rank'], Coregene_directNum), axis=1)
rank_df.to_csv('%s/BCR.gene_rank.txt'%(outdir), sep='\t', index=False)


# In[16]:


def NormalizeCoreGene(df):
    df.set_index('Gene', inplace=True)
    for col in df.columns:
        df[col] = df[col]/df[col].sum()
    
    return df


# In[17]:


gene_df = pd.read_csv('%s/BCR.V.sample.geneusage.txt'%(outdir), sep='\t')
gene_core_df = gene_df[gene_df['Gene'].isin(rank_df[rank_df['Type'] == 'Core']['Gene'])]
gene_core_df = NormalizeCoreGene(gene_core_df)
gene_core_df.to_csv('%s/core-BCR.V.sample.geneusage.txt'%(outdir), sep='\t')