#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re,csv,glob,sys,os
import pandas as pd
import numpy as np


# In[3]:


region_bcr_aa = ['FR1']*26+['CDR1']*13+['FR2']*17+['CDR2']*10+['FR3']*39
region_bcr_nt = ['FR1']*78+['CDR1']*39+['FR2']*51+['CDR2']*30+['FR3']*117
bcr_aa_dict = dict(zip(range(1, 106), region_bcr_aa))
bcr_nt_dict = dict(zip(range(1, 316), region_bcr_nt))
bcr_region_length = {'FR1':54, 'CDR1':39, 'FR2':51, 'CDR2':30, 'FR3':117}
bcr_region_aa = {'FR1':18, 'CDR1':13, 'FR2':17, 'CDR2':10, 'FR3':39}


# In[4]:


region_tcr_aa = ['FR1']*26+['CDR1']*12+['FR2']*17+['CDR2']*10+['FR3']*39
region_tcr_nt = ['FR1']*78+['CDR1']*36+['FR2']*51+['CDR2']*30+['FR3']*117
tcr_aa_dict = dict(zip(range(1, 105), region_tcr_aa))
tcr_nt_dict = dict(zip(range(1, 313), region_tcr_nt))
tcr_region_length = {'FR1':78, 'CDR1':36, 'FR2':51, 'CDR2':30, 'FR3':117}
tcr_region_aa = {'FR1':26, 'CDR1':12, 'FR2':17, 'CDR2':10, 'FR3':39}


# In[5]:




# In[19]:


outdir = sys.argv[1]
infile = sys.argv[2]
chain = 'TCR'


# In[6]:


if chain == 'BCR':
    nt_dict = bcr_nt_dict
    region_nt_dict = bcr_region_length
    region_aa_dict = bcr_region_aa
elif chain == 'TCR':
    nt_dict = tcr_nt_dict
    region_nt_dict = tcr_region_length
    region_aa_dict = tcr_region_aa
else:
    nt_dict, region_dict = {}, {}

df = pd.read_csv(infile, sep='\t')
df['Region'] = df.apply(lambda x:nt_dict[x['New_pos2']], axis=1)


data = []
SNP_total, SAP_total = 0, 0
for region, region_df in df.groupby(['Region']):
    region_nt = region_nt_dict[region]
    region_aa = region_aa_dict[region]
    SNP_num = region_df.shape[0]
    SAP_num = region_df[region_df['Judge'] == 'SAP'].shape[0]
    SNP_norm = SNP_num / region_nt
    SAP_norm = SAP_num / region_aa
    data.append([region, SNP_norm, 'SNPs', SNP_num])
    data.append([region, SAP_norm, 'SAPs', SAP_num])
    SNP_total += SNP_num
    SAP_total += SAP_num


# In[22]:


result_df = pd.DataFrame(data, columns = ['Region', 'Number', 'Group', 'nTotal'])
result_df2 = result_df.pivot('Region', 'Group', 'Number')
result_df2.to_csv('%s/FigS2E.SNP-SAP.%s.normalize.length.txt'%(outdir, chain), sep='\t')
