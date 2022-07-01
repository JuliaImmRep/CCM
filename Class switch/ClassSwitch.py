#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re,csv,glob,sys,os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from scipy import stats
from collections import Counter
from upsetplot import from_memberships
from upsetplot import plot
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rcParams.update({'font.size': 8})
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams["font.family"] = 'Arial'


# In[2]:


group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],               'Alone':['Liver', 'Gallbladder', 'Parotid_gland'],               'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],               'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],               'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
new_dict = {}
for k,v in group_dict.items():
    for n in v:
        new_dict[n] = k


# In[3]:


def JudgePublic(samples):
    if 'CCM1' in samples and 'CCM2' in samples and 'CCM3' in samples:
        return "123_Public"
    elif 'CCM1' in samples and 'CCM2' in samples:
        return "12_Public"
    elif 'CCM1' in samples and 'CCM3' in samples:
        return "13_Public"
    elif 'CCM2' in samples and 'CCM3' in samples:
        return "23_Public"
    else:
        return "Individual"


# In[4]:


def SplitSampleList(sample_list, donor):
    new_list = [i for i in sample_list.split('|') if donor in i]
    return '|'.join(new_list)


# In[5]:


def extractIsotype(sample_list):
    list2 = sorted(list(set([i.split('-')[-1] for i in sample_list.split('|')])))
    return "|".join(list2)


# In[6]:


def JudgeInterTissueShare(samples):
    tissue_list = sorted(list(set([i.split('-')[1] for i in samples.split('|')])))
    group_list = sorted(list(set([new_dict[i] for i in tissue_list])))
    flag = 'Inter_tissue'
    if len(tissue_list) == 1:
        flag = 'Private'
    return pd.Series([flag,'|'.join(tissue_list), '|'.join(group_list)])


# In[8]:


def CountUpsetData2(count, isotype):
    key_list, value_list = [], []
    for i in count.keys():
        if i != isotype:
            key_list.append(i.split('|'))
            value_list.append(count[i])
    mylist = from_memberships(key_list, data=value_list)
    return mylist


# In[9]:


def CountClassSwitch(df_count, isotype):
    count_all = Counter([])
    for donor in ['CCM1', 'CCM2', 'CCM3']:
        donor_df = df_count[df_count['sample_list'].str.contains(donor)]
        donor_df['New_list'] = donor_df.apply(lambda x:SplitSampleList(x['sample_list'], donor), axis=1)
        donor_df['Isotype_list'] = donor_df.apply(lambda x:extractIsotype(x['New_list']), axis=1)
        IgD_df = donor_df[donor_df['Isotype_list'].str.contains(isotype)]
        mycount = Counter(IgD_df['Isotype_list'])
        count_all = count_all + mycount
    return count_all


# In[10]:


def CountClassSwitchIGMD(df_count):
    count_all = Counter([])
    for donor in ['CCM1', 'CCM2', 'CCM3']:
        donor_df = df_count[df_count['sample_list'].str.contains(donor)]
        donor_df['sample_list'] = donor_df['sample_list'].str.replace('IgM', 'IgMD')
        donor_df['sample_list'] = donor_df['sample_list'].str.replace('IgD', 'IgMD')
        donor_df['New_list'] = donor_df.apply(lambda x:SplitSampleList(x['sample_list'], donor), axis=1)
        donor_df['Isotype_list'] = donor_df.apply(lambda x:extractIsotype(x['New_list']), axis=1)
        IgD_df = donor_df[donor_df['Isotype_list'].str.contains('IgMD')]
        mycount = Counter(IgD_df['Isotype_list'])
        count_all = count_all + mycount
    return count_all


# In[11]:


def CountClassSwitchGastricIntestinalIGMD(df_count):
    count_all = Counter([])
    for donor in ['CCM1', 'CCM2', 'CCM3']:
        donor_df = df_count[df_count['sample_list'].str.contains(donor)]
        donor_df['sample_list'] = donor_df['sample_list'].str.replace('IgM', 'IgMD')
        donor_df['sample_list'] = donor_df['sample_list'].str.replace('IgD', 'IgMD')
        donor_df['New_list'] = donor_df.apply(lambda x:SplitSampleList(x['sample_list'], donor), axis=1)
        donor_df['Isotype_list'] = donor_df.apply(lambda x:extractIsotype(x['New_list']), axis=1)
        IgD_df = donor_df[donor_df['Isotype_list'].str.contains('IgMD')]
        IgD_df = IgD_df[IgD_df['Groups'].str.contains('intestine|Gastric')]
        mycount = Counter(IgD_df['Isotype_list'])
        count_all = count_all + mycount
    return count_all


# In[15]:


outdir = sys.argv[1]
infile = sys.argv[2]
df = pd.read_csv(infile, sep='\t')


# In[16]:


df['Judge_Public'] = df.apply(lambda x:JudgePublic(x['sample_list']), axis=1)
public_df = df[df['Judge_Public'] != 'Individual']


# In[17]:


public_df[['Judge_InterTissue', 'Tissues', 'Groups']] = public_df.apply(lambda x:JudgeInterTissueShare(x['sample_list']), axis=1)


# In[20]:


df['Judge_Public'] = df.apply(lambda x:JudgePublic(x['sample_list']), axis=1)
df[['Judge_InterTissue', 'Tissues', 'Groups']] = df.apply(lambda x:JudgeInterTissueShare(x['sample_list']), axis=1)


# In[22]:


public_df = df[df['Judge_Public'] != 'Individual']
individual_df = df[df['Judge_Public'] == 'Individual']


# In[23]:


### IgMD:
public_MD_count = CountClassSwitchIGMD(public_df)
individual_MD_count = CountClassSwitchIGMD(individual_df)


# In[24]:


perc_public_MD = 100 - public_MD_count['IgMD']*100/sum(public_MD_count.values())
perc_individual_MD = 100 - individual_MD_count['IgMD']*100/sum(individual_MD_count.values())


# In[26]:


out = open('%s/Fig6E.result.txt'%(outdir), 'w')
out.write('Public\t%s\n'%(perc_public_MD))
out.write('Intra-individual\t%s\n'%(perc_individual_MD))
out.close()


# In[27]:


public_MD_list = CountUpsetData2(public_MD_count, 'IgMD')
individual_MD_list = CountUpsetData2(individual_MD_count, 'IgMD')


# In[28]:


fig = plt.figure(figsize=(10, 3))
plot(individual_MD_list, sort_by='cardinality', fig=fig, element_size=20)
plt.savefig('%s/FigS6D.All-IgMD.Individual.Class-Switch.png'%(outdir), dpi=1200, bbox_inches='tight')
plt.show()
plt.close()


# In[35]:


fig = plt.figure(figsize=(8, 3))
plot(public_MD_list, sort_by='cardinality', fig=fig, element_size=20)
plt.savefig('%s/Fig6F.All-IgMD.Public.Class-Switch.png'%(outdir), dpi=1200, bbox_inches='tight')
plt.show()
plt.close()


# In[30]:


gas_IgMD_count = CountClassSwitchGastricIntestinalIGMD(public_df)


# In[31]:


gas_MD_list = CountUpsetData2(gas_IgMD_count, 'IgMD')


# In[32]:


public_df2 = public_df[public_df['Groups'].str.contains('intestine|Gastric')]


# In[34]:


fig = plt.figure(figsize=(8, 3))
plot(gas_MD_list, sort_by='cardinality', fig=fig, element_size=17)
plt.savefig('%s/Fig6G.GasIntes-IgMD.Public.Class-Switch.png'%(outdir), dpi=1200, bbox_inches='tight')
plt.show()
plt.close()

