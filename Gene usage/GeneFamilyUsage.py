#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib
matplotlib.use('Agg')
import re,csv,glob,sys,os
from Bio import SeqIO
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import font_manager as fm, rcParams
import itertools as it
from scipy import stats
from scipy . spatial . distance import pdist
sns.set_style({'xtick.direction':'out', 'ytick.direction':'out'})
plt.rcParams.update({'font.size': 10})
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams["font.family"] = 'Arial'


# In[3]:


Family_color = {'IGHV1':'#e78ac3', 'IGHV2':'#8da0cb', 'IGHV3':'#fc8d62', 'IGHV4':'#66c2a5', 'IGHV5':'#ffd92f', 'IGHV6':'#a6d854', 'IGHV7':'#e5c494'}
Family_color_list = ['#e78ac3', '#8da0cb', '#fc8d62', '#66c2a5', '#ffd92f', '#a6d854', '#e5c494']
Isotype_color = {'IgM':'#7f8084', 'IgD':'#b89354', 'IgG':'#f0563f', 'IgA':'#23abd8', 'IgE':'#f8a275'}
Group_color = {'Immune related':'#FC7871', 'Alone':'#98c599', 'Gastric':'#A1AEBF', 'Small intestine':'#7a9eea', 'Large intestine':'#126f84'}


# In[4]:


group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],               'Alone':['Liver', 'Gallbladder', 'Parotid_gland'],               'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],               'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],               'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
new_dict = {}
for k,v in group_dict.items():
    for n in v:
        new_dict[n] = k


# In[5]:


outdir = sys.argv[1]
gene_file = sys.argv[2]
gene_name = gene_file.split('/')[-1].split('.')[1]
gene_name


# In[6]:


gene_df = pd.read_csv(gene_file, sep='\t')
gene_df.columns = ['Gene']+gene_df.columns.tolist()[1:]
gene_df['Family'] = gene_df['Gene'].str.split('-|S', expand=True)[0]


# In[7]:


family_groups = gene_df.groupby(['Family'])
data_list = []
for family, group in family_groups:
    new_group = group.sum()[1:-1].tolist()+[family]
    data_list.append(new_group)

mydf = pd.DataFrame(data_list, columns = gene_df.columns.tolist()[1:])
mydf.set_index('Family', inplace=True)


# In[8]:


data2 = mydf.unstack().to_frame().reset_index()
data2.columns = ['Sample', 'Family', 'Usage']
data2[['Donor', 'Tissue', 'Isotype']] = data2['Sample'].str.split('-', expand=True)
data2['Group'] = data2.apply(lambda x:new_dict[x['Tissue']], axis=1)


# In[9]:


def CreatPair(time_list, group_list):
    group_pair = [i for i in it.combinations(group_list, 2)]
    tmp_list = []
    for i in time_list:
        for n in group_pair:
            mycp1 = (i, n[0])
            mycp2 = (i, n[1])
            tmp_list.append((mycp1, mycp2))
    
    return tmp_list


# In[10]:


mypair = CreatPair(data2['Family'].unique().tolist(), data2['Isotype'].unique().tolist())


# In[12]:


plt.figure(figsize=(3, 1.5))
hue_order = ['IgM', 'IgG', 'IgA']
ax = sns.boxplot(x='Family', y='Usage', hue='Isotype', data=data2, linewidth=0.8, fliersize=1, hue_order=hue_order, palette = Isotype_color, saturation=1)
plt.legend(loc=[1,0], frameon=False)
#plt.xticks(rotation=45)
ax.tick_params(width=1)
w = 1
for linelabel in ['bottom', 'top', 'left', 'right']:
    ax.spines[linelabel].set_linewidth(w)
plt.xlabel('')
plt.ylabel('')
plt.ylim(-0.1, 1.2)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.get_legend().remove()
plt.savefig('%s/Fig3E.Family.usage.%s.isotype.png'%(outdir, gene_name), bbox_inches='tight', dpi=1200)
#plt.show()


# In[13]:


def CalT_test_ind(list1, list2):
    result = stats.ttest_ind(list1, list2)
    return result[-1]


# In[14]:


def fdr_correct(p,t):
    p2 = p.copy()
    idx = ~np.isnan(p2)
    p2[idx] = p2[idx] * idx.sum()/(p2[idx].argsort() + 1)
    p2[p2>t] = np.nan
    return p2


# In[15]:


def returnGroup(p_value):
    if p_value >= 0.05:
        return 0
    elif p_value < 0.001:
        return 3
    elif p_value < 0.01:
        return 2
    elif p_value < 0.05:
        return 1
    else:
        return -1


# In[16]:


data = []
for family, family_df in data2.groupby(['Family']):
    group_pair = [i for i in it.combinations(family_df['Isotype'].unique(), 2)]
    for n in group_pair:
        list1 = family_df[family_df['Isotype'] == n[0]]['Usage']
        list2 = family_df[family_df['Isotype'] == n[1]]['Usage']
        p_value = CalT_test_ind(list1, list2)
        data.append([family, n[0], n[1], p_value])
data_p = pd.DataFrame(data, columns = ['Family', 'Group1', 'Group2', 'P_value'])
data_p['FDR_p'] = fdr_correct(data_p['P_value'], 0.05)
data_p['FDR_group'] = data_p.apply(lambda x:returnGroup(x['FDR_p']), axis=1)
data_p.to_csv('%s/Fig3E.Isotype.t-test.txt'%(outdir), sep='\t', index=False)


# In[18]:


plt.figure(figsize=(3, 1.5))
hue_order = list(Group_color.keys())
ax = sns.boxplot(x='Family', y='Usage', hue='Group', data=data2, linewidth=0.5, fliersize=1, palette = Group_color, hue_order=hue_order, saturation=1)
plt.legend(loc=[1,0], frameon=False)
#plt.xticks(rotation=45)
ax.tick_params(width=1)
w = 1
for linelabel in ['bottom', 'top', 'left', 'right']:
    ax.spines[linelabel].set_linewidth(w)
plt.xlabel('')
plt.ylabel('')
plt.ylim(-0.1, 1.2)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.get_legend().remove()
plt.savefig('%s/Fig3F.Family.usage.%s.tissue-group.png'%(outdir, gene_name), bbox_inches='tight', dpi=1200)


# In[19]:


data = []
for family, family_df in data2.groupby(['Family']):
    group_pair = [i for i in it.combinations(family_df['Group'].unique(), 2)]
    for n in group_pair:
        list1 = family_df[family_df['Group'] == n[0]]['Usage']
        list2 = family_df[family_df['Group'] == n[1]]['Usage']
        p_value = CalT_test_ind(list1, list2)
        data.append([family, n[0], n[1], p_value])
data_p = pd.DataFrame(data, columns = ['Family', 'Group1', 'Group2', 'P_value'])
data_p['FDR_p'] = fdr_correct(data_p['P_value'], 0.05)
data_p['FDR_group'] = data_p.apply(lambda x:returnGroup(x['FDR_p']), axis=1)
data_p.to_csv('%s/Fig3F.Tissue-group.t-test.txt'%(outdir), sep='\t', index=False)
