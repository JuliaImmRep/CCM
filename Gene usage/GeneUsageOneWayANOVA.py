#!/usr/bin/env python
# coding: utf-8

# In[31]:


import matplotlib
matplotlib.use('Agg')
import re,csv,glob,sys,os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import font_manager as fm, rcParams
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import pandas as pd
import numpy as np
from scipy import stats
from scipy . spatial . distance import pdist
sns.set_style('ticks')
plt.rcParams.update({'font.size': 10})
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams["font.family"] = 'Arial'


# In[2]:


group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],               'Alone':['Liver', 'Gallbladder', 'Parotid_gland'],               'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],               'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],               'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
new_dict = {}
for k,v in group_dict.items():
    for n in v:
        new_dict[n] = k
Group_color = {'Immune related':'#FC7871', 'Alone':'#98c599', 'Gastric':'#A1AEBF', 'Small intestine':'#7a9eea', 'Large intestine':'#126f84'}


# In[14]:


outdir = sys.argv[1]
result_file = sys.argv[2]
result_df = pd.read_csv(result_file, sep='\t')
result_df.set_index('Gene', inplace=True)
result_core_df = result_df.unstack().to_frame().reset_index()
result_core_df.columns = ['Sample', 'Gene', 'Usage']


# In[15]:


result_core_df[['Donor', 'Tissue', 'Isotype']] = result_core_df['Sample'].str.split('-', expand=True)
result_core_df['Group'] = result_core_df.apply(lambda x:new_dict[x['Tissue']], axis=1)


# In[16]:


gene_groups = result_core_df.groupby(['Gene'])


# In[19]:


tissue_list, group_list, donor_list, isotype_list = [], [], [], []
for gene, group in gene_groups:
    mean_usage = group['Usage'].mean()
    if mean_usage > 0.01:
        model3 = ols('Usage~C(Donor)', data=group).fit()
        try:
            anova_table3 = anova_lm(model3)
            p3 = anova_table3.iat[0, 4]
            #if p3 < 0.05:
            donor_list.append([gene, p3])
        except:
            pass


# In[21]:


def JudgeDegree(p_value):
    if p_value < 0.001:
        return 3
    elif p_value < 0.01:
        return 2
    elif p_value < 0.05:
        return 1
    else:
        return 0


# In[28]:


gene_name = 'V'
df4 = pd.DataFrame(donor_list, columns = ['Gene', 'P_value'])
df4['Label'] = 'Individual'
all_df = df4
all_df = all_df[all_df['P_value'] < 0.05]
all_df['Gene'] = all_df['Gene'].str.replace('IGHV', '')
all_result_df = all_df.pivot('Gene', 'Label', 'P_value')
all_result_df['Individual2'] = all_result_df.apply(lambda x:JudgeDegree(x['Individual']), axis=1)
all_result_df.to_csv('%s/Fig3C.BCR.%s.ANOVA.table.txt'%(outdir, gene_name), sep='\t')


# In[29]:


all_result_df = pd.read_csv('%s/Fig3C.BCR.%s.ANOVA.table.txt'%(outdir, gene_name), sep='\t')
all_result_df.set_index('Gene', inplace=True)
all_result_df = all_result_df[['Individual2']]


# In[30]:


plot_df = all_result_df.T
plt.figure(figsize=(15, 0.35))
sns.heatmap(plot_df, cmap='Reds',linecolor='gray', linewidth=0.1, vmin=0, vmax=3)
plt.savefig('%s/Fig3C.BCR.ANOVA.%s.Gene_differ.heatmap.pdf'%(outdir, gene_name), bbox_inches='tight')
