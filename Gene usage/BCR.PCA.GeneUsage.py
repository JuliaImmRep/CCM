#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re,csv,glob,sys,os
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


Tissue_color = {'PBMC': '#cc0000', 'BM': '#d63333', 'Spleen': '#e06666', 'Thymus': '#eb9999',
                 'Tonsils': '#f5cccc', 'SLN': '#ff9966', 'ALN': '#ffad85', 'HLN': '#ffc2a3', 'MLN': '#ffd6c2',
                 'ILN': '#ffebe0', 'Liver': '#66cc99', 'Gallbladder': '#99ddbb', 'Parotid_gland': '#cceedd', 
                 'Gastric_cardia': '#00ccff', 'Gastric_body': '#55ddff', 'Gastric_pylorus': '#aaeeff',
                 'Duodenum': '#3399cc', 'Jejunum': '#77bbdd', 'Ileum': '#bbddee', 'Cecum': '#0066ff', 
                 'Ascending_colon': '#247cff', 'Transverse_colon': '#4992ff', 'Descending_colon': '#6da8ff', 
                 'Sigmoid_colon': '#92bdff', 'Rectum': '#b6d3ff'}
Donor_color = {'CCM1':'#92cfbc', 'CCM2':'#f19f60', 'CCM3':'#26b3e9'}
Isotype_color = {'IgM':'#7f8084', 'IgD':'#b89354', 'IgG':'#f0563f', 'IgA':'#23abd8', 'IgE':'#f8a275'}
Isotype_color2 = {'IgM':'#7f8084', 'IgD':'#b89354', 'IgG':'#ce8e90', 'IgA':'#81aad8', 'IgE':'#f8a275'}


# In[3]:


group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],               'Alone':['Liver', 'Gallbladder', 'Parotid_gland'],               'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],               'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],               'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
new_dict = {}
for k,v in group_dict.items():
    for n in v:
        new_dict[n] = k
Group_color = {'Immune related':'#FC7871', 'Alone':'#98c599', 'Gastric':'#A1AEBF', 'Small intestine':'#7a9eea', 'Large intestine':'#126f84'}


# In[4]:


outdir = sys.argv[1]
gene_file = sys.argv[2]


# In[5]:


gene_df = pd.read_csv(gene_file, sep='\t')
gene_df.set_index('Gene', inplace=True)


# In[14]:


gene_df2 = gene_df.T.fillna(0)


# In[15]:


X = gene_df2.values
y = gene_df2.index.tolist()


# In[8]:


pca = PCA(n_components=2)
X_p = pca.fit(X).transform(X)


# In[9]:


print (pca.explained_variance_ratio_)
print (pca.explained_variance_)


# In[16]:


plot_df = pd.DataFrame(X_p)
plot_df.columns = ['PC1', 'PC2']
plot_df['Sample'] = y
plot_df[['Donor', 'Tissue', 'Isotype']] = plot_df['Sample'].str.split('-', expand=True)
plot_df['Tissue_Group'] = plot_df.apply(lambda x:new_dict[x['Tissue']], axis=1)


# In[12]:


def PlotPCA(flag, plot_df):
    gene = 'V'
    if flag == 'Isotype':
            color = Isotype_color
    elif flag == 'Tissue':
        color = Tissue_color
    else:
        color = Group_color

    plt.figure(figsize=(2.5, 1.8))
    ax = sns.scatterplot(x='PC1', y='PC2', hue=flag, data=plot_df, palette = color, style='Donor')
    ax.tick_params(width=1)
    w = 1
    for linelabel in ['bottom', 'top', 'left', 'right']:
        ax.spines[linelabel].set_linewidth(w)
    plt.xlabel('')
    plt.ylabel('')
    ax.get_legend().remove()
    plt.savefig('%s/FigS3B-D.PCA.BCR.%s.%s.usage.png'%(outdir, flag, gene), dpi=1200, bbox_inches='tight')
    plt.show()
    plt.close()

    return 0


# In[13]:


PlotPCA('Isotype', plot_df)
PlotPCA('Tissue', plot_df)
PlotPCA('Tissue_Group', plot_df)
