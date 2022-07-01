#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib
matplotlib.use('Agg')
import re,csv,glob,sys,os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib import font_manager as fm, rcParams
import pandas as pd
import numpy as np
from scipy import stats
from scipy . spatial . distance import pdist

sns.set_style('ticks')
plt.rcParams.update({'font.size': 12})
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams["font.family"] = 'Arial'



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
Group_color = {'Immune related':'#FC7871', 'Alone':'#98c599', 'Gastric':'#A1AEBF', 'Small intestine':'#7a9eea', 'Large intestine':'#126f84'}
gene_sort = list(Tissue_color.keys())




group_dict = {'Immune related':['PBMC', 'BM', 'Spleen', 'Thymus', 'Tonsils', 'SLN', 'ALN', 'HLN', 'MLN', 'ILN'],               'Alone':['Liver', 'Gallbladder', 'Parotid_gland'],               'Gastric':['Gastric_cardia', 'Gastric_body', 'Gastric_pylorus'],               'Small intestine':['Duodenum', 'Jejunum', 'Ileum'],               'Large intestine':['Cecum','Ascending_colon','Transverse_colon', 'Descending_colon','Sigmoid_colon', 'Rectum']}
new_dict = {}
for k,v in group_dict.items():
    for n in v:
        new_dict[n] = k
Group_color = {'Immune related':'#FC7871', 'Alone':'#98c599', 'Gastric':'#A1AEBF', 'Small intestine':'#7a9eea', 'Large intestine':'#126f84'}




def GetCorrelation(df, outdir, GeneName):
    mydf = df.T.fillna(0)
    mydict = {}
    sample_list = []
    for index, row in mydf.iterrows():
        mydict[index] = np.array(row)
        sample_list.append(index)
        
    data = []
    for k1 in sorted(sample_list):
        mydata = []
        mydata.append(k1)
        for k2 in sorted(sample_list):
            x1, x2 = mydict[k1], mydict[k2]
            X = np.vstack([x1, x2])
            idx = np.argwhere(np.all(X[..., :] == 0, axis=0))
            X2 = np.delete(X, idx, axis=1)
            dis = pdist(X2)
            mydata.append(dis[0])
        data.append(mydata)
    dis_df = pd.DataFrame(data, columns = ['Sample']+sorted(sample_list))
    dis_df.to_csv('%s/Fig3A.Dis.%s.txt'%(outdir, GeneName), sep='\t', index=False)
    PlotDistance(dis_df, outdir, GeneName)
    
    
    return 0




def PlotDistance(df, outdir, GeneName):
    mydf = df
    mydf['Donor'] = mydf['Sample'].str.split('-', expand=True)[0]
    mydf['Tissue'] = mydf['Sample'].str.split('-', expand=True)[1]
    mydf['Isotype'] = mydf['Sample'].str.split('-', expand=True)[2]
    mydf['Group'] = mydf.apply(lambda x:new_dict[x['Tissue']], axis=1)
    #print (mydf.head(3))
    mydf.set_index('Sample', inplace=True)
    isotype = mydf.pop('Isotype')
    donor = mydf.pop('Donor')
    tissue = mydf.pop('Tissue')
    group = mydf.pop('Group')
    
    Isotype_col = Isotype_color
    Donor_col = Donor_color
    Tissue_col = Tissue_color
    Group_col = Group_color
    
    row_isotype_colors = isotype.map(Isotype_col)
    row_donor_colors = donor.map(Donor_col)
    row_tissue_colors = tissue.map(Tissue_col)
    row_group_colors = group.map(Group_col)
    row_colors = pd.concat([row_donor_colors, row_isotype_colors, row_group_colors, row_tissue_colors],axis=1)
    
    mydf = mydf.fillna(-1)
    for m in ['ward']:
        plt.figure()
        g = sns.clustermap(mydf, row_colors=row_colors, cmap='Spectral', col_cluster=True, figsize=(11, 10), method=m, xticklabels=False, yticklabels=False, vmax=0.3)
        #g = sns.clustermap(mydf, row_colors=row_colors, cmap='Spectral', col_cluster=True, figsize=(15, 14), method=m)
        g.ax_col_dendrogram.set_ylim(0, g.ax_col_dendrogram.get_ylim()[1]*3)
        g.ax_row_dendrogram.set_xlim(g.ax_row_dendrogram.get_xlim()[0]*2, 0)
        g.cax.set_position([0.01, 0.35, 0.03, 0.3])
        for line in g.ax_row_dendrogram.collections:
            line.set_linewidth(1.5)
        for line in g.ax_col_dendrogram.collections:
            line.set_linewidth(1.5)

        plt.savefig('%s/Fig3A.Dis.%s.heatmap.pdf'%(outdir, GeneName), dpi=600, bbox_inches='tight')
        plt.show()
        plt.close()
    
    
    return 0




outdir = sys.argv[1]
gene_file = sys.argv[2]
gene_name = sys.argv[3]




gene_df = pd.read_csv(gene_file,sep='\t', index_col=0)
GetCorrelation(gene_df, outdir, gene_name)
