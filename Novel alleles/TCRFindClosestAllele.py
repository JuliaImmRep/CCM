#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re,csv,glob,sys,os
from Bio import SeqIO
import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from collections import Counter


# In[2]:


def ExtractPosition(myseq, myref):
    mydata, myresult = [], []
    ref = myref.replace('-', 'N')
    gapnum = 0
    gappos = []
    gap_char_pos, true_char_pos, corr_pos = [], [], []
    for match in re.finditer('N+', ref):
        gapnum += match.end()-match.start()
        true_char_pos = match.end()
        mytrue_char_pos = true_char_pos-gapnum
        corr_pos.append([true_char_pos, gapnum, mytrue_char_pos])
        flag = 'I'
        insert_char = myseq[match.start(): match.end()]
        myresult.append([mytrue_char_pos, flag, myref[mytrue_char_pos], insert_char])
        for n in range(match.end()-match.start()):
            gappos.append(match.start()+n)
    ref2 = ref.replace('N', '')
    
    aa, bb = list(myref), list(myseq)
    for ind, pos in enumerate(gappos):
        mypos = pos - ind
        del aa[mypos]
        del bb[mypos]
    ref2 = ("".join(aa))
    seq2 = ("".join(bb))
    #print (ref2, seq2)
    zip_list = zip(ref2, seq2)
    for k1, v1 in enumerate(zip_list):
        if v1[0] != v1[1]:
            if v1[1] == '-':
                flag = 'D'
            else:
                flag = 'S'
            myresult.append([k1, flag, v1[0],v1[1]])
    
    return myresult


# In[3]:


# Cal PSA
def GetPairwise(seq, ref, allele_id, gene_id):
    aligndata = pairwise2.align.localms(seq, ref, 2,-1,-5,-1)[0]
    myseq = aligndata[0]
    myref = aligndata[1]
    myresult = ExtractPosition(myseq, myref)
    myresult_df = pd.DataFrame(myresult, columns = ['Pos', 'Type', 'Ref', 'Seq'])
    myresult_df['Allele'] = allele_id
    myresult_df['Gene'] = gene_id
    myresult_df.sort_values('Pos', inplace=True)
    return myresult_df


# In[17]:


outdir = sys.argv[1]
ref_file = sys.argv[2]
novel_file = sys.argv[3]


# In[18]:


ref_dict = {}
for rec in SeqIO.parse(ref_file, 'fasta'):
    ref_dict[rec.id] = str(rec.seq)


# In[19]:


novel_df = pd.read_csv(novel_file, sep='\t')
novel_dict = dict(zip(novel_df['Alleles'], novel_df['Seq']))


# In[20]:


result_list = []
for k1, v1 in novel_dict.items():
    for k2, v2 in ref_dict.items():
        pairwise_df = GetPairwise(v1, v2, k1, k2)
        result_list.append(pairwise_df)


# In[21]:


final_df = pd.concat(result_list)


# In[22]:


final_df2 = final_df


# In[33]:


num_df = final_df2.groupby(['Allele', 'Gene']).size().to_frame().reset_index()
num_df.columns = ['Allele', 'Gene', 'Number']


# In[34]:


closest_data = []
for allele, allele_df in num_df.groupby(['Allele']):
    allele_df.sort_values('Number', inplace=True)
    closest_data.append([allele, allele_df['Gene'].tolist()[0], allele_df['Number'].tolist()[0]])


# In[35]:


closest_df = pd.DataFrame(closest_data, columns = ['Allele', 'Reported', 'Mismatch'])
closest_df['Chain'] = 'TCR'
closest_df.to_csv('%s/FigS2B.TCR.closest.txt'%(outdir), sep='\t', index=False)
