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


outdir2 = sys.argv[1]
novel_file = sys.argv[2]
nogap_ref = sys.argv[3]
gap_ref = sys.argv[4]
chain = 'TCR'
blast_file = sys.argv[5]


blast_dict = {}
for rec in csv.reader(open(blast_file, 'r'), delimiter='\t'):
    if rec[1] == 'Reported':
        continue
    else:
        blast_dict[rec[0]] = rec[1]


# In[5]:


ref_dict = {}
for rec in SeqIO.parse(nogap_ref, 'fasta'):
    ref_dict[rec.id] = str(rec.seq)


# In[6]:


novel_df = pd.read_csv(novel_file, sep='\t')
novel_df['Reported'] = novel_df.apply(lambda x:blast_dict[x['Alleles']], axis=1)
novel_df['Ref'] = novel_df.apply(lambda x:ref_dict[x['Reported']], axis=1)
novel_df.head(3)


# In[7]:


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


# In[8]:


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


# In[9]:


df_list = []
for index,row in novel_df.iterrows():
    mydf = GetPairwise(row['Seq'], row['Ref'], row['Alleles'], row['Reported']) #for TCR
    df_list.append(mydf)
final_df = pd.concat(df_list)


# In[10]:


gene_dict = {}
for rec2 in SeqIO.parse(gap_ref,  'fasta'):
    gene_dict[rec2.id] = str(rec2.seq)


# In[11]:


def CorrPos(gap_seq):
    gappos, pos_result = [], []
    seq1 = gap_seq.replace('.', 'N')
    gapnum = 0
    
    for match in re.finditer('N+', seq1):
        gapnum += match.end()-match.start()
        true_char_pos = match.end()
        for n in range(match.end()-match.start()):
            gappos.append(match.start()+n)
            
    i = 0
    for ind1, item1 in enumerate(gap_seq):
        if ind1 in gappos:
            continue
        else:
            pos_result.append([item1, ind1, i])
            i += 1
    
    return pos_result


# In[12]:


myfinal_list = []
for label, group in final_df.groupby(['Allele', 'Gene']):
    allele_id = label[0]
    gene_id = label[-1]
    pos_result = CorrPos(gene_dict[gene_id])
    pos_df = pd.DataFrame(pos_result, columns = ['Char', 'gap_pos', 'remove_pos'])
    pos_dict = dict(zip(pos_df['remove_pos'], pos_df['gap_pos']))
    group['New_pos'] = group.apply(lambda x:pos_dict.get(x['Pos'], 'End'), axis=1)
    myfinal_list.append(group)
myfinal_df = pd.concat(myfinal_list)


# In[16]:


myfinal_df = myfinal_df[myfinal_df['New_pos'] != 'End']#for TCR
myfinal_df['New_pos'] = myfinal_df['New_pos'].astype(int)#for TCR
myfinal_df['New_pos2'] = myfinal_df['New_pos']+1#for TCR


# In[17]:


myfinal_df = myfinal_df[myfinal_df['New_pos2'] < 313]
plot_df = myfinal_df
plot_df2 = plot_df[plot_df['Type'] == 'S']


# In[23]:


figS2b_df = plot_df2.groupby(['Allele']).size().to_frame().reset_index()
figS2b_df.columns = ['Allele', 'SNP']
figS2b_df.sort_values('SNP', inplace=True, ascending=False)
figS2b_df.to_csv('%s/FigS2B.TCR.SNPs.txt'%(outdir2), sep='\t', index=False)


# In[31]:


ref_aa_dict = {}
for rec in SeqIO.parse(nogap_ref, 'fasta'):
    ref_aa_dict[rec.id] = str(rec.seq.translate())

novel_seq_dict = dict(zip(novel_df['Alleles'], novel_df['Seq']))
novel_aa_dict = {}
for k,v in novel_seq_dict.items():
    novel_aa_dict[k] = str(Seq(v).translate())


# In[32]:


def CalAApos(pos):
    codon_index = (int(pos-1) // 3)+1
    codon_pos = int(pos-1) % 3
    return pd.Series([codon_index, codon_pos])


# In[33]:


plot_df2[['AA_pos', 'Codon_pos']] = plot_df2.apply(lambda x:CalAApos(x['New_pos2']), axis=1)
plot_df2[['Old_AA_pos', 'Old_Codon_pos']] = plot_df2.apply(lambda x:CalAApos(x['Pos']+1), axis=1)
plot_df2['Novel_aa'] = plot_df2.apply(lambda x:novel_aa_dict[x['Allele']][x['Old_AA_pos']-1], axis=1)
plot_df2['Ref_aa'] = plot_df2.apply(lambda x:ref_aa_dict[x['Gene']][x['Old_AA_pos']-1], axis=1)


# In[34]:


def JudgeSame(a1, a2):
    if a1 == a2:
        return "No"
    else:
        return "SAP"


# In[43]:


plot_df2['Judge'] = plot_df2.apply(lambda x:JudgeSame(x['Novel_aa'], x['Ref_aa']), axis=1)
plot_df2.to_csv('%s/Summary.%s.SNP.SAP.result.txt'%(outdir2, chain), sep='\t', index=False)


# In[41]:


SNP_df = plot_df2['New_pos2'].value_counts().to_frame().reset_index()
SNP_df.columns = ['Pos', 'nSNP']
SNP_dict = dict(zip(SNP_df['Pos'], SNP_df['nSNP']))
SNP_data = []
for n in range(1, 313):
    if n in SNP_dict.keys():
        SNP_data.append([n, SNP_dict[n]])
    else:
        SNP_data.append([n, 0])
SNP_result_df = pd.DataFrame(SNP_data, columns = ['Position', 'nSNP'])
SNP_result_df.to_csv('%s/FigS3C.%s.PlotSNP.txt'%(outdir2, chain), sep='\t', index=False)