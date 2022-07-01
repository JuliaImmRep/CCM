import pandas as pd
from scipy import stats
import numpy as np

f1 = pd.read_csv('example-BCR-SHM-only-isotype-%sclone.txt' % 'all',
               sep = '\t')
f1['Group'] = 'All'
f2 = pd.read_csv('example-BCR-SHM-only-isotype-%sclone.txt' % 'public',
               sep = '\t')
f2['Group'] = 'Public'

f = pd.concat([f1,f2])

pData = []
isoList = ['IgM','IgD','IgG','IgA','IgE']
for i1 in range(5):
    d1 = list(f[(f['Isotype'] == isoList[i1]) & (f['Group'] == 'All')]['SHM'])
    d2 = list(f[(f['Isotype'] == isoList[i1]) & (f['Group'] == 'Public')]['SHM'])
    levene = stats.levene(d1,d2)[1]
    if levene > 0.05:
        pvalue = stats.ttest_ind(d1,d2, equal_var = True)[1]
    else:
        pvalue = stats.ttest_ind(d1,d2, equal_var = False)[1]
    if pvalue <= 0.05:
        D1 = np.var(d1)
        D2 = np.var(d2)
        Cohensd = (np.mean(d1)-np.mean(d2))/((((len(d1)-1)*D1 + (len(d2)-1)*D2)/(len(d1)+len(d2)))**0.5)
        pData.append([isoList[i1], len(d1), len(d2), pvalue, round(abs(Cohensd), 3)])
df = pd.DataFrame(pData, columns = ['Isotype','nAll','nPublic','P value', "Effect size (Cohen's d)"])
df.to_csv('example-Figure6B-all-public-clone-effect-size.csv', index = False)
