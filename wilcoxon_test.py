#!~/miniconda2/envs/py3/bin/python
import pandas as pd
import scipy.stats as ss

pd.options.mode.chained_assignment = None
df = pd.read_csv('least_4_depth_10.txt',sep=',')


COL = df.columns.tolist()
Pvalue = []
met = []
for i in range(3,len(COL),4):
    met.append(COL[i])
    normal=[s for s in met if 'N' in s]
    normal.sort()
    tumor=[s for s in met if 'T' in s]
    tumor.sort()
for index, row in df.iterrows():
    N = row[normal].dropna().tolist()
    T = row[tumor].dropna().tolist()      
    p = ss.ranksums(N, T)[1].round(3)
    Pvalue.append(p)
df['p_value'] = Pvalue

df.to_csv('least_4_depth_10_WilcoxonTest.txt',index = None,na_rep='NA')

df1 = df[df['p_value']<=0.05]
df1.to_csv('least_4_depth_10_WilcoxonTest_filter.txt',index = None,na_rep='NA')


