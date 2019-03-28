#!~/miniconda2/envs/py3/bin/python
import pandas as pd
import glob

# get all the file normal and tumor
normal = glob.glob("*N_prefilter.txt")
normal.sort()
tumor = glob.glob("*T_prefilter.txt")
tumor.sort()

# paired sample merge and filter sites which not simultaneously appear in normal and tumor
for i in range(0,len(normal)):
    n1 = pd.read_csv(normal[i],sep=',')
    t1 = pd.read_csv(tumor[i],sep=',')
    df = pd.merge(n1,t1,how='inner')
    prefix = tumor[i].split('_')[0] + '_' + normal[i].split('_')[0] + '.txt'
    df.to_csv(prefix,index = None,sep='\t')

