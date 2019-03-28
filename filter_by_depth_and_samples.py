#!~/miniconda2/envs/py3/bin/python

import pandas as pd
import argparse
import scipy.stats as ss
import time

'''
The script does the following with the input file:
    
1.All CPG sites were measured at least in the specified number of samples, both in normal and tumor samples
2.At each site, the average sequencing depth of all samples is at least the specified value
3.wilcoxon test p_value <= 0.05

'''

parser = argparse.ArgumentParser(description="purpose:filter this sites by average depth and samples number")
parser.add_argument('-i','--input',required = True,
                    help='input file[txt format]')
parser.add_argument('-s','--samples',required = True,type=int,
                    help='the numbers of the sample')
parser.add_argument('-outpre','--outprefix',required = True,
                    help='Prefix for output files.')
parser.add_argument('-d',"--depth",type=int,default=10,
                    help="The average depth at raw_sample_merge do not less than this number[10]")
parser.add_argument('-t',"--threshold", type=int,default=5,
                    help="The number of samples contained the same cpg site is not less than this threshold[5]")
args = parser.parse_args()

pd.options.mode.chained_assignment = None
start =time.process_time()
df = pd.read_csv(args.input,sep=',',iterator=True)
loop = True
chunkSize = 1000000
chunks = []
DepThreshold = args.depth
SampleThreshold = args.threshold


############## define 3 functions to filter
def depth_caculate(df,depth):
    COL = df.columns.tolist()
    depth_col = [s for s in COL if 'Depth' in s]
    dff = df[df[depth_col].mean(axis=1,skipna=True)>depth]
    return(dff)

def wilcoxon_test(df):
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
    return(df)

def sample_threshold(df,threshold):
    COL = df.columns.tolist()
    normal=[s for s in COL if 'N_Methy' in s]
    tumor=[s for s in COL if 'T_Methy' in s]
    N_NA = df[normal].isnull().sum(axis=1).tolist()
    T_NA = df[tumor].isnull().sum(axis=1).tolist()
    c={"N_NA" : N_NA,"T_NA" : T_NA}
    nana=pd.DataFrame(c)
    df1 = pd.concat([df,nana],axis=1)
    t1 = len(normal) - threshold
    t2 = len(tumor) - threshold
    df2 = df1[(df1['N_NA'] <=t1) & (df1['T_NA'] <=t2)]
    df3 = df2.drop(['N_NA', 'T_NA'], axis=1)
    return(df3)

############### main function
while loop:
    try:
        chunk = df.get_chunk(chunkSize)
        chunk = chunk[ ~ chunk['Chr'].str.contains('_')]
        chunk_1 = sample_threshold(chunk,SampleThreshold)
        chunk_2 = depth_caculate(chunk_1,DepThreshold)
        chunk_3 = wilcoxon_test(chunk_2)
        print('already handled one chunk!')
        chunks.append(chunk_3)
        print('============================')
    except StopIteration:
        loop = False
        print("Iteration is stopped.")

df1 = pd.concat(chunks, ignore_index=True,axis=0)
end =time.process_time()
print('all chunks has been concat.')
print('total Running time: %s Seconds'%(end-start))
print('================================')
prefix = args.outprefix + '.txt'
df1.to_csv(prefix,index = None,na_rep='NA',chunksize=20000)
prefix_1 = args.outprefix + '_filter_by_pvalue.txt'
df2 = df1[df1['p_value']<=0.05]
df2.to_csv(prefix_1,index = None,na_rep='NA',chunksize=20000)
    
##################### get simplify dataframe
COL = df1.columns.tolist()
normal=[s for s in COL if 'N_Methy' in s]
normal.sort()
title = ['End','Start','Chr'] 
for item in title:
    normal.insert(0,item)
tumor=[s for s in COL if 'T_Methy' in s]
tumor.sort()
tumor.append('p_value')

prefix_2 = args.outprefix + '_simplify.txt'
N = df1[normal]
T = df1[tumor]
df3 = pd.concat([N,T],axis=1)   
df3.to_csv(prefix_2,index = None,na_rep='NA',chunksize=10000,sep='\t')

prefix_3 = args.outprefix + '_filter_by_pvalue_simplify.txt'
N = df2[normal]
T = df2[tumor]
df4 = pd.concat([N,T],axis=1)   
df4.to_csv(prefix_3,index = None,na_rep='NA',chunksize=10000,sep='\t')
