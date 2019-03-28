#!~/miniconda2/envs/py3/bin/python
import pandas as pd
import os

# get all the file at the path
def file_name(file_dir):
    L=[]
    S=[]
    for root, dirs,files in os.walk(file_dir):
        for file in files:
            if file.find('.met.txt') != -1:
                L.append(file.split('.')[0])
                S.append(os.path.join(root, file))
    return(L,S)
sample,files = file_name('./')

df1 = pd.read_csv(files[0],sep='\t')

for i in range(1,len(sample)):
    df2 = pd.read_csv(files[i],sep='\t')
    df1 = pd.merge(df1,df2,how='outer')
    
df1.to_csv('sample_merge.txt',index = None,na_rep='NA',chunksize=10000,compression='gzip')
