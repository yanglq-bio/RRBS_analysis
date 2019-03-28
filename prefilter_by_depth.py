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

for i in range(0,len(sample)):
    df = pd.read_csv(files[i],sep='\t')
    cols = sample[i]+'_Depth'
    df1 = df[df[cols]>=5]
    df2=df1[ ~ df1['Chr'].str.contains('_')]
    prefix = sample[i]+'_prefilter.txt'
    df2.to_csv(prefix,index = None)
    