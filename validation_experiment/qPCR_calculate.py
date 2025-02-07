# -*- coding: utf-8 -*-
"""
Created on Fri May 24 19:46:13 2024

@author: lenovo
"""

import pandas as pd
import numpy as np



########Samples

samples = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\\qPCR\\WDSUB1_samples.txt' , header = 0 , usecols=[0])
genes = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\\qPCR\\WDSUB1_samples.txt' , header = 0 , usecols=[1])



samples = list(samples.dropna()['samples'])
genes = list(genes.dropna()['genes'])
genes = [str(x) for x in genes]


#######qPCR data

data = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\qPCR\\WDSUB1.txt' , header=0 , skiprows=1)

row_to_move = data.loc[90:].copy()

data = data.drop([42 , 43 , 44 , 45 , 46 , 47 , 90 , 91 , 92 , 93 , 94 , 95])

data = pd.concat([data.iloc[:54] , row_to_move , data.iloc[54:]]).reset_index(drop=True)


#####Remove NaN

data['Cp'].fillna(100, inplace=True)


#####Add Sample names


tmp_sample = []
for x in range(len(samples)):
    tmp_sample += [samples[x]] * 30


data['samples'] = tmp_sample




#####Add Gene names

tmp_genes = []
for x in range(len(samples)):
    for y in range(len(genes)):
        tmp_genes += [str(genes[y])] * 3
        


data['genes'] = tmp_genes



#####

ref_gene = 'ACTB'

ref_average = {}
for s in samples:
    tmp_s = data[data['samples'] == s]
    tmp_g = tmp_s[tmp_s['genes'] == ref_gene]
    ref_average[s] = np.mean(list(tmp_g['Cp']))
    
    
diff = np.zeros((len(genes) , len(samples) * 3))

for i in range(len(samples)):
    s = samples[i]
    ave = ref_average[samples[i]]
    for j in range(len(genes)):
        g = genes[j]
        tmp1 = data[data['samples'] == s]
        tmp2 = tmp1[tmp1['genes'] == genes[j]]
        diff[j , i * 3] =  tmp2['Cp'].iloc[0] - ave
        diff[j , i * 3 + 1] =  tmp2['Cp'].iloc[1] - ave
        diff[j , i * 3 + 2] =  tmp2['Cp'].iloc[2] - ave




result = 2**-diff





