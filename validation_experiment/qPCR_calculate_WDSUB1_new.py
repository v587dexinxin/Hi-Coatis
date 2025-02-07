# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 10:17:55 2024

@author: lenovo
"""

import pandas as pd
import numpy as np



samples = ['WT' , 'dCase' , 'CK' , 'WDSUB1-1' , 'WDSUB1-2']
genes = ['ACTB' , 'TUBB' , 'WDSUB1' , 'BAZ2B_1' , 'CD302_2' , 'CFAP210_3' , 'GLS_4' , 'LY75-CD302_5' , 'MARCHF7_6' , 'PKP4_7' , 'STAT1_8' , 'TANC1_9']
genes_raw = ['ATCB' , 'TUBB' , 'WDSUB1' , 'WD1' , 'WD2' , 'WD3' , 'WD4' , 'WD5' , 'WD6' , 'WD7' , 'WD8' , 'WD9']
#



data1 = pd.read_excel('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\\qPCR\\new\\2024-07-24 210621_WDSUB_new.xls' , sheet_name='Results' , skiprows=46)
data2 = pd.read_excel('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\\qPCR\\new\\2024-07-24 222931_WDSUB1_PARD6B.xls' , sheet_name='Results' , skiprows=46)

data1.dropna(subset=['CT'], inplace=True)
data2.dropna(subset=['CT'], inplace=True)

data = pd.concat([data1 , data2])



ref_gene = 'ATCB'


ref_average = {}
for s in samples:
    tmp_s = data[data['Sample Name'] == s]
    tmp_g = tmp_s[tmp_s['Target Name'] == ref_gene]
    ref_average[s] = np.mean(list(tmp_g['CT']))
    
    
diff = np.zeros((len(genes) , len(samples) * 2))

for i in range(len(samples)):
    s = samples[i]
    ave = ref_average[samples[i]]
    for j in range(len(genes_raw)):
        g = genes[j]
        tmp1 = data[data['Sample Name'] == s]
        tmp2 = tmp1[tmp1['Target Name'] == genes_raw[j]]
        diff[j , i * 2] =  tmp2['CT'].iloc[0] - ave
        diff[j , i * 2 + 1] =  tmp2['CT'].iloc[1] - ave



result = 2**-diff















































