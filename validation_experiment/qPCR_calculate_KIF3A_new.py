# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 11:00:21 2024

@author: lenovo
"""

import pandas as pd
import numpy as np



samples = ['WT' , 'dCas' , 'KIF1' , 'KIF2']
genes = ['ACTB' , 'TUBB' , 'KIF3A' , 'AFF4' , 'AP3S1' , 'CCNI2' , 'RAD50' , 'SEPTIN8']
genes_raw = ['ACTB' , 'TUBB' , 'KIF3A' , '1' , '2' , '3' , '4' , '5']
#



data = pd.read_excel('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\\qPCR\\new\\2024-09-12 124003_KIF3A.xls' , sheet_name='Results' , skiprows=46)

data.dropna(subset=['CT'], inplace=True)




ref_gene = 'ACTB'


ref_average = {}
for s in samples:
    tmp_s = data[data['Sample Name'] == s]
    tmp_g = tmp_s[tmp_s['Target Name'] == ref_gene]
    ref_average[s] = np.mean(list(tmp_g['CT']))
    
    
diff = np.zeros((len(genes) , len(samples) * 3))

for i in range(len(samples)):
    s = samples[i]
    ave = ref_average[samples[i]]
    for j in range(len(genes_raw)):
        g = genes[j]
        tmp1 = data[data['Sample Name'] == s]
        tmp2 = tmp1[tmp1['Target Name'] == genes_raw[j]]
        for k in range(3):
            if tmp2['CT'].iloc[k] == 'Undetermined':
                pass
            else:
                diff[j , i * 3 + k] =  tmp2['CT'].iloc[k] - ave
        


result = 2**-diff

