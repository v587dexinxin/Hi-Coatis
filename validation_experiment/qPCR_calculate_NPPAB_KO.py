# -*- coding: utf-8 -*-
"""
Created on Sat Sep 14 18:52:31 2024

@author: lenovo
"""

import pandas as pd
import numpy as np



samples = ['WT' , 'B2' , '2' , '9' , '13']
genes = ['ACTB' , 'NPPA' , 'NPPA-AS1_1' , 'NPPA-AS1_2' , 'NPPB' , 'SWAP70']
genes_raw = ['ACTB' , 'NPPA' , 'NPPA-AS1_1' , 'NPPA-AS1_2' , 'NPPB' , 'SWAP70']
#



data1 = pd.read_excel('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2_敲除\\qPCR\\2024-09-14 165256_NPPAB_KO.xls' , sheet_name='Results' , skiprows=46)

data1.dropna(subset=['CT'], inplace=True)

data = data1.copy()


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
        diff[j , i * 3] =  tmp2['CT'].iloc[0] - ave
        diff[j , i * 3 + 1] =  tmp2['CT'].iloc[1] - ave
        diff[j , i * 3 + 2] =  tmp2['CT'].iloc[2] - ave



result = 2**-diff















































