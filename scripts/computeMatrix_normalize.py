# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 13:14:39 2023

@author: lenovo
"""

import numpy as np
import pandas as pd


####-------------------------TSS_TES-----------------------------------

data = pd.read_csv('TSS_TES_plotMatrix' , header = None , sep = '\t' , skiprows=1)

data = data.fillna(0)


    
data_genes = data.iloc[:,:6]
sample1 = data.iloc[:,6:1206]
sample2 = data.iloc[:,1206:2406]
sample3 = data.iloc[:,2406:3606]
sample4 = data.iloc[:,3606:4806]


average1 = sample1.mean().mean()
average2 = sample2.mean().mean()
average3 = sample3.mean().mean()
average4 = sample4.mean().mean()


sample2 = sample2 * (100 / average2)
sample3 = sample3 * (100 / average3)
sample4 = sample4 * (100 / average4)


data_new = pd.concat([data_genes , sample1 , sample2 , sample3 , sample4] , axis = 1)
data_new.to_csv('TSS_TES_plotMatrix_norm_1' , header=None , index = None , sep = '\t')




####-------------------------union-----------------------------------

data = pd.read_csv('K562_0.1FA_union' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)


    
data_genes = data.iloc[:,:6]
sample1 = data.iloc[:,6:406]
sample2 = data.iloc[:,406:806]
sample3 = data.iloc[:,806:1206]
sample4 = data.iloc[:,1206:1606]


average1 = sample1.mean().mean()
average2 = sample2.mean().mean()
average3 = sample3.mean().mean()
average4 = sample4.mean().mean()


sample2 = sample2 * (average1 / average2)
sample3 = sample3 * (average1 / average3)
sample4 = sample4 * (average1 / average4)


data_new = pd.concat([data_genes , sample1 , sample2 , sample3 , sample4] , axis = 1)
data_new.to_csv('K562_0.1FA_union_normrpc' , header=None , index = None , sep = '\t')


####-------------------------union-----------------------------------

data = pd.read_csv('K562_0.1FA_VS_ChIP_union' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)


    
data_genes = data.iloc[:,:6]
sample1 = data.iloc[:,6:406]
sample2 = data.iloc[:,406:806]
sample3 = data.iloc[:,806:1206]
sample4 = data.iloc[:,1206:1606]
sample5 = data.iloc[:,1606:2006]
sample6 = data.iloc[:,2006:2406]
sample7 = data.iloc[:,2406:2806]
sample8 = data.iloc[:,2806:3206]
sample9 = data.iloc[:,3206:3606]



average1 = sample1.mean().mean()
average2 = sample2.mean().mean()
average3 = sample3.mean().mean()
average4 = sample4.mean().mean()
average5 = sample5.mean().mean()
average6 = sample6.mean().mean()
average7 = sample7.mean().mean()
average8 = sample8.mean().mean()
average9 = sample9.mean().mean()


sample2 = sample2 * (average1 / average2)
sample3 = sample3 * (average1 / average3)
sample4 = sample4 * (average1 / average4)
sample5 = sample5 * (average1 / average5)
sample6 = sample6 * (average1 / average6)
sample7 = sample7 * (average1 / average7)
sample8 = sample8 * (average1 / average8)
sample9 = sample9 * (average1 / average9)



data_new = pd.concat([data_genes , sample1 , sample2 , sample3 , sample4 , sample5 , sample6 , sample7 , sample8 , sample9 ] , axis = 1)
data_new.to_csv('K562_0.1FA_VS_ChIP_union_normrpc_1' , header=None , index = None , sep = '\t')




####-------------------------HiRPC_classify-----------------------------------


import pandas as pd


data = pd.read_csv('TSS_TES_plotMatrix' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)


n =  len(data.columns) // 1200

data_genes = data.iloc[:,:6]


data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * 1200 + 6) : (i + 1) * 1200 + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    


data_new.to_csv('TSS_TES_plotMatrix_norm_1' , header=None , index = None , sep = '\t')




####-------------------------K562_lowcells-----------------------------------

data = pd.read_csv('/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig7A_K562_lowcells_peaks_heatmap/new/K562_lowcells_union_1M' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)

bins_num = 400

n =  len(data.columns) // bins_num

data_genes = data.iloc[:,:6]



data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    
    
    
data_new.to_csv('/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig7A_K562_lowcells_peaks_heatmap/new/K562_lowcells_union_1M_norm_1' , header=None , index = None , sep = '\t')




####-------------------------K562_lowcells_Z-score-----------------------------------

# data = pd.read_csv('/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig7A_K562_lowcells_peaks_heatmap/new/K562_lowcells_union_1M' , header = None , sep = '\t' , skiprows=1)


# data = data.fillna(0)

# data_genes = data.iloc[:,:6]



# data_new = pd.DataFrame([])
# data_new = pd.concat([data_new , data_genes])
# sample = data.iloc[: , 6:]
# average = np.mean(sample)
# std = np.std(sample)
# z_scores = (sample - average) / std
# data_new =  pd.concat([data_new , z_scores] , axis=1)

    
    
# data_new.to_csv('/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig7A_K562_lowcells_peaks_heatmap/new/Z-score/K562_lowcells_union_1M_norm_z-score_1' , header=None , index = None , sep = '\t')



####-------------------------K562_lowcells_min-max-----------------------------------


# from sklearn.preprocessing import MinMaxScaler
# import numpy as np

# data = pd.read_csv('/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig7A_K562_lowcells_peaks_heatmap/new/K562_lowcells_union_1M' , header = None , sep = '\t' , skiprows=1)


# data = data.fillna(0)

# bins_num = 400

# n =  len(data.columns) // bins_num

# data_genes = data.iloc[:,:6]



# data_new = pd.DataFrame([])
# data_new = pd.concat([data_new , data_genes])
# for i in range(n):
#     sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
#     scaler = MinMaxScaler(feature_range=(0, 10))
#     normalized_data = scaler.fit_transform(sample)
#     normalized_data = pd.DataFrame(normalized_data)

#     data_new =  pd.concat([data_new , normalized_data] , axis=1)
    
    
    
# data_new.to_csv('/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig7A_K562_lowcells_peaks_heatmap/new/min_max/K562_lowcells_union_1M_norm_min_max_1' , header=None , index = None , sep = '\t')






####-------------------------K562_ChIP_union-----------------------------------

data = pd.read_csv('/scratch/2024-12-09/bio-shenw/Ljniu/K562/plots/Fig2G_RPC_VS_ChIPs/new/new_1/HMM_classify/K562_0.1FA_VS_ChIP_union_classify' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)

bins_num = 400

n =  len(data.columns) // bins_num

data_genes = data.iloc[:,:6]



data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    
    
    
data_new.to_csv('/scratch/2024-12-09/bio-shenw/Ljniu/K562/plots/Fig2G_RPC_VS_ChIPs/new/new_1/HMM_classify/K562_0.1FA_VS_ChIP_union_classify_norm' , header=None , index = None , sep = '\t')






####-------------------------K562_Repetitive_CNV---------------------------------

data = pd.read_csv('/scratch/2024-12-30/bio-shenw/Ljniu/K562/plots/plots_New_number_bylxx_xjs/one-dimensional-peaks_Anno/Peaks_Anno_bed/Repetitive_CNV_peaks/K562_Repetitive_CNV_union' , header = None , sep = '\t' , skiprows=1)


data = data.fillna(0)

bins_num = 400

n =  len(data.columns) // bins_num

data_genes = data.iloc[:,:6]



data_new = pd.DataFrame([])
data_new = pd.concat([data_new , data_genes])
for i in range(n):
    sample = data.iloc[: , (i * bins_num + 6) : (i + 1) * bins_num + 6]
    average = sample.mean().mean()
    sample = sample * (100 / average)
    data_new =  pd.concat([data_new , sample] , axis=1)
    
    
    
data_new.to_csv('/scratch/2024-12-30/bio-shenw/Ljniu/K562/plots/plots_New_number_bylxx_xjs/one-dimensional-peaks_Anno/Peaks_Anno_bed/Repetitive_CNV_peaks/K562_Repetitive_CNV_union_norm' , header=None , index = None , sep = '\t')














