# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 15:21:55 2024

@author: lenovo
"""

from heapq import merge
from itertools import count, islice
# from contextlib2 import ExitStack
from matplotlib.backends.backend_pdf import PdfPages
from random import random
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# from palettable.colorbrewer.qualitative import Dark2_8
import os, sys, re, time, subprocess, multiprocessing, gc, bisect, math
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
from itertools import islice
import pyBigWig



    
def Get_bw_sig(bw_file , g , start , end):
    '''
    '''
    bw = pyBigWig.open(bw_file, 'r')
    values = bw.values(g, start, end)
    return values



def Gene_classification(signals , RNA_data):
    ave = {}
    for g in chrom:
        print (g)
        tmp = signals.values(g , 0 , signals.chroms()[g])
        ave[g] = np.mean(tmp)
        
    cl1 = pd.DataFrame([]) 
    cl2 = pd.DataFrame([])
    cl3 = pd.DataFrame([])
    cl4  = pd.DataFrame([])
    
    
    for i in RNA_data.index:
        if i % 1000 == 0:
            print (i)
        gene_name = RNA_data.loc[i]['Gene_Name']
        g = RNA_data.loc[i]['Chr']
        start = RNA_data.loc[i]['Start']
        end = RNA_data.loc[i]['End']
        strand = RNA_data.loc[i]['Strand']
        if end - start <= 400:
            print (gene_name)
            continue
        if strand == '+':
            pro_s = start - 200
            pro_e = start + 400
            body_s = start + 401
            body_e = end 
        else:
            pro_s = end - 400
            pro_e = end + 200
            body_s = start 
            body_e = end - 401        
            
        v_pro = signals.values(g , pro_s , pro_e)
        v_body = signals.values(g , body_s , body_e)
        
        ratio1 = np.mean(v_pro) / ave[g]
        ratio2 = np.mean(v_body) / ave[g]
        
        if (ratio1 > 3) and (ratio2 > 1):
            cl1 = pd.concat([cl1 , RNA_data.loc[[i]]])
        elif (ratio1 > 3) and (ratio2 <= 1):
            cl2 = pd.concat([cl2 , RNA_data.loc[[i]]])
        elif (ratio1 <= 3) and (ratio2 > 1):
            cl3 = pd.concat([cl3 , RNA_data.loc[[i]]])
        elif (ratio1 <= 3) and (ratio2 <= 1):
            cl4 = pd.concat([cl4 , RNA_data.loc[[i]]])
            
    cl1['score'] = np.zeros(len(cl1))
    cl2['score'] = np.zeros(len(cl2))
    cl3['score'] = np.zeros(len(cl3))
    cl4['score'] = np.zeros(len(cl4))
    
    return (cl1 , cl2 , cl3 , cl4)



def Box_plot_4cellline(data , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':1},
            medianprops={'color':'dodgerblue','linewidth':1},
            capprops={'color':'dodgerblue','linewidth':1},
            whiskerprops={'color':'dodgerblue','linewidth':1})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4 ])
    ax.set_xticklabels(['C1' , 'C2' , 'C3' , 'C4' ] , fontsize = 10)
    ax.set_ylabel('FPKM' , fontsize = 20)
    ax.set_xlabel('K562_HiRPC_peaks_classify')
    ax.set_xlim((0.5 , 4.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

####K562
RNA_data = pd.read_csv('/scratch/2024-08-26/bio-shenw/Ljniu/K562/RNA-seq/FPKM/union_all_FPKM.csv' , header = 0)
RNA_data['K562_FPKM'] = (RNA_data['K562_WT_R1_FPKM'] + RNA_data['K562_WT_R2_FPKM']) / 2
# RNA_data = RNA_data[RNA_data['K562_FPKM'] >= 0.5]
RNA_data = RNA_data[RNA_data['Chr'].isin(chrom)]
RNA_data = RNA_data[RNA_data['End'] - RNA_data['Start'] > 401]

signals = pyBigWig.open('/scratch/2025-02-01/bio-shenw/Ljniu/K562/K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw' , 'r')

cl1 , cl2 , cl3 , cl4 = Gene_classification(signals , RNA_data)
cl1[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2025-02-01/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/classify1_all.bed' , header=None , index = None , sep = '\t')
cl2[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2025-02-01/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/classify2_all.bed' , header=None , index = None , sep = '\t')
cl3[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2025-02-01/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/classify3_all.bed' , header=None , index = None , sep = '\t')
cl4[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2025-02-01/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/classify4_all.bed' , header=None , index = None , sep = '\t')

fig = Box_plot_4cellline([cl1['K562_FPKM'] , cl2['K562_FPKM'] , cl3['K562_FPKM'] , cl4['K562_FPKM']] , 0 , 50)
run_Plot(fig , '/scratch/2024-08-26/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/K562_HiRPC_signal_classify_RNA_intensity.pdf')


####HCT116
RNA_data = pd.read_csv('/scratch/2024-08-26/bio-shenw/Ljniu/K562/RNA-seq/FPKM/union_all_FPKM.csv' , header = 0)
RNA_data['HCT116_FPKM'] = (RNA_data['HCT116_WT_R1_FPKM'] + RNA_data['HCT116_WT_R2_FPKM']) / 2
# RNA_data = RNA_data[RNA_data['K562_FPKM'] >= 0.5]
RNA_data = RNA_data[RNA_data['Chr'].isin(chrom)]
RNA_data = RNA_data[RNA_data['End'] - RNA_data['Start'] > 401]

signals = pyBigWig.open('/scratch/2024-08-26/bio-shenw/Ljniu/HCT116/hg38/one-dimensional_1/mapping/all_reps/HCT116_0.1FA_merged2_RPKM_10bp.bw' , 'r')

cl1 , cl2 , cl3 , cl4 = Gene_classification(signals , RNA_data)

cl1[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-08-26/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_HCT116/classify1_all.bed' , header=None , index = None , sep = '\t')
cl2[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-08-26/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_HCT116/classify2_all.bed' , header=None , index = None , sep = '\t')
cl3[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-08-26/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_HCT116/classify3_all.bed' , header=None , index = None , sep = '\t')
cl4[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-08-26/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_HCT116/classify4_all.bed' , header=None , index = None , sep = '\t')

fig = Box_plot_4cellline([cl1['HCT116_FPKM'] , cl2['HCT116_FPKM'] , cl3['HCT116_FPKM'] , cl4['HCT116_FPKM']] , 0 , 40)
run_Plot(fig , '/scratch/2024-08-26/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_HCT116/HCT116_HiRPC_signal_classify_RNA_intensity.pdf')



####K562 Hemin
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
RNA_data = pd.read_csv('/scratch/2024-09-02/bio-shenw/Ljniu/K562/RNA-seq/FPKM/union_all_FPKM_K562_Hemin.csv' , header = 0)
RNA_data['Hemin_FPKM'] = (RNA_data['K562_H_R1_FPKM'] + RNA_data['K562_H_R2_FPKM']) / 2
# RNA_data = RNA_data[RNA_data['K562_FPKM'] >= 0.5]
RNA_data = RNA_data[RNA_data['Chr'].isin(chrom)]
RNA_data = RNA_data[RNA_data['End'] - RNA_data['Start'] > 401]

signals = pyBigWig.open('/scratch/2024-09-02/bio-shenw/Ljniu/K562/K562_H_New/one-dimensional/mapping/bam2/all_reps/K562_Hemin_New_allreps_RPKM_10bp.bw' , 'r')

cl1 , cl2 , cl3 , cl4 = Gene_classification(signals , RNA_data)

cl1[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-09-02/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_Hemin/classify1_all.bed' , header=None , index = None , sep = '\t')
cl2[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-09-02/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_Hemin/classify2_all.bed' , header=None , index = None , sep = '\t')
cl3[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-09-02/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_Hemin/classify3_all.bed' , header=None , index = None , sep = '\t')
cl4[['Chr' , 'Start' , 'End' , 'Gene_Name' , 'score' , 'Strand']].to_csv('/scratch/2024-09-02/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_Hemin/classify4_all.bed' , header=None , index = None , sep = '\t')


fig = Box_plot_4cellline([cl1['Hemin_FPKM'] , cl2['Hemin_FPKM'] , cl3['Hemin_FPKM'] , cl4['Hemin_FPKM']] , 0 , 40)
run_Plot(fig , '/scratch/2024-09-02/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_Hemin/K562_Hemin_HiRPC_signal_classify_RNA_intensity.pdf')







