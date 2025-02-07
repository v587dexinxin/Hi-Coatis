# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 10:02:12 2024

@author: lenovo
"""


import numpy as np
import pandas as pd
from itertools import islice
import os


def Load_peaks(file , peaks_type):
    sz = os.path.getsize(file)
    if sz != 0:
        peaks = pd.read_table(file , header = None)
        if peaks_type == 'narrow':
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue' , 'lengtn']
        else:
            peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signal' , 'pvalue' , 'qvalue']
    else:
        peaks = []
    
    return peaks



    
def merge_intervals(intervals):
    if not intervals:
        return []

    # 先按照区间的开始时间排序
    intervals.sort(key=lambda x: x[0])

    merged = [intervals[0]]
    for current in intervals:
        last = merged[-1]
        if current[0] <= last[1]:  # 如果当前区间的开始时间小于等于最后一个合并区间的结束时间
            # 合并区间
            merged[-1] = (last[0], max(last[1], current[1]))
        else:
            # 不重叠的区间，直接添加到结果列表中
            merged.append(current)
    
    return merged





#####################Load data######################

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


RNA = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)

RNA['HCT116_FPKM'] = (RNA['HCT116_WT_R1_FPKM'] + RNA['HCT116_WT_R2_FPKM']) / 2

RNA = RNA.sort_values(by = ['HCT116_FPKM'] , ascending=False)

RNA = RNA[RNA['Chr'] != 'chrM']




loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\loops\\HCT116_merged6_0.1FA.hg38_loops_one_anchor_binding_union_peaks.bedpe' , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'pvalue']


peaks_116 = Load_peaks('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak', 'narrow')



##################calculation ###################

promoter_genes = []
for i in RNA.index:
    g = RNA.loc[i]['Chr']
    strand = RNA.loc[i]['Strand']
    name = RNA.loc[i]['Gene_Name']
    fpkm = RNA.loc[i]['HCT116_FPKM']
    if strand == '+':
        start = RNA.loc[i]['Start'] - 2000
        end = RNA.loc[i]['End'] 
    elif strand == '-':
        start = RNA.loc[i]['Start']
        end = RNA.loc[i]['End'] + 2000
    else:
        print (i)
    promoter_genes.append((g , start , end , name , fpkm))
    
promoter_genes = pd.DataFrame(promoter_genes)
promoter_genes.columns = ['chr' , 'start' , 'end' , 'gene_name' , 'FPKM']




peaks_name = ['HCT116_R2_q0.05_peak_30906' , 'HCT116_R2_q0.05_peak_33685' , 'HCT116_R2_q0.05_peak_33688' , 
              'HCT116_R2_q0.05_peak_45512' , 'HCT116_R2_q0.05_peak_50480' , 'HCT116_R2_q0.05_peak_51595' , 
              'HCT116_R2_q0.05_peak_8151' , 'HCT116_R2_q0.05_peak_13234' , 'HCT116_R2_q0.05_peak_21636' , 
              'HCT116_R2_q0.05_peak_26807' , 'HCT116_R2_q0.05_peak_27136']




##########在这里改变peak的选择，可输出该peak通过loop靶向的基因##################


peak_name = peaks_name[0]


##############

n = 1
for peak_name in peaks_name:    
    overlap = pd.DataFrame([])
    selected_genes = pd.DataFrame([])
    peak = peaks_116[peaks_116['name'] == peak_name]
    g = peak.iloc[0]['chr']
    start = peak.iloc[0]['start']
    end = peak.iloc[0]['end']
    tmp_loops = loops[loops['chr1'] == g]
    mask1 = (tmp_loops['start1'] <= end) & (tmp_loops['end1'] >= start)
    mask2 = (tmp_loops['start2'] <= end) & (tmp_loops['end2'] >= start)
    overlap1 = tmp_loops[mask1]
    overlap2 = tmp_loops[mask2]
    # print (overlap3)
    if (len(overlap1) != 0):
        overlap = pd.concat([overlap , overlap1])
    if (len(overlap2) != 0):
        overlap = pd.concat([overlap , overlap2])    
    overlap = overlap.drop_duplicates()
    overlap.to_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks_related_genes_fasta\\116_selected_intergenic_peaks_related_loops\\' + peak_name + '_IF.bedpe' , header=None , index = None , sep = '\t')
    
    
    for i in overlap.index:
        g = overlap.loc[i]['chr1']
        start1 = overlap.loc[i]['start1']
        end1 = overlap.loc[i]['end1']
        start2 = overlap.loc[i]['start2']
        end2 = overlap.loc[i]['end2']
        tmp_promoter = promoter_genes[promoter_genes['chr'] == g]
        mask1 = (tmp_promoter['start'] <= end1) & (tmp_promoter['end'] >= start1)
        mask2 = (tmp_promoter['start'] <= end2) & (tmp_promoter['end'] >= start2)
        overlap1 = tmp_promoter[mask1]
        overlap2 = tmp_promoter[mask2]
        if len(overlap1) != 0:
            # print (overlap1)
            selected_genes = pd.concat([selected_genes , overlap1])
        if len(overlap2) != 0:
            # print (overlap2)
            selected_genes = pd.concat([selected_genes , overlap2])
    if len(selected_genes) != 0:
        selected_genes = selected_genes.drop_duplicates()        
        selected_genes = selected_genes.sort_values(by = 'FPKM' , ascending=False)        
        # selected_genes = selected_genes[selected_genes['FPKM'] >= 10]
        # gene_names = list(selected_genes['gene_name'])
        # selected_peak_genes[peak_name] = gene_names
        
        
        print (selected_genes)
        selected_genes.to_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks_related_genes_fasta\\116_selected_intergenic_peaks_related_loops\\116_selected_intergenic_peaks_related_genes\\KO' + str(n) + '_' + peak_name + 'related_genes.txt' ,  index = None , sep = '\t')
    n += 1
        
    
    
    
    
    
    
    
    
    























