# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 10:02:12 2024

@author: lenovo
"""

import pandas as pd


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


RNA = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)

RNA['HCT116_FPKM'] = (RNA['HCT116_WT_R1_FPKM'] + RNA['HCT116_WT_R2_FPKM']) / 2

RNA = RNA.sort_values(by = ['HCT116_FPKM'] , ascending=False)

RNA = RNA[RNA['Chr'] != 'chrM']




loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\loops\\HCT116_merged6_0.1FA.hg38_loops_two_anchors_binding_common_peaks.bedpe' , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2']


peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak' , header = None , usecols=(0 , 1 , 2 , 3))
peaks.columns = ['chr' , 'start' , 'end' , 'peak_name']


promoter = []
for i in RNA.index:
    g = RNA.loc[i]['Chr']
    strand = RNA.loc[i]['Strand']
    name = RNA.loc[i]['Gene_Name']
    if strand == '+':
        start = RNA.loc[i]['Start'] - 2000
        end = RNA.loc[i]['Start'] 
    elif strand == '-':
        start = RNA.loc[i]['End']
        end = RNA.loc[i]['End'] + 2000
    else:
        print (i)
    promoter.append((g , start , end , name))
    
promoter = pd.DataFrame(promoter)
promoter.columns = ['chr' , 'start' , 'end' , 'gene_name']


peaks_intergenic = []
n = 0
for g in chrom:
     tmp_pro = promoter[promoter['chr'] == g]
     tmp_peaks = peaks[peaks['chr'] == g]
     for i in tmp_peaks.index:
         start = tmp_peaks.loc[i]['start']
         end = tmp_peaks.loc[i]['end']
         name = tmp_peaks.loc[i]['peak_name']
         mask = (tmp_pro['start'] <= end) & (tmp_pro['end'] >= start)
         overlap = tmp_pro[mask]
         if overlap.size != 0:
             pass
             n += 1
         else:
             # print (i)
             peaks_intergenic.append((g , start , end , name))
         
peaks_intergenic = pd.DataFrame(peaks_intergenic)         
peaks_intergenic.columns = ['chr' , 'start' , 'end' , 'peak_name']


loops_intergenic = pd.DataFrame([])
for g in chrom:
    print(g)
    tmp_loops = loops[loops['chr1'] == g]
    tmp_peaks = peaks_intergenic[peaks_intergenic['chr'] == g]
    for i in tmp_loops.index:
        start1 = tmp_loops.loc[i]['start1']
        end1 = tmp_loops.loc[i]['end1']
        start2 = tmp_loops.loc[i]['start2']
        end2 = tmp_loops.loc[i]['end2']
        overlap1 = tmp_peaks[(tmp_peaks['start'] <= end1) & (tmp_peaks['end'] >= start1)]
        overlap2 = tmp_peaks[(tmp_peaks['start'] <= end2) & (tmp_peaks['end'] >= start2)]
        if len(overlap1) > 1 or len(overlap2) > 1:
            loops_intergenic = pd.concat([loops_intergenic , tmp_loops.loc[i:i]] , axis = 0)
        

loops_intergenic = loops_intergenic.drop_duplicates()





peaks_name = ['HCT116_R2_q0.05_peak_704']

peaks_name = ['HCT116_R2_q0.05_peak_2366' , 'HCT116_R2_q0.05_peak_2367']

peaks_name = ['HCT116_R2_q0.05_peak_3125']

overlap = pd.DataFrame([])

for peak_name in peaks_name:
    peak = peaks_intergenic[peaks_intergenic['peak_name'] == peak_name]
    g = peak.iloc[0]['chr']
    start = peak.iloc[0]['start']
    end = peak.iloc[0]['end']
    tmp_loops = loops[loops['chr1'] == g]
    mask1 = (tmp_loops['start1'] <= end) & (tmp_loops['end1'] >= start)
    mask2 = (tmp_loops['start2'] <= end) & (tmp_loops['end2'] >= start)
    overlap1 = tmp_loops[mask1]
    overlap2 = tmp_loops[mask2]
    if (len(overlap1) != 0):
        overlap = pd.concat([overlap , overlap1])
    if (len(overlap2) != 0):
        overlap = pd.concat([overlap , overlap2])    
    overlap = overlap.drop_duplicates()
    overlap.to_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\' + peak_name + '.bedpe' , header=None , index = None , sep = '\t')
    
 
for i in overlap.index:
    g = overlap.loc[i]['chr1']
    start1 = overlap.loc[i]['start1']
    end1 = overlap.loc[i]['end1']
    start2 = overlap.loc[i]['start2']
    end2 = overlap.loc[i]['end2']
    tmp_promoter = promoter[promoter['chr'] == g]
    mask1 = (tmp_promoter['start'] <= end1) & (tmp_promoter['end'] >= start1)
    mask2 = (tmp_promoter['start'] <= end2) & (tmp_promoter['end'] >= start2)
    




    
gene_names = ['ACTG1' , 'RPLP0']

gene_names = ['RPLP0']




overlap = pd.DataFrame([])
for gene_name in gene_names:
    gene = RNA[RNA['Gene_Name'] == gene_name]
    g = gene.iloc[0]['Chr']
    start = gene.iloc[0]['Start']
    end = gene.iloc[0]['End']
    strand = gene.iloc[0]['Strand']
    if strand == '+':
        pro_s = start - 2000
        pro_e = start + 1000
    else:
        pro_s = end - 1000
        pro_e = end + 2000
    tmp_loops = loops[loops['chr1'] == g]
    mask1 = (tmp_loops['start1'] <= pro_e) & (tmp_loops['end1'] >= pro_s)
    mask2 = (tmp_loops['start2'] <= pro_e) & (tmp_loops['end2'] >= pro_s)
    overlap1 = tmp_loops[mask1]
    overlap2 = tmp_loops[mask2]
    if (len(overlap1) != 0):
        overlap = pd.concat([overlap , overlap1])
    if (len(overlap2) != 0):
        overlap = pd.concat([overlap , overlap2])
    overlap.to_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\' + gene_name + '.bedpe' , header=None , index = None , sep = '\t')
    
    
    
















