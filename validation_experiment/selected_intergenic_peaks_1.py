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


def Get_SNP_peaks(SNPs , peaks):
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    n = 0; tmp = pd.DataFrame()
    for g in chrom:
        tmp_snps = SNPs[SNPs['CHR_ID'] == g.lstrip('chr')]
        tmp_peaks = peaks[peaks['chr'] == g]
        
        for i in tmp_snps.index:
            pos = tmp_snps.loc[i]['CHR_POS']
            
            
            overlap = tmp_peaks[(tmp_peaks['start'] <= pos + 100) & (tmp_peaks['end'] >= pos - 100)]
            if len(overlap) != 0:
                n += 1
                tmp = pd.concat([tmp , overlap])
    
    
    tmp = tmp.drop_duplicates()
    return tmp


def read_genome(filename):
    file_genome = open(filename)
    dict_genome = {}
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>' :
            chrs = (line.split('>')[1]).split()[0]
            dict_genome[chrs] = []
        
        else :
            dict_genome[chrs].extend(lists)
    return dict_genome


def write_pos_to_seq(CDS , outfil):
    
    out = open(outfil , 'w')
    
    for i in CDS:
        out.writelines('>' + i + '_CDS' + '\n')
        for j in range(len(CDS[i]) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(CDS[i][n1 : n2]).upper() + '\n')
    out.close()
    
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






genome_h38 = read_genome('H:\\work\\literature_data\\genome\\hg38\\hg38.fa')

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


RNA = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\FPKM\\union_all_FPKM.csv' , header = 0)

RNA['HCT116_FPKM'] = (RNA['HCT116_WT_R1_FPKM'] + RNA['HCT116_WT_R2_FPKM']) / 2

RNA = RNA.sort_values(by = ['HCT116_FPKM'] , ascending=False)

RNA = RNA[RNA['Chr'] != 'chrM']




loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\loops\\HCT116_merged6_0.1FA.hg38_loops_one_anchor_binding_loops_IF.bedpe' , header = None , usecols= (0 , 1 , 2 , 3 , 4 , 5))
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2']


SNPs = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\CAD_related_SNPs_LD0.99_all_risk_allel_sort_seqname.csv' , header = 0)

peaks_huvec = Load_peaks('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\one-dimensional\\HUVEC_control_ls_onedimensional_q0.05_union2_peaks.narrowPeak', 'narrow')
peaks_116 = Load_peaks('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak', 'narrow')



huvec_snp_peaks = Get_SNP_peaks(SNPs , peaks_huvec)
hct116_snp_peaks = Get_SNP_peaks(SNPs , peaks_116)

union_snp_peaks = pd.concat([huvec_snp_peaks , hct116_snp_peaks])
union_snp_peaks = union_snp_peaks.drop_duplicates()
union_snp_peaks = union_snp_peaks.sort_values(by = ['chr' , 'start' , 'end'])



CVD_genes = pd.read_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Genes associated with cardiovascular disease.csv' , header = 0)
CVD_genes = CVD_genes.values.flatten().tolist()
CVD_genes = pd.Series(CVD_genes).dropna().tolist()


promoter = []
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
    promoter.append((g , start , end , name , fpkm))
    
promoter = pd.DataFrame(promoter)
promoter.columns = ['chr' , 'start' , 'end' , 'gene_name' , 'FPKM']




peaks_name = list(union_snp_peaks['name'])

peaks_name = ['HCT116_R2_q0.05_peak_14557']
peaks_name = ['HCT116_R2_q0.05_peak_9402']
peaks_name = ['control-2_S18_q0.05_peak_17977']
peaks_name = ['HCT116_R2_q0.05_peak_21354']
peaks_name = ['HCT116_R2_q0.05_peak_23077']
peaks_name = ['HCT116_R2_q0.05_peak_24379']
peaks_name = ['HCT116_R2_q0.05_peak_28484']
peaks_name = ['HCT116_R2_q0.05_peak_28682']
peaks_name = ['HCT116_R2_q0.05_peak_28703']
peaks_name = ['HCT116_R2_q0.05_peak_33276']
peaks_name = ['ls-2_S19_p0.005_peak_27930']
peaks_name = ['HCT116_R2_q0.05_peak_36729']
peaks_name = ['HCT116_R2_q0.05_peak_38139']
peaks_name = ['HCT116_R2_q0.05_peak_40278']
peaks_name = ['HCT116_R2_q0.05_peak_50057']
peaks_name = ['control-2_S18_q0.05_peak_49453']
peaks_name = ['HCT116_R2_q0.05_peak_56432']
peaks_name = ['HCT116_R2_q0.05_peak_57564']
peaks_name = ['HCT116_R2_q0.05_peak_57881']
peaks_name = ['HCT116_R2_q0.05_peak_2394']
peaks_name = ['HCT116_R2_q0.05_peak_2395']



# selected_peak_genes = {}

for peak_name in peaks_name:
    overlap = pd.DataFrame([])
    selected_genes = pd.DataFrame([])
    peak = union_snp_peaks[union_snp_peaks['name'] == peak_name]
    g = peak.iloc[0]['chr']
    start = peak.iloc[0]['start']
    end = peak.iloc[0]['end']
    tmp_loops = loops[loops['chr1'] == g]
    tmp_snps = SNPs[SNPs['CHR_ID'] == g.lstrip('chr')]
    mask1 = (tmp_loops['start1'] <= end) & (tmp_loops['end1'] >= start)
    mask2 = (tmp_loops['start2'] <= end) & (tmp_loops['end2'] >= start)
    mask3 = (tmp_snps['CHR_POS'] >= start - 100) & (tmp_snps['CHR_POS'] <= end + 100)
    overlap1 = tmp_loops[mask1]
    overlap2 = tmp_loops[mask2]
    overlap3 = tmp_snps[mask3]
    # print (overlap3)
    if (len(overlap1) != 0):
        overlap = pd.concat([overlap , overlap1])
    if (len(overlap2) != 0):
        overlap = pd.concat([overlap , overlap2])    
    overlap = overlap.drop_duplicates()
    # overlap.to_csv('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\HiRPC\\SNP_peaks_loops\\' + peak_name + '.bedpe' , header=None , index = None , sep = '\t')
    
    
    for i in overlap.index:
        g = overlap.loc[i]['chr1']
        start1 = overlap.loc[i]['start1']
        end1 = overlap.loc[i]['end1']
        start2 = overlap.loc[i]['start2']
        end2 = overlap.loc[i]['end2']
        tmp_promoter = promoter[promoter['chr'] == g]
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
        
        
        # print (selected_genes)
        
        
        for i in list(selected_genes['gene_name']):
            if i in CVD_genes:
                print (peak_name , i)



















    
    
peaks_name = ['HCT116_R2_q0.05_peak_30906' , 'HCT116_R2_q0.05_peak_33685' , 'HCT116_R2_q0.05_peak_33688' , 'HCT116_R2_q0.05_peak_45512' ,
              'HCT116_R2_q0.05_peak_50480' , 'HCT116_R2_q0.05_peak_51595' , 'HCT116_R2_q0.05_peak_8151' , 'HCT116_R2_q0.05_peak_13234' , 
              'HCT116_R2_q0.05_peak_21636' , 'HCT116_R2_q0.05_peak_26807' , 'HCT116_R2_q0.05_peak_27136']


# peaks_name = ['HCT116_R2_q0.05_peak_526']



out = open('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks.fasta' , 'w')
# out = open('H:\\work\\Postdoctoral\\GWAS疾病位点检测\\results\\CAD\\first_6000\\Confirmation_Experiment\\NPPA_NPPB_peak2.fasta' , 'w')
for peak_name in peaks_name:
    peak = peaks[peaks['peak_name'] == peak_name]
    g = peak.iloc[0]['chr']
    start = peak.iloc[0]['start'] 
    end = peak.iloc[0]['end'] 
    peaks_seq =  genome_h38[g][start:end]
    out.writelines('>' + g + ': ' + str(start) + '_' + str(end) + '_' + peak_name + '\n')
    for i in range(len(peaks_seq) // 70 + 1):
        n1 = i * 70
        n2 = (i + 1) * 70
        out.writelines(''.join(peaks_seq[n1 : n2]).upper() + '\n')
    out.writelines('\n')

out.close()



selected_peak_genes = {}

for peak_name in peaks_name:
    overlap = pd.DataFrame([])
    selected_genes = pd.DataFrame([])
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
    # overlap.to_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\' + peak_name + '.bedpe' , header=None , index = None , sep = '\t')
    
 
    for i in overlap.index:
        g = overlap.loc[i]['chr1']
        start1 = overlap.loc[i]['start1']
        end1 = overlap.loc[i]['end1']
        start2 = overlap.loc[i]['start2']
        end2 = overlap.loc[i]['end2']
        tmp_promoter = promoter[promoter['chr'] == g]
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

    selected_genes = selected_genes.drop_duplicates()        
    selected_genes = selected_genes.sort_values(by = 'FPKM' , ascending=False)        
    selected_genes = selected_genes[selected_genes['FPKM'] >= 10]
    gene_names = list(selected_genes['gene_name'])
    selected_peak_genes[peak_name] = gene_names
    
    
    # print (selected_genes)
    CDS = {}
    for name in gene_names:
        CDS[name] = []
        if '_' in name:
            n = name.split('_')[1]
        else:
            n = name
        gene = gtf[gtf['gene_name'] == n]
        gene = pd.DataFrame(gene)
        gene = gene.sort_values(by = ['start' , 'end'])
        gene = gene.drop_duplicates(keep = 'first')
        if len(gene) == 0:
            continue
        g = gene.iloc[0]['chr']
        intervals = list(zip(gene['start'], gene['end']))
        merged_intervals = merge_intervals(intervals)
        merged_df = pd.DataFrame(merged_intervals, columns=['start', 'end'])
        for i in range(len(merged_df) - 1):
            start = merged_df.iloc[i]['start']
            end = merged_df.iloc[i]['end']
            if merged_df.iloc[i]['end'] <= merged_df.iloc[i + 1]['start']:
                CDS[name] += genome_h38[g][start - 1:end]
            else:
                print(i , name)
                break
            
        start = merged_df.iloc[-1]['start'] 
        end = merged_df.iloc[-1]['end']
        
        CDS[name] += genome_h38[g][start - 1:end]
                

    write_pos_to_seq(CDS , 'H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\' + peak_name  + '_related_genes_CDS.fa')














