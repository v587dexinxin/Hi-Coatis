# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 19:50:39 2024

@author: lenovo
"""


import pandas as pd


TSS = pd.read_table('H:\\work\\literature_data\\genome\\hg38\\UCSC\\TSS_+-1bp_gene_name.bed' , header = None)
TSS.columns = ['chr' , 'start' , 'end' , 'gene_name' , 'strand']
loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\loops\\HCT116_merged6_0.1FA.hg38_loops_one_anchor_binding_union_peaks.bedpe' , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'Pvalue']
genes = pd.read_table('H:\\work\\literature_data\\genome\\hg38\\UCSC\\hg38.ncbiRefSeq.bed' , header = None)
genes.columns = ['chr' , 'start' , 'end' , 'gene_name' , 'name' , 'strand']

peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak')

RNA = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\RNA-seq\\FPKM\\HCT116-WT-R1_FPKM.txt' , header = 0)


gene_name = ['KIF3A' , 'FALEC' , 'TLCD3B' , 'SWT1' , 'PARD6B' , 'RP4' , 'WDSUB1' , 'PTK6' , 'SPATA1' , 'DMXL2']
gene_name = ['BRCA1']


interact_loops = []
for name in gene_name:
    if name == 'RP4':
        print (name)
    else:
        gene = TSS[TSS['gene_name'] == name]
        tss = gene.iloc[0]['start'] + 1
        strand = gene.iloc[0]['strand']
        g = gene.iloc[0]['chr']
        tmp_loop = loops[loops['chr1'] == g]
        if strand == '+':
            start = tss - 2000
            end = tss + 1000
        else:
            start = tss - 1000
            end = tss + 2000
        
        mask1 = (tmp_loop['start1'] <= end) & (tmp_loop['end1'] >= start)
        mask2 = (tmp_loop['start2'] <= end) & (tmp_loop['end2'] >= start)
        overlap1 = tmp_loop[mask1]
        overlap2 = tmp_loop[mask2]
        if len(overlap1) > 0:
            for i in overlap1.index:
                interact_loops.append(tuple(list(overlap1.loc[i]) + [name]))
        if len(overlap2) > 0:
            for i in overlap2.index:
                interact_loops.append(tuple(list(overlap2.loc[i]) + [name]))
                
                
            
                
                
interact_loops = pd.DataFrame(interact_loops)
interact_loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'Pvalue' , 'gene_name']

interact_loops = interact_loops.drop_duplicates()                


genes_1 = genes[genes['strand'] == '+']
genes_2 = genes[genes['strand'] == '-']

genes_1['start'] = genes_1['start'] - 2000  
genes_2['end'] = genes_2['end'] + 2000            
            
        
genes_new = pd.concat((genes_1 , genes_2))



interact_loops_new = []

for i in interact_loops.index:
    g = interact_loops.loc[i]['chr1']
    start1 = interact_loops.loc[i]['start1']
    end1 = interact_loops.loc[i]['end1']
    start2 = interact_loops.loc[i]['start2']
    end2 = interact_loops.loc[i]['end2']
    IF = interact_loops.loc[i]['IF']
    pvalue = interact_loops.loc[i]['Pvalue']
    tmp_genes = genes_new[genes_new['chr'] == g]
    mask1 = (tmp_genes['start'] <= end1) & (tmp_genes['end'] >= start1)
    mask2 = (tmp_genes['start'] <= end2) & (tmp_genes['end'] >= start2)
    overlap1 = tmp_genes[mask1]
    overlap2 = tmp_genes[mask2]
    if (len(overlap1) > 0) & (len(overlap2) > 0):
        g1 = overlap1.iloc[0]['gene_name']
        g2 = overlap2.iloc[0]['gene_name']
        # if (g1 not in gene_name) and (g2 not in gene_name):
        #     continue
        interact_loops_new.append((g , start1 , end1 , g , start2 , end2 , IF , pvalue , g1 , g2))
    else:
        print (i)
        
    
interact_loops_new = pd.DataFrame(interact_loops_new)
interact_loops_new.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'IF' , 'pvalue' , 'gene_name1' , 'gene_name2']


related_genes = {}

for name in gene_name:
    if name == 'RP4':
        related_genes[name] = []
        continue
    a = interact_loops_new[(interact_loops_new['gene_name1'] == name )| (interact_loops_new['gene_name2'] == name)]
    tmp_a = []
    for n in list(set(a['gene_name1'])) + list(set(a['gene_name2'])):
        rna = RNA[RNA['Gene Name'] == n]
        tmp_a.append((n , rna.iloc[0]['Gene Name'] , rna.iloc[0]['FPKM']))   
    tmp_a = pd.DataFrame(tmp_a)
    tmp_a = tmp_a.sort_values(by = (2))
    
    related_genes[name] = tmp_a.iloc[-12:,0]
    a.iloc[:,:8].to_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\慢病毒包装\\loops\\IF_bedpe\\' + name + '_IF.bedpe' , header = None , index = None , sep = '\t')









