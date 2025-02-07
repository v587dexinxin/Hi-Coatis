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




def Volcano_Gene_Plot(fil,title, s1 , s2 , gene_names , qmax , fmax , fc = 0.5 , q = 0.01):
    """
    """
    f = open(fil,'r')
    Gene = []
    for line in islice(f,1,None):
        line= line.strip().split(',')
        if line[6] == 'NA' or line[10] == 'NA' or (line[1] not in chrom) or (float(line[-1]) < 2 and float(line[-2]) < 2 and float(line[-3]) < 2 and float(line[-4]) < 2):
            continue
        else:
            Gene.append((line[0],line[6],-math.log10(float(line[10]) + 10**-323)))
    f.close()
    
    Gene_type = np.dtype({'names':['Gene_name' , 'FC' , 'q'],
                          'formats':['U64' , np.float64,np.float64]})
    
    Gene = np.array(Gene,dtype = Gene_type)
    # P = np.log2(Gene['M'].sum() / Gene['P'].sum())    
    
    NC_bound = fc 
    si_bound = -fc
    
    NC_mask = (Gene['FC']> NC_bound) & (Gene['q'] > -math.log10(q))
    NC_Genes = Gene[NC_mask]
    
    si_mask = (Gene['FC'] < si_bound) & (Gene['q'] > -math.log10(q))
    si_Genes = Gene[si_mask]
    
    Non_Genes = Gene[~(NC_mask | si_mask)] 
    
    
            
    fig,ax = plt.subplots(1)
    ax.scatter(NC_Genes['FC'],NC_Genes['q'], s= 10, c = 'red')
    ax.scatter(si_Genes['FC'],si_Genes['q'], s= 10, c= 'blue')
    ax.scatter(Non_Genes['FC'],Non_Genes['q'], s= 10, c = 'gray')
    ax.plot([-fmax,fmax],[-math.log10(q),-math.log10(q)], ls = '--', c = 'black', lw = 1.0)
    ax.plot([0,0],[-2,qmax], ls = '--', c = 'black', lw = 1.0)
    ax.plot([NC_bound,NC_bound],[-2,qmax], ls = '--', c = 'red', lw = 1.0)
    ax.plot([si_bound,si_bound],[-2,qmax],ls = '--',c = 'blue', lw = 1.0)
    
    ax.set_xticks([-8,-4,NC_bound,0,si_bound,4,8])
    ax.set_xticklabels(['-8','-4',str(NC_bound),'0',str(si_bound),'4','8'])
    ax.set_xlabel('log2FoldChange', size = 15)
    ax.set_ylabel('-log10(q-value)',size = 15)
    ax.set_ylim(-2,qmax)
    ax.set_xlim(-fmax,fmax)
    ax.text(fmax - 5,qmax - 10 ,s1 + '_up_genes : %d' % len(NC_Genes))
    ax.text(-fmax + 1 ,qmax - 10,s2 + '_up_genes : %d' % len(si_Genes))
    for n in gene_names:
        if n in Non_Genes['Gene_name']:
            continue
        gene = Gene[Gene['Gene_name'] == n]
        if len(gene) == 0:
            continue
        if gene['Gene_name'] in Non_Genes['Gene_name']:
            continue
        x = gene[0][1]
        y = gene[0][2]
        ax.scatter(x , y, s= 10, edgecolors='black', facecolors='none', marker='o')
        ax.text(x , y , n)
        
    ax.set_title(title)
    
    return fig


def Get_gene_start(genes):
    a = genes[genes['Strand'] == '+']
    b = genes[genes['Strand'] == '-']
    a['gene_start'] = a['Start']
    b['gene_start'] = b['End']
    tmp = pd.concat([a , b] , axis = 0)
    tmp = tmp.sort_values(by = ['Chr' , 'gene_start'])
    tmp['gene_start_end'] = tmp['gene_start'] + 1
    genes_start = tmp[['Chr' , 'gene_start' , 'gene_start_end' , 'Gene_Name']]
    return (genes_start)


def Get_diff_genes(RNA , fc , q):
    '''
    '''
    
    wt_RNA = RNA[(RNA['Chr'].isin(chrom)) & (RNA['log2FoldChange'] > fc) & (RNA['padj'] < q)]
    hemin_RNA = RNA[(RNA['Chr'].isin(chrom)) & (RNA['log2FoldChange'] < -fc) & (RNA['padj'] < q)]
    stable_RNA = RNA[(RNA['Chr'].isin(chrom)) & ~(((RNA['log2FoldChange'] > fc) & (RNA['padj'] < q)) | ((RNA['log2FoldChange'] < -fc) & (RNA['padj'] < q)))]
    
    return wt_RNA , hemin_RNA , stable_RNA




chrom = ['chr' + str(i) for i in range(1 , 23)] + ['chrX']




#####WT_VS_KO1

genes = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks_related_genes_fasta\\116_selected_intergenic_peaks_related_loops\\116_selected_intergenic_peaks_related_genes\\KO1_HCT116_R2_q0.05_peak_30906related_genes.txt' , header = 0)
gene_names = list(genes['gene_name'])


pp1 = PdfPages('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_WT_VS_KO1_DEGs_scatter_q0.05_fc0.5_1.pdf')
fig = Volcano_Gene_Plot('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_DEGs_WT_VS_KO1.csv','WT_VS_KO1_DEGs', 'WT' , 'KO1' , gene_names , 40 , 8 , fc = 0.5, q = 0.05)

pp1.savefig(fig)
pp1.close() 






#####WT_VS_KO4

genes = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks_related_genes_fasta\\116_selected_intergenic_peaks_related_loops\\116_selected_intergenic_peaks_related_genes\\KO4_HCT116_R2_q0.05_peak_45512related_genes.txt' , header = 0)
gene_names = list(genes['gene_name'])


pp1 = PdfPages('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_WT_VS_KO4_DEGs_scatter_q0.05_fc0.5_1.pdf')
fig = Volcano_Gene_Plot('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_DEGs_WT_VS_KO4.csv','WT_VS_KO4_DEGs', 'WT' , 'KO4' , gene_names , 40 , 8 , fc = 0.5, q = 0.05)

pp1.savefig(fig)
pp1.close() 








#####WT_VS_KO5

genes = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks_related_genes_fasta\\116_selected_intergenic_peaks_related_loops\\116_selected_intergenic_peaks_related_genes\\KO5_HCT116_R2_q0.05_peak_50480related_genes.txt' , header = 0)
gene_names = list(genes['gene_name'])


pp1 = PdfPages('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_WT_VS_KO5_DEGs_scatter_q0.05_fc0.5_1.pdf')
fig = Volcano_Gene_Plot('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_DEGs_WT_VS_KO5.csv','WT_VS_KO5_DEGs', 'WT' , 'KO5' , gene_names , 150 , 9 , fc = 0.5, q = 0.05)

pp1.savefig(fig)
pp1.close() 





#####WT_VS_KO7

genes = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\selected_intergenic_peaks_related_genes_fasta\\116_selected_intergenic_peaks_related_loops\\116_selected_intergenic_peaks_related_genes\\KO7_HCT116_R2_q0.05_peak_8151related_genes.txt' , header = 0)
gene_names = list(genes['gene_name'])


pp1 = PdfPages('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_WT_VS_KO7_DEGs_scatter_q0.05_fc0.5_1.pdf')
fig = Volcano_Gene_Plot('H:\\work\\niulongjian\\HiRPC_processed_data\\验证实验\\基因间区peak敲除\\fasta\\RNA-seq\\DEGs\\HCT116_DEGs_WT_VS_KO7.csv','WT_VS_KO7_DEGs', 'WT' , 'KO7' , gene_names , 60 , 6 , fc = 0.5, q = 0.05)

pp1.savefig(fig)
pp1.close() 



