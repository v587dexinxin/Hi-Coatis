# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 15:11:14 2021

@author: xxli
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




def Volcano_Gene_Plot(fil,title, s1 , s2 , fc = 0.5 , q = 0.01):
    """
    """
    f = open(fil,'r')
    Gene = []
    for line in islice(f,1,None):
        line= line.strip().split(',')
        if line[6] == 'NA' or line[10] == 'NA' or (line[1] not in chrom):
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
    ax.plot([-15,15],[-math.log10(q),-math.log10(q)], ls = '--', c = 'black', lw = 1.0)
    ax.plot([0,0],[-2,200], ls = '--', c = 'black', lw = 1.0)
    ax.plot([NC_bound,NC_bound],[-2,200], ls = '--', c = 'red', lw = 1.0)
    ax.plot([si_bound,si_bound],[-2,200],ls = '--',c = 'blue', lw = 1.0)
    
    ax.set_xticks([-8,-4,NC_bound,0,si_bound,4,8])
    ax.set_xticklabels(['-8','-4',str(NC_bound),'0',str(si_bound),'4','8'])
    ax.set_xlabel('log2FoldChange', size = 15)
    ax.set_ylabel('-log10(q-value)',size = 15)
    ax.set_ylim(-2,200)
    ax.set_xlim(-15,15)
    ax.text(3,150,s1 + '_up_genes : %d' % len(NC_Genes))
    ax.text(-14,150,s2 + '_down_genes : %d' % len(si_Genes))
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


def Get_diff_genes(genes , fc , q):
    '''
    '''
    
    wt_RNA = RNA[(RNA['Chr'].isin(chrom)) & (RNA['log2FoldChange'] > fc) & (RNA['padj'] < q)]
    hemin_RNA = RNA[(RNA['Chr'].isin(chrom)) & (RNA['log2FoldChange'] < -fc) & (RNA['padj'] < q)]
    stable_RNA = RNA[(RNA['Chr'].isin(chrom)) & ~(((RNA['log2FoldChange'] > fc) & (RNA['padj'] < q)) | ((RNA['log2FoldChange'] < -fc) & (RNA['padj'] < q)))]
    
    return wt_RNA , hemin_RNA , stable_RNA




chrom = ['chr' + str(i) for i in range(1 , 23)] + ['chrX']


###########K562_HCT116#################

pp1 = PdfPages('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\DEGs\\K562_HCT116_RNA_DEGs_scatter_q0.001_fc0.5.pdf')


fig = Volcano_Gene_Plot('H:\\work\\niulongjian\\HiRPC_processed_data\\K562_HCT116_RNA-seq\\DEGs\\DEGs_K562_VS_HCT116.csv' , 'K562_VS_HCT116_DEGs' , 'K562' , 'HCT116' , fc = 0.5, q = 0.001)

pp1.savefig(fig)

pp1.close() 


###########K562_Hemin#################
pp1 = PdfPages('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_RNA-seq\\DEGs\\K562_Hemin_RNA_DEGs_scatter_q0.001_fc0.5.pdf')


fig = Volcano_Gene_Plot('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_RNA-seq\\DEGs\\K562_DEGs_WT_VS_Hemin.csv','K562_VS_Hemin_DEGs', 'K562' , 'Hemin' , fc = 0.5, q = 0.001)

pp1.savefig(fig)

pp1.close() 
