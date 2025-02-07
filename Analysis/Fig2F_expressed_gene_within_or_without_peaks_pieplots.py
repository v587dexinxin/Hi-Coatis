# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 21:22:08 2024

@author: lenovo
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    


#####################################K562


RNA = pd.read_csv('H:\\work\\literature_data\\K562\\total_RNA_seq\\FPKM\\union_all_FPKM.csv' , header=0)


strand1 = RNA[RNA['Strand'] == '+']
strand2 = RNA[RNA['Strand'] == '-']

strand1['pstart'] = strand1['Start'] - 2000
strand1['pend'] = strand1['End'] + 3000
strand2['pstart'] = strand2['Start'] - 3000
strand2['pend'] = strand2['End'] + 2000


RNA_new = pd.concat([strand1 , strand2])

peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak' , header = None)
peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signalValue' , 'pvalue' , 'qvalue' , 'peak']
# peaks = peaks[peaks['score'] >= np.percentile(peaks['score'] , 70)]

# peaks = pd.read_table('H:\\work\\literature_data\\K562\\ChIP-seq\\hg38\\peaks\\K562_H3K27ac_hg38_ENCFF153YGA.bed' , header = None)
# peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signalValue' , 'pvalue' , 'qvalue' , 'peak']


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


peaks_gene = 0
peaks_exp_gene = 0
for g in chrom:
    print (g)
    tmp_gene = RNA_new[RNA_new['Chr'] == g]
    tmp_peaks = peaks[peaks['chr'] == g]
    for i in tmp_gene.index:
        start = tmp_gene.loc[i]['pstart']
        end = tmp_gene.loc[i]['pend']
        fpkm = tmp_gene.loc[i]['FPKM']
        mask = (tmp_peaks['start'] <= end) &  (tmp_peaks['end'] >= start)
        overlap = tmp_peaks[mask]
        if len(overlap) > 0:
            peaks_gene += 1
            if fpkm > 2:
                peaks_exp_gene += 1
    

ratio_562 = peaks_exp_gene / len(RNA_new[RNA_new['FPKM'] > 2])




#####################################HCT116


RNA = pd.read_csv('H:\\work\\literature_data\\HCT116\\hg38\\RNA\\total_RNA_seq\\FPKM\\union_all_FPKM.csv' , header=0)


strand1 = RNA[RNA['Strand'] == '+']
strand2 = RNA[RNA['Strand'] == '-']

strand1['pstart'] = strand1['Start'] - 2000
strand1['pend'] = strand1['End'] + 3000
strand2['pstart'] = strand2['Start'] - 3000
strand2['pend'] = strand2['End'] + 2000


RNA_new = pd.concat([strand1 , strand2])

peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\HCT116\\HCT116_HiRPC_0.1FA\\one-dimensional\\peaks\\HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak' , header = None)
peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score' , 'strand' , 'signalValue' , 'pvalue' , 'qvalue' , 'peak']
# peaks = peaks[peaks['score'] >= np.percentile(peaks['score'] , 70)]

chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


peaks_gene = 0
peaks_exp_gene = 0
for g in chrom:
    print (g)
    tmp_gene = RNA_new[RNA_new['Chr'] == g]
    tmp_peaks = peaks[peaks['chr'] == g]
    for i in tmp_gene.index:
        start = tmp_gene.loc[i]['pstart']
        end = tmp_gene.loc[i]['pend']
        fpkm = tmp_gene.loc[i]['FPKM']
        mask = (tmp_peaks['start'] <= end) &  (tmp_peaks['end'] >= start)
        overlap = tmp_peaks[mask]
        if len(overlap) > 0:
            peaks_gene += 1
            if fpkm > 2:
                peaks_exp_gene += 1
    

ratio_116 = peaks_exp_gene / len(RNA_new[RNA_new['FPKM'] > 2])




# left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
# size_axes = [left, bottom, width, height]
# fig = plt.figure(figsize = (12, 12))
# ax = fig.add_axes(size_axes)
# ax.bar([1 , 2] , [1 , 1] ,  color = 'darkblue')
# ax.bar([1 , 2] , [ratio_562 , ratio_116] , color = 'gold' , label = 'expressed gene with peaks')
# ax.legend()

# ax.set_xticks([1 , 2]) 
# ax.set_xticklabels(['K562' , 'HCT116'])
# ax.set_ylabel('percentage')




left, bottom, width, height = 0.2, 0.2, 0.4, 0.4
size_axes1 = [left, bottom, width, height]
size_axes2 = [left + 0.25, bottom, width, height]

fig = plt.figure(figsize = (12, 8))
ax1 = fig.add_axes(size_axes1)
ax1.pie([ratio_562 , 1-ratio_562] ,  
       colors = ['salmon' , 'lightblue'] , autopct='%1.1f%%', startangle=90 , 
       textprops={'fontsize': 20, 'fontweight': 'bold'})

ax1.set_title('K562', fontsize=30, fontweight='bold')

ax2 = fig.add_axes(size_axes2)
wedges, texts, autotexts = ax2.pie([ratio_116 , 1-ratio_116] , labels=None,  
       colors = ['salmon' , 'lightblue'] , autopct='%1.1f%%', startangle=90 , 
       textprops={'fontsize': 20, 'fontweight': 'bold'})

labels = ['expressed gene with peaks' , 'without_peaks'] 

ax2.legend(wedges, labels, loc="lower center", bbox_to_anchor=(0.5, -0.3))

ax2.set_title('HCT116', fontsize=30, fontweight='bold')

run_Plot(fig , 'H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\Fig2F_K562_HCT116_expressed_gene_with_peaks\\K562_HCT116_expressed_gene_with_VS_without_peaks_pie_plots.pdf')

