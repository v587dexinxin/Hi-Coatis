# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 12:59:37 2024

@author: lenovo
"""

from __future__ import division
import numpy as np
import pandas as pd
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
import pyBigWig




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
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

    
    
def Get_bw_sig(bw_file , g , start , end):
    '''
    '''
    bw = pyBigWig.open(bw_file, 'r')
    values = bw.values(g, start, end)
    return values




def peaks_sig(peaks , sig_fil):
    '''
    '''
    peaks_sig = []
    
    bw = pyBigWig.open(sig_fil, 'r')
    
    
    for g in chrom:
        print (g)
        tmp_peaks = peaks[peaks['chr'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            sig = bw.values(g, start, end)
            peaks_sig.append(np.mean(sig))
            
    return peaks_sig
    
            
    
    
    

def Box_plot_4cellline(data , vmin , vmax , title):                
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
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['random' , 'compartmentA' , 'compartmentB' , 'boundary'] , fontsize = 10)
    ax.set_ylabel('Signal intensity' , fontsize = 20)
    ax.set_xlabel(title)
    ax.set_xlim((0.5 , 4.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig



chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


compartmentA = pd.read_table('/scratch/2024-03-25/bio-shenw/literature/K562/in_situHiC/Encode/K562_insitu_HiC_compartmentA_hg38_ENCFF749YOA.bed' , header=None)
compartmentA.columns = ['chr' , 'start' , 'end']
compartmentB = pd.read_table('/scratch/2024-03-25/bio-shenw/literature/K562/in_situHiC/Encode/K562_insitu_HiC_compartmentB_hg38_ENCFF749YOA.bed' , header=None)
compartmentB.columns = ['chr' , 'start' , 'end']
boundary = pd.read_table('/scratch/2024-03-25/bio-shenw/literature/K562/in_situHiC/Encode/K562_insitu_HiC_boundary_hg38_ENCFF271SAF.bed' , header=None)
boundary.columns = ['chr' , 'start' , 'end']

random = []

hg38 = pd.read_table('/scratch/2024-03-25/bio-shenw/ref/Human/hg38/hg38.chrom.size' , header = None)
hg38.columns = ['chr' , 'length']



for g in chrom:
    tmp = hg38[hg38['chr'] == g]
    length = int(tmp['length'])
    max_ = length // 5000 - 1
    for i in range(max_):
        start = i * 5000
        end = (i + 1) * 5000
        random.append((g , start , end))
        
random = pd.DataFrame(random , columns=boundary.columns)

        




RPC_signals = '/scratch/2024-03-25/bio-shenw/Ljniu/K562/K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw'


compartmentA_sig = peaks_sig(compartmentA , RPC_signals)
compartmentB_sig = peaks_sig(compartmentB , RPC_signals)
boundary_sig = peaks_sig(boundary , RPC_signals)
random_sig = peaks_sig(random , RPC_signals)



data = [random_sig , compartmentA_sig , compartmentB_sig , boundary_sig ]


fig = Box_plot_4cellline(data , 1 , 100 , 'K562 peaks signals')


run_Plot(fig , '/scratch/2024-03-25/bio-shenw/Ljniu/K562/plots/Fig2C_RPC_peaks_around_TADs_compartment/K562_0.1FA_compartment_boundary.pdf')



