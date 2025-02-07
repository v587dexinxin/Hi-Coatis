# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 01:28:07 2024

@author: lenovo
"""

from __future__ import division
import numpy as np 
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
import scipy
from scipy.stats import ranksums
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap


chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


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


def peak_intensity(bedfile , peaks):
    data = pd.read_table(bedfile , header = None)
    data.columns = ['chr' , 'start' , 'end' , 'score']
    intensity = []
    for g in chrom:
        print (g)
        tmp_peaks = peaks[peaks['chr'] == g]
        tmp_sig = data[data['chr'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            mask = (tmp_sig['start'] <= end) & (tmp_sig['end'] >= start)
            overlap = tmp_sig[mask]
            if len(overlap) != 0:
                intensity.append(np.array(overlap['score']).mean())
    return intensity







def Box_plot_4cellline(data , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})


    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['50K' , '25K' , '5K' , '1K' ] , fontsize = 20)
    ax.set_ylabel('Peak intensity' , fontsize = 30)
    ax.set_xlim((0.5 , 4.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig



def Box_plot_5cellline(data , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 , 
            boxprops={'color': 'darkred','linewidth':2},
            medianprops={'color':'darkred','linewidth':2},
            capprops={'color':'darkred','linewidth':2},
            whiskerprops={'color':'darkred','linewidth':2})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'dodgerblue','linewidth':2},
            medianprops={'color':'dodgerblue','linewidth':2},
            capprops={'color':'dodgerblue','linewidth':2},
            whiskerprops={'color':'dodgerblue','linewidth':2})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkgreen','linewidth':2},
            medianprops={'color':'darkgreen','linewidth':2},
            capprops={'color':'darkgreen','linewidth':2},
            whiskerprops={'color':'darkgreen','linewidth':2})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkorange','linewidth':2},
            medianprops={'color':'darkorange','linewidth':2},
            capprops={'color':'darkorange','linewidth':2},
            whiskerprops={'color':'darkorange','linewidth':2})

    ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'deeppink','linewidth':2},
            medianprops={'color':'deeppink','linewidth':2},
            capprops={'color':'deeppink','linewidth':2},
            whiskerprops={'color':'deeppink','linewidth':2})
    # d1 = np.round(wilcoxon(data[0] , data[1])[1] , 5)
    # d2 = np.round(wilcoxon(data[0] , data[2])[1] , 5)
    # d3 = np.round(wilcoxon(data[1] , data[2])[1] , 5)
    
    
    # d1 = np.round(scipy.stats.ranksums(data[0] , data[1])[1] , 5)
    # d2 = np.round(scipy.stats.ranksums(data[0] , data[2])[1] , 5)
    # d3 = np.round(scipy.stats.ranksums(data[1] , data[2])[1] , 5)

    
    ax.set_xticks([1 , 2 , 3 , 4 , 5])
    ax.set_xticklabels(['1M' , '50K' , '25K' , '5K' , '1K' ] , fontsize = 20)
    ax.set_ylabel('Peak intensity' , fontsize = 30)
    ax.set_xlim((0.5 , 5.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))
    
    return fig






##-------------K562_1M_peaks-------

from sklearn.preprocessing import MinMaxScaler
import numpy as np

def norm(data , val):
    average = np.array(data).mean()
    data_new = np.array(data) * val /average
    return data_new

def norm_min_max(data , val):
    
    data = np.array(data)
    min_val = np.min(data)
    max_val = np.max(data)
    normalized_data = (data - min_val) / (max_val - min_val)
    return(normalized_data)


def norm_zscore(data , val):
    
    data = np.array(data)
    
    mean = np.mean(data)
    std = np.std(data)
    normalized_data = (data - mean) / std
    return(normalized_data)


    

peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak', header = None , usecols=(0 , 1 , 2))

peaks.columns = ['chr' , 'start' , 'end']



K562_01_FA = peak_intensity('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\signals\\all_reps\\K562_0.1FA_allreps_RPKM_10bp.bedgraph' , peaks)
intensity_50k = peak_intensity('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_lowcells\\one-dimensional\\signals\\all_reps\\bedgraph\\K562_lowcells_50k_allreps_RPKM_10bp.bedgraph' , peaks)
intensity_25k = peak_intensity('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_lowcells\\one-dimensional\\signals\\all_reps\\bedgraph\\K562_lowcells_25k_allreps_RPKM_10bp.bedgraph' , peaks)
intensity_5k = peak_intensity('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_lowcells\\one-dimensional\\signals\\all_reps\\bedgraph\\K562_lowcells_5k_allreps_RPKM_10bp.bedgraph' , peaks)
intensity_1k = peak_intensity('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_lowcells\\one-dimensional\\signals\\all_reps\\bedgraph\\K562_lowcells_1k_allreps_RPKM_10bp.bedgraph' , peaks)


K562_01_FA_new = norm(K562_01_FA, 100)
intensity_50k_new = norm(intensity_50k, 100)
intensity_25k_new = norm(intensity_25k, 100)
intensity_5k_new = norm(intensity_5k, 100)
intensity_1k_new = norm(intensity_1k, 100)


K562_01_FA_new = norm_min_max(K562_01_FA, 100)
intensity_50k_new = norm_min_max(intensity_50k, 100)
intensity_25k_new = norm_min_max(intensity_25k, 100)
intensity_5k_new = norm_min_max(intensity_5k, 100)
intensity_1k_new = norm_min_max(intensity_1k, 100)


K562_01_FA_new = norm_zscore(K562_01_FA, 100)
intensity_50k_new = norm_zscore(intensity_50k, 100)
intensity_25k_new = norm_zscore(intensity_25k, 100)
intensity_5k_new = norm_zscore(intensity_5k, 100)
intensity_1k_new = norm_zscore(intensity_1k, 100)


K562_01_FA_new = np.log2(K562_01_FA)
intensity_50k_new = np.log2(intensity_50k)
intensity_25k_new = np.log2(intensity_25k)
intensity_5k_new = np.log2(intensity_5k)
intensity_1k_new = np.log2(intensity_1k)





data = [K562_01_FA_new , intensity_50k_new , intensity_25k_new , intensity_5k_new , intensity_1k_new]
 



fig = Box_plot_5cellline(data , 2 , 11)


run_Plot(fig , 'H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\Fig7A_lowcells_peak_intensity\\K562_lowcells_peak_intensity_boxplot_2_log2.pdf')







