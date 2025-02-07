# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 15:25:35 2024

@author: lenovo
"""


import pandas as pd
import cooler
import numpy as np
import pyBigWig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages




def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in chrom:
        matrix = HiC_Lib.matrix(balance=True).fetch(g)
        matrix[np.isnan(matrix)] = 0
        Lib_new[g] = matrix
    return Lib_new


def Get_bw_sig(bw_file , g , start , end):
    '''
    '''
    bw = pyBigWig.open(bw_file, 'r')
    values = bw.values(g, start, end)
    return values




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
    ax.boxplot(data[2] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'darkred','linewidth':1},
            medianprops={'color':'darkred','linewidth':1},
            capprops={'color':'darkred','linewidth':1},
            whiskerprops={'color':'darkred','linewidth':1})
    ax.boxplot(data[3] , positions=[5] , showfliers=False, widths = 0.7 ,
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


    ax.set_xticks([1 , 2 , 3 , 4 , 5])
    ax.set_xticklabels(['left_+' , 'left_-' , '' , 'right_+' , 'right_-'] , fontsize = 10)
    ax.set_ylabel('CTCF Signal intensity' , fontsize = 20)
    ax.set_xlabel(title)
    ax.set_xlim((0.5 , 5.5))
    # ax.set_title(cl + ',TAD_numbers:' + str(len(tads[cl])))
    ax.set_ylim((vmin , vmax))

    return fig



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    




chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
R = 25000



##Load CTCF motif data

ctcf = pd.read_table('/scratch/2024-05-06/bio-shenw/ref/memes/fimo/hg38/v11/one_by_one/CTCF_fimo/fimo.tsv' , header = 0)
ctcf = ctcf.drop(ctcf.tail(3).index)


##Load CTCF peaks

ctcf_peaks = pd.read_table('/scratch/2024-05-06/bio-shenw/literature/K562/ChIP_seq/hg38/peaks/K562_CTCF_hg38_ENCFF901CBP.bed' , header=None , usecols=(0 , 1 , 2 , 3 , 4))
ctcf_peaks.columns = ['chr' , 'start' , 'end' , 'name' , 'score']


##Load CTCF data

ctcf_bw = pyBigWig.open('/scratch/2024-05-06/bio-shenw/literature/K562/ChIP_seq/hg38/K562_CTCF_hg38_ENCFF682MFJ.bigWig', 'r')


##Load matrix data
RPC = cooler.Cooler("/scratch/2024-05-06/bio-shenw/Ljniu/K562/matrix/cool/K562_RPC_0.1FA/Combined_562_0.1FA_merged5_.hg38.mapq_30.1000_balanced.mcool::/resolutions/25000")
RPC_m = Get_nan_zero_Matrix(RPC)


 
##Load Strip_TADs_data
left_strip = pd.read_table('/scratch/2024-05-06/bio-shenw/Ljniu/K562/K562_0.1_FA/TADs/TADs_25K_balanced/Strip_TADs/Strip_10K/K562_0.1FA_left_strip_10K.bedpe' , header = None)
right_strip = pd.read_table('/scratch/2024-05-06/bio-shenw/Ljniu/K562/K562_0.1_FA/TADs/TADs_25K_balanced/Strip_TADs/Strip_10K/K562_0.1FA_right_strip_10K.bedpe' , header = None)

left_strip.columns = ['chr' , 'pos1' , 'pos2' , 'chr2' , 'pos3' , 'pos4']
right_strip.columns = ['chr' , 'pos1' , 'pos2' , 'chr2' , 'pos3' , 'pos4']

strips = {'left' : left_strip,
                   'right' : right_strip}




###direction peaks


n = 0 ; direction_peaks = {}
for s in strips:
    strip_1 = strips[s]
    direction_peaks[s] = []
    for g in chrom:
        tmp_ctcf = ctcf[ctcf['sequence_name'] == g]
        tmp_peaks = ctcf_peaks[ctcf_peaks['chr'] == g]
        tmp_strip_1 = strip_1[strip_1['chr'] == g]
        for i in tmp_strip_1.index:
            start = tmp_strip_1.loc[i]['pos1'] 
            end = tmp_strip_1.loc[i]['pos1'] + 50000
            mask = (tmp_peaks['start'] <= end) & (tmp_peaks['end'] >= start)
            overlap = tmp_peaks[mask]
            if len(overlap) != 0:
                n += 1
                for j in overlap.index:
                    direction_peaks[s].append((g , overlap.loc[j]['start'] , overlap.loc[j]['end'] , overlap.loc[j]['score']))
                    

                
for k in direction_peaks:
    v = pd.DataFrame(direction_peaks[k])
    v.columns = ['chr' , 'start' , 'end' , 'score']
    direction_peaks[k] = v
     
     

    
###direction motifs

motifs = {}       
for s in direction_peaks:
    d_peaks = direction_peaks[s]
    motifs[s] = []
    for g in chrom:
        tmp_d_peaks = d_peaks[d_peaks['chr'] == g]
        tmp_ctcf = ctcf[ctcf['sequence_name'] == g]
        for i in tmp_d_peaks.index:
            start = tmp_d_peaks.loc[i]['start']
            end = tmp_d_peaks.loc[i]['end']
            mask = ((start <= tmp_ctcf['stop']) & (end >= tmp_ctcf['start']))
            overlap = tmp_ctcf[mask]
            if len(overlap) != 0:
                for j in overlap.index:
                    ctcf_s = int(overlap.loc[j]['start'])
                    ctcf_e = int(overlap.loc[j]['stop'])
                    strand = overlap.loc[j]['strand']
                    motifs[s].append((g , ctcf_s , ctcf_e , strand , 1 , '.' , 1 , 1 , 1 , 1))
    
          
                    
for k in motifs:
    v = pd.DataFrame(motifs[k])
    v.columns = ['chr' , 'start' , 'end' , 'strand' , 'score' , '6' , '7' , '8' , '9' , '10']
    motifs[k] = v
    



x =  ['left_strip' , 'right_strip']   
y1 = [1 , 1]
y2 = [len(motifs['left'][motifs['left']['strand'] == '+']) / len(motifs['left']) , 
      len(motifs['right'][motifs['right']['strand'] == '+']) / len(motifs['right'])]
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.bar(x , y1 , color = 'gold')
ax.bar(x , y2 , color = 'midnightblue')

run_Plot(fig , '/scratch/2024-05-06/bio-shenw/Ljniu/K562/plots/Fig3D_strip_TADs_VS_CTCF_direction/Strip_TADs_enrichment_CTCF_motif_Direction_percentage_barplot.pdf')




                    
                    
                  
            
###direction signals boxplot


signals = {}


for s in motifs:
    signals[s] = {'+':[] , '-':[]}
    motif = motifs[s]
    for g in chrom:
        tmp_peaks = ctcf_peaks[ctcf_peaks['chr'] == g]
        tmp_motif = motif[motif['chr'] == g]
        for i in tmp_motif.index:
            start = tmp_motif.loc[i]['start']
            end = tmp_motif.loc[i]['end']
            strand = tmp_motif.loc[i]['strand']
            mask = ((start <= tmp_peaks['end']) & (end >= tmp_peaks['start']))
            overlap = tmp_peaks[mask].iloc[0]
            bw = ctcf_bw.values(g , overlap['start'] , overlap['end'])
            signals[s][strand].append(np.mean(bw))

fig = Box_plot_4cellline([signals['left']['+'] , signals['left']['-'] , signals['right']['+'] , signals['right']['-']] , 0 , 120 , 'Strip TADs VS CTCF motif Direction')
run_Plot(fig , '/scratch/2024-05-06/bio-shenw/Ljniu/K562/plots/Fig3D_strip_TADs_VS_CTCF_direction/Strip_TADs_VS_CTCF_motif_Direction_boxplot.pdf')



