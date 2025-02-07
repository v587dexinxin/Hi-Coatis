# -*- coding: utf-8 -*-
"""
Created on Thu May  9 19:14:01 2024

@author: lenovo
"""

from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import sys
import os
import pandas as pd
import cooler
from skimage.transform import resize
 
 
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                           ['#FFFFFF' , '#CD0000'])
my_cmap.set_bad('#2672a1')


def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = True, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , length = 5 ,labelsize = 30  )

def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)







def acquireSingleIns(matrix_data_chr,bound,left_right,category):
    ins=0
    start_site=0;end_site=matrix_data_chr.shape[0]
    if ((bound-left_right<start_site)|(end_site<bound+left_right)):        
        return ins    
    aa=matrix_data_chr[bound-left_right:bound-0,bound+0:bound+left_right]
    b1=[[matrix_data_chr[i,j] for i in range(bound-left_right,bound-0) if j>i] 
            for j in range(bound-left_right,bound-0)]
    b2=[[matrix_data_chr[i,j] for i in range(bound+0,bound+left_right) if j>i] 
            for j in range(bound+0,bound+left_right)]
    
    aa_zero=sum([sum(np.array(item)==0) for item in aa])
    b1_zero=sum([sum(np.array(item)==0) for item in b1])
    b2_zero=sum([sum(np.array(item)==0) for item in b2])
    if aa_zero+b1_zero+b2_zero>=left_right:
        return ins    
    aa_sum=sum([sum(item) for item in aa])
    b1_sum=sum([sum(item) for item in b1])
    b2_sum=sum([sum(item) for item in b2])
    if aa_sum>0:
        if(category=='divide'):
            ins=np.log2((aa_sum+b1_sum+b2_sum)/float(aa_sum))
        elif(category=='average'):
            ins=aa_sum/float(left_right)/left_right
        else:
            print('the calc type went wrong')
    return ins




def Get_nan_zero_Matrix(HiC_Lib):
    '''
    '''
    Lib_new = {}
    for g in chrom:
        matrix = HiC_Lib.matrix(balance=True).fetch(g)
        matrix[np.isnan(matrix)] = 0
        Lib_new[g] = matrix
    return Lib_new




def normalize_TAD_interaction(Lib , Tad , length = 40):
    # IS = np.zeros(length + 20)
    resize_matrix = np.zeros((length , length))
    n = 0
    for i in Tad.index:
        chro = Tad.loc[i]['chr']
        start = Tad.loc[i]['start'] // R
        end = Tad.loc[i]['end'] // R
        d = end - start
        if end - start < 20:
            continue
        startHiC = start - d
        endHiC = end + d
        if startHiC < 0:
            continue

        matrix = Lib[chro][startHiC:endHiC , startHiC:endHiC]
        resize_m = resize(matrix , (length + 20 , length + 20), order = 3)
        resize_matrix += resize_m[10:50,10:50] 
        n += 1
        # for bound in range(length + 20):
        #     Is = acquireSingleIns(resize_m , bound  , 10 , 'average')
        #     IS[bound] += Is

    resize_matrix = resize_matrix / n
    # IS = IS / n
    # IS = acquireSingleIns(resize_matrix , 10  , 10 , 'divide')
    print (n)
    return resize_matrix
    
    
    
    

def normalize_TAD_interaction_list(Lib , Tad , length = 40):
    IS = []
    resize_matrix = np.zeros((length + 20 , length + 20))
    n = 0
    for i in Tad.index:
        chro = Tad.loc[i]['chr']
        start = Tad.loc[i]['start'] // R 
        end = Tad.loc[i]['end'] // R
        d = end - start 
        if end - start < 20:
            continue
        startHiC = start - d
        endHiC = end + d
        if startHiC < 0:
            continue
        
        matrix = Lib[chro][startHiC:endHiC , startHiC:endHiC]
        resize_m = resize(matrix , (length + 20 , length + 20), order = 3)
        resize_matrix += resize_m
        n += 1
        

            
    
    resize_matrix = resize_matrix / n
    
    for bound in range(length + 20):
        Is = acquireSingleIns(resize_matrix , bound  , 10 , 'divide')
        IS.append(Is)
        
    resize_matrix = resize_matrix[10:50,10:50] 
    # IS = IS / n
    # IS = acquireSingleIns(resize_matrix , 10  , 10 , 'divide')
    print (n , IS[10:50])
    return resize_matrix , IS[10:50]
    
    
    

    
def Heatmap_Plot(matrix ,  c):
    # matrix = matrix / matrix.sum()
    vmax = np.percentile(matrix , 85)
    vmin = np.percentile(matrix , 3)
    
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                    extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, vmin = vmin , origin = 'lower')
    cxlim = ax1.get_xlim()
    cylim = ax1.get_ylim()
    ## Ticks and Labels
    ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
    labels = ['-1/2TAD' , 'Boundary' , '' , 'Boundary' , '+1/2TAD']
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)
    ax1.set_yticks(ticks)
    ax1.set_yticklabels(labels, rotation = 'horizontal')
    # ax1.set_xlabel(c , fontsize = 30 )
    ax1.set_xlabel(c + '_IS: ' + str(IS), fontsize = 30 )
    
    ax1.set_xlim(cxlim)
    ax1.set_ylim(cylim)                    
    ## Colorbar
    ax2 = fig.add_axes([Left + 0.6 , HB - 0.06 , 0.1 , 0.035])
    fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [vmin ,vmax])
    return fig



def Insulation_Plot(data):
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    ax1.plot(np.arange(len(data[0])) , data[0] , label = 'RPC')
    ax1.plot(np.arange(len(data[1])) , data[1] , label = 'PolII')
    ax1.plot(np.arange(len(data[2])) , data[2] , label = 'HiChIP')
    
    
    ## Ticks and Labels
#    ticks = list(np.linspace(0 , len(matrix) , 7).astype(float))
#    labels = ['-TAD' , '-1/2TAD' , 'Boundary' , '' , 'Boundary' , '+1/2TAD' , '+TAD']
#    ax1.set_xticks(ticks)
#    ax1.set_xticklabels(labels)
#    ax1.set_yticks(ticks)
#    ax1.set_yticklabels(labels, rotation = 'horizontal')
    ax1.set_ylabel('Insulation Score' , fontsize = 20 )
    ax1.legend()
    

    return fig



def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()



chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
R = 10000
 
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6
 

#HiC Data Process
RPC = cooler.Cooler("/scratch/2024-05-06/bio-shenw/Ljniu/K562/matrix/cool/K562_RPC_0.1FA/Combined_562_0.1FA_merged5_.hg38.mapq_30.1000_balanced.mcool::/resolutions/10000")
RPC_m = Get_nan_zero_Matrix(RPC)
PolII = cooler.Cooler("/scratch/2024-05-06/bio-shenw/Ljniu/K562/matrix/cool/K562_RNAPII_chiaPET/K562_RNAPII_ChIA-PET_hg38_ENCFF582UXM_balanced.mcool::/resolutions/10000")
PolII_m = Get_nan_zero_Matrix(PolII)
hichip = cooler.Cooler("/scratch/2024-05-06/bio-shenw/literature/K562/HiChIP/hg38/workspace/results/coolers_library_group/K562_HiChIP_H3K27ac.hg38.mapq_30.1000.mcool::/resolutions/10000")
hichip_m =  Get_nan_zero_Matrix(hichip)
HiC_Data = {'RPC':RPC_m ,
            'RNAPolII':PolII_m , 
            'HiChIP':hichip_m}

 

##Load Strip_TADs_data
left_strip = pd.read_table('/scratch/2024-05-06/bio-shenw/Ljniu/K562/K562_0.1_FA/TADs/TADs_25K_balanced/Strip_TADs/Strip_10K/K562_0.1FA_left_strip_10K.bedpe' , header = None)
right_strip = pd.read_table('/scratch/2024-05-06/bio-shenw/Ljniu/K562/K562_0.1_FA/TADs/TADs_25K_balanced/Strip_TADs/Strip_10K/K562_0.1FA_right_strip_10K.bedpe' , header = None)

left_strip.columns = ['chr' , 'pos1' , 'pos2' , 'chr2' , 'pos3' , 'pos4']
right_strip.columns = ['chr' , 'pos1' , 'pos2' , 'chr2' , 'pos3' , 'pos4']



strip_tads = {'left' : left_strip,
              'right' : right_strip}


for s in strip_tads:
    strip_1 = strip_tads[s]
    s_new = []
    for i in strip_1.index:
        g = strip_1.loc[i]['chr']
        if s == 'left':
            start = strip_1.loc[i]['pos1']
            end = strip_1.loc[i]['pos4']
        else:
            start = strip_1.loc[i]['pos3']
            end = strip_1.loc[i]['pos1']
        s_new.append((g , start , end))
    tmp = pd.DataFrame(s_new)
    tmp.columns = ['chr' , 'start' , 'end']
    strip_tads[s] = tmp
    



strip_side = 'left'
 
cells = ['RPC' , 'RNAPolII' , 'HiChIP']
# IS_data = []
for c in cells:
    HiC_Lib = HiC_Data[c]
    matrix = normalize_TAD_interaction(HiC_Data[c] , strip_tads[strip_side] , 40)
    if strip_side == 'left':
        IS = acquireSingleIns(matrix , 10  , 10 , 'divide')
    elif strip_side == 'right':
        IS = acquireSingleIns(matrix , 30  , 10 , 'divide')
    print (IS)
    fig = Heatmap_Plot(matrix , c)
    run_Plot(fig, '/scratch/2024-05-06/bio-shenw/Ljniu/K562/plots/Fig3B_strip_Loop_heatmapplot/K562_' + c + '_' + strip_side + '_TADs_interaction_10K_divide.pdf')
    






strip_side = 'right'
 
cells = ['RPC' , 'RNAPolII' , 'HiChIP']
IS_data = []
for c in cells:
    HiC_Lib = HiC_Data[c]
    matrix , IS = normalize_TAD_interaction_list(HiC_Data[c] , strip_tads[strip_side] , 40)
    IS_data.append(IS)
fig1 = Insulation_Plot(IS_data)
run_Plot(fig1, '/scratch/2024-05-06/bio-shenw/Ljniu/K562/plots/Fig3B_strip_Loop_heatmapplot/K562_TADs_' + strip_side + '_IS_lineplots.pdf')


   
   
   
   
   
   
   
   
