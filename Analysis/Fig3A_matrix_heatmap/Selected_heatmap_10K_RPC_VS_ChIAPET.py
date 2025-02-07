# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 22:23:16 2023

@author: lenovo
"""

from __future__ import division
import numpy as np
import cooler
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF' ,'#CD0000'])
my_cmap.set_bad('#D3D3D3')




def nan_to_zero(matrix):
    '''
    '''
    nanmask = np.isnan(matrix)
    matrix[nanmask] = 0
    return matrix

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'M'])
    else:
        return ''.join([str(i_part), 'M'])
    
    

####matrix#######
c_rpc = cooler.Cooler("/scratch/2024-10-21/bio-shenw/Ljniu/K562/matrix/cool/K562_RPC_0.1FA/Combined_562_0.1FA_merged5_.hg38.mapq_30.1000.mcool::/resolutions/10000")
c_chiapet = cooler.Cooler("/scratch/2024-10-21/bio-shenw/Ljniu/K562/matrix/cool/K562_RNAPII_chiaPET/K562_RNAPII_ChIA-PET_hg38_ENCFF582UXM.mcool::/resolutions/10000")


rpc_matrix = c_rpc.matrix(balance=False).fetch('chr11')
chiapet_matrix = c_chiapet.matrix(balance=False).fetch('chr11')

rpc_matrix = nan_to_zero(rpc_matrix)
chiapet_matrix = nan_to_zero(chiapet_matrix)



upper_triangle_matrix = np.triu(chiapet_matrix , k = 1)
upper_triangle_matrix = upper_triangle_matrix / upper_triangle_matrix.sum()
lower_triangle_matrix = np.tril(rpc_matrix , k = -1)
lower_triangle_matrix = lower_triangle_matrix / lower_triangle_matrix.sum()

matrix_0 = upper_triangle_matrix + lower_triangle_matrix

np.fill_diagonal(matrix_0, 0)




#####Out_Files#####
OutFolder = '/scratch/2024-10-21/bio-shenw/Ljniu/K562/plots/Fig1_K562_0.1FA_hic_heatmap/RPC_VS_ChIAPET'
OutFil = 'Selected_K562_0.1FA_VS_RNAPII_chiapet_heatmap_10K_chr11_55.2_65M.pdf'
pp = PdfPages(os.path.join(OutFolder , OutFil))


####interval#####


interval = [('chr11' , 55200000 , 65000000)]
interval_1 = [('chr11' , 62000000 , 63000000)]


####Plot####
size = (12, 12)   
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
R = 10000

startHiC = interval[0][1] // R 
endHiC = interval[0][2] // R 

matrix = matrix_0[startHiC:endHiC , startHiC:endHiC]

ticks = list(np.linspace(0 , matrix.shape[1] , 2).astype(float))
pos = [((startHiC + t) * R) for t in ticks]
labels = [properU(p) for p in pos]
nonzero = matrix[np.nonzero(matrix)]
vmax = np.percentile(nonzero, 75)
i_start = interval_1[0][1] // R - startHiC
i_end = interval_1[0][2] // R - startHiC

fig = plt.figure(figsize = size)
ax = fig.add_axes([Left  , HB , width , HH])
sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none', vmax = vmax,extent = (0, matrix.shape[1], matrix.shape[0] , 0) , origin = 'upper')
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xticklabels(labels , fontsize=15)
ax.set_yticklabels(labels , fontsize=15)
ax.set_xlabel('RPC_Chr11' , fontsize=20)
ax.set_ylabel('RNAPolII_ChIAPET_Chr11' , fontsize=20)
ax.yaxis.set_label_position("right")
# ax.set_xlim(0, matrix.shape[1])
# ax.set_ylim(0, matrix.shape[0])
ax.set_title('K562_0.1FA_VS_RNAPolII_ChIAPET_10K_chr11_55.2_65M' , fontsize=25)

# ##Selected_interval

ax.plot([i_start , i_start] , [i_end , i_start] , c = 'black')
ax.plot([i_end , i_start] , [i_start , i_start] , c = 'black')
ax.plot([i_end , i_end] , [i_end , i_start] , c = 'black')
ax.plot([i_end , i_start] , [i_end , i_end] , c = 'black')


## Colorbar
ax = fig.add_axes([Left + 0.5 , HB - 0.1 , 0.1 , 0.035])
cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
cbar.set_ticks([0 , vmax])




pp.savefig(fig)
pp.close()
