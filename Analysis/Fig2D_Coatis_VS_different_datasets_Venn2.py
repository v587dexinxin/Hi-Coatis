from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2, venn2_circles
import pandas as pd
import os
import pandas as pd 



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

def Common_peaks(peaks1 , peaks2):
    n = 0 ; common = []
    chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
    
    for g in chrom:
        tmp1 = peaks1[peaks1['chr'] == g]
        tmp2 = peaks2[peaks2['chr'] == g]
        for i in tmp1.index:
            start = tmp1.loc[i]['start']
            end = tmp1.loc[i]['end']
            mask = (tmp2['start'] <= end) & (tmp2['end'] >= start)
            overlap = tmp2[mask]
            if overlap.size != 0:
                n += 1
                c = tuple(pd.concat([tmp1.loc[i] , overlap.iloc[0]] , axis = 0))
                common.append(c)
    print (n)
    common = pd.DataFrame(common)
    
    return common

def plot_venn2(n1 , n2 , n3 , title , out):
    fig = plt.figure(figsize = (10, 10))
    venn2(subsets=(n1 , n2 , n3), set_labels=('K562_0.1FA_Peaks', 'K562_' + i))
    for text_obj in fig.findobj(matplotlib.text.Text):
        text_obj.set_fontsize(20)
        
    plt.title(title , fontsize = 20)
    run_Plot(fig , out)
                


    
###------------------files----------------------

chip_path = 'H:\\work\\literature_data\\K562\\ChIP-seq\\hg38\\peaks'
pro_path = 'H:\\work\\literature_data\\K562\\PRO-seq\\hg38\\peaks'
Rchip_path = 'H:\\work\\literature_data\\K562\\R_chip\\macs2_remove_input'
polII_path = 'H:\\work\\literature_data\\K562\\ChIP-seq\\hg38\\processed_byself\peaks'
peak_path = 'H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks'


chip_files = os.listdir(chip_path)
chip_files = [x for x in chip_files if '.bed' in x]   


##----------------------------------------------------

pro_peaks = Load_peaks(os.path.join(pro_path , 'K562_PRO-seq_SRR8137173_q0.05_peaks.narrowPeak'), 'narrow')
Rchip_peaks = Load_peaks(os.path.join(Rchip_path , 'K562_D210N_V5ChIP_union_q0.05_peaks.narrowPeak'), 'narrow')
polII_peaks = Load_peaks(os.path.join(polII_path , 'K562_POLR2A_ChIP_union2_q0.05_peaks.narrowPeak'), 'narrow')
peaks_01FA = Load_peaks(os.path.join(peak_path , 'K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak') , 'narrow')
RNA = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_RNA-seq\\FPKM\\union_all_FPKM_K562_Hemin.csv' , header = 0)
RNA['K562_WT_FPKM'] = (RNA['K562_WT_R1_FPKM'] + RNA['K562_WT_R2_FPKM']) / 2

RNA = RNA.drop_duplicates(subset = ['Gene_Name'] , keep = 'first')
expressed_rna = RNA[RNA['K562_WT_FPKM'] >= 2]





chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

 
Peaks_all = {}
for f in chip_files:
    sample = f.split('_')[1]
    peaks = Load_peaks(os.path.join(chip_path , f) , 'narrow')
    Peaks_all[sample] = peaks
    
Peaks_all['PRO'] = pro_peaks
Peaks_all['RChIP'] = Rchip_peaks
Peaks_all['PolII'] = polII_peaks


Commons = {}
for i in Peaks_all.keys():
    print (i)
    Commons[i] = Common_peaks(peaks_01FA , Peaks_all[i])
    n1 = len(peaks_01FA) - len(Commons[i])
    n2 = len(Peaks_all[i]) - len(Commons[i])
    n3 = len(Commons[i])
    plot_venn2(n1 , n2 , n3 , '0.1FA_Peaks_VS_' + i + '_peaks_overlap' ,'H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\Fig2C_RPC_PRO_RChIP_POl2A_overlap_venn2\\reps_union\\' + i + '_peaks_overlap_venn2.pdf')
    




##--------------genes----------------

genes = []
for i in expressed_rna.index:
    g = expressed_rna.loc[i]['Chr']
    strand = expressed_rna.loc[i]['Strand']
    start = expressed_rna.loc[i]['Start']
    end = expressed_rna.loc[i]['End']
    fpkm = expressed_rna.loc[i]['K562_WT_FPKM']
    if strand == '+':
        start = start - 2000
    else:
        end = end + 2000
    genes.append((g , start , end , fpkm))
    
    
genes = pd.DataFrame(genes)
genes.columns = ['chr' , 'start' , 'end' , 'fpkm']


n = 0 ; common = []
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']

for g in chrom:
    tmp1 = peaks_01FA[peaks_01FA['chr'] == g]
    tmp2 = genes[genes['chr'] == g]
    for i in tmp2.index:
        start = tmp2.loc[i]['start']
        end = tmp2.loc[i]['end']
        mask = (tmp1['start'] <= end) & (tmp1['end'] >= start)
        overlap = tmp1[mask]
        if overlap.size != 0:
            n += 1
            
        
print (n)

n1 = len(peaks_01FA) - n
n2 = len(genes) - n
n3 = n

plot_venn2(n1 , n2 , n3 , '0.1FA_Peaks_VS_expressed_genes_overlap' ,'H:\\work\\niulongjian\\K562\\hg38\\0.1FA\\plots\\expressed_genes_VS_0.1FA_peaks_overlap_venn2.pdf')



#######################eRNA#####################################
peaks1 = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks1.columns = ['chr' , 'start' , 'end']

peaks_distal = pd.read_table('H:\\work\\literature_data\\K562\\eRNA\\Distal_K562_f1df8e8f-c5a5-4731-9135-2359852b41bd.bed' , usecols = (0 , 1 , 2) , header = None)
peaks_proximal = pd.read_table('H:\\work\\literature_data\\K562\\eRNA\\Proximal_K562_f1df8e8f-c5a5-4731-9135-2359852b41bd.bed' , usecols = (0 , 1 , 2) , header = None)



peaks_erna = pd.concat([peaks_distal , peaks_proximal])
peaks_erna.columns = ['chr' , 'start' , 'end']


peaks_erna = peaks_erna.reset_index(drop = True)



chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']


n = 0
for g in chrom:
    print (g)
    tmp_peaks = peaks1[peaks1['chr'] == g]
    tmp_erna = peaks_erna[peaks_erna['chr'] == g]
    for i in tmp_peaks.index:
        start = tmp_peaks.loc[i]['start']
        end = tmp_peaks.loc[i]['end']
        mask = (start <= tmp_erna['end']) & (end >= tmp_erna['start'])
        overlap = tmp_erna[mask]
        if len(overlap) > 0:
            n += 1
            
            
            
plot_venn2(len(peaks1) - n , len(peaks_erna) - n , n , 'K562 Coatis VS eRNA' , 'H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\K562_HiCoatis_VS_eRNA\\K562_Coatis_VS_eRNA_Venn2_1.pdf')





