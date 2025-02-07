import pandas as pd
import matplotlib.pyplot as plt


def calculate_overlap(row):
    if row['chr'] != g:
        return 0  # 不同染色体无重叠
    overlap_start = max(row['start'], start)
    overlap_end = min(row['end'], end)
    overlap = max(0, overlap_end - overlap_start)  # 确保不为负
    return overlap



#######Load files
peaks = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak' , usecols = (0 , 1 , 2 , 4) , header = None)
peaks.columns = ['chr' , 'start' , 'end' , 'score']
peaks_pro = pd.read_table('H:\\work\\literature_data\\K562\\PRO-seq\\hg38\\peaks\\K562_PRO-seq_SRR8137173_q0.01_peaks.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks_pro.columns = ['chr' , 'start' , 'end']
peaks_rchip = pd.read_table('H:\\work\\literature_data\\K562\\R_chip\\macs2_remove_input\\K562_D210N_V5ChIP_union_q0.05_peaks.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks_rchip.columns = ['chr' , 'start' , 'end']
peaks_polII = pd.read_table('H:\\work\\literature_data\\K562\\ChIP-seq\\hg38\\peaks\\K562_PoLR2A_hg38_ENCFF286QTF.bed' , usecols = (0 , 1 , 2) , header = None)
peaks_polII.columns = ['chr' , 'start' , 'end']

Peaks = {'Coatis':peaks , 'PRO' : peaks_pro , 'RChIP':peaks_rchip , 'PolII':peaks_polII } 

hmm = pd.read_table('H:\\work\\literature_data\\chromHMM\\HiCoatis_VS_Capture_HiC_VS_ChIAPET_loops\\调控元件注释\\wgEncodeBroadHmmK562HMM.hg19_to_hg38.bed' , header = None , usecols=(0 , 1 , 2 , 3))
hmm.columns = ['chr' , 'start' , 'end' , 'Anno']


#######Peaks Annotation
chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
Anno = {}

for k in Peaks.keys():
    out = open('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\Peaks_Anno\\' + k + '_peaks_Anno.bed' , 'w')
    n = 0; Anno[k] = []
    for g in chrom:
        print (g)
        tmp_peaks = Peaks[k][Peaks[k]['chr'] == g]
        tmp_hmm = hmm[hmm['chr'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            mask = (tmp_hmm['start'] <= end) & (tmp_hmm['end'] >= start)
            overlap = tmp_hmm[mask]
            if len(overlap) != 0:
                n += 1
                length = overlap.apply(calculate_overlap, axis=1)
                anno = overlap.loc[length.idxmax()]['Anno']
                Anno[k].append(anno)
                out.writelines('\t'.join([g , str(start) , str(end) , anno]) + '\n')
                
                
            else:
                out.writelines('\t'.join([g , str(start) , str(end) , 'NA']) + '\n')
                pass
                
    out.close()
    


#########barplots
matrix = pd.DataFrame([])

for k in ['Coatis' , 'PRO' , 'RChIP' , 'PolII']:
    df = pd.DataFrame({'names': Anno[k]})
    a = df.value_counts()
    a['4_5_Strong_Enhancer'] = a['4_Strong_Enhancer'] + a['5_Strong_Enhancer']
    a['6_7_Weak Enhancer'] = a['6_Weak_Enhancer'] + a['7_Weak_Enhancer']
    a['14_15_Repetitive/CNV'] = a['14_Repetitive/CNV'] + a['15_Repetitive/CNV']
    
    a_df = pd.DataFrame(a.values, index=a.index, columns=[k])
    matrix = pd.concat((matrix , a_df) , axis = 1)
    
matrix = matrix.drop(['4_Strong_Enhancer' , '5_Strong_Enhancer' , '6_Weak_Enhancer' , '7_Weak_Enhancer' , '14_Repetitive/CNV' , '15_Repetitive/CNV'])


colors = plt.cm.tab20.colors[:15]
    
matrix = matrix.T
matrix.plot(kind="bar", stacked=True , color=colors)
    

