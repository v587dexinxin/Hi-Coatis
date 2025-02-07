# -*- coding: utf-8 -*-
"""
Created on Sun Oct 20 11:22:58 2024

@author: lenovo
"""

import pandas as pd




def read_genome(filename):
    file_genome = open(filename)
    dict_genome = {}
    for line in file_genome:
        line = line.strip('\n')
        lists = list(line)
        if len(lists) > 1 and lists[0] == '>' :
            chrs = (line.split('>')[1]).split()[0]
            dict_genome[chrs] = []
        
        else :
            dict_genome[chrs].extend(lists)
    return dict_genome





def write_pos_to_seq(peaks , flanking , outfil):
    
    out = open(outfil , 'w')
    
    for i in peaks:
        g = i[0]
        start = i[1]
        end = i[2]
        f_start = start - flanking
        f_end = end + flanking
        seq = genome_h38[g][f_start:f_end]
        peak_name = i[3]
        gene_name = i[4]
        out.writelines('>' + g + ': ' + str(start) + '_' + str(end) + '_' +  peak_name + '_' + gene_name + '\n')
        for j in range(len(seq) // 70 + 1):
            n1 = j * 70
            n2 = (j + 1) * 70
            out.writelines(''.join(seq[n1 : n2]).upper() + '\n')
    out.close()
    
    
    



genome_h38 = read_genome('/scratch/2024-11-18/bio-shenw/ref/Human/hg38/hg38.fa')

# genes expression and the promoter regions
gene_exp_file = "/scratch/2024-11-18/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_HCT116/HCT116_union_gene_classify.txt"
gene_exp_df = pd.read_csv(gene_exp_file, header=0, sep="\t", index_col=None)
gene_exp_df["gID"] = gene_exp_df.index.values
gene_exp_df

c2 = gene_exp_df[gene_exp_df['Level'] == 'c2']
c2 = c2.sort_values(by = 'mean_FPKM')


peaks = pd.read_table('/scratch/2024-11-18/bio-shenw/Ljniu/HCT116/hg38/one-dimensional_1/mapping/common_union/HCT116_0.1FA_onedimensional_q0.05_union2_peaks.narrowPeak' , header = None)
peaks.columns = ['chr' , 'start' , 'end' , 'peak_name' , 'score' , 'strand' , 'sigvalue' , 'pvalue' , 'qvalue' , 'peak']



selected_peaks = ['HCT116_R2_q0.05_peak_32802' , 'HCT116_R2_q0.05_peak_5705' , 'HCT116_R2_q0.05_peak_49672' , 
                  'HCT116_R2_q0.05_peak_20715' , 'HCT116_R2_q0.05_peak_59050' , 'HCT116_R2_q0.05_peak_42657' , 
                  'HCT116_R2_q0.05_peak_9799' , 'HCT116_R2_q0.05_peak_50104' , 'HCT116_R2_q0.05_peak_50118' ,
                  'HCT116_R2_q0.05_peak_11494' , 'HCT116_R2_q0.05_peak_32928' , 'HCT116_R2_q0.05_peak_14440' ,
                  'HCT116_R2_q0.05_peak_25924' , 'HCT116_R2_q0.05_peak_43203' , 'HCT116_R2_q0.05_peak_55364' , 
                  'HCT116_R2_q0.05_peak_25796' , 'HCT116_R2_q0.05_peak_52543' , 'HCT116_R2_q0.05_peak_2578' , 
                  'HCT116_R2_q0.05_peak_16615' , 'HCT116_R2_q0.05_peak_16526' , 'HCT116_R2_q0.05_peak_11592' , 
                  'HCT116_R2_q0.05_peak_12526' , 'HCT116_R2_q0.05_peak_36656' , 'HCT116_R2_q0.05_peak_36664' , 
                  'HCT116_R2_q0.05_peak_37046' , 'HCT116_R2_q0.05_peak_31880' , 'HCT116_R2_q0.05_peak_32269' , 
                  'HCT116_R2_q0.05_peak_35972' , 'HCT116_R2_q0.05_peak_35182' , 'HCT116_R2_q0.05_peak_4114']




selected_peaks_1 = ['HCT116_R2_q0.05_peak_38366' , 'HCT116_R2_q0.05_peak_52845' , 'HCT116_R2_q0.05_peak_28265' , 
                  'HCT116_R2_q0.05_peak_32647' , 'HCT116_R2_q0.05_peak_30690' , 'HCT116_R2_q0.05_peak_31503' , 
                  'HCT116_R2_q0.05_peak_50707' , 'HCT116_R2_q0.05_peak_3730' , 'HCT116_R2_q0.05_peak_31987' ,
                  'HCT116_R2_q0.05_peak_32610' , 'HCT116_R2_q0.05_peak_16549' , 'HCT116_R2_q0.05_peak_34291' ,
                  'HCT116_R2_q0.05_peak_26100' , 'HCT116_R2_q0.05_peak_26415' , 'HCT116_R2_q0.05_peak_8530' , 
                  'HCT116_R2_q0.05_peak_43038' , 'HCT116_R2_q0.05_peak_26503' , 'HCT116_R2_q0.05_peak_4718' , 
                  'HCT116_R2_q0.05_peak_13773' , 'HCT116_R2_q0.05_peak_56192' , 'HCT116_R2_q0.05_peak_44020' , 
                  'HCT116_R2_q0.05_peak_5973' , 'HCT116_R2_q0.05_peak_2539' , 'HCT116_R2_q0.05_peak_20401' , 
                  'HCT116_R2_q0.05_peak_40261' , 'HCT116_R2_q0.05_peak_46909' , 'HCT116_R2_q0.05_peak_59396' , 
                  'HCT116_R2_q0.05_peak_53229' , 'HCT116_R2_q0.05_peak_48554' , 'HCT116_R2_q0.05_peak_49840' , 
                  'HCT116_R2_q0.05_peak_3485' , 'HCT116_R2_q0.05_peak_31381' , 'HCT116_R2_q0.05_peak_49117' , 
                  'HCT116_R2_q0.05_peak_42279' , 'HCT116_R2_q0.05_peak_51664' , 'HCT116_R2_q0.05_peak_4141' , 
                  'HCT116_R2_q0.05_peak_8809' , 'HCT116_R2_q0.05_peak_40945' , 'HCT116_R2_q0.05_peak_30851' ,
                  'HCT116_R2_q0.05_peak_59300' , 'HCT116_R2_q0.05_peak_56740' , 'HCT116_R2_q0.05_peak_34544' ,
                  'HCT116_R2_q0.05_peak_22880' , 'HCT116_R2_q0.05_peak_57020' , 'HCT116_R2_q0.05_peak_32742' , 
                  'HCT116_R2_q0.05_peak_39211' , 'HCT116_R2_q0.05_peak_17086' , 'HCT116_R2_q0.05_peak_44452' , 
                  'HCT116_R2_q0.05_peak_56743' , 'HCT116_R2_q0.05_peak_30099' , 'HCT116_R2_q0.05_peak_42295' , 
                  'HCT116_R2_q0.05_peak_43232' , 'HCT116_R2_q0.05_peak_23469' , 'HCT116_R2_q0.05_peak_54432' , 
                  'HCT116_R2_q0.05_peak_15132' , 'HCT116_R2_q0.05_peak_28254' , 'HCT116_R2_q0.05_peak_44488' , 
                  'HCT116_R2_q0.05_peak_15271' , 'HCT116_R2_q0.05_peak_32908' , 'HCT116_R2_q0.05_peak_6308']




selected_pos = []
for peak_name in selected_peaks_1:
    peak = peaks[peaks['peak_name'] == peak_name]
    g = peak.iloc[0]['chr']
    start = peak.iloc[0]['start']
    end = peak.iloc[0]['end']
    tmp_gene = c2[c2['Chr'] == g]
    mask = (tmp_gene['pstart'] <= end) & (tmp_gene['pend'] >= start)
    gene = tmp_gene[mask]
    gene_name = gene.iloc[0]['Gene_Name']
    selected_pos.append((g , start , end , peak_name , gene_name))
    
    
    
selected_pos
    
    
    
write_pos_to_seq(selected_pos , 1000 , '/scratch/2024-11-18/bio-shenw/Ljniu/K562/Confirmation_Experiment/classify_HCT116_cl2/classify_HCT116_classify2_gene_peaks+-1K_1.fasta')
write_pos_to_seq(selected_pos , 500 , '/scratch/2024-11-18/bio-shenw/Ljniu/K562/Confirmation_Experiment/classify_HCT116_cl2/classify_HCT116_classify2_gene_peaks+-0.5K_1.fasta')









