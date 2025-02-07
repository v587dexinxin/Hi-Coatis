# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 22:02:25 2024

@author: lenovo
"""

import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt





filename = "/scratch/2024-12-16/bio-shenw/Ljniu/K562/K562_0.1_FA/loops/ChIA-pipline_8K+/one_anchors_binding_loops/K562_merge5_0.1FA.hg38_peaks_one_anchors_binding_loops.bedpe"
hicoatis_contact_df = pd.read_csv(filename, sep="\t", names=["chr1", "s1", "e1", "chr2","s2", "e2", "count","-log10qval"])

filename = "/scratch/2024-12-16/bio-shenw/literature/K562/ChIA-PET/loops/K562_RNAPII_ChIAPET_hg38_ENCFF511QFN_one_anchor_binding_loops.bedpe"
chiapet_contact_df = pd.read_csv(filename, sep="\t", names=["chr1", "s1", "e1", "chr2","s2", "e2"])


filename = "/scratch/2024-12-16/bio-shenw/literature/K562/Capture_HiC/K562/loops/K562_Capture_hic.hg38_peaks_one_anchors_binding_loops.bedpe"
capture_contact_df = pd.read_csv(filename, sep="\t", names=["chr1", "s1", "e1", "chr2","s2", "e2", "count","-log10qval"])



chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']
flanking = 5000

hicoatis_contact_df['s1_f'] = hicoatis_contact_df['s1'] - flanking
hicoatis_contact_df['e1_f'] = hicoatis_contact_df['e1'] + flanking
hicoatis_contact_df['s2_f'] = hicoatis_contact_df['s2'] - flanking
hicoatis_contact_df['e2_f'] = hicoatis_contact_df['e2'] + flanking



chiapet_contact_df['s1_f'] = chiapet_contact_df['s1'] - flanking
chiapet_contact_df['e1_f'] = chiapet_contact_df['e1'] + flanking
chiapet_contact_df['s2_f'] = chiapet_contact_df['s2'] - flanking
chiapet_contact_df['e2_f'] = chiapet_contact_df['e2'] + flanking




capture_contact_df['s1_f'] = capture_contact_df['s1'] - flanking
capture_contact_df['e1_f'] = capture_contact_df['e1'] + flanking
capture_contact_df['s2_f'] = capture_contact_df['s2'] - flanking
capture_contact_df['e2_f'] = capture_contact_df['e2'] + flanking






def loops_overlap(loops1 , loops2 , flanking):
    common = []
    for g in chrom:
        print (g)
        tmp1 = loops1[loops1['chr1'] == g]
        tmp2 = loops2[loops2['chr1'] == g]
        for i in tmp1.index:
            start1 = tmp1.loc[i]['s1_f']
            end1 = tmp1.loc[i]['e1_f']
            start2 = tmp1.loc[i]['s2_f']
            end2 = tmp1.loc[i]['e2_f']
            mask1 = (tmp2['s1_f'] <= end1) & (tmp2['e1_f'] >= start1)
            overlap1 = tmp2[mask1]
            if len(overlap1) != 0:
                mask2 = (overlap1['s2_f'] <= end2) & (overlap1['e2_f'] >= start2)
                overlap2 = overlap1[mask2]
                if len(overlap2) != 0:
                    common.append((g , start1 + flanking , end1 - flanking , g , start2 + flanking , end2 - flanking))
    common = pd.DataFrame(common , columns=["chr1", "s1", "e1", "chr2","s2", "e2"])
    return common


overlap12 = loops_overlap(chiapet_contact_df , hicoatis_contact_df , flanking)
overlap12.to_csv('/scratch/2024-12-16/bio-shenw/Ljniu/K562/plots/plots_New_number_bylxx_xjs/loops_overlap/K562_ChIAPET_VS_HiCoatis_loops_common.bedpe' , header = None , index = None , sep = '\t')

overlap13 = loops_overlap(chiapet_contact_df , capture_contact_df , flanking)
overlap13.to_csv('/scratch/2024-12-16/bio-shenw/Ljniu/K562/plots/plots_New_number_bylxx_xjs/loops_overlap/K562_ChIAPET_VS_Capture_loops_common.bedpe' , header = None , index = None , sep = '\t')

overlap23 = loops_overlap(hicoatis_contact_df , capture_contact_df , flanking)
overlap23.to_csv('/scratch/2024-12-16/bio-shenw/Ljniu/K562/plots/plots_New_number_bylxx_xjs/loops_overlap/K562_HiCoatis_VS_Capture_loops_common.bedpe' , header = None , index = None , sep = '\t')


overlap23['s1_f'] = overlap23['s1'] - flanking
overlap23['e1_f'] = overlap23['e1'] + flanking
overlap23['s2_f'] = overlap23['s2'] - flanking
overlap23['e2_f'] = overlap23['e2'] + flanking








                        
overlap123 = loops_overlap(chiapet_contact_df , overlap23 , flanking)                  
overlap123.to_csv('/scratch/2024-12-16/bio-shenw/Ljniu/K562/plots/plots_New_number_bylxx_xjs/loops_overlap/K562_ChIAPET_VS_HiCoatis_VS_Capture_loops_common3.bedpe' , header = None , index = None , sep = '\t')
                   

#####################windows Venn3 plots

filename = "H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\loops\\K562_merge5_0.1FA.hg38_peaks_one_anchors_binding_loops.bedpe"
hicoatis_contact_df = pd.read_csv(filename, sep="\t", names=["chr1", "s1", "e1", "chr2","s2", "e2", "count","-log10qval"])

filename = "H:\\work\\literature_data\\K562\\ChIA-PET\\hg38\\K562_RNAPII_ChIAPET_hg38_ENCFF511QFN_one_anchor_binding_loops.bedpe"
chiapet_contact_df = pd.read_csv(filename, sep="\t", names=["chr1", "s1", "e1", "chr2","s2", "e2"])


filename = "H:\\work\\literature_data\\K562\\Capture Hi-C\\loops\\K562_Capture_hic.hg38_peaks_one_anchors_binding_loops.bedpe"
capture_contact_df = pd.read_csv(filename, sep="\t", names=["chr1", "s1", "e1", "chr2","s2", "e2", "count","-log10qval"])



overlap12 = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\K562_loops_Venn3\\K562_ChIAPET_VS_HiCoatis_loops_common.bedpe' , header = None , sep = '\t')
overlap13 = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\K562_loops_Venn3\\K562_ChIAPET_VS_Capture_loops_common.bedpe' , header = None , sep = '\t')
overlap23 = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\K562_loops_Venn3\\K562_HiCoatis_VS_Capture_loops_common.bedpe' , header = None , sep = '\t')
overlap123 = pd.read_csv('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\K562_loops_Venn3\\K562_ChIAPET_VS_HiCoatis_VS_Capture_loops_common3.bedpe' , header = None , sep = '\t')




overlap = {'100' : len(hicoatis_contact_df) - len(overlap12) - len(overlap23) + len(overlap123),
           '010' : len(chiapet_contact_df) - len(overlap12) - len(overlap13) + len(overlap123),
           '001' : len(capture_contact_df) - len(overlap13) - len(overlap23) + len(overlap123),
           '110' : len(overlap12) - len(overlap123),
           '101' : len(overlap23) - len(overlap123),
           '011' : len(overlap13) - len(overlap123),
           '111' : len(overlap123)}



venn = venn3(subsets=overlap , set_labels=('HiCoatis', 'ChIAPET', 'Capture HiC'))




