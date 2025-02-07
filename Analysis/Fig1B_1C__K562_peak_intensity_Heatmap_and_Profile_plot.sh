####K562 Heatmap of signal intensity under different crosslinking conditions.
##Union of different peaks.
cat ../../../K562_0.1_FA/one-dimensional/macs2/common_union/K562_0.1FA_onedimensional_q0.05_union_peaks.bed ../../../K562_no_FA/data_new/mapping/macs2/common_union/K562_noFA_union_q0.05_peaks.bed ../../../K562_1FA/one-dimensional/Combined_562_1FA_q0.05_peaks.bed ../../../K562_RnAT/one-dimensional_new/mapping/macs2/common_union/K562_RnAT_union_q0.05_peaks.bed > union_K562_0.1FA_onFA_1FA_RnAT_all.bed
sort -k1,1 -k2,2n union_K562_0.1FA_onFA_1FA_RnAT_all.bed > union_K562_0.1FA_onFA_1FA_RnAT_all_sorted.bed
bedtools merge -i union_K562_0.1FA_onFA_1FA_RnAT_all_sorted.bed > union_K562_0.1FA_onFA_1FA_RnAT_merge_sorted.bed

##Heatmap and Profile plot.
computeMatrix reference-point -S ../../K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw ../../K562_no_FA/data_new/mapping/bam2/merged_2/K562_noFA_allreps_RPKM_10bp.bw ../../K562_1FA/one-dimensional/Combined_562_1FA_RPKM_10bp.bw ../../K562_RnAT/one-dimensional_new/mapping/bam2/merged_2/K562_RnAT_allreps_RPKM_10bp.bw -R ./union_peaks/union_K562_0.1FA_onFA_1FA_RnAT_merge_sorted.bed --referencePoint center -a 2000 -b 2000 -out K562_0.1FA_noFA_1FA_RnAT_union.gz --skipZeros --missingDataAsZero --numberOfProcessors 20
plotHeatmap -m K562_0.1FA_noFA_1FA_RnAT_union.gz -out K562_0.1FA_noFA_1FA_RnAT_union.pdf --heatmapHeight 15 --refPointLabel peaks_center  --regionsLabel union_peaks --samplesLabel 0.1FA noFA 1FA RnAT --colorList white,red


####K562 Heatmap of signal intensity under different cells number.
##The normalization method for peak intensity across different cell numbers can be found in the scripts folder.
computeMatrix reference-point -S ../../../K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_50k_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_25k_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_5k_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_1k_allreps_RPKM_10bp.bw -R ../../../K562_0.1_FA/one-dimensional/macs2/common_union/K562_0.1FA_onedimensional_q0.05_union_peaks.bed --referencePoint center -a 2000 -b 2000 -out K562_lowcells_union_1M.gz --skipZeros --missingDataAsZero
plotHeatmap -m K562_lowcells_union_1M_norm.gz -out K562_lowcells_union_1M_peaks_heatmap_norm.pdf --heatmapHeight 15 --refPointLabel peaks_center  --regionsLabel 1M_peaks --samplesLabel 1M 50k 25k 5k 1k --colorList white,red


#####K562 Heatmap of signal intensity under different datasets.
##The normalization method for peak intensity across different cell numbers can be found in the scripts folder.
computeMatrix reference-point -S ../../../K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw ../../../../../literature/K562/PRO_seq/hg38/mapping/K562_PRO_seq_SRR8137173_10bp.bw ../../../../../literature/K562/R_chip/mapping/bam2/K562_D210N_V5ChIP_SRR5379780_Rep1_RPKM_10bp.bw ../../../../../literature/K562/ChIP_seq/hg38/K562_POLR2A_hg38_ENCFF124WLE.bigWig -R ../../../K562_0.1_FA/one-dimensional/macs2/common_union/K562_0.1FA_onedimensional_q0.05_common_peaks.narrowPeak --referencePoint center -a 2000 -b 2000 -out K562_0.1FA_union.gz --skipZeros --missingDataAsZero
plotHeatmap -m K562_0.1FA_union_normrpc.gz -out K562_0.1FA_PRO_RChip_RNAP2A_heatmap_normrpc.pdf --heatmapHeight 15 --refPointLabel peaks_center  --regionsLabel union_peaks --samplesLabel 0.1FA PRO-seq RChIP PolII --colorList white,red





