############The gene classification criteria and associated gene expression information can be found in the Gene_classification_based_on_Hi-Coatis_signal_intensity.py script in the scripts folder.

#######K562
computeMatrix scale-regions --regionsFileName classify1.bed classify2.bed classify3.bed classify4.bed --scoreFileName /scratch/2024-07-29/bio-shenw/Ljniu/K562/K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw --outFileNameMatrix TSS_TES_outFileNameMatrix --regionBodyLength 6000 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --startLabel TSS --endLabel TES --skipZeros --numberOfProcessors 20 --outFileName TSS_TES_plotMatrix.gz
plotHeatmap -m TSS_TES_plotMatrix.gz -out K562_HIRPC_classify.pdf --heatmapHeight 15 --startLabel TSS --endLabel TES --plotFileFormat pdf  --regionsLabel c1 c2 c3 c4 --colorList white,red


######K562 Hemin
computeMatrix scale-regions --regionsFileName classify1_all.bed classify2_all.bed classify3_all.bed classify4_all.bed --scoreFileName /scratch/2024-09-09/bio-shenw/Ljniu/K562/K562_H_New/one-dimensional/mapping/bam2/all_reps/K562_Hemin_New_allreps_RPKM_10bp.bw --outFileNameMatrix TSS_TES_outFileNameMatrix --regionBodyLength 6000 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --startLabel TSS --endLabel TES --skipZeros --numberOfProcessors 20 --outFileName TSS_TES_plotMatrix.gz
plotHeatmap -m TSS_TES_plotMatrix.gz -out K562_Hemin_HIRPC_classify.pdf --heatmapHeight 15 --startLabel TSS --endLabel TES --plotFileFormat pdf  --regionsLabel c1 c2 c3 c4 --colorList white,red
