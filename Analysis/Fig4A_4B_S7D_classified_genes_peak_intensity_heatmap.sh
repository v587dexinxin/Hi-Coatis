############

#######K562
computeMatrix scale-regions --regionsFileName classify1.bed classify2.bed classify3.bed classify4.bed --scoreFileName /scratch/2024-07-29/bio-shenw/Ljniu/K562/K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw --outFileNameMatrix TSS_TES_outFileNameMatrix --regionBodyLength 6000 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --startLabel TSS --endLabel TES --skipZeros --numberOfProcessors 20 --outFileName TSS_TES_plotMatrix.gz
plotHeatmap -m TSS_TES_plotMatrix.gz -out K562_HIRPC_classify.pdf --heatmapHeight 15 --startLabel TSS --endLabel TES --plotFileFormat pdf  --regionsLabel c1 c2 c3 c4 --colorList white,red


######K562 Hemin
