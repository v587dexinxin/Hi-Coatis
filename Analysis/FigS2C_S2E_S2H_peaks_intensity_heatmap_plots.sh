
##########
computeMatrix scale-regions --regionsFileName /scratch/2023-12-11/bio-shenw/ref/Human/annotation/hg38.ncbiRefSeq/hg38.ncbiRefSeq.bed --scoreFileName ../../K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw ../../../../literature/K562/PRO_seq/hg38/mapping/K562_PRO_seq_SRR8137173_10bp.bw ../../../../literature/K562/R_chip/mapping/bam2/K562_D210N_V5ChIP_SRR5379780_Rep1_RPKM_10bp.bw ../../../../literature/K562/ChIP_seq/hg38/K562_POLR2A_hg38_ENCFF124WLE.bigWig --outFileNameMatrix TSS_TES_outFileNameMatrix --regionBodyLength 6000 --beforeRegionStartLength 3000 --afterRegionStartLength 3000 --startLabel TSS --endLabel TES --skipZeros --numberOfProcessors 20 --outFileName TSS_TES_plotMatrix.gz



plotHeatmap -m TSS_TES_plotMatrix.gz -out K562_0.1FA_PRO_RChip_RNAP2A_TSS_TES.pdf --heatmapHeight 15 --startLabel TSS --endLabel TES --plotFileFormat pdf  --regionsLabel all_genes --samplesLabel 0.1FA PRO-seq RChIP Pol2A --colorList white,red 


plotHeatmap -m TSS_TES_plotMatrix_norm_1.gz -out K562_0.1FA_PRO_RChip_RNAP2A_TSS_TES_norm_1.pdf --heatmapHeight 15 --startLabel TSS --endLabel TES --plotFileFormat pdf  --regionsLabel all_genes --samplesLabel 0.1FA PRO-seq RChIP Pol2A --colorList white,red --zMax 300 --yMax 400 --averageTypeSummaryPlot mean
