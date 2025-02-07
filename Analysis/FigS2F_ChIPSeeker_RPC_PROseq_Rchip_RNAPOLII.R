library(ChIPseeker)
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("org.Mm.eg.db")
library("GenomicFeatures")

spompe <- makeTxDbFromGFF('H:/work/literature_data/genome/hg38/genecode/gencode.v40.chr_patch_hapl_scaff.annotation.gtf')
files <- list(K562_0.1FA = c('H:/work/niulongjian/HiRPC_processed_data/K562/K562_HiRPC_0.1FA/one-dimensional/peaks/K562_0.1FA_onedimensional_q0.05_union_peaks.bed') ,
              PRO_seq = c('H:/work/literature_data/K562/PRO-seq/hg38/peaks/K562_PRO-seq_SRR8137173_q0.05_peaks.bed') , 
              RChIP = c('H:/work/literature_data/K562/R_chip/macs2_remove_input/K562_D210N_V5ChIP_union_q0.05_peaks.bed') , 
              RNAPoLII = c('H:/work/literature_data/K562/ChIP-seq/hg38/processed_byself/peaks/K562_POLR2A_ChIP_union2_q0.05_peaks.bed'))


##plotAnnoBar
peakAnnoList <- lapply(files , annotatePeak , TxDb = spompe , tssRegion=c(-2000, 2000) , overlap = "all" , addFlankGeneInfo = TRUE, flankDistance = 100000,verbose = FALSE)


pdf('H:/work/niulongjian/HiRPC_processed_data/plots/plots_New_number_bylxx_xjs/Fig2D_HiRPC_proseq_rchip_rnapolII_peaksAnno/K562_peakAnno_score.pdf')
plotAnnoBar(peakAnnoList)
dev.off()


coatis <- as.data.frame(peakAnnoList$K562_0.1FA)[, 1:6]
PRO_seq <- as.data.frame(peakAnnoList$PRO_seq)[, 1:6]
RChIP <- as.data.frame(peakAnnoList$RChIP)[, 1:6]
RNAPoLII <- as.data.frame(peakAnnoList$RNAPoLII)[, 1:6]


write.table(coatis, file = "H:/work/niulongjian/HiRPC_processed_data/plots/plots_New_number_bylxx_xjs/Fig2D_HiRPC_proseq_rchip_rnapolII_peaksAnno/K562_Coatis_peaks_Anno.bed", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(PRO_seq, file = "H:/work/niulongjian/HiRPC_processed_data/plots/plots_New_number_bylxx_xjs/Fig2D_HiRPC_proseq_rchip_rnapolII_peaksAnno/K562_PRO_seq_peaks_Anno.bed", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(RChIP, file = "H:/work/niulongjian/HiRPC_processed_data/plots/plots_New_number_bylxx_xjs/Fig2D_HiRPC_proseq_rchip_rnapolII_peaksAnno/K562_RChIP_peaks_Anno.bed", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(RNAPoLII, file = "H:/work/niulongjian/HiRPC_processed_data/plots/plots_New_number_bylxx_xjs/Fig2D_HiRPC_proseq_rchip_rnapolII_peaksAnno/K562_RNAPoLII_peaks_Anno.bed", sep = "\t", row.names = FALSE, quote = FALSE)



