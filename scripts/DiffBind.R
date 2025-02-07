####K562_VS_Hemin
library(DiffBind)
dbObj <- dba(sampleSheet="H:/work/niulongjian/HiRPC_processed_data/K562/K562_HiRPC_Hemin_new/peaks/DiffBind/Samplesheet_WT_Hemin_HiRPC.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("H:/work/niulongjian/HiRPC_processed_data/K562/K562_HiRPC_Hemin_new/peaks/DiffBind/WT_VS_Hemin_new_PCA.pdf", pointsize=10)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

#dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="H:/work/niulongjian/HiRPC_processed_data/K562/K562_HiRPC_Hemin_new/peaks/DiffBind/HiRPC_WT_vs_Hemin_new_deseq2.csv", sep=",", quote=F,row.names = FALSE)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="H:/work/niulongjian/HiRPC_processed_data/K562/K562_HiRPC_Hemin_new/peaks/DiffBind/HiRPC_WT_vs_Hemin_new_edgeR.csv", sep=",", quote=F, row.names = FALSE)


####K562_VS_Hemin_H3K27ac
library(DiffBind)
dbObj <- dba(sampleSheet="H:/work/niulongjian/HiRPC_processed_data/K562/K562_ChIP-seq/peaks/DiffBind/Samplesheet_WT_Hemin_H3K27ac.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("H:/work/niulongjian/HiRPC_processed_data/K562/K562_ChIP-seq/peaks/DiffBind/WT_VS_Hemin_H3K27ac_PCA.pdf", pointsize=10)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

#dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="H:/work/niulongjian/HiRPC_processed_data/K562/K562_ChIP-seq/peaks/DiffBind/H3K27ac_WT_vs_Hemin_deseq2.csv", sep=",", quote=F,row.names = FALSE)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="H:/work/niulongjian/HiRPC_processed_data/K562/K562_ChIP-seq/peaks/DiffBind/H3K27ac_WT_vs_Hemin_edgeR.csv", sep=",", quote=F, row.names = FALSE)





####K562_VS_HCT116

library(DiffBind)
dbObj <- dba(sampleSheet="H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_HiRPC_0.1FA/DiffBind/K562_VS_HCT116.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
pdf("H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_HiRPC_0.1FA/DiffBind/K562_VS_HCT116_PCA.pdf", pointsize=10)
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()

#dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
out <- as.data.frame(comp1.deseq)
write.table(out, file="H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_HiRPC_0.1FA/DiffBind/HiRPC_K562_vs_HCT116_deseq2.csv", sep=",", quote=F,row.names = FALSE)
out <- as.data.frame(comp1.edgeR)
write.table(out, file="H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_HiRPC_0.1FA/DiffBind/HiRPC_K562_vs_HCT116_edgeR.csv", sep=",", quote=F, row.names = FALSE)







