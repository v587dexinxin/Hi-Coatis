library(DESeq2)
Data <- read.table("H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_RNA-seq/readscount/union_all_reads_count.csv", header=T, row.names=1, sep=",")
Data <- cbind(Gene_Name=rownames(Data), as.data.frame(Data))
FPKM <- read.table("H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_RNA-seq/FPKM/union_all_FPKM.csv", header=T, row.names=1, sep=",")
FPKM <- FPKM[c('K562_WT_R1_FPKM','K562_WT_R2_FPKM','HCT116_WT_R1_FPKM','HCT116_WT_R2_FPKM')]
FPKM <- cbind(Gene_Name=rownames(FPKM), as.data.frame(FPKM))
sample <- read.table("H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_RNA-seq/readscount/K562_VS_HCT116.csv", header=T, row.names=1, com='', quote='', check.names=F, sep=",", colClasses="factor")
data <- Data[c('K562_WT_R1_Count','K562_WT_R2_Count','HCT116_WT_R1_Count','HCT116_WT_R2_Count')]
data <- data[rowSums(data)>2,]

ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data,
                                            colData = sample,  design= ~ conditions)

dds <- DESeq(ddsFullCountTable)

rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)


sampleA <- 'K562'
sampleB <- 'HCT116'


contrastV <- c("conditions", sampleA, sampleB)
res <- results(dds,  contrast=contrastV)

baseA <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleA]

if (is.vector(baseA)){
  baseMeanA <- as.data.frame(baseA)
} else {
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- sampleA
head(baseMeanA)


baseB <- counts(dds, normalized=TRUE)[, colData(dds)$conditions == sampleB]
if (is.vector(baseB)){
  baseMeanB <- as.data.frame(baseB)
} else {
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- sampleB
head(baseMeanB)


res <- cbind(baseMeanA, baseMeanB, as.data.frame(res))
res <- cbind(Gene_Name=rownames(res), as.data.frame(res))
res$padj[is.na(res$padj)] <- 1

head(res)
res <- merge(Data,res,by="Gene_Name")
res <- merge(res,FPKM,by="Gene_Name")
res <- res[, c('Gene_Name' , 'Chr' , 'Strand' , 'Start' , 'End' , 'baseMean' , 'log2FoldChange' , 'lfcSE' , 'stat' , 'pvalue' , 'padj' , 'K562_WT_R1_FPKM','K562_WT_R2_FPKM','HCT116_WT_R1_FPKM','HCT116_WT_R2_FPKM')]
res <- res[order(res$padj),]


write.table(res, file='H:/work/niulongjian/HiRPC_processed_data/K562_HCT116_RNA-seq/DEGs/K562_DEGs_WT_VS_Hemin.csv', sep=",", quote=F, row.names=F)





