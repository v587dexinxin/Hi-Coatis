# Hi-Coatis
we developed Hi-Coatis (High-throughput Capture of Actively Transcribed Region-Interacting Sequences), a novel chromatin conformation capture technology that seamlessly integrates the detection of active transcription signals with three-dimensional chromatin interaction studies. Hi-Coatis operates without antibodies or probes, enabling low-input cells experiments with high resolution and robustness, capturing more than 93% of expressed genes and over 60,000 regulatory loci. 

This is the code used to make a pipeline for analysing Hi-Coatis data. In this work, all data generated in this study have been deposited to the Genome Sequence Archive in BGI Data Center (https://bigd.big.ac.cn/gsa/) with the GSA accession number HRA009396. We use the human hg38 genome as reference, and we obtained the public Hi-C, Chip-seq, ATAC-seq and RNA-seq datasets of  K562 and HCT116 cell lines from ENDCODE database.

# Software
Trimmomatic (v0.39)
HISAT2 (v2.1.1)
SAMtools (v1.3.1)
HTSeq (v2.0.2)
DESeq2 (v1.34.0)
Bowtie2 (v2.2.5)
Deeptools (v3.5)
MACS2 (v2.2.8)
Distiller-nf mapping pipeline (v0.3.4)
Bwa mem (v0.7.17-r1188)
Pairtools (v0.3.0)
Juicer Tools (v1.22.01)
Cooler (0.10.2)
HiCExplorer (v2.2.1.1)
Stripenn (v1.1.65)
CrossMap (v0.7.0)
ChIPSeeker (v1.30.3)
FIMO (v5.5.3)
ClusterProfiler (v0.1.dev14)
HiCPeaks (v0.3.7)
