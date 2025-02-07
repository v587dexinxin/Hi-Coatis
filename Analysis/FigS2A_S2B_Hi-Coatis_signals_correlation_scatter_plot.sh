##########K562 Coatis signals 2 replicates
multiBigwigSummary bins -b ../../K562_0.1_FA/one-dimensional/mapping/GZ23089190-562w1_1-562w1_1_RPKM_10bp.bw ../../K562_0.1_FA/one-dimensional/mapping/GZ23089191-562w2_2-562w2_2_RPKM_10bp.bw --labels K562_RPC_R1 K562_RPC_R2 -out K562_RPC_scores_per_bin.npz --outRawCounts K562_RPC_scores_per_bin.tab
plotCorrelation -in K562_RPC_scores_per_bin.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of K562 peaks 2 replicates" --whatToPlot scatterplot  --plotNumbers -o K562_SpearmanCorr_2replicates.pdf  --outFileCorMatrix K562_PearsonCorr.tab --log1p --removeOutliers


##########K562 Hi-Coatis one dimension signal VS Coatis signals
multiBigwigSummary bins -b ../../K562_0.1FA_allreps_RPKM_10bp.bw ../../../../old/one-dimensional/mapping/bam2/merged_2/K562_0.1FA_subset3000w_allreps_RPKM_10bp.bw  --binSize 50000 --labels RPC HiRPC  -out K562_scores_per_bin.npz --outRawCounts K562_scores_per_bin.tab
plotCorrelation -in K562_scores_per_bin.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of K562 peaks" --whatToPlot scatterplot  --plotNumbers -o K562_SpearmanCorr.pdf  --outFileCorMatrix K562_PearsonCorr.tab --log1p --removeOutliers
