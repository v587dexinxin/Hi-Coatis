multiBigwigSummary bins -b ../../../K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_1k_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_5k_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_25k_allreps_RPKM_10bp.bw ../../../K562_lowcells/one-dimensional/all_reps/signals/K562_lowcells_50k_allreps_RPKM_10bp.bw --labels 1M 1k 5k 25k 50k  -out K562_lowcells_scores_per_bin.npz --outRawCounts K562_lowcells_scores_per_bin.tab
plotCorrelation -in K562_lowcells_scores_per_bin.npz --corMethod spearman --skipZeros --plotTitle "Spearman Correlation of K562 lowcells peaks" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o K562_lowcells_SpearmanCorr.pdf  --outFileCorMatrix K562_lowcells_SpearmanCorr.tab --zMin 0.8 --zMax 1
