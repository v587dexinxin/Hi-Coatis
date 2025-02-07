#####FRiP calculation command line.
KAS-Analyzer FRiP -o K562_Coatis_PRO_RChIP_PolII_FRiP -p peaks_files.txt -l labels.txt -k bed.txt

##peaks_files.txt
/scratch/2024-12-16/bio-shenw/Ljniu/K562/K562_0.1_FA/one-dimensional/macs2/common_union/K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak
/scratch/2024-12-16/bio-shenw/literature/K562/PRO_seq/hg38/mapping/K562_PRO_seq_SRR8137173_q0.05_peaks.narrowPeak
/scratch/2024-12-16/bio-shenw/literature/K562/R_chip/mapping/macs2_remove_input/common_union/K562_D210N_V5ChIP_union_q0.05_peaks.narrowPeak
/scratch/2024-12-16/bio-shenw/literature/K562/ChIP_seq/hg38/processed_byself/mapping/macs2/common_union/K562_POLR2A_ChIP_union2_q0.05_peaks.narrowPeak
##bed.txt
/scratch/2024-12-16/bio-shenw/Ljniu/K562/K562_0.1_FA/one-dimensional/mapping/all_reps/K562_0.1FA_allreps.sorted.q20.rmdup.bed
/scratch/2024-12-16/bio-shenw/literature/K562/PRO_seq/hg38/mapping/K562_PRO_seq_SRR8137173.sorted.q20.bed
/scratch/2024-12-16/bio-shenw/literature/K562/R_chip/mapping/bam2/all_reps/K562_RChIP_all_reps.bed
/scratch/2024-12-16/bio-shenw/literature/K562/ChIP_seq/hg38/processed_byself/mapping/bam2/all_reps/K562_POLR2A_all_reps.bed
##labels.txt
K562_Coatis
K562_PRO_seq
K562_RChIP
K562_PolII


######K562_Coatis_PRO_RChIP_PolII_FRiP_bubble_plots python scripts
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

FRiP = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\one-dimensional-FRIP\\K562_Coatis_PRO_RChIP_PolII_FRiP_FRiP.txt' , header=0)
FRiP = FRiP[FRiP['Types'] == 'Inside']
peaks_coatis = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\K562_0.1FA_onedimensional_q0.05_union_peaks.narrowPeak' , usecols = (0 , 1 , 2 , 4) , header = None)
peaks_coatis.columns = ['chr' , 'start' , 'end' , 'score']
peaks_pro = pd.read_table('H:\\work\\literature_data\\K562\\PRO-seq\\hg38\\peaks\\K562_PRO-seq_SRR8137173_q0.01_peaks.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks_pro.columns = ['chr' , 'start' , 'end']
peaks_rchip = pd.read_table('H:\\work\\literature_data\\K562\\R_chip\\macs2_remove_input\\K562_D210N_V5ChIP_union_q0.05_peaks.narrowPeak' , usecols = (0 , 1 , 2) , header = None)
peaks_rchip.columns = ['chr' , 'start' , 'end']
peaks_polII = pd.read_table('H:\\work\\literature_data\\K562\\ChIP-seq\\hg38\\peaks\\K562_PoLR2A_hg38_ENCFF286QTF.bed' , usecols = (0 , 1 , 2) , header = None)
peaks_polII.columns = ['chr' , 'start' , 'end']


# 数据
methods = list(FRiP['labels'])
frip_values = list(FRiP['Percentage'])
peaks = [len(peaks_coatis), len(peaks_pro) , len(peaks_rchip) , len(peaks_polII)]

# 颜色和大小
colors = plt.cm.Blues(np.interp(frip_values, (min(frip_values), max(frip_values)), (0.2, 1)))
sizes = [frip * 20 for frip in frip_values]  # 大小与frip_values成正比

# 创建图形
fig, ax = plt.subplots(figsize=(8, 5))

# 绘制散点图
for i in [3 , 2 , 1 , 0]:
    method = methods[i]
    ax.hlines(y=method, xmin=0, xmax=peaks[i], color='black', linestyles='solid', linewidth=1)
    ax.scatter(peaks[i], methods[i], s=sizes[i], color=colors[i], label=f"FRiP: {frip_values[i]:.2f}", edgecolors="black", zorder=2)
    
# # 添加颜色条
# sm = plt.cm.ScalarMappable(cmap="Blues", norm=plt.Normalize(vmin=min(frip_values), vmax=max(frip_values)))
# sm.set_array([])
# cb = fig.colorbar(sm, ax=ax, orientation='horizontal', label='FRiP Values', shrink=0.6, pad=0.2)

# 添加图例和标题
legend = ax.legend(title="Methods", loc="lower right", bbox_to_anchor=(1.55, 0), handletextpad=0.5, borderpad=0.5, labelspacing=1.5)
ax.set_title("Peak Numbers by Method")
ax.set_xlabel("Peaks")
ax.set_xlim(0, max(peaks) + 10000)
ax.set_ylabel("Methods")
ax.set_ylim(-0.5, 3.5)

# 调整布局
plt.tight_layout()
plt.show()


