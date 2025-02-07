# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 21:12:32 2024

@author: lenovo
"""



import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


data = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\one-dimensional\\peaks\\Peaks_Anno\\Coatis_peaks_Anno.bed' , header = None)
data.columns = ['chr' , 'start' , 'end' , 'Anno']

loops = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\K562\\K562_HiRPC_0.1FA\\loops\\K562_merge5_0.1FA.hg38_peaks_one_anchors_binding_loops.bedpe' , header = None)
loops.columns = ['chr1' , 'start1' , 'end1' , 'chr2' , 'start2' , 'end2' , 'count' , 'pvalue']

gene_exp_df = pd.read_table('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\HiRPC_classify\\classify_K562\\K562_union_gene_classify.txt' , header = 0)



classify = {}

for i in data.index:
    g = data.loc[i]['chr']
    start = data.loc[i]['start']
    end = data.loc[i]['end']
    anno = data.loc[i]['Anno']
    if  not pd.isna(anno):
        if anno not in classify:
            classify[anno] = []
            classify[anno].append((g , start , end , anno))
        else:
            classify[anno].append((g , start , end , anno))
    else:
        pass
classify['4_5_Strong_Enhancer'] = classify['4_Strong_Enhancer'] + classify['5_Strong_Enhancer']
classify['6_7_Weak_Enhancer'] = classify['6_Weak_Enhancer'] + classify['7_Weak_Enhancer']
classify['14_15_Repetitive/CNV'] = classify['14_Repetitive/CNV'] + classify['15_Repetitive/CNV']
    

for k in classify.keys():
    classify[k] = pd.DataFrame(classify[k] , columns=['chr' , 'start' , 'end' , 'Anno'])
    
    
        
keys = ['1_Active_Promoter' , '2_Weak_Promoter' , '3_Poised_Promoter' , '4_5_Strong_Enhancer' , '6_7_Weak_Enhancer',
        '8_Insulator' , '9_Txn_Transition' , '10_Txn_Elongation' , '11_Weak_Txn' , '12_Repressed' , '13_Heterochrom/lo' ,
        '14_15_Repetitive/CNV']




chrom = ['chr' + str(x) for x in range(1 , 23)] + ['chrX']      



#######Genes Classify

expressed_genes = gene_exp_df[gene_exp_df['mean_FPKM'] >= 2]

low_fpkm = np.percentile(expressed_genes['mean_FPKM'] , 33.3)
high_fpkm = np.percentile(expressed_genes['mean_FPKM'] , 66.67)

noexp_gene = gene_exp_df[gene_exp_df['mean_FPKM'] < 2]
low_gene = gene_exp_df[(gene_exp_df['mean_FPKM'] >= 2) & (gene_exp_df['mean_FPKM'] < low_fpkm)]
middle_gene = gene_exp_df[(gene_exp_df['mean_FPKM'] >= low_fpkm) & (gene_exp_df['mean_FPKM'] < high_fpkm)]
high_gene = gene_exp_df[(gene_exp_df['mean_FPKM'] >= high_fpkm)]


all_gene_name = list(gene_exp_df['Gene_Name'])
noexp_gene_name = list(noexp_gene['Gene_Name'])
low_gene_name = list(low_gene['Gene_Name'])
middle_gene_name = list(middle_gene['Gene_Name'])
high_gene_name = list(high_gene['Gene_Name'])


Genes_all = {'all' : gene_exp_df , 'no_expressed' : noexp_gene , 'low_expressed' : low_gene , 'middle_expressed' : middle_gene , 'high_expressed' : high_gene}
Genes_name = {'all' : all_gene_name , 'no_expressed' : noexp_gene_name , 'low_expressed' : low_gene_name , 'middle_expressed' : middle_gene_name , 'high_expressed' : high_gene_name}




##############loops number boxplot
    
loops_num = {}    
for k in keys:
    print (k)
    loops_num[k] = []
    
    
    for g in chrom:
        tmp_peaks = classify[k][classify[k]['chr'] == g]
        tmp_loops = loops[loops['chr1'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            mask1 = (start <= tmp_loops['end1']) & (end >= tmp_loops['start1'])
            mask2 = (start <= tmp_loops['end2']) & (end >= tmp_loops['start2'])
            overlap1 = tmp_loops[mask1]
            overlap2 = tmp_loops[mask2]
            n = len(overlap1) + len(overlap2)
            loops_num[k].append((n))
            
            


df = pd.DataFrame(dict([(key, pd.Series(values)) for key, values in loops_num.items()]))

# 将 DataFrame 转换为长格式
df_melted = df.melt(var_name="Keys", value_name="Values")

# 绘制箱线图
plt.figure(figsize=(12, 6))
sns.boxplot(x="Keys", y="Values", data=df_melted, palette="Set3" , showfliers=False)

# 添加标题和标签
plt.title("Different classify peaks related loop numbers", fontsize=14)
plt.ylabel("Loop numbers", fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()

# 显示图表
plt.show()



#####################loops Anno related respective level genes
loops_num = {}    
for k in keys:
    print (k)
    loops_num[k] = []
    
    
    for g in chrom:
        tmp_peaks = classify[k][classify[k]['chr'] == g]
        tmp_loops = loops[loops['chr1'] == g]
        tmp_gene = gene_exp_df[gene_exp_df['Chr'] == g]
        for i in tmp_peaks.index:
            start = tmp_peaks.loc[i]['start']
            end = tmp_peaks.loc[i]['end']
            mask1 = (start <= tmp_loops['end1']) & (end >= tmp_loops['start1'])
            mask2 = (start <= tmp_loops['end2']) & (end >= tmp_loops['start2'])
            overlap1 = tmp_loops[mask1]
            overlap2 = tmp_loops[mask2]
            r_gene = [] ; l_gene = []
            if len(overlap1) != 0:
                for j in overlap1.index:
                    l_s = overlap1.loc[j]['start2']
                    l_e = overlap1.loc[j]['end2']
                    mask3 = (l_s <= tmp_gene['pend']) & (l_e >= tmp_gene['pstart'])
                    overlap3 = tmp_gene[mask3]
                    r_gene.extend(list(overlap3['Gene_Name']))
            if len(overlap2) != 0:  
                for j in overlap2.index:
                    l_s = overlap2.loc[j]['start1']
                    l_e = overlap2.loc[j]['end1']
                    mask4 = (l_s <= tmp_gene['pend']) & (l_e >= tmp_gene['pstart'])
                    overlap4 = tmp_gene[mask4]
                    l_gene.extend(list(overlap4['Gene_Name']))
            genes = list(set(l_gene + r_gene))
            loops_num[k].append(genes)
                
                    


###############boxplots        
            
# 原始字典数据
genes_num = {}        
        
for i in ['all', 'no_expressed', 'low_expressed', 'middle_expressed', 'high_expressed']:
    g_name = Genes_name[i]
    genes_num[i] = {}
    for j in loops_num:
        genes_num[i][j] = []
        tmp_num = loops_num[j]
        for k in tmp_num:
            n = 0
            for l in k:
                if l in g_name:
                    n += 1
            if n != 0:
                genes_num[i][j].append(n)




# 转换为适合绘图的长格式 DataFrame
df = pd.DataFrame.from_dict(genes_num, orient="index").stack().reset_index()
df.columns = ["Group", "Category", "Value"]
df = df.explode("Value")  # 将列表展开为单独的行
df["Value"] = df["Value"].astype(float)  # 确保值为数值类型

# 绘制箱线图
plt.figure(figsize=(8, 6))
sns.boxplot(x="Category", y="Value", hue="Group", data=df, palette="Set2" , showfliers=False)
plt.title("Boxplot for Dictionary Data", fontsize=16)
plt.xticks(rotation=45)
plt.xlabel("Category", fontsize=12)
plt.ylabel("Value", fontsize=12)
plt.legend(title="Group")
plt.tight_layout()
plt.show()


plt.savefig('')



                    
            
############sns barplots
matrix = []
for level in ['no_expressed' , 'low_expressed' , 'middle_expressed' , 'high_expressed']:
    matrix.append([])
    o_num = genes_num[level]
    for anno in keys:
        sum_genes = sum(o_num[anno])
        matrix[-1].append(sum_genes)

        
matrix = pd.DataFrame(matrix)
matrix.columns = keys        
matrix.index = ['no_expressed' , 'low_expressed' , 'middle_expressed' , 'high_expressed']      


for c in ['no_expressed' , 'low_expressed' , 'middle_expressed' , 'high_expressed']:
    matrix.loc[c] /= len(Genes_name[c])
    
for c in keys:
    matrix[c] /= len(classify[c])
    
    



# 转换为适合绘图的长格式
df_long = matrix.reset_index().melt(id_vars="index", var_name="Category", value_name="Value")

# 绘制柱状图
plt.figure(figsize=(10, 6))
sns.barplot(data=df_long, x="Category", y="Value", hue="index", palette="viridis")
plt.title("Barplot Categorized by Expression Levels", fontsize=16)
plt.xlabel("HMM Chromatin State Categories", fontsize=12)
plt.ylabel("Average loops number", fontsize=12)
plt.legend(title="Expression Levels")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()


        
plt.savefig('H:\\work\\niulongjian\\HiRPC_processed_data\\plots\\plots_New_number_bylxx_xjs\\K562_loops_Anno\\Anno_peaks_related_average_loops_number_barplots.pdf')
    
    
    









