
'''
get only one gene from the gtf file
'''
#import gtfparse
import gffutils
import pandas as pd
import gzip
import json
import pybedtools
from pybedtools import BedTool
import numpy as np
import os


# In[296]:


import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import seaborn as sns
from matplotlib import rcParams 
rcParams['pdf.fonttype'] = 42 # True font
rcParams['font.size'] =  8  
rcParams['grid.linewidth'] =  0.5 
rcParams['lines.color'] = 'b' 
rcParams['lines.linewidth'] = 1 
rcParams['lines.markersize'] = 3
rcParams['lines.markeredgewidth'] = 0 # set Marker with no edgelines
rcParams['axes.linewidth'] = 0.5
rcParams['axes.titlesize'] = 12
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['legend.title_fontsize'] = 8



outfolder = '/scratch/2024-11-12/bio-shenw/Ljniu/K562/plots/HiRPC_classify/plots/lxx_new/K562/'


'''
#读入文件格式: chr1, s1, e1, chr2, s2, e2, count, pval
第一步筛选文件中chrom1和chrom2列与anchor_region["chrom"]都相同的数据行
'''

filename = "/scratch/2024-11-25/bio-shenw/Ljniu/K562/K562_0.1_FA/loops/ChIA-pipline_8K+/one_anchors_binding_loops/K562_merge5_0.1FA.hg38_peaks_one_anchors_binding_loops.bedpe"
if ".gz" in filename:
    fileID = gzip.open(filename, "rt")
else:
    fileID = open(filename, "r")

contact_df = pd.read_csv(fileID, sep="\t", names=["chr1", "s1", "e1",
                                                               "chr2","s2", "e2", "count","-log10qval"])
## set contact ID
contact_df["ID"] = contact_df.index.values
contact_df 


# genes expression and the promoter regions
gene_exp_file = "/scratch/2024-11-25/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/K562_union_gene_classify.txt"
gene_exp_df = pd.read_csv(gene_exp_file, header=0, sep="\t", index_col=None)
gene_exp_df["gID"] = gene_exp_df.index.values
gene_exp_df



Genenamedict = gene_exp_df["Gene_Name"].to_dict()



gene_exp_df["Level"].value_counts()




gene_ratio = gene_exp_df["Level"].value_counts()/gene_exp_df["Level"].value_counts().sum()
gene_ratio 




# Enhancer and super enhancers
TEfile = "/scratch/2024-11-25/bio-shenw/literature/K562/Enhancer/SEdb_K562_01_0039_TE_hg38.bed"
TE_df = pd.read_csv(TEfile, sep="\t", header=0, index_col=None )
#print(TE_df.head())
## Super enhancer elements
SEfile = "/scratch/2024-11-25/bio-shenw/literature/K562/Enhancer/SEdb_K562_01_0039_SE_ele_hg38.bed"
SE_df = pd.read_csv(SEfile, sep="\t", header=0, index_col=None )
SE_df["ele_id"] = SE_df.index.values + 1
SE_df.loc[:, "ele_id"]  = SE_df["ele_id"].apply(lambda x: f"SE_ELE_01_0039{x:05d}")
#print(SE_df.head())
# Merge two enhancer dataframe
adf = TE_df.iloc[:, [0,1,2,3]]
adf.columns = ["chrom", "start", "end", "name"]
bdf = SE_df.iloc[:, [0,1,2,5]]
bdf.columns = ["chrom", "start", "end", "name"]

Enhancer_df = pd.concat([adf, bdf], ignore_index=True).sort_values(by=["chrom", "start"], ignore_index=True)
print(f"TE:{len(adf)}; SE_ele:{len(bdf)}")
del(adf, bdf)

Enhancer_df


# In[305]:


# region overlaps
def overlap_length(row):
    '''
    Calculate the overlap length
    '''
    if row["chrom"] != row["pchrom"]:
        result = 0 # Interchromosomal
    else:
        result = min(row['end'], row["pend"]) - max(row["start"], row["pstart"])
    return abs(result)

def Overlap_records(c_df, p_df):
    '''
    Find overlap of contact_df in both anchors with promoter and enhancers
    '''
    # overlaps
    # 将pandas dataframe转为 pybedtools BedTool 对象
    p_bed = pybedtools.BedTool.from_dataframe( p_df )
    c_bed = pybedtools.BedTool.from_dataframe( c_df ) # contact interaction
    # 找到 contact_bed, one of the anchor, overlap with elements promoter_bed overlap 的记录
    overlap_bed = c_bed.intersect( p_bed, wa=True, wb=True)
    overlap_df = overlap_bed.to_dataframe()
    pybedtools.cleanup()
    colnames = ["chrom", "start", "end", "ID" , "pchrom", "pstart", "pend", "gID"]
    overlap_df.columns = colnames
    overlap_df["intlen"] = overlap_df.apply(overlap_length, axis=1)
    overlap_df = overlap_df.sort_values(by=["ID", "intlen"], ascending=False).drop_duplicates(subset=["ID"], keep="first")
    overlap_df = overlap_df.loc[:, colnames]
    return(overlap_df)

## promoter and contact anchors
A1_cols = ["chr1", "s1", "e1", "ID"]
A2_cols = ["chr1", "s2", "e2", "ID"]
# promoter and contact anchor1, anchor2
A1_promoter_df = Overlap_records(contact_df.loc[:, A1_cols], gene_exp_df.iloc[:, [7,8,9,4]] )
A2_promoter_df =  Overlap_records(contact_df.loc[:, A2_cols], gene_exp_df.iloc[:, [7,8,9,4]] )
A2_promoter_df.head()


# In[306]:


# ### Assign anchor to overlap promoter
contact_df["A1p"] = ""
IDs = A1_promoter_df["ID"].values
gIDs = A1_promoter_df["gID"].values
contact_df.loc[IDs, "A1p"] = gIDs

contact_df["A2p"] = ""
IDs = A2_promoter_df["ID"].values
gIDs = A2_promoter_df["gID"].values
contact_df.loc[IDs, "A2p"] = gIDs
contact_df


# In[307]:


## Enhancer and contact anchors
A1_cols = ["chr1", "s1", "e1", "ID"]
A2_cols = ["chr1", "s2", "e2", "ID"]
# Enhancer and contact anchor1, anchor2
A1_enhancer_df = Overlap_records(contact_df.loc[:, A1_cols], Enhancer_df )
A2_enhancer_df =  Overlap_records(contact_df.loc[:, A2_cols], Enhancer_df )
A2_enhancer_df.head()


# In[308]:


# ### Assign anchor to overlap enhancer
contact_df["A1e"] = ""
IDs = A1_enhancer_df["ID"].values
gIDs = A1_enhancer_df["gID"].values
contact_df.loc[IDs, "A1e"] = gIDs

contact_df["A2e"] = ""
IDs = A2_enhancer_df["ID"].values
gIDs = A2_enhancer_df["gID"].values
contact_df.loc[IDs, "A2e"] = gIDs
contact_df


# In[309]:


P = (contact_df["A1p"] != "") & (contact_df["A2e"] != "")
contact_df.loc[P, :]


# In[310]:


# Contact type: PP, PE, EE, other
## Promoter-Promoter PP
contact_df["PP"] = 0
P1 = (contact_df["A1p"] != "")  & (contact_df["A2p"] != "")
contact_df.loc[P1, "PP"] = 1

## Promoter-Enhancer PE
contact_df["PE"] = 0
P2 = ( (contact_df["A1p"] != "")  & (contact_df["A2e"] != "") ) | ( (contact_df["A2p"] != "")  & (contact_df["A1e"] != "") )
contact_df.loc[P2, "PE"] = 1

## Enhancer-Enhancer EE
contact_df["EE"] = 0
P3 = (contact_df["A1e"] != "")  & (contact_df["A2e"] != "")
contact_df.loc[P3, "EE"] = 1

contact_df["other"] = 1
contact_df.loc[(P1|P2|P3), "other"] = 0
contact_df


# In[311]:


contact_df.loc[contact_df["other"]==0, :]


# In[312]:


# contact distance
contact_df["dist"] = (contact_df["s1"] + contact_df["e1"])/2  - (contact_df["s2"] + contact_df["e2"])/2 
contact_df["dist"] = abs( contact_df["dist"].astype("int").values )


# In[313]:


contact_df


# In[314]:


# Plot venn
from venn import venn
from matplotlib_venn import venn3, venn3_circles
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    

            
c1_gID = list(gene_exp_df[gene_exp_df['Level'] == 'c1']['Gene_Name'])
c2_gID = list(gene_exp_df[gene_exp_df['Level'] == 'c2']['Gene_Name'])
c3_gID = list(gene_exp_df[gene_exp_df['Level'] == 'c3']['Gene_Name'])
c4_gID = list(gene_exp_df[gene_exp_df['Level'] == 'c4']['Gene_Name'])



classify = { 'c1' : contact_df[contact_df['A1p'].isin(c1_gID) | contact_df['A2p'].isin(c1_gID)] ,
             'c2' : contact_df[contact_df['A1p'].isin(c2_gID) | contact_df['A2p'].isin(c2_gID)] ,
             'c3' : contact_df[contact_df['A1p'].isin(c3_gID) | contact_df['A2p'].isin(c3_gID)] ,
             'c4' : contact_df[contact_df['A1p'].isin(c4_gID) | contact_df['A2p'].isin(c4_gID)] , 
             'total_concat' : contact_df , 
             'total_classify' : pd.concat([contact_df[contact_df['A1p'].isin(c1_gID) | contact_df['A2p'].isin(c1_gID)] , 
                                           contact_df[contact_df['A1p'].isin(c1_gID) | contact_df['A2p'].isin(c1_gID)] , 
                                           contact_df[contact_df['A1p'].isin(c3_gID) | contact_df['A2p'].isin(c3_gID)] , 
                                           contact_df[contact_df['A1p'].isin(c4_gID) | contact_df['A2p'].isin(c4_gID)]]) }
            


for c in ['c1' , 'c2' , 'c3' , 'c4' , 'total_concat'  , 'total_classify' ]:
    tmp_contact_df = classify[c]
    print(f"Significant contacts:{len(tmp_contact_df) }")
    
    # Only one anchor is Promoter the other anchor unkonwn
    P1 = (tmp_contact_df["A1p"] != "")&(tmp_contact_df["A2p"] == "")&(tmp_contact_df["A2e"] == "")  
    P2 =  (tmp_contact_df["A2p"] != "")&(tmp_contact_df["A1p"] == "")&(tmp_contact_df["A1e"] == "")   
    print(f"Promoter-other contacts: {sum(P1|P2)} ({100*sum(P1|P2)/len(tmp_contact_df):.1f}%)")
    
    # Only one anchor is Enhancer the other anchor unkonwn
    P1 = (tmp_contact_df["A1e"] != "")&(tmp_contact_df["A2e"] == "")&(tmp_contact_df["A2p"] == "")  
    P2 =  (tmp_contact_df["A2e"] != "")&(tmp_contact_df["A1e"] == "")&(tmp_contact_df["A1p"] == "")   
    print(f"Enhancer-other contacts:{sum(P1|P2)} ({100*sum(P1|P2)/len(tmp_contact_df):.1f}%)")
    
    P = tmp_contact_df["other"]==0
    print(f"Promoters/Enhancers in both anchors contacts: {sum(P)} {100*sum(P)/len(tmp_contact_df):.1f}%")
    
    # None Promoter and Enhancers contacts
    Po = (tmp_contact_df["A1p"] == "")&(tmp_contact_df["A2p"] == "")&(tmp_contact_df["A1e"] == "")&(tmp_contact_df["A2e"] == "")  
    print(f"None Promoter_enhancer contacts:{sum(Po)} {100*sum(Po)/len(tmp_contact_df):.1f}%")


    df = tmp_contact_df.loc[P, ["PP", "PE", "EE"]]
    # pe, pp, ee sets ids
    pe_set = set(df[df["PE"] == 1].index)
    pp_set = set(df[df["PP"] == 1].index)
    ee_set = set(df[df["EE"] == 1].index)
    ## venn plot
    fig = plt.figure(figsize = (10, 10))
    venn3([pe_set , pp_set , ee_set] , ('PE' , 'PP' , 'EE'))
    ax = plt.gca()
    for text in ax.texts:
        text.set_fontsize(20)  # 调整字体大小，例如 14
    run_Plot(fig, os.path.join( outfolder , c + '_Venn3.pdf'))


    
    
#####write loops information     
contact_df.replace('', 'no_values', inplace=True)    
contact_df.to_csv('/scratch/2024-10-14/bio-shenw/Ljniu/K562/plots/HiRPC_classify/classify_K562/loops_classify/K562_loops_classify.txt' , header = True , index = None , sep = '\t')    
    


# P-P subtype
def ExpressionCombination(row):
    '''
    P-P's expression level 
    '''
    gene1, gene2 = row["A1p"], row["A2p"]
    l1, l2 = genelevel_dict[gene1], genelevel_dict[gene2]
    l12 = sorted([l1,l2])
    ltype = f"{l12[0]}-{l12[1]}"
    return(ltype)

genelevel_dict = gene_exp_df.set_index('Gene_Name')['Level'].to_dict()
gene_ratio = (gene_exp_df["Level"].value_counts()/len(gene_exp_df)).to_dict()

PP_df = contact_df.loc[contact_df["PP"]==1, ["ID","count", "-log10qval", "dist", "A1p", "A2p"] ]
PP_df["PP_type"] = PP_df.apply(ExpressionCombination, axis=1)
PP_df


# In[317]:


PP_enrich =  PP_df["PP_type"].value_counts()
PP_enrich = pd.DataFrame(PP_enrich).reset_index()
PP_enrich.columns = ["PP_Type", "obs_contacts"]
PP_enrich["obs_ratio"] = PP_enrich["obs_contacts"].values/PP_enrich["obs_contacts"].sum()

def Exp_ratio(pp_type):
    [l1, l2] = pp_type.split("-")
    if l1 == l2:
        eratio = gene_ratio[l1]*gene_ratio[l2]
    if l1 != l2:
        eratio = gene_ratio[l1]*gene_ratio[l2] + gene_ratio[l2]*gene_ratio[l1]
    return(eratio)
PP_enrich["exp_ratio"] = PP_enrich["PP_Type"].apply(lambda x: Exp_ratio(x) )
PP_enrich["log2(Obs/Exp)"] = np.log2( PP_enrich["obs_ratio"].values/ PP_enrich["exp_ratio"].values )
PP_enrich


# In[318]:


# PP type order
orderlist = ["c1-c1","c2-c2", "c1-c3", "c1-c2", "c2-c3", "c2-c4" , "c1-c4", "c3-c3", "c3-c4", "c4-c4"]
PP_enrich["PP_Type"] = pd.Categorical(PP_enrich["PP_Type"], categories=orderlist, ordered=True)
PP_enrich = PP_enrich.sort_values("PP_Type")
PP_enrich


# In[319]:


# 绘制P-P contacts distribution plot
fig = plt.figure(figsize = (10, 10))
plt.pie(PP_enrich["obs_contacts"].values, labels=PP_enrich["PP_Type"].values, autopct='%1.1f%%')
plt.title('Promoter-Promoter interactions percentages')
run_Plot(fig, os.path.join( outfolder , 'Promoter_Promoter_interactions_percentages.pdf'))


# In[385]:


# 绘制P-P contacts distribution plot
fig = plt.figure(figsize = (10, 10))
plt.pie(PP_enrich["exp_ratio"].values, labels=PP_enrich["PP_Type"].values, autopct='%1.1f%%')
plt.title('Promoter-Promoter interactions (Expected)')
run_Plot(fig,  os.path.join( outfolder , 'Promoter-Promoter interactions (Expected).pdf'))



# In[321]:


# P-P enrichment bar plot
fig = plt.figure(figsize = (10, 10))
plt.bar(PP_enrich["PP_Type"].values, PP_enrich["log2(Obs/Exp)"].values)

# 添加标题和标签
plt.xlabel('Promoter-Promote interactions')
plt.ylabel('log2(Obs/Exp)')
plt.title('Enrichment of Promoter-Promoter interactions')
plt.xticks(PP_enrich["PP_Type"].values, fontsize=8, rotation=45, color='k')

run_Plot(fig,  os.path.join( outfolder , 'Enrichment of Promoter-Promoter interactions.pdf'))



# In[324]:


plotdata = PP_df.reset_index(drop=True).copy()
orderlist = ["c1-c1","c2-c2", "c1-c3", "c1-c2", "c2-c3", "c2-c4" , "c1-c4", "c3-c3", "c3-c4", "c4-c4"]


plotdata["PP_type"] = pd.Categorical(plotdata["PP_type"], categories=orderlist, ordered=True)
plotdata = plotdata.sort_values("PP_type")
plotdata.loc[:,"log10dist"] = np.log10( plotdata["dist"].values)
plotdata.loc[:,"log10(pval+1)"] = np.log10( plotdata["-log10qval"].values + 1)



count = plotdata['PP_type'].value_counts()
count = pd.DataFrame(count)


ratio_dict = []
for i in ["c1-c1","c2-c2", "c1-c3", "c1-c2", "c2-c3", "c2-c4" , "c1-c4", "c3-c3", "c3-c4", "c4-c4"]:
    loop_num = count.loc[i]['PP_type']
    tmp = plotdata[plotdata['PP_type'] == i]
    gene_num = len(set(tmp['A1p']) | set(tmp['A2p']))
    print (i , loop_num / gene_num)
    ratio_dict.append((i , loop_num / gene_num))

ratio_dict = pd.DataFrame(ratio_dict , columns=['PP_type' , 'ratio'])

## subplot of contact distance and contact strength
fig, axs = plt.subplots(2, 1,figsize=(8,4))
axs = axs.flatten() 
### contact distance
sns.boxplot(x="PP_type", y="log10dist", data=plotdata, ax=axs[0])

### contact strength logq
sns.barplot(x="PP_type", y="ratio", data=ratio_dict, ax=axs[1])

# ### contact strength logq
# sns.barplot(x="PP_type", y="count", data=plotdata, ax=axs[2])

run_Plot(fig,  os.path.join( outfolder , 'Contact distance and average contact number.pdf'))








# Promoter interactions
PIdict = {}
for n, row in PP_df.iterrows():
    p1, p2 = row["A1p"], row["A2p"]
    if p1 not in PIdict:
        PIdict[p1] = []
    PIdict[p1].append(p2)
    if p2 not in PIdict:
        PIdict[p2] = []
    PIdict[p2].append(p1)


# In[329]:


Promoters = []
Levelist = []
Interact_Pnum = []
High_IP = []
Mid_IP = []
Low_IP = []
No_IP = []

for pkey, plist in PIdict.items():
    Promoters.append(pkey)
    Levelist.append(  genelevel_dict[pkey] )
    plist = list( set(plist) )
    Interact_Pnum.append( len(plist) )
    IP_levels = [ genelevel_dict[p]  for p in plist ]
    High_IP.append( IP_levels.count("c1") )
    Mid_IP.append( IP_levels.count("c2") )
    Low_IP.append( IP_levels.count("c3") )
    No_IP.append( IP_levels.count("c4") )

Promoter_Int_summary = pd.DataFrame({"Promoter":Promoters,
                                     "Level": Levelist,
                                     "P_num": Interact_Pnum,
                                     "C1_num":High_IP,
                                    "C2_num":Mid_IP,
                                    "C3_num":Low_IP,
                                    "C4_num":No_IP})

orderlist = ["c1","c2", "c3", "c4"]
Promoter_Int_summary["Level"] = pd.Categorical(Promoter_Int_summary["Level"], categories=orderlist, ordered=True)
Promoter_Int_summary= Promoter_Int_summary.sort_values("Level")
Promoter_Int_summary



## subplot of contact distance and contact strength
fig = plt.figure(figsize=(4,4))
### contact distance
sns.boxplot(x="Level", y="P_num", data=Promoter_Int_summary)
plt.xlabel("Genes of different classify")
plt.ylabel("Contact promoters")
run_Plot(fig,  os.path.join( outfolder , 'Genes of different classify VS Contact promoters.pdf'))



# In[332]:


# PEI
## enhancer

### ELE and SE ID dict
ELE_SE_dict = SE_df.loc[:, ["ele_id", "se_id"]].set_index("ele_id", drop=True)["se_id"].to_dict()

PE_df = contact_df.loc[contact_df["PE"]==1, :].copy()

# Promoter-Enhancer interactions
PEdict = {}
EPdict = {}
PEpair = {}
for n, row in PE_df.iterrows():
    p1, p2 = row["A1p"], row["A2p"]
    e1, e2 = row["A1e"], row["A2e"]
    # SE_ELE id to SE id
    if "SE" in e1:
        e1 = ELE_SE_dict[e1]
    if "SE" in e2:
        e2 = ELE_SE_dict[e2]
        
    if p1 != "" and e2 !="":
        if p1 not in PEdict:
            PEdict[p1] = []
            PEpair[(p1, e2)] = 1
        PEdict[p1].append(e2)
        if e2 not in EPdict:
            EPdict[e2] = []
        EPdict[e2].append(p1)
        
    if p2 != "" and e1 !="":
        if p2 not in PEdict:
            PEdict[p2] = []
            PEpair[(p2, e1)] = 1
        PEdict[p2].append(e1)
        if e1 not in EPdict:
            EPdict[e1] = []
        EPdict[e1].append(p2)
PEdict


# In[333]:


for key, Elist in PEdict.items():
    if key in PIdict:
        Elist


# In[334]:


PEpairsummary = {}
for p_e in PEpair.keys():
    (pkey, ekey) = p_e
    level = genelevel_dict[pkey]
    if "SE" in ekey:
        Etype = "SE"
    else:
        Etype = "TE"
    PEtype = f"{level}-{Etype}"
    if PEtype not in PEpairsummary:
        PEpairsummary[PEtype] = 0
    PEpairsummary[PEtype] += 1
    
PEpairsummary = pd.DataFrame(PEpairsummary, index=["count"]).T.reset_index(drop=False)
PEpairsummary.columns = ["PEtype", "Intcount"]
PEpairsummary[["Plevel","Etype"]] = PEpairsummary["PEtype"].str.split("-", expand=True)
PEpairsummary = PEpairsummary.sort_values(by=["Etype", "Plevel"]).reset_index(drop=True)

'''
gene ratio
c1 0.2544186145436757      
c2 0.08165131760852458
c3 0.2249824531552417
c4 0.438947614692558 
'''

P = PEpairsummary["Etype"] == "SE"
PEpairsummary.loc[P, "Obs_ratio"] = PEpairsummary.loc[P, "Intcount"].values / PEpairsummary.loc[P, "Intcount"].sum()
PEpairsummary.loc[P, "Exp_ratio"] = [0.2544186, 0.0816513, 0.2249825,  0.4389476]

P = PEpairsummary["Etype"] == "TE"
PEpairsummary.loc[P, "Obs_ratio"] = PEpairsummary.loc[P, "Intcount"].values / PEpairsummary.loc[P, "Intcount"].sum()
PEpairsummary.loc[P, "Exp_ratio"] = [0.2544186, 0.0816513, 0.2249825,  0.4389476]

PEpairsummary["log2(Obs/Exp)"] = np.log2( PEpairsummary["Obs_ratio"].values/ PEpairsummary["Exp_ratio"].values )
PEpairsummary


# In[335]:


# P-E enrichment bar plot
fig = plt.figure(figsize=(6,4))
plt.bar(PEpairsummary["PEtype"].values, PEpairsummary["log2(Obs/Exp)"].values)

# 添加标题和标签
plt.xlabel('Promoter-Enhancer interaction')
plt.ylabel('log2(Obs/Exp)')
plt.title('Enrichment of Promoter-Enhancer interactions')
plt.xticks(PEpairsummary["PEtype"].values, fontsize=8, rotation=45, color='k')
plt.ylim([-5.5, 5.5])
# 显示图表
run_Plot(fig,  os.path.join( outfolder , 'Enrichment of Promoter-Enhancer interactions.pdf'))



# In[336]:


# P-E summary
Promoters = []
Levelist = []
Interact_Enum = []
SEs = []
TEs = []
for pkey, Elist in PEdict.items():
    Promoters.append(pkey)
    Levelist.append(  genelevel_dict[pkey] )
    Elist = list( set(Elist) )
    Interact_Enum.append( len(Elist) )
    SEs.append(  len( list(filter(lambda x: "SE" in x, Elist) ) )  )
    TEs.append(  len( list(filter(lambda x: "TE" in x, Elist) ) )  )
    
PE_Int_summary = pd.DataFrame({"Promoter":Promoters,
                                     "Level": Levelist,
                                     "E_num": Interact_Enum,
                                     "SE_num":SEs,
                                    "TE_num":TEs})

orderlist = ["c1","c2", "c3", "c4"]
PE_Int_summary["Level"] = pd.Categorical(PE_Int_summary["Level"], categories=orderlist, ordered=True)
PE_Int_summary= PE_Int_summary.sort_values("Level")
PE_Int_summary


# In[356]:


# exfile = "Promoter-Enhancer_interaction_summary.csv"
# PE_Int_summary.to_csv(exfile, header=True, sep="\t")


# In[338]:


##  Promoters-Enahncer
fig = fig = plt.figure(figsize=(4,4))
### contact distance
sns.boxplot(x="Level", y="E_num", data=PE_Int_summary, showfliers=False)
plt.xlabel("Gene classify")
plt.ylabel("The number of enhancers")
#plt.ylim([0,30])
run_Plot(fig,  os.path.join( outfolder , 'Gene Expression VS The number of enhancers.pdf'))




# Promoter Ehancer 互作的堆叠直方图
plotdf = PE_Int_summary.copy()
plotdf.loc[plotdf["SE_num"]>=4,"P_num"] = 4
# 将E_num列进行分桶，划分为0-1, 1-2, ..., 5-6等6个区间
plotdf['E_num_bin'] = pd.cut(plotdf['SE_num'], bins=np.arange(0, 5, 1) )
# 根据E_num_bin和Level的分组，计算每个组的ID数
grouped = plotdf.groupby(['E_num_bin', 'Level']).size().reset_index(name='count')

# 计算每个组的百分比
grouped['percent'] = 0
for l in ["c1","c2", "c3", "c4"]:
    P = grouped["Level"]==l
    lsum = grouped.loc[P, "count"].sum()
    grouped.loc[P, 'percent'] = 100*grouped.loc[P, 'count'] / lsum 

# 绘制堆叠直方图
fig = plt.figure(figsize=(8, 4))
sns.barplot(x='E_num_bin', y='percent', hue='Level', data=grouped, ci=None)

# 添加y轴标签和标题
plt.ylabel('Promoters (%)')
plt.xlabel("Number of Super-enhancers") 
plt.title('Multiple Super-enhancer interaction')
plt.xticks(rotation=90)  # 因x轴标签可能过多，进行旋转显示
run_Plot(fig,  os.path.join( outfolder , 'Number of Super-enhancers VS Multiple Super-enhancer interaction.pdf'))



# In[340]:


# Promoter Ehancer 互作的堆叠直方图
plotdf = PE_Int_summary.copy()
plotdf.loc[plotdf["TE_num"]>=4,"P_num"] = 4
# 将E_num列进行分桶，划分为0-1, 1-2, ..., 5-6等6个区间
plotdf['E_num_bin'] = pd.cut(plotdf['TE_num'], bins=np.arange(0, 5, 1) )
# 根据E_num_bin和Level的分组，计算每个组的ID数
grouped = plotdf.groupby(['E_num_bin', 'Level']).size().reset_index(name='count')

# 计算每个组的百分比
grouped['percent'] = 0
for l in ["c1","c2", "c3", "c4"]:
    P = grouped["Level"]==l
    lsum = grouped.loc[P, "count"].sum()
    grouped.loc[P, 'percent'] = 100*grouped.loc[P, 'count'] / lsum 

# 绘制堆叠直方图
fig = plt.figure(figsize=(8, 4))
sns.barplot(x='E_num_bin', y='percent', hue='Level', data=grouped, ci=None)

# 添加y轴标签和标题
plt.ylabel('Promoters (%)')
plt.xlabel("Number of Typical-enhancers") 
plt.title('Multiple Typical-enhancer interaction')
plt.xticks(rotation=90)  # 因x轴标签可能过多，进行旋转显示
run_Plot(fig,  os.path.join( outfolder , 'Number of Typical-enhancers VS Multiple Typical-enhancer interaction.pdf'))



# In[341]:


# Promoter Ehancer 互作的堆叠直方图
plotdf = PE_Int_summary.copy()
plotdf.loc[plotdf["E_num"]>=6,"P_num"] = 6
# 将E_num列进行分桶，划分为0-1, 1-2, ..., 5-6等6个区间
plotdf['E_num_bin'] = pd.cut(plotdf['E_num'], bins=np.arange(0, 7, 1) )
# 根据E_num_bin和Level的分组，计算每个组的ID数
grouped = plotdf.groupby(['E_num_bin', 'Level']).size().reset_index(name='count')

# 计算每个组的百分比
grouped['percent'] = 0
for l in ["c1","c2", "c3", "c4"]:
    P = grouped["Level"]==l
    lsum = grouped.loc[P, "count"].sum()
    grouped.loc[P, 'percent'] = 100*grouped.loc[P, 'count'] / lsum 

# 绘制堆叠直方图
fig = plt.figure(figsize=(5, 4))
sns.barplot(x='E_num_bin', y='percent', hue='Level', data=grouped, ci=None)

# 添加y轴标签和标题
plt.ylabel('Promoters (%)')
plt.xlabel("Number of enhancers") 
plt.title('Multiple enhancer interaction')
plt.xticks(rotation=90)  # 因x轴标签可能过多，进行旋转显示
run_Plot(fig,  os.path.join( outfolder , 'Number of enhancers VS Multiple enhancer interaction.pdf'))



# In[342]:


## E-P summary
Enhancers = []
EType = []
Interact_Pnum = []
Pgenes = []
High_num, Mid_num, Low_num, No_num = [],[],[],[]

for Ekey, plist in EPdict.items():
    Enhancers.append(Ekey)
    if "SE" in Ekey:
        EType.append(  "SE" )
    else:
        EType.append(  "TE" )
    plist = list( set(plist) )
    genelist =  [ Genenamedict[pID] for pID in plist  ]
    Pgenes.append( ",".join(genelist) )
    Interact_Pnum.append( len(plist) )
    
    levellist = [ genelevel_dict[pID] for pID in plist ]
    High_num.append( levellist.count("c1") )
    Mid_num.append( levellist.count("c2") )
    Low_num.append( levellist.count("c3") )
    No_num.append( levellist.count("c4") )
    
EP_Int_summary = pd.DataFrame({"Enhancers":Enhancers,
                                "EType": EType,
                                 "P_num": Interact_Pnum,
                                "P_genes": Pgenes,
                                "c1_num": High_num,
                                "c2_num": Mid_num,
                                "c3_num": Low_num,
                                "c4_num": No_num })
EP_Int_summary




exfile = os.path.join(outfolder , 'K562_TE_SE_intPromoters.summary.csv')
EP_Int_summary.to_csv(exfile, header=True, sep="\t",index=False)





import seaborn as sns

## Enahncer - Promoters
fig = plt.figure(figsize=(4,4))
### contact distance
sns.boxplot(x="EType", y="P_num", data=EP_Int_summary, showfliers=False)
plt.xlabel("Mulitple Promoters interactions")
plt.xlabel("Enhancer types")
plt.ylabel("The number of promoters")
#plt.ylim([0,30])
run_Plot(fig,  os.path.join( outfolder , 'Enhancer types VS The number of promoters.pdf'))


## Enahncer - Promoters

fig = plt.figure(figsize=(4,4))
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]
ax = fig.add_axes(size_axes)
### contact distance
te_data = EP_Int_summary[EP_Int_summary['EType'] == 'TE']
se_data = EP_Int_summary[EP_Int_summary['EType'] == 'SE']
xlabel=['c1_SE' , 'c2_SE' , 'c3_SE' , 'c4_SE' , 'c1_TE' , 'c2_TE' , 'c3_TE' , 'c4_TE']
y = [se_data['c1_num'] , se_data['c2_num'] , se_data['c3_num'] , se_data['c4_num'] , te_data['c1_num'] , te_data['c2_num'] , te_data['c3_num'] , te_data['c4_num']]
 

ax.boxplot([y[0] , y[4]] , positions=[1 , 6] , showfliers=False, widths = 0.7 , 
        boxprops={'color': 'darkred','linewidth':1},
        medianprops={'color':'darkred','linewidth':1},
        capprops={'color':'darkred','linewidth':1},
        whiskerprops={'color':'darkred','linewidth':1})
ax.boxplot([y[1] , y[5]] , positions=[2 , 7] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'dodgerblue','linewidth':1},
        medianprops={'color':'dodgerblue','linewidth':1},
        capprops={'color':'dodgerblue','linewidth':1},
        whiskerprops={'color':'dodgerblue','linewidth':1})
ax.boxplot([y[2] , y[6]] , positions=[3 , 8] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'green','linewidth':1},
        medianprops={'color':'green','linewidth':1},
        capprops={'color':'green','linewidth':1},
        whiskerprops={'color':'green','linewidth':1})
ax.boxplot([y[3] , y[7]] , positions=[4 , 9] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'deeppink','linewidth':1},
        medianprops={'color':'deeppink','linewidth':1},
        capprops={'color':'deeppink','linewidth':1},
        whiskerprops={'color':'deeppink','linewidth':1})

ax.set_xticks([1 , 2 , 3 , 4 , 6 , 7 , 8 , 9])
ax.set_xticklabels(xlabel , fontsize = 5)
ax.set_xlim((0.5 , 9.5))    
    

ax.set_xlabel("Enhancer types")
ax.set_ylabel("The number of promoters")
#plt.ylim([0,30])
run_Plot(fig,  os.path.join( outfolder , 'Enhancer types VS The number of promoters.pdf'))


## Promoters - Enhancers


PE_classify = {}
for i in PEdict.keys():
    level = genelevel_dict[i]
    k_s = level + '_SE'
    k_t = level + '_TE'
    if k_s not in PE_classify.keys():
        PE_classify[k_s] = []
    else:
        pass
    if k_t not in PE_classify.keys():
        PE_classify[k_t] = []
    else:
        pass
    s = 0 ; t = 0
    for j in PEdict[i]:
        if 'SE' in j:
            s += 1
        elif 'TE' in j:
            t += 1
        else:
            print (j)
    if s != 0:
        PE_classify[k_s].append(s)
    if t != 0:
        PE_classify[k_t].append(t)




fig = plt.figure(figsize=(4,4))
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]
ax = fig.add_axes(size_axes)
### contact distance
# te_data = EP_Int_summary[EP_Int_summary['EType'] == 'TE']
# se_data = EP_Int_summary[EP_Int_summary['EType'] == 'SE']
xlabel=['c1_SE' , 'c2_SE' , 'c3_SE' , 'c4_SE' , 'c1_TE' , 'c2_TE' , 'c3_TE' , 'c4_TE']
# y = [se_data['c1_num'] , se_data['c2_num'] , se_data['c3_num'] , se_data['c4_num'] , te_data['c1_num'] , te_data['c2_num'] , te_data['c3_num'] , te_data['c4_num']]
y = [PE_classify['c1_SE'] , PE_classify['c2_SE'] , PE_classify['c3_SE'] , PE_classify['c4_SE'] , PE_classify['c1_TE'] , PE_classify['c2_TE'] , PE_classify['c3_TE'] , PE_classify['c4_TE']]

ax.boxplot([y[0] , y[4]] , positions=[1 , 6] , showfliers=False, widths = 0.7 , 
        boxprops={'color': 'darkred','linewidth':1},
        medianprops={'color':'darkred','linewidth':1},
        capprops={'color':'darkred','linewidth':1},
        whiskerprops={'color':'darkred','linewidth':1})
ax.boxplot([y[1] , y[5]] , positions=[2 , 7] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'dodgerblue','linewidth':1},
        medianprops={'color':'dodgerblue','linewidth':1},
        capprops={'color':'dodgerblue','linewidth':1},
        whiskerprops={'color':'dodgerblue','linewidth':1})
ax.boxplot([y[2] , y[6]] , positions=[3 , 8] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'green','linewidth':1},
        medianprops={'color':'green','linewidth':1},
        capprops={'color':'green','linewidth':1},
        whiskerprops={'color':'green','linewidth':1})
ax.boxplot([y[3] , y[7]] , positions=[4 , 9] , showfliers=False, widths = 0.7 ,
        boxprops={'color': 'deeppink','linewidth':1},
        medianprops={'color':'deeppink','linewidth':1},
        capprops={'color':'deeppink','linewidth':1},
        whiskerprops={'color':'deeppink','linewidth':1})

ax.set_xticks([1 , 2 , 3 , 4 , 6 , 7 , 8 , 9])
ax.set_xticklabels(xlabel , fontsize = 5)
ax.set_xlim((0.5 , 9.5))    
    

ax.set_xlabel("Enhancer types")
ax.set_ylabel("The number of Enhancers")
#plt.ylim([0,30])
run_Plot(fig,  os.path.join( outfolder , 'Enhancer types VS The number of enhancers.pdf'))


