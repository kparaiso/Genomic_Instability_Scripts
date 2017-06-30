# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 23:32:58 2016

@author: Wesley
"""

import pandas as pd
from os import walk, listdir


if __name__ == '__main__':
    
    wd = "D:/Graeber Lab/Genomic_instability/Winter_break/LOH_8p/Jan_10/"   
    gsea_dir = [x[0] for x in walk(wd)]
    del gsea_dir[0::2]
    
    # grab the gene sets we are interested in from c2
    GeneSetList = []
    c = pd.read_excel("D:/Graeber Lab/Genomic_instability/Winter_break/Uveal/across_tumors_ICNA_v_8ploss.xlsx")
    GeneSetList += c["GENE_SET"].tolist()
    GeneSetList = [x.rstrip() for x in GeneSetList]
    
    # get the tumor group names and their corresponding measures (i.e. CoreTumor1/2, ICNA/Bkpt rvalues)
    tumor_group = ['-'.join((i.split("/")[-1]).split(".")[0].split('_')[2:5]) for i in gsea_dir]
#    tumor_group = list(set(tumor_group))
    df = pd.DataFrame(index=GeneSetList, columns=tumor_group)
    
    # put the corresponding gene set rank in each tumor group into a dataframe
    for i in gsea_dir:
        tumor = '-'.join((i.split("/")[-1]).split(".")[0].split('_')[2:5])
        files = listdir(i)
        na_pos = "" 
        na_neg = ""
        for f in files:
            if ("gsea_report_for_YES" in f) and (".xls" in f):
                na_pos = f
            if ("gsea_report_for_NO" in f) and (".xls" in f):
                na_neg = f
        df1 = pd.read_table(i+'/'+na_pos, index_col=0)
        df2 = pd.read_table(i+'/'+na_neg, index_col=0)
        df_ = pd.concat([df1, df2])
        for g in GeneSetList:
           try:
               if (pd.isnull(df.loc[g, tumor])):
                   df.set_value(g, tumor, df_.loc[g, "NES"])
           except:
               continue
#        print(df)
          
    df = df.astype(float)
#    df = df.reindex_axis(df.mean(axis=1).sort_values(ascending=False).index, axis=0, )
    df.to_excel('genomic_instability_BRCA_SKCM_new.xlsx')
    
    # draw the heatmap
    import seaborn as sns
    import matplotlib.pyplot as plt
    fig = plt.figure()
    fig.set_size_inches(100, 100)
    sns.heatmap(df, annot=True, fmt=".2f",annot_kws={"size":30})
    plt.xticks(rotation=0, fontsize=40)
    plt.yticks(rotation=0, fontsize=40)
#    #plt.rc('font',size=18)
    plt.tight_layout()
    plt.savefig("genomic_instability_BRCA_SKCM_new.png", dpi=200, bbox_inches='tight')
