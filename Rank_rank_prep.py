 # -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 10:20:33 2016

@author: Wesley
"""

import pandas as pd
from os.path import join
from os import walk, listdir
from collections import defaultdict


def generate_rrho_file(input_file, wd):
    rf = pd.read_table(input_file, index_col=0)
    columns_ = rf.columns.tolist()
    for i in range(len(columns_)):
        for j in range(i+1, len(columns_)):
            rrho_file = pd.DataFrame({"Unigene":"N/A", "Gene":rf.index.tolist(), columns_[j]:rf[columns_[j]], columns_[i]:rf[columns_[i]]})
            j_rank = rrho_file[columns_[j]].rank(ascending=False).astype(int)
            i_rank = rrho_file[columns_[i]].rank(ascending=False).astype(int)
            rrho_file.insert(2, columns_[j]+"_rank", j_rank)
            rrho_file.insert(3, columns_[i]+"_rank", i_rank)
            rrho_file = rrho_file[["Unigene", "Gene", columns_[j]+"_rank", columns_[i]+"_rank", columns_[j], columns_[i]]]
            output_fo = join(wd, "RRHO_"+columns_[i]+"-"+columns_[j]+".txt")
            rrho_file.to_csv(output_fo, sep="\t",index=False)
            print("RRHO_"+columns_[i]+"-"+columns_[j]+".txt", "generated")


def generate_rrho_file_2(input_f1, input_f2, wd):
    rf1 = pd.read_table(input_f1, index_col=0)
    rf2 = pd.read_table(input_f2, index_col=0)
    cnv1 = input_f1.split("/")[-1].split("_")[0]
    cnv2 = input_f2.split("/")[-1].split("_")[0]
    overlap_index = rf1.index.intersection(rf2.index)
    rf = pd.DataFrame({"Unigene":"N/A", "Gene":overlap_index.tolist(), 
                      cnv1:rf1.loc[overlap_index.tolist(), "log2FoldChange"], 
                      cnv2:rf2.loc[overlap_index.tolist(), "log2FoldChange"]})
#    rf.index = overlap_index    
#    rf[cnv1] = rf1.loc[overlap_index.tolist(), "log2FoldChange"]
#    rf[cnv2] = rf2.loc[overlap_index.tolist(), "log2FoldChange"]
    cnv1_rank = rf[cnv1].rank(ascending=False).astype(int)
    cnv2_rank = rf[cnv2].rank(ascending=False).astype(int)
    rf[cnv1+"_rank"] = cnv1_rank
    rf[cnv2+"_rank"] = cnv2_rank
    rf = rf[["Unigene", "Gene", cnv1+"_rank", cnv2+"_rank", cnv1, cnv2]]
    output_fo = join(wd, "RRHO_"+cnv1+"-"+cnv2+".txt")
    rf.to_csv(output_fo, sep="\t", index=False)


def generate_rrho_file_GSEA(input_f1, input_f2, wd1, wd2, wd):
    # rf1 = pd.read_table(join(wd1, input_f1), index_col=0)
    # rf2 = pd.read_table(join(wd2, input_f2), index_col=0)
    rf1 = input_f1
    rf2 = input_f2
    cnv1 = '-'.join(wd1.split("/")[-1].split("_")[2:5]).split('.')[0]
    cnv2 = '-'.join(wd2.split("/")[-1].split("_")[2:5]).split('.')[0]

    overlap_index = rf1.index.intersection(rf2.index)
#    print(rf2.loc[overlap_index.tolist(), "NES"])
    rf = pd.DataFrame({"Unigene":"N/A", "GeneSets":overlap_index.tolist(), 
                      cnv1: rf1.loc[overlap_index.tolist(), "NES"], 
                      cnv2: rf2.loc[overlap_index.tolist(), "NES"]})
#    print(rf2.loc[overlap_index.tolist(), "NES"])
    cnv1_rank = rf[cnv1].rank(ascending=False).astype(int)
    cnv2_rank = rf[cnv2].rank(ascending=False).astype(int)
    rf[cnv1+"_rank"] = cnv1_rank
    rf[cnv2+"_rank"] = cnv2_rank
    rf = rf[["Unigene", "GeneSets", cnv1+"_rank", cnv2+"_rank", cnv1, cnv2]]
    
    cnv1 = cnv1.split('-')
    cnv2 = cnv2.split('-')
    assert cnv1[0] == cnv2[0] or cnv1[1] == cnv2[1]
    if cnv1[0] == cnv2[0]:
        cnv1 = '-'.join(cnv1)
        cnv2 = '-'.join(cnv2[1:])
    else:
        cnv1 = '-'.join(cnv1)
        cnv2 = cnv2[0]

    output_fo = join(wd, "RRHO_"+cnv1+"_"+cnv2+".txt")
    rf.to_csv(output_fo, sep="\t", index=False)
    print(output_fo)
    return rf
    

def go_through_GSEA(wd):
    gsea_dir = [x[1] for x in walk(wd)][0]
    for i in range(len(gsea_dir)):
        na_pos = ''
        na_neg = ''
        files_i = listdir(join(wd,gsea_dir[i]))
        wd_i = wd + gsea_dir[i]
        for f in files_i:
            if ("gsea_report_for_YES" in f) and (".xls" in f):
                na_pos = f
            if ("gsea_report_for_NO" in f) and (".xls" in f):
                na_neg = f
        df1 = pd.read_table(join(wd_i, na_pos), index_col=0)
        df2 = pd.read_table(join(wd_i, na_neg), index_col=0)
        df_i = pd.concat([df1, df2])
        for j in range(i+1, len(gsea_dir)):
            files_j = listdir(join(wd,gsea_dir[j]))
            wd_j = wd + gsea_dir[j]
            for f in files_j:
                if ("gsea_report_for_YES" in f) and (".xls" in f):
                    na_pos = f
                if ("gsea_report_for_NO" in f) and (".xls" in f):
                    na_neg = f
            df1 = pd.read_table(join(wd_j, na_pos), index_col=0)
            df2 = pd.read_table(join(wd_j, na_neg), index_col=0)
            df_j = pd.concat([df1, df2])
            cnv1 = gsea_dir[i].split("_")[2:5]
            cnv2 = gsea_dir[j].split("_")[2:5]
            if (cnv1[0] == cnv2[0] or cnv1[1] == cnv2[1]):
                generate_rrho_file_GSEA(df_i, df_j, wd_i, wd_j, wd)
            else:
                continue


def geneset_overlap_cancers(wd):
    gsea_dir = [x[1] for x in walk(wd)][0]
    all_cancer_cnv = pd.DataFrame()
    for i in range(len(gsea_dir)):
        cancer_cnv = '-'.join((gsea_dir[i].split('.')[0]).split('_')[2:])
        na_pos = ''
        na_neg = ''
        files_i = listdir(join(wd, gsea_dir[i]))
        wd_i = wd + gsea_dir[i]
        for f in files_i:
            if ("gsea_report_for_YES" in f or "gsea_report_for_na_pos" in f) and (".xls" in f):
                na_pos = f
            if ("gsea_report_for_NO" in f or "gsea_report_for_na_neg" in f) and (".xls" in f):
                na_neg = f
        df1 = pd.read_table(join(wd_i, na_pos), index_col=0)
        df2 = pd.read_table(join(wd_i, na_neg), index_col=0)
        df_i = pd.concat([df1, df2])
        print(df_i['NES'])
        if (len(all_cancer_cnv.index)):
            all_cancer_cnv[cancer_cnv] = df_i['NES']
        else:
            all_cancer_cnv = pd.DataFrame(index=df_i.index)  
            all_cancer_cnv[cancer_cnv] = df_i['NES']
    
    all_cancer_cnv.to_excel(wd+"All_cancers_CNV_NES_Score__.xlsx")
    

if __name__ == "__main__":
#    wd = "C:/Users/wesle/OneDrive/College/Graeber Lab/Genomic_instability/LOH_8p/Correlations"
    # wd = "C:/Users/wesle/OneDrive/College/Graeber Lab/Genomic_instability/LOH_8p/Correlations/DGE_stages/"
    # filename = "GID.txt"
#    input_f = join(wd, filename)
#    generate_rrho_file(input_f, wd)
#    rank_rank_files = [x[2] for x in walk(wd)][0]
#    for i in range(len(rank_rank_files)-2):
#        for j in range(i+1, len(rank_rank_files)-2):
#            if ("YES" in rank_rank_files[i] and "YES" in rank_rank_files[j]):
#                file1 = join(wd, rank_rank_files[i])
#                print(file1)
#                file2 = join(wd, rank_rank_files[j])
#                generate_rrho_file_2(file1, file2, wd)

    wd = "C:/Users/wesle/OneDrive/College/Graeber Lab/Genomic_instability/Winter_2017/Thres_0.2_GSEA/"
    # go_through_GSEA(wd)
    geneset_overlap_cancers(wd)