 # -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 10:20:33 2016

@author: Wesley
"""

import pandas as pd
from os.path import join
from os import walk

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
    
if __name__ == "__main__":
#    wd = "C:/Users/wesle/OneDrive/College/Graeber Lab/Genomic_instability/LOH_8p/Correlations"
    wd = "C:/Users/wesle/OneDrive/College/Graeber Lab/Genomic_instability/LOH_8p/Correlations/DGE_stages/"
#    filename = "GID.txt"
#    input_f = join(wd, filename)
#    generate_rrho_file(input_f, wd)
    rank_rank_files = [x[2] for x in walk(wd)][0]
    for i in range(len(rank_rank_files)-2):
        for j in range(i+1, len(rank_rank_files)-2):
            if ("YES" in rank_rank_files[i] and "YES" in rank_rank_files[j]):
                file1 = join(wd, rank_rank_files[i])
                print(file1)
                file2 = join(wd, rank_rank_files[j])
                generate_rrho_file_2(file1, file2, wd)
