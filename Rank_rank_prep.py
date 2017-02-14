 # -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 10:20:33 2016

@author: Wesley
"""

import pandas as pd
from os.path import join

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

if __name__ == "__main__":
    wd = "C:/Users/wesle/OneDrive/College/Graeber Lab/Genomic_instability/LOH_8p/Correlations"
    filename = "GID.txt"
    input_f = join(wd, filename)
    generate_rrho_file(input_f, wd)
