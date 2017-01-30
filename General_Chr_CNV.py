"""
Analysis for BRCA and SKCM chromosome 8p loss data
"""
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
#from collections import OrderedDict
from sklearn.decomposition import PCA
#import collections
#import filter_functions as ff

# remove the duplicate columns
def uniquify(df_columns):
    seen = set()

    for item in df_columns:
        fudge = 1
        newitem = item

        while newitem in seen:
            fudge += 1
            newitem = "{}_{}".format(item, fudge)

        yield newitem
        seen.add(newitem)
    return list(seen)

# filter the normal samples from the BRCA and SKCM data
class Aneuploidy:

    def __init__(self, cancerType, RSEM_Gene_data, SNP_data, chromosome, arm, cond="loss"):
        self.cancer = cancerType
        self.rsem = RSEM_Gene_data
        self.snp = SNP_data
        self.snp_patients = pd.DataFrame()
        self.altered_chr = []
        self.normal_chr = []
        self.cond = cond
        self.chr = chromosome
        self.arm = arm
        self.samples_target = pd.DataFrame()
        self.chr_category = pd.DataFrame()
    
    # remove the normal samples from the segment file
    def remove_normal_samples(self):
        index_ = set(self.snp.index.tolist())
        normal_sample = list(filter(lambda i : i[0].split("-")[3] == ("10A"or"10B"or"11A"or"11B"or"12A"or"12B"or"13A"or"13B"or"14A"or"14B"), index_))       
        self.snp_patients = self.snp.drop(normal_sample, axis=0)
        print("patients' samples", len(set(self.snp_patients.index.tolist())))
        return self.snp_patients
        
    # return two lists of samples
    # first lists with targeted chromosomal CNV, second list without such CNV
    # chr stands for the targeted chromosome we are looking AttributeError
    # threshold stands for the location of probes we are cutting off
    # homozygous_deletion paper, defines threshold of one-copy loss (hemizygous deletion) as <-threshold
    def chr_CNV(self, threshold, threshold_start=0, threshold_end=0):
        self.snp_patients = self.snp_patients.sortlevel()
        index_ = set(filter(lambda x: x[1] == self.chr, self.snp_patients.index.tolist()))
        print("length_index_"+str(self.chr)+": ", len(index_))
        if self.arm == "":
            for i in index_:
                segment_lengths = np.array(self.snp_patients.loc[i].End-self.snp_patients.loc[i].Start)
                segment_means = np.array(self.snp_patients.loc[i].Segment_Mean)
                check = np.dot(segment_lengths,segment_means)/sum(segment_lengths)
                if self.cond == "loss":
                    if check < -threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold<check<threshold:
                        self.normal_chr.append(i[0])
                elif self.cond == "gain":
                    if check > threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold<check<threshold:
                        self.normal_chr.append(i[0])
        elif self.arm == "p":
            for i in index_:
                segment_lengths = np.array(self.snp_patients.loc[i][self.snp_patients.loc[i].End<threshold_end].End - 
                                           self.snp_patients.loc[i][self.snp_patients.loc[i].End<threshold_end].Start)
                segment_means = np.array(self.snp_patients.loc[i][self.snp_patients.loc[i].End<threshold_end].Segment_Mean)
                check = np.dot(segment_lengths,segment_means)/sum(segment_lengths)              
                if self.cond == "loss":
                    if check < -threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold<check<threshold:
                        self.normal_chr.append(i[0])
                elif self.cond == "gain":
                    if check > threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold<check<threshold:
                        self.normal_chr.append(i[0])
        elif self.arm == "q":
            for i in index_:
                segment_lengths = np.array(self.snp_patients.loc[i][self.snp_patients.loc[i].Start>threshold_start].End - 
                                           self.snp_patients.loc[i][self.snp_patients.loc[i].Start>threshold_start].Start)
                segment_means = np.array(self.snp_patients.loc[i][self.snp_patients.loc[i].Start>threshold_start].Segment_Mean)
                check = np.dot(segment_lengths,segment_means)/sum(segment_lengths)
                if self.cond == "loss":
                    if check < -threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold<check<threshold:
                        self.normal_chr.append(i[0])
                elif self.cond == "gain":
                    if check > threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold<check<threshold:
                        self.normal_chr.append(i[0])
        
        print(self.cancer+"_chromosome_"+str(self.chr)+self.arm+self.cond+"altered samples: ", len(self.altered_chr), '\n', 
              self.cancer+"_chromosome_"+"normal samples: ", len(self.normal_chr))
        return self.altered_chr, self.normal_chr
    
    def set_samples_altered(self, indexCol):
        samples = []
        try:
            self.rsem = self.rsem.reset_index().drop_duplicates(subset = indexCol, keep='last').set_index(indexCol)
        except:
            print("Wrong input for index column name")
            
        # remove duplicate columns and rows
        rsem_cols = self.rsem.columns.tolist()
        rsem_cols = ["-".join(x.split("-")[0:4]) for x in rsem_cols]
        self.rsem.columns = uniquify(rsem_cols)
        
        self.altered_chr = ["-".join(x.split("-")[0:4]) for x in self.altered_chr]
        self.altered_chr = list(filter(lambda i : i.split("_")[0] in self.altered_chr, self.rsem.columns))
#        self.altered_chr = [x[:-1] for x in self.altered_chr]
        self.normal_chr = ["-".join(x.split("-")[0:4]) for x in self.normal_chr]
        self.normal_chr = list(filter(lambda i : i.split("_")[0] in self.normal_chr, self.rsem.columns))
#        self.no_altered_chr = [x[:-1] for x in self.no_altered_chr]
        
#        self.altered_chr = list(set(self.altered_chr))
#        self.no_altered_chr = list(set(self.no_altered_chr))
        
#        print("altered_chr samples: ", self.altered_chr, '\n', "No_altered_chr samples: ", self.no_altered_chr)
        CNV_samples = self.altered_chr + self.normal_chr
        print(self.cancer+"_chromosome_"+str(self.chr)+self.arm+self.cond+"altered samples length: ", len(self.altered_chr))
        print(self.cancer+"_chromosome_"+"normal samples length: ", len(self.normal_chr))
#        rsem_cols = [x[:-1] for x in rsem_cols]

        for i in CNV_samples:
            if i in self.rsem.columns.tolist():
                samples.append(i)
        self.samples_target = self.rsem[samples]
#        self.samples_8p.sort_index(axis=1,inplace=True)
        
        # filter keratin and immunome genes
        immune_keratin_genes = pd.read_excel("/home/rshen/genomic_instability/keratin_immune_etc.xlsx",sheetname="Table S4")
        immune_keratin_gene_list = immune_keratin_genes["Gene"].tolist()
        for i in immune_keratin_gene_list:
            try:
                self.samples_8p = self.samples_8p.drop([i], axis=0)
            except:
                continue        
        #        print(self.samples_8p)
        return self.samples_target
#                
    def set_category(self, condition):
        self.chr_category = pd.DataFrame({"Sample_ID": self.samples_target.columns.tolist()})
        self.chr_category.set_index("Sample_ID",inplace=True)
        for i in self.chr_category.index.tolist():
            if (i in self.altered_chr):
                self.chr_category.set_value(i, "chr"+str(self.chr)+self.arm+"_CNV", "YES")
            elif (i in self.normal_chr):
                self.chr_category.set_value(i, "chr"+str(self.chr)+self.arm+"_CNV", "NO")
            else:
                print(i, "Error! Sample doesn't relate to chromosome "+str(self.chr)+self.arm)
        self.chr_category.sort_values(axis=0, by="chr"+str(self.chr)+self.arm+"_CNV", inplace=True)
        self.samples_target = self.samples_target[self.chr_category.index.tolist()]
        
    # put normalized or raw_counts in condition
    def output_(self, condition):
        self.chr_category.to_csv(self.cancer+"_"+str(self.chr)+self.arm+"_"+self.cond+"_category_" + condition + ".txt", sep="\t")
        print(self.cancer+"_"+str(self.chr)+self.arm+"_"+self.cond+"_category","finish output")
        self.samples_target.to_csv(self.cancer+"_"+str(self.chr)+self.arm+"_"+self.cond+"_sameples_" + condition + ".txt", sep="\t")
        print(self.cancer+"_"+str(self.chr)+self.arm+"_"+self.cond+"_samples","finish output")

    def PCA_plot(self):
        pca = PCA(n_components=4, whiten=True)
        transf = pca.fit_transform(self.samples_target.T)
        variance_ratio = pca.explained_variance_ratio_
        loadings = pca.components_
        fig_sample = plt.gcf()
        fig_sample.set_size_inches(14,14)
        #print(transf)
        plt.xlabel("PC1: " + str(round(variance_ratio[0],3)))
        plt.ylabel("PC2: " + str(round(variance_ratio[1],3)))
        colName = self.samples_target.columns.tolist()
        for n in range(len(colName)):
            if (colName[n] in self.altered_chr):
                altered, = plt.plot(transf[n,0],transf[n,1], markersize=8, color='red', alpha=1, label="chromosome"+str(self.chr)+self.arm+"_"+self.cond)
            if (colName[n] in self.normal_chr):
                normal, = plt.plot(transf[n,0],transf[n,1], markersize=8, color='blue', alpha=1, label='normal_samples')
        plt.legend(loc='best', scatterpoints=1, handles=[altered, normal])
        plt.title("PCA_plot_" + self.cancer + "_"  + str(self.chr)+self.arm+"_"+self.cond)
        fig_sample.savefig("PCA_" + self.cancer + '_' + str(self.chr)+self.arm+"_"+self.cond + "_Jan_22.png", dpi=100)
        PCA_loadings = pd.DataFrame(loadings, index=["PC1", "PC2","PC3","PC4"], columns=self.samples_target.index.tolist())
    #        print(self.cancer, "loadings", loadings)
        PCA_loadings.to_csv(self.cancer + "_" +str(self.chr)+self.arm+"_"+self.cond +"_PCA_loadings_Jan_22.txt", sep="\t")
        
    #    Incrementalpca = PCA(n_components=4, whiten=True)
    #    transf = pca.fit_transform(self.samples_target.T)
    #    variance_ratio = pca.explained_variance_ratio_
    #    loadings = pca.components_
    #    fig_sample = plt.gcf()
    #    fig_sample.set_size_inches(14,14)
    #    plt.xlabel("PC1: " + str(round(variance_ratio[0],3)))
    #    plt.ylabel("PC2: " + str(round(variance_ratio[1],3)))
    #    colName = self.samples_target.columns.tolist()
    #     for n in range(len(colName)):
    #         if (colName[n] in self.altered_chr):
    #             print("altered plot")
    #             altered, = plt.plot(transf[n,0],transf[n,1], markersize=8, color='red', alpha=1, label="chromosome"+str(self.chr)+self.arm+"_"+self.cond)
    #         if (colName[n] in self.normal_chr):
    #             print("normal plot")
    #             normal, = plt.plot(transf[n,0],transf[n,1], markersize=8, color='blue', alpha=1, label='normal_samples')
    #     plt.legend(loc='best', scatterpoints=1, handles=[altered, normal])
    #     plt.title("PCA_plot_" + self.cancer + "_"  + str(self.chr)+self.arm+"_"+self.cond)
    #     fig_sample.savefig("PCA_" + self.cancer + '_' + str(self.chr)+self.arm+"_"+self.cond + "_Jan_22.png", dpi=100)
    #     PCA_loadings = pd.DataFrame(loadings, index=["PC1", "PC2","PC3","PC4"], columns=self.samples_target.index.tolist())
    # #        print(self.cancer, "loadings", loadings)
    #     PCA_loadings.to_csv(self.cancer + "_" +str(self.chr)+self.arm+"_"+self.cond +"_PCA_loadings_Jan_22.txt", sep="\t")

def GNI(tumor, chr, arm, var, seg, RNA_, CNV_cutoff, chr_arm_cutoff):
    aneuploidy = Aneuploidy(tumor, RNA_, seg, chr, arm, var)
    aneuploidy.remove_normal_samples()
    if (arm == "p"):
        aneuploidy.chr_CNV(threshold=CNV_cutoff, threshold_end=chr_arm_cutoff)
    elif (arm == "q"):
        aneuploidy.chr_CNV(threshold=CNV_cutoff, threshold_start=chr_arm_cutoff)
    else:
        aneuploidy.chr_CNV(threshold=CNV_cutoff)

    print(tumor,chr, arm, "Grouping done.")
    aneuploidy.set_samples_altered("GeneSymbol")
    print(tumor, chr, arm, "samples filtering done.")
    aneuploidy.set_category("normalized")
    print(tumor, chr, arm, "GSEA preparation done.")
    #aneuploidy.PCA_plot()
    aneuploidy.output_("normalized")
    print("Output done.")

if __name__ == '__main__':
    # Investigate the CNV in chromosome 1q gain, 3 loss, 6p gain, 6q loss, 8p loss, 8q gain, 9p loss, 18q loss
    chr_alter_dict = {"loss": [(3,''),(6,'q'),(8,'p'),(9,'p'),(18,'q')],"gain":[(1,'q'),(6,'p'),(8,'q')]}

    # each chromosome arm's cutoff value in segment files
    chr_arm_cufoff = {(3,''):0,(6,'q'):6.0E7,(8,'p'):4.5E7,(9,'p'):5.0E7,(6,'p'):6.0E7,(8,'q'):5.0E7,(1,'q'):1.5E8,(18,'q'):1.5E7}

    # read in the segment file and RNA data
    BRCA_ = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/BRCA__CNV.seg.txt", index_col=[0,1])
    BRCA_RNA = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/BRCA_normalized_results_simplified.txt", index_col=0)

    SKCM_ = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/SKCM__CNV.seg.txt", index_col=[0,1])
    SKCM_RNA = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/SKCM_normalized_results_simplified.txt", index_col=0)

    UVM_ = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/UVM__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Aug-04-2015.seg.txt", index_col=[0,1])
    UVM_RNA = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/UVM_normalized_results_processed_No_keratin_immune.txt", index_col=0)

    # LOH case, CNV_threshold 0.5
    # BRCA case
    
    for variation in chr_alter_dict.keys():
        for chr_arm in chr_alter_dict[variation]:
            GNI("BRCA", chr_arm[0], chr_arm[1], variation, BRCA_, BRCA_RNA, 0.5, chr_arm_cufoff[chr_arm])
    
    for variation in chr_alter_dict.keys():
        for chr_arm in chr_alter_dict[variation]:
            GNI("SKCM", chr_arm[0], chr_arm[1], variation, SKCM_, SKCM_RNA, 0.5, chr_arm_cufoff[chr_arm])

    for variation in chr_alter_dict.keys():
        for chr_arm in chr_alter_dict[variation]:
            GNI("UVM", chr_arm[0], chr_arm[1], variation, UVM_, UVM_RNA, 0.5, chr_arm_cufoff[chr_arm])

    # BRCA_gainof1q = Aneuploidy("BRCA", BRCA_RNA, BRCA_, 1, "q", "gain")
    # BRCA_gainof1q.remove_normal_samples()
    # BRCA_gainof1q.chr_CNV(threshold_start=3.6E7)
    # print("BRCA Grouping done.")
    # BRCA_gainof1q.set_samples_altered("GeneSymbol")
    # print("BRCA samples filtering done.")
    # BRCA_gainof1q.set_category("normalized")
    # print("BRCA GSEA preparation done.")
    # BRCA_gainof1q.PCA_plot()
    # BRCA_gainof1q.output_("normalized")

    # SKCM_lossOf18q = Aneuploidy("SKCM", SKCM_RNA, SKCM_, 18, "q", "loss")
    # SKCM_lossOf18q.remove_normal_samples()
    # SKCM_lossOf18q.chr_CNV(threshold_start=1.4E7)
    # print("SKCM Grouping done.")
    # SKCM_lossOf18q.set_samples_altered("Genes")
    # print("SKCM samples filtering done.")
    # SKCM_lossOf18q.set_category("normalized")
    # print("SKCM GSEA preparation done.")
    # SKCM_lossOf18q.PCA_plot()
    # SKCM_lossOf18q.output_("normalized")

