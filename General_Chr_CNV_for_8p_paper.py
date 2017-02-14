"""
Analysis for BRCA and SKCM chromosome 8p loss data
"""
import matplotlib
import pandas as pd
import numpy as np

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from collections import defaultdict
#import seaborn as sns
from scipy import stats

# import collections

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

def hinton(matrix, max_weight=None, ax=None):
    """Draw Hinton diagram for visualizing a weight matrix."""
    ax = ax if ax is not None else plt.gca()

    if not max_weight:
        max_weight = 2**np.ceil(np.log(np.abs(matrix).max())/np.log(2))

    ax.patch.set_facecolor('lightgray')
    ax.set_aspect('equal', 'box')
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    for (x, y), w in np.ndenumerate(matrix):
        color = 'red' if w > 0 else 'blue'
        size = np.sqrt(np.abs(w))
        rect = plt.Rectangle([x - size / 2, y - size / 2], size, size,
                             facecolor=color, edgecolor=color)
        ax.add_patch(rect)

    nticks = matrix.shape[0]
    ax.xaxis.tick_top()
    ax.set_xticks(range(nticks))
    ax.set_xticklabels(list(matrix.columns), rotation=90)
    ax.set_yticks(range(nticks))
    ax.set_yticklabels(matrix.columns)
    ax.grid(False)

    ax.autoscale_view()
    ax.invert_yaxis()

# filter the normal samples from the BRCA and SKCM data
class Aneuploidy:
    def __init__(self, cancerType, RSEM_Gene_data, SNP_data, chromosome, arm, cond="loss"):
        self.cancer = cancerType
        self.rsem = RSEM_Gene_data
        # self.snp = SNP_data
        self.snp_patients = SNP_data
        self.altered_chr = []
        self.normal_chr = []
        self.cond = cond
        self.chr = chromosome
        self.arm = arm
        self.samples_target = pd.DataFrame()
        self.chr_category = pd.DataFrame()
        self.instability_scores = defaultdict(float)
        self.Instability_score_samples = pd.DataFrame()

    # remove the normal samples from the segment file
    def remove_normal_samples(self):
        index_ = set(self.snp.index.tolist())
        normal_sample = list(filter(lambda i: i[0].split("-")[3] == ("10A" or "10B" or "11A" or "11B" or "12A" or "12B" or "13A" or "13B" or "14A" or "14B"), index_))
        self.snp_patients = self.snp.drop(normal_sample, axis=0)
        print(self.cancer, "patients' samples", len(self.snp_patients.index.tolist()))
        # filter keratin and immunome genes
        immune_keratin_genes = pd.read_excel("/home/rshen/genomic_instability/keratin_immune_etc.xlsx", sheetname="Table S4")
        immune_keratin_gene_list = immune_keratin_genes["Gene"].tolist()
        # immune_keratin_gene_list_ = list(filter(lambda i: i in self.rsem.index.tolist(), immune_keratin_gene_list))
        # self.rsem.drop(immune_keratin_gene_list_, axis=0, inplace=True)
        # print("immune_keratin_genes removed.")
        self.snp_patients.to_csv(self.cancer + str(self.chr) + self.arm + "_" + self.cond + "_patients_sameples_ONLY_" + ".txt", sep="\t")
        return self.snp_patients

    # return two lists of samples
    # first lists with targeted chromosomal CNV, second list without such CNV
    # chr stands for the targeted chromosome we are looking AttributeError
    # threshold stands for the location of probes we are cutting off
    # homozygous_deletion paper, defines threshold of one-copy loss (hemizygous deletion) as <-threshold
    def chr_cnv(self, threshold, threshold_start=0, threshold_end=0):
        self.snp_patients = self.snp_patients.sortlevel()
        index_ = set(filter(lambda x: x[1] == self.chr, self.snp_patients.index.tolist()))
        print("length_index_chr" + str(self.chr) + ": ", len(index_))
        if self.arm == "":
            for i in index_:
                # segment_ends = np.array(self.snp_patients.loc[i].End)
                # segment_starts = np.array(self.snp_patients.loc[i].Start)
                segment_lengths = np.array(self.snp_patients.loc[i].End - self.snp_patients.loc[i].Start)
                segment_means = np.array(self.snp_patients.loc[i].Segment_Mean)
                check = np.dot(segment_lengths, segment_means) / sum(segment_lengths)
                # use the stringent approach
                if self.cond == "loss":
                    if check < -threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold < check < threshold:
                        self.normal_chr.append(i[0])
                elif self.cond == "gain":
                    if check > threshold:
                        self.altered_chr.append(i[0])
                    elif -threshold < check < threshold:
                        self.normal_chr.append(i[0])
        elif self.arm == "p":
            # specific way for grouping to replicate the 8p LOH paper
            for i in index_:
                segment_ends = np.array(self.snp_patients.loc[i].End)
                segment_starts = np.array(self.snp_patients.loc[i].Start)
                segment_means = np.array(self.snp_patients.loc[i].Segment_Mean)
                # select the probes testing the interested chromosomal arm
                segment_end = np.array(list(filter(lambda x: x < threshold_end, segment_ends)) + [threshold_end])
                segment_start = segment_starts[:len(segment_end)]
                segment_lengths = segment_end - segment_start
                segment_length = np.array(list(filter(lambda x: x > 0, segment_lengths)))
                segment_mean = segment_means[:len(segment_length)]
                dict_ = dict(zip(segment_mean, segment_length))
                # check = np.dot(segment_length, segment_mean)
                # if (i[0] == "TCGA-A8-A08B-01A-11D-A011-01"):
                #     print(segment_end, segment_start)
                if self.cond == "loss":
                    keys = list(filter(lambda x: x < -threshold, dict_.keys()))
                    segments = [dict_[key] for key in keys]
                    if sum(segments) >= (threshold_end - threshold_start) / 2:
                        self.altered_chr.append(i[0])
                    # Note: problem here!!! doesn't distinguish gain from normal
                    else:
                        self.normal_chr.append(i[0])
                elif self.cond == "gain":
                    keys = list(filter(lambda x: x > threshold, dict_.keys()))
                    segments = [dict_[key] for key in keys]
                    if sum(segments) >= (threshold_end - threshold_start) / 2:
                        self.altered_chr.append(i[0])
                    # Note: problem here!!! doesn't distinguish gain from normal
                    else:
                        self.normal_chr.append(i[0])
        elif self.arm == "q":
            for i in index_:
                segment_ends = np.array(self.snp_patients.loc[i].End)
                segment_starts = np.array(self.snp_patients.loc[i].Start)
                segment_means = np.array(self.snp_patients.loc[i].Segment_Mean)
                segment_start = np.array(list(filter(lambda x: x > threshold_start, segment_starts)) + [threshold_start])
                segment_end = segment_ends[::-1][:len(segment_start)]
                segment_lengths = segment_end - segment_start
                segment_length = np.array(list(filter(lambda x: x > 0, segment_lengths)))
                segment_mean = segment_means[::-1][:len(segment_length)]
                dict_ = dict(zip(segment_mean, segment_length))
                # check = np.dot(segment_lengths,segment_means)/sum(segment_lengths)
                if self.cond == "loss":
                    keys = list(filter(lambda x: x < -threshold, dict_.keys()))
                    segments = [dict_[key] for key in keys]
                    if sum(segments) >= (threshold_end - threshold_start) / 2:
                        self.altered_chr.append(i[0])
                    # Note: problem here!!! doesn't distinguish gain from normal
                    else:
                        self.normal_chr.append(i[0])
                elif self.cond == "gain":
                    keys = list(filter(lambda x: x > threshold, dict_.keys()))
                    segments = [dict_[key] for key in keys]
                    if sum(segments) >= (threshold_end - threshold_start) / 2:
                        self.altered_chr.append(i[0])
                    # Note: problem here!!! doesn't distinguish gain from normal
                    else:
                        self.normal_chr.append(i[0])
        print(self.cancer + "_chromosome_" + str(self.chr) + self.arm + self.cond + "altered samples: ",
              len(self.altered_chr), '\n', self.cancer + "_chromosome_" + "normal samples: ", len(self.normal_chr))
        return self.altered_chr, self.normal_chr

    def calculate_Instability_score(self):
        sample_index = set([x[0] for x in self.snp_patients.index.tolist()])
        print("Number of patient samples to calculate instability score:", len(sample_index))
        for i in sample_index:
            segments_lens = np.array(self.snp_patients.loc[i].End - self.snp_patients.loc[i].Start)
            segments_means = np.array(abs(self.snp_patients.loc[i].Segment_Mean))
            instability_score = np.dot(segments_lens, segments_means)
            skimmed_index = "-".join(i.split("-")[0:4])
            self.instability_scores[skimmed_index] = instability_score
        rsem_cols = self.rsem.columns.tolist()
        rsem_cols = ["-".join(x.split("-")[0:4]) for x in rsem_cols]
        self.rsem.columns = uniquify(rsem_cols)
        self.instability_scores = {k: v for k, v in self.instability_scores.items() if k in self.rsem.columns}
        self.Iscore = pd.DataFrame.from_dict(self.instability_scores, orient='index')
        self.Iscore.columns = ["instability_score"]
        self.Iscore.sort_values(axis=0, by="instability_score", ascending=False, inplace=True)
        self.Instability_score_samples = self.rsem[self.Iscore.index.tolist()]

    def set_samples_altered(self, indexCol):
        samples = []
        try:
            self.rsem = self.rsem.reset_index().drop_duplicates(subset=indexCol, keep='last').set_index(indexCol)
        except:
            print("Wrong input for index column name")

        # remove duplicate columns and rows
        rsem_cols = self.rsem.columns.tolist()
        rsem_cols = ["-".join(x.split("-")[0:4]) for x in rsem_cols]
        self.rsem.columns = uniquify(rsem_cols)

        self.altered_chr = ["-".join(x.split("-")[0:4]) for x in self.altered_chr]
        self.altered_chr = list(filter(lambda i: i.split("_")[0] in self.altered_chr, self.rsem.columns))
        #        self.altered_chr = [x[:-1] for x in self.altered_chr]
        self.normal_chr = ["-".join(x.split("-")[0:4]) for x in self.normal_chr]
        self.normal_chr = list(filter(lambda i: i.split("_")[0] in self.normal_chr, self.rsem.columns))
        #        self.no_altered_chr = [x[:-1] for x in self.no_altered_chr]

        # check if the loss of 8p samples by computation consistent with eyeballing
        #loss_8p = pd.read_excel("/home/rshen/genomic_instability/chromosome8p/LOH_8p_paper/BRCA_group_by_8p_loss.xlsx", index_col=0)
        #loss_8p_samples = loss_8p.loc[loss_8p["8p_status"] == "L"].index.tolist()
        #overlaps = list(filter(lambda i: i.split("_")[0] in loss_8p_samples, self.altered_chr))
        #print("Overlaps between eyeballing and computation: ", len(overlaps))
        #print("False positive: ", list(filter(lambda i: i.split("_")[0] in loss_8p_samples, self.normal_chr)))

        cnv_samples = self.altered_chr + self.normal_chr
        print(self.cancer + "_chromosome_" + str(self.chr) + self.arm + self.cond + "altered samples length: ",
              len(self.altered_chr))
        print(self.cancer + "_chromosome_" + "normal samples length: ", len(self.normal_chr))
        #        rsem_cols = [x[:-1] for x in rsem_cols]

        for i in cnv_samples:
            if i in self.rsem.columns.tolist():
                samples.append(i)
        self.samples_target = self.rsem[samples]

    def set_category(self, condition):
        self.chr_category = pd.DataFrame({"Sample_ID": self.samples_target.columns.tolist()})
        self.chr_category.set_index("Sample_ID", inplace=True)
        for i in self.chr_category.index.tolist():
            if (i in self.altered_chr):
                self.chr_category.set_value(i, "chr" + str(self.chr) + self.arm + "_CNV", "YES")
            elif (i in self.normal_chr):
                self.chr_category.set_value(i, "chr" + str(self.chr) + self.arm + "_CNV", "NO")
            else:
                print(i, "Error! Sample doesn't relate to chromosome " + str(self.chr) + self.arm)
        self.chr_category.sort_values(axis=0, by="chr" + str(self.chr) + self.arm + "_CNV", inplace=True)
        self.samples_target = self.samples_target[self.chr_category.index.tolist()]

    # put normalized or raw_counts in condition
    def output_(self, condition):
        self.chr_category.to_csv(
            self.cancer + str(self.chr) + self.arm + "_" + self.cond + "_category_0.2" + condition + ".txt", sep="\t")
        self.samples_target.to_csv(
            self.cancer + str(self.chr) + self.arm + "_" + self.cond + "_sameples_0.2" + condition + ".txt", sep="\t")
        # self.Iscore.to_csv(self.cancer+"_Instability_Score_" + ".txt", sep="\t")
        # self.Instability_score_samples.to_csv(self.cancer+"_Instability_Score_samples" + ".txt", sep="\t")

    def pca_plot(self):
        pca = PCA(n_components=4, whiten=True)
        transf = pca.fit_transform(self.samples_target.T)
        variance_ratio = pca.explained_variance_ratio_
        loadings = pca.components_
        fig_sample = plt.gcf()
        fig_sample.set_size_inches(14, 14)
        # print(transf)
        plt.xlabel("PC1: " + str(round(variance_ratio[0], 3)))
        plt.ylabel("PC2: " + str(round(variance_ratio[1], 3)))
        colName = self.samples_target.columns.tolist()
        for n in range(len(colName)):
            if (colName[n] in self.altered_chr):
                # print("altered plot")
                altered, = plt.plot(transf[n, 0], transf[n, 1], marker='o', markersize=8, color='red', alpha=1,
                                    label="chromosome" + str(self.chr) + self.arm + "_" + self.cond)
            if (colName[n] in self.normal_chr):
                # print("normal plot")
                normal, = plt.plot(transf[n, 0], transf[n, 1], marker='o', markersize=8, color='blue', alpha=1,
                                   label='normal_samples')
        plt.legend(loc='best', scatterpoints=1, handles=[altered, normal])
        plt.title("PCA_plot_" + self.cancer + "_" + str(self.chr) + self.arm + "_" + self.cond)
        fig_sample.savefig("PCA_" + self.cancer + '_' + str(self.chr) + self.arm + "_" + self.cond + "_Feb_12.png", dpi=100)
        PCA_loadings = pd.DataFrame(loadings, index=["PC1", "PC2", "PC3", "PC4"],
                                    columns=self.samples_target.index.tolist())
        #        print(self.cancer, "loadings", loadings)
        PCA_loadings.to_csv(self.cancer + "_" + str(self.chr) + self.arm + "_" + self.cond + "_PCA_loadings_Feb_12.txt", sep="\t")

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

def GNI(tumor, chr, arm, var, seg, RNA_, CNV_cutoff, start, end):
    aneuploidy = Aneuploidy(tumor, RNA_, seg, chr, arm, var)
    aneuploidy.remove_normal_samples
    if (arm == "p"):
        aneuploidy.chr_cnv(threshold=CNV_cutoff, threshold_start=start, threshold_end=end)
    elif (arm == "q"):
        aneuploidy.chr_cnv(threshold=CNV_cutoff, threshold_start=start, threshold_end=end)
    else:
        aneuploidy.chr_cnv(threshold=CNV_cutoff)
    print(tumor, chr, arm, "Grouping done.")
    aneuploidy.set_samples_altered("GeneSymbol")
    print(tumor, chr, arm, "samples filtering done.")
    aneuploidy.set_category("normalized")
    print(tumor, chr, arm, "GSEA preparation done.")
    aneuploidy.pca_plot()
    aneuploidy.output_("normalized")
    return aneuploidy

if __name__ == '__main__':
    # Investigate the CNV in chromosome 1q gain, 3 loss, 6p gain, 6q loss, 8p loss, 8q gain, 9p loss, 18q loss
    """
    chr_alter_dict = {"loss": [(3, ''), (6, 'q'), (8, 'p'), (9, 'p'), (18, 'q')],
                      "gain": [(1, 'q'), (6, 'p'), (8, 'q'),(5, 'q')]}

    # each chromosome arm's cutoff value in segment files
    chr_arm_cufoff = {(3, ''): (0, 0), (6, 'q'): (6.0E7, 1.7E8), (8, 'p'): (2E7, 4.5E7), (9, 'p'): (5.0E7, 1.4E8), (6, 'p'): (1E6, 6.0E7), 
                      (8, 'q'): (4.8E7, 1.5E8), (1, 'q'): (1.3E8, 2.5E8), (5, 'q'): (5E7, 1.8E8),  (18, 'q'): (1.9E7, 7.7E7)}

    genomic_instability_df = pd.DataFrame()

    BRCA_RNA = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/BRCA_normalized_results_simplified.txt", index_col=0)
    BRCA_ = pd.read_table("./BRCA8p_loss_patients_sameples_ONLY_.txt", index_col=[0, 1])
    # LOH case, CNV_threshold 0.2
    # BRCA case
    for variation in chr_alter_dict.keys():
        threshold = 0.3
        if variation == 'gain':
            threshold = 0.58 
        for chr_arm in chr_alter_dict[variation]:
            aneuploidy = GNI("BRCA", chr_arm[0], chr_arm[1], variation, BRCA_, BRCA_RNA, threshold, start=chr_arm_cufoff[chr_arm][0], end=chr_arm_cufoff[chr_arm][1])
            samples = aneuploidy.samples_target
            altered_chr = aneuploidy.altered_chr
            altered_samples = samples[altered_chr]
            genomic_instability_df['{}_{}_{}'.format(chr_arm[0],chr_arm[1],variation)] = altered_samples.mean(axis=1)
            print(genomic_instability_df)
    """
    # genomic_instability_df.to_csv("GID.txt",sep='\t')
    genomic_instability_df = pd.read_table("GID.txt",index_col=0)
    # correlation approach doesn't work here!
    corrmat = genomic_instability_df.corr()
    corrmat.to_csv("corr_only.txt",sep='\t')
    print("corrmat finished")
    #fig1 = plt.figure()
    #sns.heatmap(corrmat, vmax=1., square=False).xaxis.ticktop()
    #plt.savefig("genomic_instability_corr.png", dpi=200, bbox_inches='tight')
    
    fig2 = plt.figure()
    hinton(corrmat)
    plt.savefig("genomic_instability_corr_hinton.png", dpi=200, bbox_inches='tight')
    
    correlation_p_value = pd.DataFrame()
    #correlation_p_value.set_index(['p-value','corr'])
    for i in (genomic_instability_df.columns):
        for j in (genomic_instability_df.columns):
            correlation_p_value['{}_{}'.format(i, j)] = stats.pearsonr(genomic_instability_df[i], genomic_instability_df[j])
    
    correlation_p_value.to_csv("corr_p_val.txt",sep='\t')
    print("pearson finished")
    # read in the segment file and RNA data
    # BRCA_ = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/BRCA__CNV.seg.txt", index_col=[0,1])

    # SKCM_ = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/SKCM__CNV.seg.txt", index_col=[0,1])
    # SKCM_RNA = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/SKCM_normalized_results_simplified.txt", index_col=0)

    # UVM_ = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/UVM__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Aug-04-2015.seg.txt", index_col=[0,1])
    # UVM_RNA = pd.read_table("/home/rshen/genomic_instability/chromosome8p/TCGA_data/UVM_normalized_results_processed_No_keratin_immune.txt", index_col=0)

    # # calculate instability scores
    # BRCA_lossOf8p = Aneuploidy("BRCA", BRCA_RNA, BRCA_, 1, "q", "gain")
    # BRCA_lossOf8p.remove_normal_samples()
    # BRCA_lossOf8p.calculate_Instability_score()
    # BRCA_lossOf8p.output_("")
    # BRCA_lossOf8p.chr_CNV(threshold_start=3.6E7)
    # print("BRCA Grouping done.")
    # BRCA_lossOf8p.set_samples_altered("GeneSymbol")
    # print("BRCA samples filtering done.")
    # BRCA_lossOf8p.set_category("normalized")
    # print("BRCA GSEA preparation done.")
    # BRCA_lossOf8p.PCA_plot()
    # BRCA_lossOf8p.output_("normalized")

    # SKCM_lossOf18q = Aneuploidy("SKCM", SKCM_RNA, SKCM_, 18, "q", "loss")
    # SKCM_lossOf18q.remove_normal_samples()
    # SKCM_lossOf18q.calculate_Instability_score()
    # SKCM_lossOf18q.output_("")

    # UVM_lossOf18q = Aneuploidy("UVM", UVM_RNA, UVM_, 18, "q", "loss")
    # UVM_lossOf18q.remove_normal_samples()
    # UVM_lossOf18q.calculate_Instability_score()
    # UVM_lossOf18q.output_("")

    # print("end of ISCORE")
    """
    BRCA_lossOf8p = Aneuploidy("BRCA", BRCA_RNA, BRCA_, 8, "p", "loss")
    # BRCA_lossOf8p.remove_normal_samples()
    # BRCA_lossOf8p.calculate_Instability_score()
    BRCA_lossOf8p.chr_cnv(threshold=0.3, threshold_start=0, threshold_end=4.5E7)
    print("BRCA Grouping done.")
    BRCA_lossOf8p.set_samples_altered("GeneSymbol")
    print("BRCA samples filtering done.")
    BRCA_lossOf8p.set_category("normalized")
    print("BRCA GSEA preparation done.")
    BRCA_lossOf8p.pca_plot()
    BRCA_lossOf8p.output_("normalized")

    cancer_stage = pd.read_excel("/home/rshen/genomic_instability/chromosome8p/LOH_8p_paper/BRCA.clin.merged.xlsx", sheetname='Sheet2')
    cancer_stage['patient.stage_event.tnm_categories.pathologic_categories.pathologic_t'] = [i[:2] if len(i) > 2 else i
                         for i in cancer_stage['patient.stage_event.tnm_categories.pathologic_categories.pathologic_t']]
    stage_group = cancer_stage.groupby("patient.stage_event.tnm_categories.pathologic_categories.pathologic_t")
    altered_samples = ['-'.join(i.split('-')[0:3]).lower() for i in BRCA_lossOf8p.altered_chr]
    normal_samples = ['-'.join(i.split('-')[0:3]).lower() for i in BRCA_lossOf8p.normal_chr]
    for key in stage_group.groups.keys():
        target_samples = list(cancer_stage.loc[stage_group.groups[key]]['patient.bcr_patient_barcode'])
        print(key, "loss of 8p", len(set(target_samples).intersection(altered_samples))/len(target_samples))
        print(key, "normal", len(set(target_samples).intersection(normal_samples))/len(target_samples))

    from sys import exit
    exit(0)
    """


    # SKCM_lossOf18q = Aneuploidy("SKCM", SKCM_RNA, SKCM_, 18, "q", "loss")
    # SKCM_lossOf18q.remove_normal_samples()
    # SKCM_lossOfBRCA_lossOf8p.calculate_Instability_score()

    # UVM_lossOf18q = Aneuploidy("UVM", UVM_RNA, UVM_, 18, "q", "loss")
    # UVM_lossOf18q.remove_normal_samples()
    # UVM_lossOfBRCA_lossOf8p.calculate_Instability_score()
    # SKCM_lossOf18q.chr_CNV(threshold_start=1.4E7)
    # print("SKCM Grouping done.")
    # SKCM_lossOf18q.set_samples_altered("Genes")
    # print("SKCM samples filtering done.")
    # SKCM_lossOf18q.set_category("normalized")
    # print("SKCM GSEA preparation done.")
    # SKCM_lossOf18q.PCA_plot()
    # SKCM_lossOf18q.output_("normalized")
