## The script calculates the genomic instability signatures (CNV) of different cancer types
## Run the DGE analysis based on the CNV of patients

library(DESeq2)
library(made4)

# read in the segment mean file
segFile <- read.table("BRCA.snp__genome_cnv_hg18__seg.seg.txt", as.is = TRUE)
colnames(segFile) <- as.character(unlist(segFile[1,]))
segFile <- segFile[-1,]


cleanUp_SampleNames <- function(x){
  split_list <- strsplit(x, split="-")
  new_name <- paste(unlist(split_list)[1:4], collapse="-")
  new_name <- substr(new_name, 1, nchar(new_name)-1)
  return(new_name)
}

cleanUp.SampleNames <- function(x){
  split_list <- strsplit(x, split="[.]")
  new_name <- paste(unlist(split_list)[1:4], collapse="-")
  new_name <- substr(new_name, 1, nchar(new_name)-1)
  return(new_name)
}

groupByCNV <- function(segfile, chr, arm=NULL, thres, arm_start=0, arm_end=0, immune=FALSE){
  # the DNA segment mean data of the patients
  cnv_cond <- data.frame(row.names = unique(segfile$Sample))
  cnv_cond$CNV <- -Inf
  
  for (i in unique(segfile$Sample)){
    print(i)
    if (arm == "p"){
      # potential problems, need to be double checked
      # as.numeric() function vs. factor levels
      subset <- subset(segfile, Sample==i & Chromosome==chr & as.numeric(End) <= arm_end)
      head(subset)
      subset$segment_lengths <- as.numeric(subset$End) - as.numeric(subset$Start)
      head(subset)
      cnv_score <- sum(as.numeric(subset$segment_lengths) * as.numeric(subset$Segment_Mean))/sum(subset$segment_length)
      print("cnv_score: ")
      print(cnv_score)
      
      if (is.nan(cnv_score)) next
      
      if (thres > 0){
        if (cnv_score > thres){
          cnv_cond[i,] = 1
        }
        else if (cnv_score < -thres){
          cnv_cond[i,] = -1
        }
        else{
          cnv_cond[i,] = 0
        }
      }
      else{
        if (cnv_score < thres){
          cnv_cond[i,] = 1
        }
        else if (cnv_score > -thres){
          cnv_cond[i,] = -1
        }
        else{
          cnv_cond[i,] = 0
        }
      }
    }
    
    else if (arm == "q"){
      subset <- subset(segfile, Sample==i & Chromosome==chr & as.numeric(Start) >= arm_start)
      subset$segment_lengths <- as.numeric(subset$End) - as.numeric(subset$Start)
      cnv_score <- sum(as.numeric(subset$segment_lengths) * as.numeric(subset$Segment_Mean))/sum(subset@segment_length)
      
      if (is.nan(cnv_score)) next
      
      if (thres > 0){
        if (cnv_score > thres){
          cnv_cond[i,] = 1
        }
        else if (cnv_score < -thres){
          cnv_cond[i,] = -1
        }
        else{
          cnv_cond[i,] = 0
        }
      }
      else{
        if (cnv_score < thres){
          cnv_cond[i,] = 1
        }
        else if (cnv_score > -thres){
          cnv_cond[i,] = -1
        }
        else{
          cnv_cond[i,] = 0
        }
      }
    }
    
    else if (is.null(arm)){
      subset <- subset(segfile, Sample==i & Chromosome==chr)
      subset$segment_lengths <- as.numeric(subset$End) - as.numeric(subset$Start)
      cnv_score <- sum(as.numeric(subset$segment_lengths) * as.numeric(subset$Segment_Mean))/sum(subset$segment_length)
      
      if (is.nan(cnv_score)) next
      
      
      if (thres > 0){
        if (cnv_score > thres){
          cnv_cond[i,] = 1
        }
        else if (cnv_score < -thres){
          cnv_cond[i,] = -1
        }
        else{
          cnv_cond[i,] = 0
        }
      }
      else{
        if (cnv_score < thres){
          cnv_cond[i,] = 1
        }
        else if (cnv_score > -thres){
          cnv_cond[i,] = -1
        }
        else{
          cnv_cond[i,] = 0
        }
      }
    }
  }
  
  cnv_cond_df <- cnv_cond
  if (thres < 0){
    write.table(cnv_cond_df, file=paste("chr",chr,arm,"loss",".txt"), sep="\t")
  }
  else{
    write.table(cnv_cond_df, file=paste("chr",chr,arm,"gain",".txt"), sep="\t")
  }
  
  ## output the groupings
  row.names(cnv_cond) <- lapply(row.names(cnv_cond), cleanUp_SampleNames)
  cnv_cond <- subset(cnv_cond, cnv_cond$CNV > -Inf)
  return (cnv_cond)
}


DGEA <- function(rna_file, cnv_group, cnv_type){
  rna_data <- read.table(rna_data, header=T)
  # preprocessing
  # remove duplicate rows
  rna_data <- rna_data[!duplicated(rna_data$GeneSymbol),]
  row.names(rna_data) <- rna_data$GeneSymbol
  rna_data <- rna_data[,-1]
  
  # group the samples by cnv conditions
  colnames(rna_data) <- lapply(colnames(rna_data), cleanUpSampleNames)
  # remove the duplicated patients (but keep track of them, not done yet)
  # potentially useful for tracking the patient's prognosis
  rna_data <- rna_data[, !duplicated(colnames(rna_data))]
  
  # reorder the colnames
  overlap_samples <- comparelists(row.names(cnv_group), colnames(rna_data))$intersect
  rna_data_ <- rna_data[overlap_samples]
  # rna_data <- rna_data[, match(names(cnv_group), colnames(rna_data))]
  
  colData <- subset(cnv_group, rownames(cnv_group) %in% overlap_samples)
  
  rna_data <- round(rna_data)
  head(rna_data)
  
  dds <- DESeqDataSetFromMatrix(count=rna_data, colData=colData, design = ~CNV)
  dds <- dds[rowSums(counts(dds)) > 1,]
  ddsSF <- DESeq(dds)
  
  res_cnv <- results(ddsSF, constrast = c("CNV", 1, 0))
  res_cnv <- res_cnv[order(res_cnv$padj),]
  res_cnv <- na.omit(res_cnv)
  write.table(res_cnv, file=paste0("res_cnv", cnv_type), sep="\t")
  
  return(rna_data)
  
}


cnv_group <- groupByCNV(segFile, 8, "p", thres=-0.3, arm_start=2E7, arm_end=4.5E7)
rna_file <- "BRCA_genes_results_processed_raw_counts.txt"
rna_data <- read.table(rna_file, header=T)
# preprocessing
# remove duplicate rows
rna_data <- rna_data[!duplicated(rna_data$GeneSymbol),]
row.names(rna_data) <- rna_data$GeneSymbol
rna_data <- rna_data[,-1]

# group the samples by cnv conditions
colnames(rna_data) <- lapply(colnames(rna_data), cleanUp.SampleNames)
# remove the duplicated patients (but keep track of them, not done yet)
# potentially useful for tracking the patient's prognosis
rna_data <- rna_data[, !duplicated(colnames(rna_data))]

# reorder the colnames
overlap_samples <- comparelists(row.names(cnv_group), colnames(rna_data))$intersect
rna_data <- rna_data[,overlap_samples]
# colnames(rna_data) <- reorder.factor(colnames(rna_data), new.order=overlap_samples)
# rna_data <- rna_data[, match(names(cnv_group), colnames(rna_data))]

colData <- subset(cnv_group, rownames(cnv_group) %in% overlap_samples)

rna_data <- round(rna_data)
head(rna_data)

dds <- DESeqDataSetFromMatrix(count=rna_data, colData=colData, design = ~CNV)
dds <- dds[rowSums(counts(dds)) > 1,]
ddsSF <- DESeq(dds)

res_cnv <- results(ddsSF, constrast = c("CNV", 1, 0))
res_cnv <- res_cnv[order(res_cnv$padj),]
res_cnv <- na.omit(res_cnv)
write.table(res_cnv, file=paste0("res_cnv", cnv_type), sep="\t")



