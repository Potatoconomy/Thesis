######################################################################################
# Libraries -----------------------------------------------------------------------------------
rm(list=ls())
load_libraries <- function(){
  library(gplots)
  library(viridis)
  library(ggplot2)
  library(reshape2)
  library(biomaRt)
  library(matrixStats)
  library(DESeq2)
  library(RColorBrewer)
  library(pheatmap)
  library(ggpubr)
  library(edgeR)
  library(ashr)
  library(apeglm)
  library(tidyverse)
  library(ggrepel)
  library(rrcov)
  library(dendextend)
  library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
}
load_libraries()

######################################################################################
# Load in data ---------------------------------------------------------------------
load_data <- function(main_path) {
  f  <- paste(path, '.csv', sep='')
  df <- read.table(file=f, header=TRUE, sep=',')
  df <- df[order(row.names(df)),]
  return (df)
}

path <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/counts_Patrick'

human_1 <- load_data(path)
rownames(human_1) <- human_1$gene_name

drops <- c("X","gene_name")
human_1 <- human_1[ , !(names(human_1) %in% drops)]

#Gene Removal List
removal_list <- function(df){
  #Gene Removal List, mitochondrial and ribosomal
  remove_1 <- grep('^RPS', rownames(df))
  remove_2 <- grep('^MT-', rownames(df))
  remove_3 <- grep('^RPL', rownames(df))
  remove <- c(remove_1, remove_2, remove_3)
  if (length(remove) > 0) {df <- df[-remove ,]}
  return(df)
}

human_1 <- removal_list(human_1)
human_1 <- human_1 %>% drop_na()

names <- c('Healthy_Control_1','Healthy_Control_2','Healthy_Control_3',
           'Healthy_Heated_1', 'Healthy_Heated_2', 'Healthy_Heated_3',
           'Tumor_Control_1',  'Tumor_Control_2',  'Tumor_Control_3',
           'Tumor_Heated_1',   'Tumor_Heated_2',   'Tumor_Heated_3')

human_1 <- human_1[, names]

rm(list=c('path','load_data','removal_list','drops'))
######################################################################################
# Drop low count genes
human_1 <- human_1[rowSums(human_1) >= 10 ,]

col_data <- data.frame('Samples' = names)
#col_data$Replicates <-
col_data$Tissue    <- factor(c(rep( 'Healthy',6), rep('Tumor',6)))
col_data$Condition <- factor(c(rep( c(rep('Ctrl',3),rep('Heat',3)), 2))) # ccc hhh ccc hhh
col_data$Batch     <- factor(c(rep( '1', 12)))

dds1 <- DESeqDataSetFromMatrix(human_1,
                              colData = col_data,
                              design = ~ Condition + Tissue)

dds1 <- DESeq(dds1)
vst1 <- vst(dds1)

dds2 <- DESeqDataSetFromMatrix(human_1,
                               colData = col_data,
                               design = ~ Condition)

dds2 <- DESeq(dds2)
vst2 <- vst(dds2)

#res <- results(dds, contrast=c("Condition","Heat","Ctrl")) #, alpha=0.1, contrast=c('Condition','Heat','Ctrl'))
#dds <- estimateSizeFactors(dds)
#counts <- counts(dds, normalized=TRUE)

p <- plotPCA(vst1, intgroup=c("Tissue", "Condition"))
ggsave(plot = p, filename = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/PCA/Both_PCA.svg')
#plotPCA(vst2, intgroup=c("Tissue", "Condition"))





