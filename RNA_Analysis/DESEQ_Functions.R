rm(list=ls())
library(gplots)
library(ggplot2)
library(matrixStats)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(edgeR)
library(ashr)
library(apeglm)
library(tidyverse)
library(ggrepel)
library(rrcov)



# Load in data ---------------------------------------------------------------------
load_data <- function(main_path) {
  f  <- paste(main_path, '.csv', sep='')
  df <- read.table(file=f, header=TRUE,row.names=1,sep=',')
  df <- df[order(row.names(df)),]
  return (df)
}

#Gene Removal List
removal_list <- function(df){
  remove_1 <- startsWith(rownames(df),'MT-')
  df <- df[!remove_1,]
  remove_2 <- startsWith(rownames(df),'RPS')
  df <- df[!remove_2,]
}

# Plot dispersions---------------------------------------------------------------------
plot_disps <- function(file_name,data){
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  plotDispEsts(data, main=file_name)
  dev.off()
}

# MA Plots ##-------------------------------------------------------------
plot_MA <- function(file_name,data){
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  DESeq2::plotMA(data, main=file_name,ylim=c(-2,2))
  dev.off()
}


# PCA Plots ##--------------------------------------------------------------
plot_pca <- function(file_name, data){
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  DESeq2::plotPCA(data)   # dev.off() doesn't work in here for some reason
}

# Sample Distances -------------------------------------------------------------------------
plot_distancematrix <- function(file_name, data){
  sampleDists <- dist(t(assay(data)))
  sampleDistMatrix <- as.matrix(sampleDists)
  #rownames(sampleDistMatrix) <- data$condition
  colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col=colors)
  dev.off()
}

# Histograms -------------------------------------------------------------------------------------
plot_hists <- function(file_name, data){
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  hist(assay(data),breaks=20)
  dev.off()
}


# heatmap count matrix -------------------------------------------------------------------------------
plot_heatmap <- function(file_name, res_data, assay_data, n_genes){
  #res_data == lfc or res
  #assay_data == rlog or vst
  selected <- na.omit(res_data)
  selected <- subset(selected, padj<0.1)
  selected <- selected[order(selected$log2FoldChange),]
  
  selected_upreg <- tail(selected, n_genes)
  selected_dereg <- selected[1: n_genes,]
  selected <- c(rownames(selected_upreg), rownames(selected_dereg))
  df <- as.matrix(assay(assay_data))
  df <- df[selected,]
  df_z <- t(scale(t(df)))
  
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  pheatmap(df_z)
  dev.off()
}

# Volcano plots--------------------------------------------------------------------------
plot_volcano <- function(file_name, data, log_level, padj_level, xrange, ymax){
  data <- data[complete.cases(data),]
  df <- data.frame("log2FoldChange"=data$log2FoldChange,"padj"=data$padj)
  
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange > log_level & df$padj < padj_level] <- "UP"
  df$diffexpressed[df$log2FoldChange < -log_level & df$padj < padj_level] <- "DOWN"
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  df$dflabel <- NA
  
  df$dflabel[df$diffexpressed != "NO"] <- rownames(df)[df$diffexpressed != "NO"]
  
  p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed,
                           label=dflabel)) + 
    geom_point() + 
    theme_minimal() + 
    geom_text_repel() + 
    geom_vline(xintercept=c(-log_level, log_level), col="red") +
    geom_hline(yintercept=-log10(padj_level), col="red") +
    scale_color_manual(values=c("blue", "black", "red")) +
    coord_cartesian(xlim = c(-xrange, xrange), ylim = c(-1, ymax)) +
    ggtitle(file_name)
  
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  print(p)
}

plot_volcano_pvalue <- function(file_name, data, log_level, pvalue_level, xrange, ymax){
  data <- data[complete.cases(data),]
  df <- data.frame("log2FoldChange"=data$log2FoldChange,"pvalue"=data$pvalue)
  rownames(df) <- rownames(data)
  
  df$diffexpressed <- "NO"
  df$diffexpressed[df$log2FoldChange > log_level & df$pvalue < pvalue_level] <- "UP"
  df$diffexpressed[df$log2FoldChange < -log_level & df$pvalue < pvalue_level] <- "DOWN"
  
  mycolors <- c("blue", "red", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  df$dflabel <- NA
  
  df$dflabel[df$diffexpressed != "NO"] <- rownames(df)[df$diffexpressed != "NO"]
  
  p <- ggplot(data=df, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed,
                           label=dflabel)) + 
    geom_point() + 
    theme_minimal() + 
    geom_text_repel() + 
    geom_vline(xintercept=c(-log_level, log_level), col="red") +
    geom_hline(yintercept=-log10(pvalue_level), col="red") +
    scale_color_manual(values=c("blue", "black", "red")) +
    coord_cartesian(xlim = c(-xrange, xrange), ylim = c(-1, ymax)) +
    ggtitle(file_name)
  
  svg(paste(fig_save_path, file_name, '.svg', sep=""))#, 1000, 1000, pointsize=20)
  print(p)
}
#----------------------------------------------------------------------------------

