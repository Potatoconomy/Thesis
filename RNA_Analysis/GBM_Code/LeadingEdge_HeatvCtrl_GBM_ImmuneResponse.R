# Patrick Campbell
# Leading Edge GSEA Script
# März 2021

rm(list=ls())
################################################
##~ User Paths ~################################
################################################
df_dir  <- '/home/patrick/Desktop/Figures_Human_GBM/'
path_dir <- '/home/patrick/Desktop/Figures_Human_GBM/'
fig_save_path <- '/home/patrick/Desktop/Figures_Human_GBM/LeadingEdge/'
gene_pathways <- '/home/patrick/Desktop/Gene_Pathways/allgenesets.RDS'
################################################
##~ Imports ~###################################
################################################
load_libraries <- function(){
  library(msigdbr)
  library(zeallot)
  library(clusterProfiler)
  library(fgsea)
  library(org.Hs.eg.db)
  library(DOSE)
  library(DESeq2)
  library(tibble)
  library(RColorBrewer)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(GOSemSim)
  library(pheatmap)
  library(enrichplot)
  library(ggplot2)
}

load_libraries()
################################################
##~ Functions ~#################################
################################################

save_pheatmap_svg <- function(x, filename, width=20, height=20) {
  # Save better svg figures
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  svg(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

loadRData <- function(DEObject_name){
  #loads an RData file, and returns it
  # DEObject_name: name of DEObject that was saved in your directory, ie: 'lfc'
  file_name <- paste(path_dir, DEObject_name, sep='')
  load(file_name)
  DEObject <- get(ls()[!(ls() %in% c("file_name","DEObject_name"))])
  return(DEObject)
}


get_rank <- function(df, method){
  # Input is an LFC DESeq2 dataframe, function returns ranked gene list for GSEA
  if (method == 'fancy'){
    df$rank <- sign(df$log2FoldChange) * -log10(df$pvalue)  
  } else if (method == 'simple') {
    df$rank <- df$log2FoldChange
  }
  selected <- na.omit(df)
  selected <- selected[order(selected$log2FoldChange),]
  selected <- selected[, c("rank")]
  print(selected)
  return(data.frame(selected))
}

get_rank2<- function(df, method){
  # Input is an LFC DESeq2 dataframe, function returns ranked gene list for GSEA
  if (method == 'fancy'){
    df$rank <- sign(df$log2FoldChange) * -log10(df$pvalue)  
  } else if (method == 'simple') {
    df$rank <- df$log2FoldChange
  }
  selected <- na.omit(df)
  selected <- selected[order(selected$log2FoldChange),]
  selected <- selected[, c("rank")]
  print(selected)
  return(data.frame(selected))
}

get_geneList <- function(df){
  # Assume input straight from get_degs()
  l <- df[,1]
  names(l) <- rownames(df)
  l <- sort(l, decreasing = TRUE)
  return(l)
}

get_genes <- function(df, cutoff){
  # Assume input straight from get_geneList
  l <- names(df)[abs(geneList)>=cutoff]
  if (length(l) < 20) {
    print("fewer than 20 genes with this cutoff value")
    print(length(l))
  }
  return(l)
}

get_pathways <- function(){
  pathways <- readRDS(gene_pathways, character())
  pathways_list <- split(x = pathways$gene, f = pathways$ont)  # set up for fgsea
  pathways_list_ <- pathways_list[c('BP.GO_REGULATION_OF_INFLAMMATORY_RESPONSE', 
                                    'BP.GO_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                                    'HALLMARK_HYPOXIA',
                                    'BP.GO_REGULATION_OF_B_CELL_ACTIVATION')]
  return(pathways_list_)
}

get_annotation_scale <- function(LFC){
  anno_max <- ceiling(max(LFC$log2FoldChange))
  anno_min <- floor(min(LFC$log2FoldChange))
  return(list(anno_max, anno_min))
}

preparation_function1 <- function(DESeqName){
  S4       <- loadRData(DESeqName)
  degs     <- get_rank(S4, 'fancy')
  geneList <- get_geneList(degs)
  pathways_list <- get_pathways()
  c(anno_max, anno_min) %<-% get_annotation_scale(S4)
  
  fgseaRes <- fgsea(pathways = pathways_list, 
                    stats    = geneList,
                    minSize  = 15,
                    maxSize  = 250,
                    nperm=10000)
  
  #fgseaRes <- fgseaRes[fgseaRes$padj <= 0.05 ,]
  
  return(list(S4, geneList, pathways_list, fgseaRes, anno_max, anno_min))
}

interesting_pathways <- function(FGSEA) {
  # Here are some interesting pathways selected for our own data. The negative pathways were the primary negative pathways
  # in our dataset. They should be changed to match the primary negative pathways in each individual dataset.

  negative_pathways <- c('BP.GO_REGULATION_OF_INFLAMMATORY_RESPONSE', 
                         'BP.GO_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                         'HALLMARK_HYPOXIA',
                         'BP.GO_REGULATION_OF_B_CELL_ACTIVATION')
  
  combined_pathways <- c(negative_pathways)
  

  fgsea_combined <- FGSEA[FGSEA$pathway %in% combined_pathways]
  
  return(fgsea_combined)
}

make_heatmap <- function(DF, LFC, fname, anno_min=-2, anno_max=2) {
  LFC.SIG <- na.omit(LFC)
  LFC.SIG <- LFC.SIG[LFC.SIG$padj <= 0.1 , ]
  LFC.SIG <- na.omit(LFC.SIG)
  
  DF.SIG <- DF[DF$padj <= 0.15 , ]

  p_vals.pathways <- data.frame(DF.SIG$padj)
  nes <- data.frame(DF.SIG$NES)
  rownames(nes) <- DF.SIG$pathway
  rownames(p_vals.pathways) <- DF.SIG$pathway
  l <- list()
  
  for (i in 1:nrow(DF.SIG)){
    for (j in DF.SIG[i, leadingEdge])
      for (k in j){
        l <- append(x = l, values = k)
      }
  }
  
  l <- unique(l)
  l <- as.character(l)
  l <- sort(l)
  
  print(rownames(LFC.SIG))
  print(l)
  print(rownames(LFC.SIG) %in% l)
  
  LFC.SIG <- LFC.SIG[rownames(LFC.SIG) %in% l , ]
  
  LFC.SIG <- na.omit(LFC.SIG)
  p_vals.genes <- LFC.SIG[, 'padj']
  names(p_vals.genes) <- rownames(LFC.SIG)
  p_vals.genes <- as.data.frame(p_vals.genes)
  X <- data.frame(matrix(ncol = nrow(LFC.SIG), nrow = nrow(DF.SIG)))
  colnames(X) <- rownames(LFC.SIG)
  rownames(X) <- DF.SIG$pathway

  
  for (i in 1:nrow(X)){
    for (gene in colnames(X)){
      if (length(grep(gene, DF.SIG[i, leadingEdge], value = F)) != 0) {
        # if (sign(LFC.SIG[gene, "log2FoldChange"]) == 1) {a <- 1}
        # else if (sign(LFC.SIG[gene, "log2FoldChange"]) == -1) {a <- -1}
        X[i, gene] <- LFC.SIG[gene, "log2FoldChange"]
        if (X[i, gene] >=  2) {X[i, gene] <-  2}
        if (X[i, gene] <= -2) {X[i, gene] <- -2}
      }
    }
  }
  
  X_ <- X
  X_[is.na(X)] <- 0
  X[nrow(X)+1 , ] <- NA
  nes[nrow(nes)+1 ,] <- 0
  
  #a_ <- max(anno_min, anno_max)
  mat_breaks = seq(-2, 2, 0.1)
  
  #nes=round(nes, digits=1)
  #color <- rev(magma(length(mat_breaks)))
  #color <- brewer.pal(n=length(mat_breaks), name="RdBu")
  primary_colors <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdYlBu")))(length(mat_breaks))
  #anno_colors <- colorRampPalette(rev(brewer.pal(n = 9, name = "PiYG")))(length(mat_breaks))
  
  # X colnames need to be ordered
  # p_vals.genes needs to follow the X colnames order
  x_order <- sort(colnames(X))
  X <- X[, x_order]
  p_vals.genes$genes <- rownames(p_vals.genes)
  p_vals.genes <- as.data.frame(p_vals.genes[colnames(X) , "p_vals.genes"])
  rownames(p_vals.genes) <- colnames(X)

  p <- pheatmap(X,
                cluster_rows = FALSE, cluster_cols = FALSE,
                na_col = 'white',
                color = primary_colors,
                annotation_row = nes,
                #annotation_col = p_vals.genes,
                cellwidth = 10, cellheight = 10,
                breaks = mat_breaks)
  
  save_pheatmap_svg(p, paste(fig_save_path, fname, "_pheatmap.svg",sep=""))
  dev.off()
  
  # Make barplots to save onto the inkscape fig
  p_vals.genes$genes <- rownames(p_vals.genes)
  colnames(p_vals.genes) <- c("pvalue", "genes")
  q <- ggplot(data=p_vals.genes, aes(x=genes, y=pvalue)) +
    geom_bar(stat="identity") +
    theme_bw() + 
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) 
  q
  ggsave(filename = paste(fig_save_path, "genes_pval_bars_", fname, '.svg', sep=""), plot = q)
  #dev.off()
  
  #### Now do the same with PATHWAYS TODO!!!!
  p_vals.pathways$pathways <- rownames(p_vals.pathways)
  colnames(p_vals.pathways)
  p <- ggplot(data=p_vals.pathways, aes(x=pathways, y=DF.SIG$padj)) +
    geom_bar(stat="identity") +
    theme_bw() + 
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) 
  p
  ggsave(filename = paste(fig_save_path, "pathways_pval_bars_", fname, '.svg', sep=""), plot = p)
  
  return(list(X, nes, p_vals.pathways, p_vals.genes))
}
make_many(fgsea_combined)

make_many <- function(fgsea_combined){
  # Calls function: make_heatmap for all pathways
  anno_min = -2
  anno_max = 2

  c(X, nes, p_vals.pathways, p_vals.genes) %<-% make_heatmap(DF = fgsea_combined,
                                                             LFC = DESeqObject,
                                                             fname = 'combined',
                                                             anno_min = anno_min,
                                                             anno_max = anno_max)
}
################################################
##~ Call Functions ~############################
################################################
# Choose which LFC object to work with
c(DESeqObject, geneList, pathways_list, fgseaRes, anno_max, anno_min) %<-% preparation_function1('LFC')
fgsea_combined %<-% interesting_pathways(fgseaRes)
make_many(fgsea_combined)


a <- get_rank2(DESeqObject, 'fancy')
a$gene_names <- rownames(a)
a$rank <- a[1]
a <- a[c('gene_names', 'rank')]
path_genes <- unlist(pathways_list)
aa <- a[path_genes ,]
aa <- na.omit(aa)
aa <- unique(aa)
aa <- aa[order(aa$rank),]




