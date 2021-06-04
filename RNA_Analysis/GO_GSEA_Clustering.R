# Patrick Campbell
# GO GSEA Clustering Script
# März 2021

rm(list=ls())
################################################
##~ User Paths ~################################
################################################
df_dir  <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/Dataframes/'
path_dir <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/Figures_Patrick/'
fig_save_path <- '/home/patrick/Desktop/Masters Thesis Figures/Unedited/GOSemSim/'
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
    df$rank <- sign(df$log2FoldChange) * df$pvalue  
  } else if (method == 'simple') {
    df$rank <- df$log2FoldChange
  }
  selected <- na.omit(df)
  selected <- selected[order(selected$log2FoldChange),]
  selected <- selected[, c("rank")]
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

get_gseGO <- function(gene_list, GO_type, keytype, pcutoff) {
  # Citation: ClusterProfiler
  # Usage:
  #   gene_list -- output from function: get_geneList
  #   GO_type   -- 'MF', 'BP', or 'CC'
  #   keytype   -- 'SYMBOL' or 'ENTREZID'  
  #   pcutoff   -- 0.05
  org <- "org.Hs.eg.db"
  ont <- GO_type
  hs <- godata(org, ont = ont)
  gse <- gseGO(gene_list,
               OrgDb = org,
               keyType = keytype,
               ont = ont,
               pvalueCutoff = pcutoff,
               minGSSize = 1)
  
  # Split into positive and negative
  gse.POS <- gse[sign(gse$enrichmentScore) == 1  , ]
  gse.NEG <- gse[sign(gse$enrichmentScore) == -1 , ]
  go_ids.pos <- rownames(gse.POS[gse.POS$p.adjust < 0.1 , ])
  go_ids.neg <- rownames(gse.NEG[gse.NEG$p.adjust < 0.1 , ])
  
  return(list(hs, gse.POS, gse.NEG, go_ids.pos, go_ids.neg))
}

get_clusters <- function(hsGO, GOIDs, gseGO, fig_name, qcutoff){
  # Cluster GO IDs into heatmap
  # Usage:
  #   hsGO: hs output from function: get_gseGO
  #   GOIDs: go_ids.___ output from function: get_gseGO
  #   gseGO: gse.___ output from function: get_gseGO
  #   fig_name: string
  #   qcutoff: quantile cutoff for heatmap.
  
  i <- 0 
  while (i==0) {
    GOSIM <- as.data.frame(mgoSim(GOIDs, GOIDs, semData=hsGO, measure="Wang", combine=NULL))
    GOSIM.cutoff <- as.integer(quantile(rowSums(GOSIM), qcutoff))  # update this cutoff later
    GOSIM <- GOSIM[rowSums(GOSIM)>=GOSIM.cutoff, colSums(GOSIM)>=GOSIM.cutoff]  # Basic filtering 
    GOSIM.OUT <- pheatmap(GOSIM)
    
    
    height <- 0.6 * max(GOSIM.OUT$tree_col[[2]])
    print(paste('quantile cutoff', qcutoff))
    cat("Enter a new height if not satisfied with number of GOIDS on plot:\nPress [enter] to continue")
    qcutoff <- readline()
    if (qcutoff == "") {i <- 1} else {qcutoff <- as.numeric(qcutoff)}
    dev.off()
  }
  
  print(paste('height', height))
  i <- 0
  while (i == 0){
    GOSIM.OUT <- pheatmap(GOSIM)
    dev.off()
    GOSIM.sorted <- sort(cutree(GOSIM.OUT$tree_row, h=height))
    # Get clusters into a gse results data frame, gse.BP.df
    GOSIM.DF <- as.data.frame(GOSIM.sorted)
    GOSIM.DF$ID <- rownames(GOSIM.DF)
    a <- gseGO[gseGO$ID %in% rownames(GOSIM.DF) ,]
    a <- merge(a, GOSIM.DF, by.x = "ID", by.y = "ID")
    GOSIM.DF <- a
    col_annos <- data.frame(GOSIM.DF$GOSIM.sorted)
    colnames(col_annos) <- c("Cluster")
    rownames(col_annos) <- GOSIM.DF$ID
    
    q <- pheatmap(GOSIM, 
                  annotation_col = col_annos, 
                  show_rownames = FALSE, 
                  show_colnames = FALSE) #, annotation_row = col_annos)
    
    cat("Enter a new height if not satisfied with clusters:\nPress [enter] to continue")
    height <- readline()
    if (height == "") { i <- 1} else {height <- as.numeric(height)}
  }
  
  
  # SLIGHT PROBLEM HERE.  THIS CAN REMOVE ENTIRE CLUSTERS WITHOUT LE GENES
  #GOSIM.DF <- GOSIM.DF[grep("/", GOSIM.DF$core_enrichment) , ]  # KEEP ONLY DATASETS WITH LEADING EDGE
  write.csv2(GOSIM.DF, file = paste(fig_save_path, fig_name, ".csv", sep = ""))
  f <- paste(fig_save_path, fig_name, '.svg', sep='')
  save_pheatmap_svg(q, f)
  return(GOSIM.DF)
}

master_function <- function(DESeqName, save_dir){
  # Get DESeq object and ranked list.
  S4       <- loadRData(DESeqName)
  degs     <- get_rank(S4, 'simple')
  geneList <- get_geneList(degs)
  
  # RUN gseGO() from ClusterProfiler on all 3 GO Databanks
  c(hsGO.MF, gseGO.MF.POS, gseGO.MF.NEG, GOIDs.MF.POS, GOIDs.MF.NEG) %<-% get_gseGO(
    gene_list = geneList, 
    GO_type = "MF", 
    keytype  = "SYMBOL", 
    pcutoff = 1)
  
  c(hsGO.BP, gseGO.BP.POS, gseGO.BP.NEG, GOIDs.BP.POS, GOIDs.BP.NEG) %<-% get_gseGO(
    gene_list = geneList,
    GO_type = "BP",
    keytype  = "SYMBOL",
    pcutoff = 1)
  
  c(hsGO.CC, gseGO.CC.POS, gseGO.CC.NEG, GOIDs.CC.POS, GOIDs.CC.NEG) %<-% get_gseGO(
    gene_list = geneList,
    GO_type = "CC",
    keytype  = "SYMBOL",
    pcutoff = 1)
  
  # TODO: Automate qcutoff based on size of significant pathways
  # Cluster into heatmaps
  GOSIM.MF.POS.DF <- get_clusters(hsGO  = hsGO.MF,
                                  GOIDs = GOIDs.MF.POS,
                                  gseGO = gseGO.MF.POS,
                                  fig_name = paste('GOSEMSIM_MF_POS_', DESeqName, sep=''),
                                  qcutoff = 0.)  
  dev.off()
  
  
  GOSIM.MF.NEG.DF <- get_clusters(hsGO  = hsGO.MF,
                                  GOIDs = GOIDs.MF.NEG,
                                  gseGO = gseGO.MF.NEG,
                                  fig_name = paste('GOSEMSIM_MF_NEG_', DESeqName, sep=''),
                                  qcutoff = 0.)  
  dev.off()
  
  GOSIM.BP.POS.DF <- get_clusters(hsGO  = hsGO.BP,
                                  GOIDs = GOIDs.BP.POS,
                                  gseGO = gseGO.BP.POS,
                                  fig_name = paste('GOSEMSIM_BP_POS_', DESeqName, sep=''),
                                  qcutoff = 0.)  
  dev.off()
  
  
  GOSIM.BP.NEG.DF <- get_clusters(hsGO  = hsGO.BP,
                                  GOIDs = GOIDs.BP.NEG,
                                  gseGO = gseGO.BP.NEG,
                                  fig_name = paste('GOSEMSIM_BP_NEG_', DESeqName, sep=''),
                                  qcutoff = 0.) 
  dev.off()
  
  GOSIM.CC.POS.DF <- get_clusters(hsGO  = hsGO.CC,
                                  GOIDs = GOIDs.CC.POS,
                                  gseGO = gseGO.CC.POS,
                                  fig_name = paste('GOSEMSIM_CC_POS_', DESeqName, sep=''),
                                  qcutoff = 0.)  
  dev.off()
  
  
  GOSIM.CC.NEG.DF <- get_clusters(hsGO  = hsGO.CC,
                                  GOIDs = GOIDs.CC.NEG,
                                  gseGO = gseGO.CC.NEG,
                                  fig_name = paste('GOSEMSIM_CC_NEG_', DESeqName, sep=''),
                                  qcutoff = 0.)  
  dev.off()
}
################################################
##~ Call Functions ~############################
################################################
LFC_HEALTHY   <- master_function('lfc_healthy')
LFC_TUMOR     <- master_function('lfc_tumor')
LFC_NOTHEATED <- master_function('lfc_notheated')
LFC_HEATED    <- master_function('lfc_heated')



