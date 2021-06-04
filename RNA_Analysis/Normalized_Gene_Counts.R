# Normalized Gene Counts
# Patrick Campbell
# März 2021


rm(list=ls())
################################################
##~ User Paths ~################################
################################################
counts_file <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/counts_Patrick.csv'
fig_save_path <-'/home/patrick/Desktop/Masters Thesis Figures/Gene_Count_Figs/'
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
  library(edgeR)
  library(RColorBrewer)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(reshape2)
  library(circlize)
  library(GOSemSim)
  library(pheatmap)
  library(enrichplot)
  library(ggplot2)
  library(conflicted)
  conflict_prefer("select", "dplyr")
}

load_libraries()


################################################
##~ Functions ~#################################
################################################
load_data <- function(main_path) {
  # load in a counts table
  # columns: NULL gene_name sample1 sample2...
  
  # f  <- paste(path, '.csv', sep='')
  f <- main_path
  df <- read.table(file=f, header=TRUE, sep=',')
  df <- df[order(row.names(df)),]
  rownames(df) <- df$gene_name
  drops <- c("X","gene_name") # X is column name of index that R gives to unnamed index column
  df <- df[ , !(names(df) %in% drops)]
  return (df)
}

removal_list <- function(df){
  #Gene Removal List, mitochondrial and ribosomal
  remove_1 <- grep('^RPS',   rownames(df))
  remove_2 <- grep('^MT-',   rownames(df))
  remove_3 <- grep('^RPL',   rownames(df))
  remove_4 <- grep('MALAT1', rownames(df))
  remove <- c(remove_1, remove_2, remove_3, remove_4)
  
  if (length(remove) > 0) {df <- df[-remove ,]}
  df <- df %>% drop_na()
  return(df)
}

################################################
##~ Call Functions ~############################
################################################

COUNTS <- load_data(counts_file)

COUNTS <- removal_list(COUNTS)


################################################
##~ Experiment Specific Info ~##################
################################################
# COUNTS has multiple different conditions inside it. This should be changed to match individual experiments
TumHeat <- select(COUNTS, 'Tumor_Heated_1', 'Tumor_Heated_2', 'Tumor_Heated_3') 
TumCtrl <- select(COUNTS, 'Tumor_Control_1', 'Tumor_Control_2', 'Tumor_Control_3')
HltHeat <- select(COUNTS, 'Healthy_Heated_1', 'Healthy_Heated_2', 'Healthy_Heated_3')
HltCtrl <- select(COUNTS, 'Healthy_Control_1', 'Healthy_Control_2', 'Healthy_Control_3')

Tumor <- select(COUNTS, 
                'Tumor_Heated_1', 
                'Tumor_Heated_2', 
                'Tumor_Heated_3',
                'Tumor_Control_1', 
                'Tumor_Control_2', 
                'Tumor_Control_3')

Healthy <- select(COUNTS, 
                  'Healthy_Heated_1', 
                  'Healthy_Heated_2', 
                  'Healthy_Heated_3', 
                  'Healthy_Control_1', 
                  'Healthy_Control_2', 
                  'Healthy_Control_3')

Heated <-  select(COUNTS, 
                  'Tumor_Heated_1', 
                  'Tumor_Heated_2', 
                  'Tumor_Heated_3',
                  'Healthy_Heated_1', 
                  'Healthy_Heated_2', 
                  'Healthy_Heated_3')

NotHeated <- select(COUNTS,
                    'Tumor_Control_1', 
                    'Tumor_Control_2', 
                    'Tumor_Control_3', 
                    'Healthy_Control_1', 
                    'Healthy_Control_2', 
                    'Healthy_Control_3')

################################################
##~ Library Size ~##############################
################################################
group_heatctrl <- c('Heat','Heat','Heat','Ctrl','Ctrl','Ctrl')  # experiment always comes first

DGE_Healthy <-   DGEList(counts=Healthy,   group=group_heatctrl)
DGE_Tumor   <-   DGEList(counts=Tumor,     group=group_heatctrl)

keepH         <- filterByExpr(DGE_Healthy, group=group_heatctrl)
keepT         <- filterByExpr(DGE_Tumor, group=group_heatctrl)

DGE_Healthy   <- DGE_Healthy[keepH, keep.lib.sizes=FALSE]
DGE_Tumor     <- DGE_Tumor[keepT, keep.lib.sizes=FALSE]


##
group_tumorhealthy <- c('Tumor','Tumor','Tumor','Healthy','Healthy','Healthy')

DGE_Heated <-    DGEList(counts=Heated,    group=group_tumorhealthy)
DGE_NotHeated <- DGEList(counts=NotHeated, group=group_tumorhealthy)

keepHeated    <- filterByExpr(DGE_Heated, group=group_tumorhealthy)
keepNotHeated <- filterByExpr(DGE_NotHeated, group=group_tumorhealthy)

DGE_Heated    <- DGE_Heated[keepHeated, keep.lib.sizes=FALSE]
DGE_NotHeated <- DGE_NotHeated[keepNotHeated, keep.lib.sizes=FALSE]


rm(list=c('Healthy', 'Heated', 'Tumor', 'NotHeated', 'TumCtrl','HltHeat',
          'TumHeat','HltCtrl','keepH','keepT','keepHeated','keepNotHeated'))

################################################
##~ DESeq2 ~####################################
################################################
# Column Data and Conditions
colData_healthy   <- data.frame(names     = colnames(DGE_Healthy$counts),  
                                condition = DGE_Healthy$samples$group)

colData_tumor     <- data.frame(names     = colnames(DGE_Tumor$counts),    
                                condition = DGE_Tumor$samples$group)

colData_heated    <- data.frame(names     = colnames(DGE_Heated$counts),   
                                condition = DGE_Heated$samples$group)

colData_notheated <- data.frame(names     = colnames(DGE_NotHeated$counts),
                                condition = DGE_NotHeated$samples$group)

# Create DDS for each group
dds_healthy <- DESeqDataSetFromMatrix(countData   = as.matrix(DGE_Healthy$counts),
                                      colData     = colData_healthy,
                                      design      = ~ condition)

dds_tumor <- DESeqDataSetFromMatrix(countData     = as.matrix(DGE_Tumor$counts),
                                    colData       = colData_tumor,
                                    design        = ~ condition)

dds_heated <- DESeqDataSetFromMatrix(countData    = as.matrix(DGE_Heated$counts),
                                     colData      = colData_heated,
                                     design       = ~ condition)

dds_notheated <- DESeqDataSetFromMatrix(countData = as.matrix(DGE_NotHeated$counts),
                                        colData   = colData_notheated,
                                        design    = ~ condition)



dds_healthy   <- DESeq(dds_healthy)
dds_tumor     <- DESeq(dds_tumor)
dds_heated    <- DESeq(dds_heated)
dds_notheated <- DESeq(dds_notheated)

# Results
res_healthy   <- results(dds_healthy,   alpha=0.1, contrast = c('condition', 'Heat',  'Ctrl'))
res_tumor     <- results(dds_tumor,     alpha=0.1, contrast = c('condition', 'Heat',  'Ctrl'))
res_heated    <- results(dds_heated,    alpha=0.1, contrast = c('condition', 'Tumor', 'Healthy'))
res_notheated <- results(dds_notheated, alpha=0.1, contrast = c('condition', 'Tumor', 'Healthy'))

save(dds_tumor, file=paste('/media/patrick/Elements/Hyperthermia/CookedHuman1/DESeq2/','dds_tumor_1',sep=''))
save(dds_healthy, file=paste('/media/patrick/Elements/Hyperthermia/CookedHuman1/DESeq2/','dds_healthy_1',sep=''))
save(dds_heated, file=paste('/media/patrick/Elements/Hyperthermia/CookedHuman1/DESeq2/','dds_heated_1',sep=''))
save(dds_notheated, file=paste('/media/patrick/Elements/Hyperthermia/CookedHuman1/DESeq2/','dds_notheated_1',sep=''))


################################################
##~ Normalized Count Plotting ~#################
################################################
plotGeneCounts <- function(dds, genes, fig_name){
  # Plot genes boxplots for a vector of gene names
  gene <- intersect(genes, rownames(dds)) # Prevent errors
  difff <- setdiff(genes, gene)
  
  if (length(difff) > 0) {cat('Following genes were not measured:', difff, '\n')}
  
  if (length(gene) == 1){
    df <- plotCounts(dds, gene, intgroup = 'condition', returnData = TRUE)
    p <- ggplot(df, aes(x=condition, y=count)) + 
      geom_point(position=position_jitter(w=0.1,h=0)) + 
      scale_y_log10(breaks=c(25,100,400))
    print(p)
  }
  
  else if (length(gene) > 1){
    df <- data.frame()
    l <- list()
    j = 1
    
    for (i in gene){
      data <- plotCounts(dds, i, intgroup = 'condition', returnData = TRUE)
      colnames(data) <- c(i, 'condition')
      l[[j]] <- data
      j <- j + 1
    }
    df <- do.call(cbind, l)
    df.m <- melt(df, id.var = "condition")
    
    p <- ggplot(data = df.m, aes(x=variable, y=value)) + 
      geom_boxplot(aes(fill=condition)) + 
      ggtitle(fig_name)
    ggsave(filename = paste(fig_save_path, fig_name, '.svg', sep=''), plot = p)
  }
  
}


plotting_function <- function(HEATSHOCK, dds, comma = TRUE) {
  # Plots Heatshock Proteins that were measured in the experiment
  # This function probably fails if a category is of size 0. Sue me and ill fix it
  if (comma == TRUE) {heatshocks <- read.csv2(paste(main, HEATSHOCK, '.csv', sep=''), sep = ' ')}
  else {heatshocks <- read.csv(paste(main, HEATSHOCK, '.csv', sep=''), sep = ' ')}
  heatshocks %>% drop_na()  # need to test this line
  # need to turn factors into numeric vectors here... stupid fucking R
  
  rownames(heatshocks) <- heatshocks$X
  
  
  heatshocks <- heatshocks[heatshocks$padj <= 0.1 , ]
  chaperones <- 'Chaperones'
  small <- 'Small'
  kda40 <- '40kDA'
  kda70 <- '70kDA'
  kda90 <- '90kDA'
  
  kDA70 <- heatshocks[heatshocks$type == kda70 , ]
  kDA90 <- heatshocks[heatshocks$type == kda90 , ]
  kDA40 <- heatshocks[heatshocks$type == kda40 , ]
  Chaps <- heatshocks[heatshocks$type == chaperones , ]
  Small <- heatshocks[heatshocks$type == small , ]
  
  plotGeneCounts(dds, gene=rownames(kDA70), paste(HEATSHOCK, kda70))
  plotGeneCounts(dds, gene=rownames(kDA90), paste(HEATSHOCK, kda90))
  plotGeneCounts(dds, gene=rownames(kDA40), paste(HEATSHOCK, kda40))
  plotGeneCounts(dds, gene=rownames(Chaps), paste(HEATSHOCK, chaperones))
  plotGeneCounts(dds, gene=rownames(Small), paste(HEATSHOCK, small))
}

# Now i work with files from my HeatshockGenes.R file (On my github: Potatoconomy)
main <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/Figures_Patrick/HeatShockProteins/'
plotting_function(HEATSHOCK = 'lfc_healthyHEATSHOCKS', dds = dds_healthy, comma=TRUE)
plotting_function(HEATSHOCK = 'lfc_tumorHEATSHOCKS', dds = dds_tumor, comma=FALSE)
plotting_function(HEATSHOCK = 'lfc_notheatedHEATSHOCKS', dds = dds_notheated, comma=TRUE)
plotting_function(HEATSHOCK = 'lfc_heatedHEATSHOCKS', dds = dds_heated,comma=FALSE)



