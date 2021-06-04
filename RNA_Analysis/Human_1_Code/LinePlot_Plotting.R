# TPM Counts

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

# Collaps Tech Reps
hc <- rowSums(select(human_1, names[1:3]))
ht <- rowSums(select(human_1, names[4:6]))
tc <- rowSums(select(human_1, names[7:9]))
tt <- rowSums(select(human_1, names[10:12]))

# Get per million scaling factor (msf)
# RPM normalizes for seq depth
hc_msf <- sum(hc)/1000000
ht_msf <- sum(ht)/1000000
tc_msf <- sum(tc)/1000000
tt_msf <- sum(tt)/1000000

# Divide the gene counts by the per million scaling factor
hc_RPM <- hc/hc_msf
ht_RPM <- ht/ht_msf
tc_RPM <- tc/tc_msf
tt_RPM <- tt/tt_msf

# Set up genes that we want to look for
genes        <- c('HSPA1A', 'HSPA1B', 'HSPH1', 'HSPA5', 'HSPA8', 'HSP90AA1', 'HSP90AB1')
gene_lengths <- c(2404,     2517,     27416,    6530,   5742,     59008,      7723)
gene_df <- data.frame('gene_id' = genes, 'lengths' = gene_lengths)

# Isolate these genes only
# Probably an error will appear if genes were removed earlier
hc_RPM_genes <- as.data.frame(hc_RPM[names(hc_RPM) %in% genes])
hc_RPM_genes$gene_id <- rownames(hc_RPM_genes)
colnames(hc_RPM_genes) <- c('counts', 'gene_id')

ht_RPM_genes <- as.data.frame(ht_RPM[names(ht_RPM) %in% genes])
ht_RPM_genes$gene_id <- rownames(ht_RPM_genes)
colnames(ht_RPM_genes) <- c('counts', 'gene_id')

tc_RPM_genes <- as.data.frame(tc_RPM[names(tc_RPM) %in% genes])
tc_RPM_genes$gene_id <- rownames(tc_RPM_genes)
colnames(tc_RPM_genes) <- c('counts', 'gene_id')

tt_RPM_genes <- as.data.frame(tt_RPM[names(tt_RPM) %in% genes])
tt_RPM_genes$gene_id <- rownames(tt_RPM_genes)
colnames(tt_RPM_genes) <- c('counts', 'gene_id')


# Merge with gene lengths df
hc <- merge(hc_RPM_genes, gene_df)
ht <- merge(ht_RPM_genes, gene_df)
tc <- merge(tc_RPM_genes, gene_df)
tt <- merge(tt_RPM_genes, gene_df)

# Divide counts by lengths to get the RPKM gene count
hc$RPKM_hc <- hc$counts / hc$lengths * 1000
ht$RPKM_ht <- ht$counts / ht$lengths * 1000
tc$RPKM_tc <- tc$counts / tc$lengths * 1000
tt$RPKM_tt <- tt$counts / tt$lengths * 1000


df <- data.frame('gene_id' = hc$gene_id)
df <- merge(df, hc[, c('gene_id', 'RPKM_hc')])
df <- merge(df, ht[, c('gene_id', 'RPKM_ht')])
df <- merge(df, tc[, c('gene_id', 'RPKM_tc')])
df <- merge(df, tt[, c('gene_id', 'RPKM_tt')])

# Log2 transform

# Add 1 to each so no negatives
df_ <- log2(df[, c('RPKM_hc', 'RPKM_ht', 'RPKM_tc', 'RPKM_tt')] + 1)
df[, c('RPKM_hc', 'RPKM_ht', 'RPKM_tc', 'RPKM_tt')] <- df_



# Get long format dataframes for healthy and tumor

df_healthy <- melt(df[, c('gene_id', 'RPKM_hc', 'RPKM_ht')])
df_tumor <- melt(df[, c('gene_id', 'RPKM_tc', 'RPKM_tt')])


# Plotting
# Limits
# Stack Overflow
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
min <- 0
max <- roundUpNice(max(df[, c('RPKM_tc', 'RPKM_tt', 'RPKM_hc', 'RPKM_ht')]))


# Healthy
data_ends_h <- df_healthy %>% filter(variable == "RPKM_ht")

lp_healthy <- ggplot(df_healthy, aes(x = variable , y = value, group=gene_id)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  theme(legend.position = "top") +
  geom_text_repel(aes(label = gene_id), 
                  data = data_ends_h,
                  fontface ="plain", 
                  color = "black",
                  size = 3) + 
  ylim(min, max)

lp_healthy
ggsave(filename ='/home/patrick/Desktop/Masters Thesis Figures/Gene_Count_Figs/LinePlots/healthy_log.svg',
       plot = lp_healthy)
# Tumor
data_ends_t <- df_tumor %>% filter(variable == "RPKM_tt")

lp_tumor <- ggplot(df_tumor, aes(x = variable , y = value, group=gene_id)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  theme(legend.position = "top") + 
  geom_text_repel(aes(label = gene_id), 
                  data = data_ends_t,
                  fontface ="plain", 
                  color = "black", 
                  size = 3) +
  ylim(min, max)

lp_tumor
ggsave(filename ='/home/patrick/Desktop/Masters Thesis Figures/Gene_Count_Figs/LinePlots/tumor_log.svg',
       plot = lp_tumor)





##########################################################################
### DESEQ2 Normalization ~################################################
##########################################################################
rm(list=ls())
healthy_file="/media/patrick/Elements/Hyperthermia/CookedHuman1/GSEA/healthy_counts_normed_PATRICK.csv"
tumor_file="/media/patrick/Elements/Hyperthermia/CookedHuman1/GSEA/tumor_counts_normed_PATRICK.csv"

healthy_cts <- read.csv2(healthy_file, sep = "\t")
tumor_cts <- read.csv2(tumor_file, sep = "\t") 

# log2( sample_mean + 1)
hc <- log2((rowMeans(healthy_cts[, c('Healthy_Control_1','Healthy_Control_2' ,'Healthy_Control_3')]))+1)
ht <- log2((rowMeans(healthy_cts[, c('Healthy_Heated_1' ,'Healthy_Heated_2'  ,'Healthy_Heated_3')]))+1)
tc <- log2((rowMeans(tumor_cts[,   c('Tumor_Control_1'  ,'Tumor_Control_2'   ,'Tumor_Control_3')]))+1)
tt <- log2((rowMeans(tumor_cts[,   c('Tumor_Heated_1'   ,'Tumor_Heated_2'    ,'Tumor_Heated_3')]))+1)



# Set up genes that we want to look for
genes <- c('HSPA1A', 'HSPA1B', 'HSPH1', 'HSPA5', 'HSPA8', 'HSP90AA1', 'HSP90AB1', 'DNAJB4', 'DNAJC15', 'CCT5', 'CCT8')
#gene_lengths <- c(2404,     2517,     27416,    6530,   5742,     59008,      7723)
#gene_df <- data.frame('gene_id' = genes, 'lengths' = gene_lengths)

# errors if genes not all in here
hc <- hc[names(hc) %in% genes]
ht <- ht[names(ht) %in% genes]
tc <- tc[names(tc) %in% genes]
tt <- tt[names(tt) %in% genes]

h_df <- data.frame('gene' = names(hc), 'hc' = hc, 'ht' = ht)  # could later add in what type of heat shock the gene is
t_df <- data.frame('gene' = names(tc), 'tc' = tc, 'tt' = tt)

# Get long format dataframes for healthy and tumor

h_df <- melt(h_df)
t_df <- melt(t_df)


# Plotting
# Limits
# Stack Overflow
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
  if(length(x) != 1) stop("'x' must be of length 1")
  10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}
min <- 0
max <- roundUpNice(max(h_df$value, t_df$value))


# Healthy
data_ends_h <- h_df %>% filter(variable == "ht")

lp_healthy <- ggplot(h_df, aes(x = variable , y = value, group=gene)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  theme(legend.position = "top") +
  geom_text_repel(aes(label = gene), 
                  data = data_ends_h,
                  fontface ="plain", 
                  color = "black",
                  size = 3) + 
  ylim(min, max)

lp_healthy
ggsave(filename ='/home/patrick/Desktop/Masters Thesis Figures/Gene_Count_Figs/LinePlots/healthy_log.svg',
       plot = lp_healthy)
# Tumor
data_ends_t <- t_df %>% filter(variable == "tt")

lp_tumor <- ggplot(t_df, aes(x = variable , y = value, group=gene)) +
  geom_line(aes(color = variable)) +
  geom_point() +
  theme(legend.position = "top") + 
  geom_text_repel(aes(label = gene), 
                  data = data_ends_t,
                  fontface ="plain", 
                  color = "black", 
                  size = 3) +
  ylim(min, max)

lp_tumor
ggsave(filename ='/home/patrick/Desktop/Masters Thesis Figures/Gene_Count_Figs/LinePlots/tumor_log.svg',
       plot = lp_tumor)










