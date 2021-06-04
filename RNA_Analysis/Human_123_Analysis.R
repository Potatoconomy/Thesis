#'DESEQ2'

# Cooked Humans


# Libraries -----------------------------------------------------------------------------------
rm(list=ls())
library(gplots)
library(viridis)
library(ggplot2)
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
source("/home/patrick/PatrickCode/RNA_Analysis/DESEQ_Functions.R")  # Set appropriate path here

fig_save_path <- '/media/patrick/Elements/Hyperthermia/FinalFigures/'


#######################################################################################################################
## Step 1: Organize DataFrame #########################################################################################
#######################################################################################################################

path_1 <- '/media/patrick/Elements/Hyperthermia/Cooked_Human_CountMatrixes/human_1_counts'
path_2 <- '/media/patrick/Elements/Hyperthermia/Cooked_Human_CountMatrixes/human_2_counts'
path_3 <- '/media/patrick/Elements/Hyperthermia/Cooked_Human_CountMatrixes/human_3_counts'

human_1 <- load_data(path_1)
human_2 <- load_data(path_2)
human_3 <- load_data(path_3)

rownames(human_1) <- human_1$gene_name
rownames(human_2) <- human_2$gene_name
rownames(human_3) <- human_3$gene_name

drops <- c("X","gene_name")
human_1 <- human_1[ , !(names(human_1) %in% drops)]
human_2 <- human_2[ , !(names(human_2) %in% drops)]
human_3 <- human_3[ , !(names(human_3) %in% drops)]

human_1 <- removal_list(human_1)
human_2 <- removal_list(human_2)
human_3 <- removal_list(human_3)

human_1 <- human_1 %>% drop_na()
human_2 <- human_2 %>% drop_na()
human_3 <- human_3 %>% drop_na()

human_1_Ctrl <- as.data.frame(rowSums(human_1[,c(1,2,3)]))
human_2_Ctrl <- as.data.frame(rowSums(human_2[,c(1,2,3)]))
human_3_Ctrl <- as.data.frame(rowSums(human_3[,c(1,2,3)]))

human_1_Heat <- as.data.frame(rowSums(human_1[,c(4,5,6)]))
human_2_Heat <- as.data.frame(rowSums(human_2[,c(4,5,6)]))
human_3_Heat <- as.data.frame(rowSums(human_3[,c(4,5,6)]))

names(human_1_Ctrl) <- c('Human_1_Ctrl')
names(human_2_Ctrl) <- c('Human_2_Ctrl')
names(human_3_Ctrl) <- c('Human_3_Ctrl')
names(human_1_Heat) <- c('Human_1_Heat')
names(human_2_Heat) <- c('Human_2_Heat')
names(human_3_Heat) <- c('Human_3_Heat')

human_1 <- merge(human_1_Heat, human_1_Ctrl, by=0, all=FALSE)
human_2 <- merge(human_2_Heat, human_2_Ctrl, by=0, all=FALSE)
human_3 <- merge(human_3_Heat, human_3_Ctrl, by=0, all=FALSE)

rownames(human_1) <- human_1$Row.names
rownames(human_2) <- human_2$Row.names
rownames(human_3) <- human_3$Row.names

human_1 <- human_1[2:length(human_1)]
human_2 <- human_2[2:length(human_2)]
human_3 <- human_3[2:length(human_3)]

df_counts <- merge(human_1, human_2, by=0, all=FALSE)
rownames(df_counts) <- df_counts$Row.names
df_counts <- df_counts[2:length(df_counts)]
df_counts <- merge(df_counts, human_3, by=0, all=FALSE)
rownames(df_counts) <- df_counts$Row.names
df_counts <- df_counts[2:length(df_counts)]

#rm(list=c('path_1','path_2','path_3','load_data','removal_list','drops','human_1','human_2','human_3'))

# heats <- c('Heat_1_A',
#            'Heat_1_B',
#            'Heat_1_C',
#            'Heat_2_A',
#            'Heat_2_B',
#            'Heat_2_C',
#            'Heat_3_A',
#            'Heat_3_B',
#            'Heat_3_C')
# ctrls <- c('Ctrl_1_A',
#            'Ctrl_1_B',
#            'Ctrl_1_C',
#            'Ctrl_2_A',
#            'Ctrl_2_B',
#            'Ctrl_2_C',
#            'Ctrl_3_A',
#            'Ctrl_3_B',
#            'Ctrl_3_C')

#heated  <- select(df_counts, heats)
#control <- select(df_counts, ctrls)


heats <- c('Human_1_Heat', 'Human_2_Heat', 'Human_3_Heat')
ctrls <- c('Human_1_Ctrl', 'Human_2_Ctrl', 'Human_3_Ctrl')
df_counts <- df_counts[, c(heats,ctrls)]


#######################################################################################################################
## Step 2: DESeq2 #####################################################################################################
#######################################################################################################################
group_heatctrl <- c('Heat','Heat','Heat','Ctrl','Ctrl','Ctrl')

DGE <- DGEList(counts = df_counts, group = group_heatctrl)
keep <- filterByExpr(DGE, group = group_heatctrl)
DGE  <- DGE[keep, keep.lib.sizes = FALSE]

# reps <- c('Heat_1','Heat_1','Heat_1','Heat_2','Heat_2','Heat_2','Heat_3','Heat_3','Heat_3',
#          'Ctrl_1','Ctrl_1','Ctrl_1','Ctrl_2','Ctrl_2','Ctrl_2','Ctrl_3','Ctrl_3','Ctrl_3')

colData <- data.frame(names=colnames(DGE$counts),  condition=DGE$samples$group)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(DGE$counts),
                              colData = colData,
                              design = ~ condition)

#dds <- collapseReplicates(dds, reps) # already collapsed

dds <- DESeq(dds)
res   <- results(dds) #, alpha=0.1, contrast=c('condition','Heat','Ctrl'))
rlog <- rlog(dds)
vst <- vst(dds)
lfc <- lfcShrink(dds, coef='condition_Heat_vs_Ctrl', type='apeglm')   #coef=2, res=res, type = 'apeglm')
lfc_ <- na.omit(lfc)
lfc_sig <- lfc_[lfc_$padj < 0.1,]

save(lfc, file=paste('/media/patrick/Elements/Hyperthermia/Cooked_Human_123_LFC/','LFC',   sep=''))

#rm(list=c('DGE', 'keep', 'group_heatctrl', 'colData', 'reps', 'ctrls', 'heats'))

#######################################################################################################################
## Step 3: Analyze ####################################################################################################
#######################################################################################################################

## Normalized Counts ##
norm_counts <- counts(dds, normalized=TRUE)

write.table(as.data.frame(norm_counts), 
            file="/media/patrick/Elements/Hyperthermia/Normalized_Counts/CookedHuman_123_counts_normed.csv",
            sep='\t')

## # Plot dispersions------------------------------------------------------------------------
plot_disps("Dispersions", dds)
dev.off()

# MA Plots ##--------------------------------------------------------------------------------
plot_MA("MA_LFC", lfc)
dev.off()

# Volcano Plots ##---------------------------------------------------------------------------
plot_volcano('Volcano_Healthy', lfc, 1., 0.05, 4, 30)
dev.off()

# PCA Plots ##--------------------------------------------------------------
plot_pca("PCA", vst) + geom_text(aes(label=name),vjust=2)
dev.off()

# Distancematrix ##------------------------------------------------------------------
plot_distancematrix("Distance_Matrix_Healthy_vst", vst)










