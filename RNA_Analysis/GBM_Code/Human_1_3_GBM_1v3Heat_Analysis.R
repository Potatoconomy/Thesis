#'DESEQ2'

# Cooked Humans 1 and 3 -- GBM Humans


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

fig_save_path <- '/media/patrick/Elements/Hyperthermia/Figures_Human_GBM/Heat_Comparison/'


#######################################################################################################################
## Step 1: Organize DataFrame #########################################################################################
#######################################################################################################################

path_1 <- '/media/patrick/Elements/Hyperthermia/Cooked_Human_CountMatrixes/human_1_counts'
path_3 <- '/media/patrick/Elements/Hyperthermia/Cooked_Human_CountMatrixes/human_3_counts'

human_1 <- load_data(path_1)
human_3 <- load_data(path_3)

rownames(human_1) <- human_1$gene_name
rownames(human_3) <- human_3$gene_name

drops <- c("X","gene_name")
human_1 <- human_1[ , !(names(human_1) %in% drops)]
human_3 <- human_3[ , !(names(human_3) %in% drops)]

#human_1 <- removal_list(human_1)  # Drops Mt- and RPS genes
#human_3 <- removal_list(human_3)

human_1 <- human_1 %>% drop_na()
human_3 <- human_3 %>% drop_na()

human_1_Heat <- as.data.frame(human_1[,c(4,5,6)])
human_3_Heat <- as.data.frame(human_3[,c(4,5,6)])

names(human_1_Heat) <- c('Human_1_Heat_A', 'Human_1_Heat_B', 'Human_1_Heat_C')
names(human_3_Heat) <- c('Human_3_Heat_A', 'Human_3_Heat_B', 'Human_3_Heat_C')

human <- merge(human_1_Heat, human_3_Heat, by=0, all=FALSE)

rownames(human) <- human$Row.names

human <- human[2:length(human)]

human <- human[rowSums(human) > 15, ] # Filter low counts

df_counts <- human

heats <- c('Human_1_Heat_A', 'Human_1_Heat_B', 'Human_1_Heat_C')
ctrls <- c('Human_3_Heat_A', 'Human_3_Heat_B', 'Human_3_Heat_C')
df_counts <- df_counts[, c(heats,ctrls)]

rm(list=c('path_1','path_3','load_data','removal_list','drops','human_1','human_3'))
#######################################################################################################################
## Step 2: DESeq2 #####################################################################################################
#######################################################################################################################
group_heatctrl <- c('h1','h1','h1','h3','h3','h3')

DGE <- DGEList(counts = df_counts, group = group_heatctrl)
keep <- filterByExpr(DGE, group = group_heatctrl)
DGE  <- DGE[keep, keep.lib.sizes = FALSE]

# reps <- c('Heat_1','Heat_1','Heat_1','Heat_2','Heat_2','Heat_2','Heat_3','Heat_3','Heat_3',
#          'Ctrl_1','Ctrl_1','Ctrl_1','Ctrl_2','Ctrl_2','Ctrl_2','Ctrl_3','Ctrl_3','Ctrl_3')

colData <- data.frame(names=colnames(DGE$counts),  condition=DGE$samples$group, batch = c('A','A','A','B','B','B'))

dds <- DESeqDataSetFromMatrix(countData = as.matrix(DGE$counts),
                              colData = colData,
                              design = ~condition)

#dds <- collapseReplicates(dds, reps) # already collapsed

dds <- DESeq(dds)
res   <- results(dds, alpha=0.1, contrast = c('condition', 'h3', 'h1'))#, lfcThreshold = 0) #, alpha=0.1, contrast=c('condition','Heat','Ctrl'))
#rlog <- rlog(dds)
vst <- vst(dds)
lfc <- lfcShrink(dds, coef='condition_h3_vs_h1', type='apeglm')   #coef=2, res=res, type = 'apeglm')
lfc_ <- na.omit(lfc)
lfc_sig <- lfc_[lfc_$padj < 0.1,]
lfc_sig2 <- lfc_[lfc_$pvalue < 0.05,]
lfc_ <- lfc_[order(lfc_$padj), ]
plot(lfc_$padj)
lfc_ <- lfc_[order(lfc_$pvalue), ]
plot(lfc_$pvalue)
plot(lfc_$pvalue[0:1000], lfc_$log2FoldChange[0:1000])

save(lfc, file=paste(fig_save_path,'LFC',   sep=''))

plot(lfc_sig2$pvalue, lfc_sig2$log2FoldChange)
dev.off()
#rm(list=c('DGE', 'keep', 'group_heatctrl', 'colData', 'reps', 'ctrls', 'heats'))

#######################################################################################################################
## Step 3: Analyze ####################################################################################################
#######################################################################################################################

## Normalized Counts ##
norm_counts <- counts(dds, normalized=TRUE)

write.table(as.data.frame(norm_counts), 
            file=paste(fig_save_path, 'CookedHuman_13_heat_comparison_counts_normed.csv',sep=''),
            sep='\t')


norm_counts <- as.data.frame(norm_counts)
norm_counts$Human_1_Heated <- rowSums(norm_counts[1:3]/3)
norm_counts$Human_3_Heated <- rowSums(norm_counts[4:6]/3)
norm_counts <- norm_counts[c('Human_1_Heated','Human_3_Heated')]
norm_counts <- norm_counts[order(norm_counts$Human_1_Heated) , ]
norm_counts_ <- norm_counts
norm_counts_ <- norm_counts_[norm_counts_$Human_1_Heated<400 , ]
norm_counts_ <- norm_counts_[norm_counts_$Human_3_Heated<400 , ]
norm_counts_ <- norm_counts_[norm_counts_$Human_1_Heated>10 , ]
norm_counts_ <- norm_counts_[norm_counts_$Human_3_Heated>10 , ]

lfc_ordered <- lfc[order(lfc$log2FoldChange) , ]
labels <- c(rownames(lfc_ordered[1:50, ]), rownames(lfc_ordered[3110:3160 , ]))

p <- ggplot(norm_counts_, aes(x=Human_1_Heated, y=Human_3_Heated)) + geom_point() + xlim(0,400) + ylim(0,400) +
  geom_text(
    label=labels, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )

svg(filename = paste(fig_save_path, 'normedcounts.svg'))
print(p)
dev.off()

## # Plot dispersions------------------------------------------------------------------------
plot_disps("Dispersions", dds)

# MA Plots ##--------------------------------------------------------------------------------
plot_MA("MA_LFC", lfc)

plot_MA("MA_RES", res)

# Volcano Plots ##---------------------------------------------------------------------------
plot_volcano_pvalue('Volcano_GBM_lfc', lfc, 0.5, 0.05, 7, 15)
dev.off()

plot_volcano_pvalue('Volcano_GBM_res', res, 0.5, 0.05, 7, 15)
dev.off()

# PCA Plots ##--------------------------------------------------------------
plot_pca("PCA", vst) + geom_text(aes(label=name),vjust=2)
dev.off()

# Distancematrix ##------------------------------------------------------------------
plot_distancematrix("Distance_Matrix_Healthy_vst", vst)

