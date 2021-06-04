# 4 Sample Venn Diagram
# Patrick Campbell
# April 2021

rm(list=ls())
################################################
##~ User Paths ~################################
################################################
path_dir <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/DESeq2/LFC/'
fig_save_path <-'/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/'
################################################
##~ Imports ~###################################
################################################
load_libraries <- function(){
  library(DESeq2)
  #library(ggVennDiagram)
  library(ggvenn)
  library(RColorBrewer)
  library(viridis)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
}
load_libraries()
################################################
##~ Functions ~###################################
################################################
loadRData <- function(DEObject_name){
  #loads an RData file, and returns it
  # DEObject_name: name of DEObject that was saved in your directory, ie: 'lfc'
  file_name <- paste(path_dir, DEObject_name, sep='')
  load(file_name)
  DEObject <- get(ls()[!(ls() %in% c("file_name","DEObject_name"))])
  return(DEObject)
}

prepData <- function(DEObject1, DEObject2, DEObject3, DEObject4, type, padj.value, names) {
  #DEObject is the LFC or Results data frame from DESeq2
  #type is character: positive, negative, or both and determines the type of signifcance to count for the venndiagram
  
  ###############################################
  DEObject1 <- na.omit(DEObject1)
  DEObject2 <- na.omit(DEObject2)
  DEObject3 <- na.omit(DEObject3)
  DEObject4 <- na.omit(DEObject4)
  
  if (type == 'positive') {
    DEObject1 <- DEObject1[DEObject1$log2FoldChange >= 0 ,]
    DEObject2 <- DEObject2[DEObject2$log2FoldChange >= 0 ,]
    DEObject3 <- DEObject3[DEObject3$log2FoldChange >= 0 ,]
    DEObject4 <- DEObject4[DEObject4$log2FoldChange >= 0 ,]
  }
  
  else if (type == 'negative') {
    DEObject1 <- DEObject1[DEObject1$log2FoldChange <= 0 ,]
    DEObject2 <- DEObject2[DEObject2$log2FoldChange <= 0 ,]
    DEObject3 <- DEObject3[DEObject3$log2FoldChange <= 0 ,]
    DEObject4 <- DEObject4[DEObject4$log2FoldChange <= 0 ,]
  }
  
  DEObject1 <- rownames(DEObject1[DEObject1$padj < padj.value ,])
  DEObject2 <- rownames(DEObject2[DEObject2$padj < padj.value ,])
  DEObject3 <- rownames(DEObject3[DEObject3$padj < padj.value ,])
  DEObject4 <- rownames(DEObject4[DEObject4$padj < padj.value ,])
  X <- list('1' = DEObject1, '2' = DEObject2, '3' = DEObject3, '4' = DEObject4)
  names(X) <- names
  return(X)
}

returnShared <- function(X1, X2) {
  # 
  stopifnot(class(X1) == 'character' && class(X2) == 'character')
  return(intersect(X1, X2))
}

returnNotShared <- function(X1, X2) {
  # 
  stopifnot(class(X1) == 'character' && class(X2) == 'character')
  not_shared1 <- setdiff(X1, X2)
  not_shared2 <- setdiff(X2, X1)
  return(list(not_shared1, not_shared2))
}


################################################
##~ Load Data ~#################################
################################################
lfc_tumor <- loadRData('lfc_tumor')
lfc_healthy <- loadRData('lfc_healthy')
lfc_heated <- loadRData('lfc_heated')
lfc_notheated <- loadRData('lfc_notheated')

################################################
##~ Call Functions ~############################
################################################

# Both
X <- prepData(lfc_healthy, 
              lfc_tumor, 
              lfc_heated, 
              lfc_notheated, 
              'both', 
              0.1, 
              names= c('Healthy','Tumor', 'Heated','Not Heated'))


p <- ggvenn(X[1:2]) #, show_percentage = FALSE) #, show_elements = TRUE,label_sep = "\n")
print(p)
pp <- returnShared(X[[1]], X[[2]])
ppp <- returnNotShared(X[[1]], X[[2]])


healthies     <- lfc_healthy[ppp[[1]] , ]
healthy_ups   <- healthies[healthies$log2FoldChange > 0 , ]
healthy_downs <- healthies[healthies$log2FoldChange < 0 , ]
print(nrow(healthy_downs))
print(nrow(healthy_ups))

tumories     <- lfc_tumor[ppp[[2]] , ]
tumor_ups   <- tumories[tumories$log2FoldChange > 0 , ]
tumor_downs <- tumories[tumories$log2FoldChange < 0 , ]
print(nrow(tumor_downs))
print(nrow(tumor_ups))


bothies <- # not possible bc genes may not be consistent log2fc between healthy n tumor


ggsave(plot = p, filename = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Both/HealthyvTumor_VennDiagram_both.svg')
write_csv(as.data.frame(pp), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Both/HealthyvTumor_VennDiagram_both_shared.csv')
write_csv(as.data.frame(ppp[[1]]), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Both/HealthyvTumor_VennDiagram_both_healthy.csv')
write_csv(as.data.frame(ppp[[2]]), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Both/HealthyvTumor_VennDiagram_both_tumor.csv')


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
# Positive
X <- prepData(lfc_healthy, 
              lfc_tumor, 
              lfc_heated, 
              lfc_notheated, 
              'positive', 
              0.1, 
              names= c('Healthy','Tumor', 'Heated','Not Heated'))


p <- ggvenn(X[1:2]) #, show_percentage = FALSE) #, show_elements = TRUE,label_sep = "\n")
print(p)
pp <- returnShared(X[[1]], X[[2]])
ppp <- returnNotShared(X[[1]], X[[2]])
ggsave(plot = p, filename = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Positive/HealthyvTumor_VennDiagram_both.svg')
write_csv(as.data.frame(pp), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Positive/HealthyvTumor_VennDiagram_both_shared.csv')
write_csv(as.data.frame(ppp[[1]]), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Positive/HealthyvTumor_VennDiagram_both_healthy.csv')
write_csv(as.data.frame(ppp[[2]]), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Positive/HealthyvTumor_VennDiagram_both_tumor.csv')


# Negative
X <- prepData(lfc_healthy, 
              lfc_tumor, 
              lfc_heated, 
              lfc_notheated, 
              'negative', 
              0.1, 
              names= c('Healthy','Tumor', 'Heated','Not Heated'))


p <- ggvenn(X[1:2]) #, show_percentage = FALSE) #, show_elements = TRUE,label_sep = "\n")
print(p)
pp <- returnShared(X[[1]], X[[2]])
ppp <- returnNotShared(X[[1]], X[[2]])
ggsave(plot = p, filename = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Negative/HealthyvTumor_VennDiagram_both.svg')
write_csv(as.data.frame(pp), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Negative/HealthyvTumor_VennDiagram_both_shared.csv')
write_csv(as.data.frame(ppp[[1]]), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Negative/HealthyvTumor_VennDiagram_both_healthy.csv')
write_csv(as.data.frame(ppp[[2]]), 
          f = '/home/patrick/Desktop/Masters Thesis Figures/Human_1/VennDiagram/Negative/HealthyvTumor_VennDiagram_both_tumor.csv')

