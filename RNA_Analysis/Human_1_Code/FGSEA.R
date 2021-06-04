# Pretty Figures fgsea
rm(list=ls())

# Input: DESEQ2 outputs
# Run the DESEQ program first to get the datasets
#####################################################################################################
####### Load Libaries ~##############################################################################
#####################################################################################################
load_libraries <- function()  {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  #library(GO.db)
  library(enrichplot)
  library(DOSE)
  library(ggupset)
  library(zeallot)
  library(ggridges)
  library(fgsea)
  library(data.table)
  library(ggplot2)
  library(msigdbr)
  library(tidyverse)
  library(DESeq2)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)
  library(GOSemSim)
  library(pheatmap)
  }
load_libraries()

#####################################################################################################
####### Functions ~##################################################################################
#####################################################################################################
loadRData <- function(DEObject_name){
  #loads an RData file, and returns it
  # DEObject_name: name of DEObject that was saved in your directory, ie: 'lfc'
  file_name <- paste(path_dir, DEObject_name, sep='')
  load(file_name)
  DEObject <- get(ls()[!(ls() %in% c("file_name","DEObject_name"))])
  return(DEObject)
}

get_rank <- function(df, method){
  #input is a results dataframe, function returns ranked genes
  if (method == 'fancy'){
    df$rank <- sign(df$log2FoldChange) * -log10(df$pvalue)
  } else if (method == 'simple') {
    df$rank <- df$log2FoldChange
  }
  selected <- na.omit(df)
  selected <- selected[order(selected$log2FoldChange),]
  selected <- selected[, c("rank")]
  return(data.frame(selected))
}

get_geneList <- function(df){
  # assume input straight from get_degs()
  l <- df[,1]
  names(l) <- rownames(df)
  l <- sort(l, decreasing = TRUE)
  return(l)
}

get_genes <- function(df, cutoff){
  # assume input straight from get_geneList
  l <- names(df)[abs(geneList)>=cutoff]
  if (length(l) < 20) {
    print("fewer than 20 genes with this cutoff value")
    print(length(l))
  }
  return(l)
}

lister <- function(pathways, grepper) {
  x <- pathways[grep(grepper, names(pathways))]
  print(length(pathways[!names(pathways) %in% names(x)]))
  #kpathways_list <<- pathways[!names(pathways) %in% names(x)]
  return(x)
}
#####################################################################################################
path_dir <- '/media/patrick/Elements/Hyperthermia/CookedHuman1/DESeq2/LFC/'
fig_save_path <- '/home/patrick/Desktop/practicecrap/'

df_h <- loadRData('lfc_healthy')
df_t <- loadRData('lfc_tumor')

rank_method <- 'fancy'

degs_h <- get_rank(df_h, rank_method)
geneList_h <- get_geneList(degs_h)
df_h$SYMBOL <- rownames(df_h)

degs_t <- get_rank(df_t, rank_method)
geneList_t <- get_geneList(degs_t)
df_t$SYMBOL <- rownames(df_t)


rm(list=c('get_rank','get_geneList','get_genes','path_dir'))
#####################################################################################################
## fgsea
# Kevins data
f_kev <- "/home/patrick/Desktop/Gene_Pathways/allgenesets.RDS"
kpathways <- readRDS(f_kev, character())
kpathways_list <- split(x = kpathways$gene, f = kpathways$ont)  # set up for fgsea
rm(list=c('kpathways', 'f_kev'))



k.BP.GO <- lister(kpathways_list, "^BP.GO_")
k.CC.GO <- lister(kpathways_list, "^CC.GO_")
k.MF.GO <- lister(kpathways_list, "^MF.GO_")


all_gene_sets <- msigdbr(species = 'Homo sapiens', category = 'H')
all_gene_sets$entrez_gene <- as.character(all_gene_sets$entrez_gene)
msigdbr_list <- split(x = all_gene_sets$human_gene_symbol, f = all_gene_sets$gs_name)  # set up for fgsea
#msigdbr_list_chr <- as.character(msigdbr_list) 


active_list <- k.BP.GO

min_size = 20
max_size = 400
nperms   = 1000000

fgseaRes_h <- fgsea(pathways = active_list, 
                    stats    = geneList_h,
                    minSize  = min_size,
                    maxSize  = max_size,
                    nperm    = nperms)

fgseaRes_t <- fgsea(pathways = active_list, 
                    stats    = geneList_t,
                    minSize  = min_size,
                    maxSize  = max_size,
                    nperm    = nperms)

# fgseaRes_h_msig <- fgsea(pathways = msigdbr_list, 
#                     stats    = geneList_h,
#                     minSize  = 1,
#                     maxSize  = 500,
#                     nperm    = 1000000)
# 
# fgseaRes_t_msig <- fgsea(pathways = msigdbr_list, 
#                     stats    = geneList_t,
#                     minSize  = 1,
#                     maxSize  = 500,
#                     nperm    = 1000000)



fgseaRes_h_ <- fgseaRes_h[fgseaRes_h$padj <= 0.2,]
fgseaRes_t_ <- fgseaRes_t[fgseaRes_t$padj <= 0.2,]
# fgseaRes_h_msig <- fgseaRes_h_msig[fgseaRes_h_msig$pval <= 0.05,]
# fgseaRes_t_msig <- fgseaRes_t_msig[fgseaRes_t_msig$pval <= 0.05,]



fgseaRes_h_ <- fgseaRes_h_[order(fgseaRes_h_$NES) ,]
fgseaRes_t_ <- fgseaRes_t_[order(fgseaRes_t_$NES) ,]
# fgseaRes_h_msig <- fgseaRes_h_msig[order(fgseaRes_h_msig$NES) ,]
# fgseaRes_t_msig <- fgseaRes_t_msig[order(fgseaRes_t_msig$NES) ,]

# # Tidy the results:
# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES)) # order by normalized enrichment score (NES)
# fgseaResTidy
f_h <- paste(fig_save_path, 'gsea_healthy_heatresponse.svg', sep='')
f_t <- paste(fig_save_path, 'gsea_tumor_heatresponse.svg', sep='')

pathway <- 'BP.GO_CELLULAR_RESPONSE_TO_HEAT'
svg(f_h)
plotEnrichment(active_list[[pathway]],
               geneList_h) + labs(title="Heat Response")
dev.off()

svg(f_t)
plotEnrichment(active_list[[pathway]],
               geneList_t) + labs(title="Heat Response")
dev.off()



############################################################################################
####### Cluster Profiler GSEA ~#############################################################
############################################################################################
# get_gseGO <- function(gene_list, GO_type, keytype, pcutoff) {
#   # Citation: ClusterProfiler
#   # Usage:
#   #   gene_list -- output from function: get_geneList
#   #   GO_type   -- 'MF', 'BP', or 'CC'
#   #   keytype   -- 'SYMBOL' or 'ENTREZID'  
#   #   pcutoff   -- 0.05
#   org <- "org.Hs.eg.db"
#   ont <- GO_type
#   hs <- godata(org, ont = ont)
#   gse <- gseGO(gene_list,
#                OrgDb = org,
#                keyType = keytype,
#                ont = ont,
#                pvalueCutoff = pcutoff,
#                minGSSize = 1)
#   
#   # Split into positive and negative
#   # gse.POS <- gse[sign(gse$enrichmentScore) == 1  , ]
#   # gse.NEG <- gse[sign(gse$enrichmentScore) == -1 , ]
#   # go_ids.pos <- rownames(gse.POS[gse.POS$p.adjust < 0.1 , ])
#   # go_ids.neg <- rownames(gse.NEG[gse.NEG$p.adjust < 0.1 , ])
#   # 
#   #return(list(hs, gse.POS, gse.NEG, go_ids.pos, go_ids.neg, gse))
#   return(gse)
# }
# 
# 
# gse_h <- get_gseGO(gene_list = geneList_h, GO_type = 'BP', keytype = 'SYMBOL', pcutoff   = 1)
# dotplot(gse_h, showCategory=10, split=".sign") + facet_grid(.~.sign)

############################################################################################
####### Curated Pathways GSEA ~#############################################################
############################################################################################
interesting_pathways <- function(pathways) {
  # Here are some interesting pathways selected for our own data. The negative pathways were the primary negative pathways
  # in our dataset. They should be changed to match the primary negative pathways in each individual dataset.
  #metabolic_pathways <- c('HALLMARK_FATTY_ACID_METABOLISM', 
  #                        'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 
  #                        'HALLMARK_GLYCOLYSIS',
  #                        'BP.GO.AEROBIC_RESPIRATION')
  
  #miRNA_pathways <- c('GO:2000637', 'GO:2000636', 'GO:0035196', 'GO:0030918')
  
  apoptosis_pathways <- c('BP.GO_APOPTOTIC_SIGNALING_PATHWAY',
                          'BP.GO_CELL_DEATH_IN_RESPONSE_TO_OXIDATIVE_STRESS', 
                          'BP.GO_NEURON_DEATH',
                          'HALLMARK_HYPOXIA',
                          'KEGG_MAPK_SIGNALING_PATHWAY',
                          'BP.GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY')
                          #'BIOCARTA_AKT_PATHWAY') #stat3 is not really apoptosis, but can lead to t cell suppression
                           #'BIOCARTA_CASPASE_PATHWAY',
                           #'#'BIOCARTA_STAT3_PATHWAY',
  
  growth_pathways    <- c('BP.GO_CELL_GROWTH', 
                          'BP.GO_REPRODUCTION',
                          'BP.GO_POSITIVE_REGULATION_OF_VASCULATURE_DEVELOPMENT',
                          'BP.GO_EPITHELIAL_CELL_DEVELOPMENT')
  #                          'BP.GO_DEVELOPMENTAL_CELL_GROWTH', 

  heat_pathways     <-   c('BP.GO_CHAPERONE_MEDIATED_PROTEIN_FOLDING', 
                           'BP.GO_CELLULAR_RESPONSE_TO_HEAT', 
                           'BP.GO_PROTEIN_FOLDING',
                           'BP.GO_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN',
                           'BP.GO_RESPONSE_TO_TEMPERATURE_STIMULUS',
                           'REACTOME_HSF1_DEPENDENT_TRANSACTIVATION')
                          #'BP.GO_DE_NOVO_PROTEIN_FOLDING',
  
  other_pathways <- c('BP.GO_REGULATION_OF_INFLAMMATORY_RESPONSE',
                         'BP.GO_REGULATION_OF_LEUKOCYTE_PROLIFERATION',
                         'BP.GO_MYELOID_LEUKOCYTE_MIGRATION',
                         'BP.GO_LYMPHOCYTE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                         'BP.GO_GRANULOCYTE_MIGRATION',
                         'BP.GO_DNA_REPLICATION',
                      'BP.GO_T_CELL_ACTIVATION')
  
  combined_pathways <- c(apoptosis_pathways, growth_pathways, heat_pathways, other_pathways)
  
  pathways <- pathways[names(pathways) %in% combined_pathways]
  return(pathways)
}

pathways <- interesting_pathways(kpathways_list)
active_list <- pathways
min_size = 1
max_size = 1000
nperms   = 1000000

fgseaRes_h_interesting <- fgsea(pathways = active_list, 
                    stats    = geneList_h,
                    minSize  = min_size,
                    maxSize  = max_size,
                    nperm    = nperms)

fgseaRes_t_interesting <- fgsea(pathways = active_list, 
                    stats    = geneList_t,
                    minSize  = min_size,
                    maxSize  = max_size,
                    nperm    = nperms)

# fgseaRes_h <- fgseaRes_h[fgseaRes_h$pval <= 0.05,]
# fgseaRes_t <- fgseaRes_t[fgseaRes_t$pval <= 0.05,]

fgseaRes_h_interesting2 <- fgseaRes_h_interesting[order(fgseaRes_h_interesting$NES) ,]
fgseaRes_t_intersting2 <- fgseaRes_t_interesting[order(fgseaRes_t_interesting$NES) ,]







  