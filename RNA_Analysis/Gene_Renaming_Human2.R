# Gene Renaming Human 2

library("biomaRt")
library("dplyr")

make_df <- function(dir) {
  files <- list.files(path=dir, pattern="*.txt", full.names=TRUE, recursive=FALSE)
  
  read_in <- function(f) {
    df <- read.table(file=f, sep='\t', header=TRUE, stringsAsFactors = FALSE, colClasses = c(NA,"NULL","NULL","NULL","NULL","NULL",NA))
    colnames(df) <- c('gene_id', 'count')
    return(df)
  }
  
  x <- lapply(files, read_in)
  x1 <- x[1]
  
  for (i in 2:length(x)) {
    x1 <- merge(x1, x[i], by=c('gene_id'))
  }
  
  colnames(x1) <- c('gene_id', '01', '02', '03', '04', '05', '06')
  return(x1)
}

query_ensembl <- function(df, dir) {
  ensembl <- useMart("ensembl")
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  
  query <- getBM(attributes=c('ensembl_gene_id_version', 'hgnc_symbol'),
                 mart=ensembl,useCache = FALSE)
  
  idx <- match(df$gene_id, query$ensembl_gene_id_version)
  df$gene_name <- query$hgnc_symbol[idx]
  
  for (i in 1:nrow(df)) {
    if (df$gene_name[i] == "" || is.na(df$gene_name[i]) ) {
      df$gene_name[i] <- df$gene_id[i]
    }
  }
  
  df <- select(df, 'gene_name', '01', '02', '03', '04', '05', '06')
  
  
  #f_save <- paste(gsub('.{4}$', '', f),'.csv',sep='')
  f_save <- paste(dir, 'counts.csv')
  write.csv(df, f_save)
  
  return(df)
}


dir <- "/media/patrick/Elements/Hyperthermia/Human_2/Counts/"
x <- make_df(dir)
x <- query_ensembl(x, dir)


# Names for Human 1
colnames(x) <- c('gene_name', 
                 '2_Ctrl_A',
                 '2_Ctrl_B',
                 '2_Ctrl_C',
                 '2_Heat_A',
                 '2_Heat_B',
                 '2_Heat_C')
                 #'3_Ctrl_A',
                 #'3_Ctrl_B',
                 #'3_Ctrl_C',
                 #'3_Heat_A',
                 #'3_Heat_B',
                 #'3_Heat_C')

x <- x[!duplicated(x$gene_name),]

dir_2 <- dir <- "/media/patrick/Elements/Hyperthermia/Human_2/CountMatrix/"
f_save <- paste(dir_2, 'human_2_counts.csv', sep='')
write.csv(x, f_save)
