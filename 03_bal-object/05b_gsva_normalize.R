library(DESeq2)
library(vsn)


options(error = function() {traceback(3); q(1)})

BASE <- '/projects/b1038/Pulmonary/nmarkov/ssc-v2/data/31_bal-object/pseudobulk-gsva/'


sanitize_name <- function(name) {
  name <- gsub(' ', '_', name, fixed = TRUE)
  name <- gsub('*', '', name, fixed = TRUE)
  name <- gsub(';', '_and', name, fixed = TRUE)
  name <- gsub('/', '_', name, fixed = TRUE)
  return(name)
}


process_cell_type <- function(cell_type_path) {
  cell_type_dir <- dirname(dirname(cell_type_path))
  writeLines(sprintf("Starting %s", basename(cell_type_dir)))
  
  meta <- read.csv(cell_type_path, row.names = 1)
  rownames(meta) <- sub("\\_[ATCG]+$", "", rownames(meta))
  
  if (nrow(meta) < 6 || length(unique(meta$Sex)) < 2) {
    return()
  }
  
  counts <- read.csv(sub('-meta.csv', '.txt', cell_type_path), sep='\t', check.names = FALSE, row.names = 1)
  
  dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~Sex)
  dds <- DESeq(dds, fitType = "local")
  
  transformed <- assay(vst(dds))
  write.table(transformed, sprintf('%s/data/%s', cell_type_dir, 'transformed.tsv'))
}


metas <- list.files(path = BASE, full.names = TRUE, recursive = TRUE, pattern = '*-meta.csv')
for (meta_path in metas) {
  process_cell_type(meta_path)
}