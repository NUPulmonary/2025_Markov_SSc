library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(vsn)
library(R.utils)


options(error = function() {traceback(3); q(1)})

BASE <- '/projects/b1038/Pulmonary/nmarkov/ssc-v2/data/31_bal-object/pseudobulk-gsva/'
OUT_DIR <- '/projects/b1038/Pulmonary/nmarkov/ssc-v2/data/31_bal-object/15_pseudobulk-mmf/'

SAMPLES <- read.csv('/projects/b1038/Pulmonary/nmarkov/ssc-v2/code/00all-samples.csv', check.names = FALSE)
CLINICAL <- read.csv('/projects/b1038/Pulmonary/nmarkov/ssc-v2/code/00clinical-v2.csv', check.names = FALSE)


plot_lib_sizes <- function(data) {
    df <- data.frame(lib.size = sort(colSums(data)))
    ggplot(df, aes(x = reorder(rownames(df), lib.size), y = lib.size)) +
        geom_col() +
        ggtitle("Sorted library sizes") +
        theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
        xlab("")
}

process_cell_type <- function(cell_type) {
  print(sprintf("Processing %s", cell_type))
  data <- read.table(
    sprintf("%s/%s/data/%s.txt", BASE, cell_type, cell_type), 
    header = TRUE, 
    check.names = FALSE
  )
  coldata <- read.csv(
    sprintf("%s/%s/data/%s-meta.csv", BASE, cell_type, cell_type), 
    check.names = FALSE
  )
  coldata <- coldata[coldata$n_cells > 10, ]

  coldata <- merge(coldata, SAMPLES[, c('External Sample ID', 'Sample')], by = 'External Sample ID')
  coldata$study_code <- gsub('_followup', '', coldata$Sample)
  coldata <- merge(coldata, CLINICAL[, c('study_code', 'patient_on_mmf')], by = 'study_code')
  coldata$mmf <- c(True = 'MMF', False = 'Naive')[coldata$patient_on_mmf]

  # Select only SSc patients
  coldata <- coldata[coldata$Status == "SSc", ]
  
  if (count(coldata$mmf == 'MMF') < 2) {
    print("Only one or two MMF samples, skipping")
    return()
  }

  if (count(coldata$mmf == 'Naive') < 2) {
    print("Only one or two control samples, skipping")
    return()
  }

  
  dir <- sprintf("%s/%s", OUT_DIR, cell_type)
  R.utils::mkdirs(dir)
  write.csv(
    coldata,
    sprintf("%s/meta.csv", dir)
  )
  
  coldata$n_cells <- scale(coldata$n_cells)
  counts <- as.matrix(data[, coldata[["External Sample ID"]]])
  dds <- DESeqDataSetFromMatrix(counts, coldata, ~ mmf)
  
  pdf(sprintf("%s/lib-sizes.pdf", dir), width = 6, height = 4)
  plot(plot_lib_sizes(counts))
  dev.off()
  
  # Library complexity
  # TODO: explore other measures of library complexity, like dispersion, skew??
  d <- data.frame(
    n_genes = colSums(counts(dds) > 0),
    uniq_counts = apply(counts(dds), 2, function(x) {return(length(unique(x)))})
  )
  
  # Plot number of genes detected in each sample
  pdf(sprintf("%s/n_genes.pdf", dir), width = 6, height = 4)
  plot(ggplot(d, aes(x = reorder(rownames(d), n_genes), y = n_genes)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    xlab(""))
  dev.off()
  
  # Number of unique counts per library
  pdf(sprintf("%s/n_unique_counts.pdf", dir), width = 6, height = 4)
  plot(ggplot(d, aes(x = reorder(rownames(d), uniq_counts), y = uniq_counts)) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    xlab(""))
  dev.off()
  
  # Plot heatmap of library similarity
  sampleDists <- dist(t(assay(vst(dds))))
  
  pdf(sprintf("%s/lib-dist.pdf", dir), width = 8, height = 8)
  pheatmap(
    as.matrix(sampleDists), 
    cluster_cols = TRUE, 
    cluster_rows = TRUE, 
    show_colnames = FALSE,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists
  )
  dev.off()
  
  status_sex_dir <- sprintf("%s/%s/degs", OUT_DIR, cell_type)
  R.utils::mkdirs(status_sex_dir)
  
  dds <- DESeq(dds, fitType = "parametric")
  pdf(sprintf("%s/disp-parametric.pdf", status_sex_dir), width = 6, height = 4)
  plotDispEsts(dds)
  dev.off()
  
  dds <- DESeq(dds, fitType = "local")
  pdf(sprintf("%s/disp-local.pdf", status_sex_dir), width = 6, height = 4)
  plotDispEsts(dds)
  dev.off()
  
  pdf(sprintf("%s/pca.pdf", status_sex_dir), width = 6, height = 4)
  d <- plotPCA(
    vst(dds, blind = FALSE), 
    ntop = 2000,
    intgroup = c("mmf", "Sex"), 
    returnData = TRUE
  )
  percentVar <- round(100 * attr(d, "percentVar"))
  
  plot(ggplot(d, aes(
    x = PC1, 
    y = PC2, 
    color = mmf, 
    shape = Sex
  )) +
    geom_point(size = 2) + 
    geom_label_repel(aes(label = colnames(dds)), show.legend = FALSE, fill = "white", size = 2) +
    scale_shape_manual(values = c("female" = 21, "male" = 23)) +
    scale_color_manual(values = c("Naive" = "blue", "MMF" = "red")) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    coord_fixed() +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))
  )
  dev.off()
  
  degs <- as.data.frame(results(dds, contrast = c("mmf", "Naive", "MMF"), alpha = 0.05))
  degs <- degs[!is.na(degs$padj), ]
  degs$sign <- ""
  degs$sign[(degs$log2FoldChange < 0) & (degs$padj < 0.05)] <- "Downregulated"
  degs$sign[(degs$log2FoldChange > 0) & (degs$padj < 0.05)] <- "Upregulated"
  n_downreg <- sum(degs$sign == "Downregulated")
  degs$sign[degs$sign == "Downregulated"] <- sprintf("Up in MMF (%d)", n_downreg)
  n_upreg <- sum(degs$sign == "Upregulated")
  degs$sign[degs$sign == "Upregulated"] <- sprintf("Up in Naive (%d)", n_upreg)
  degs <- degs[order(degs$sign),]
  
  pdf(sprintf("%s/ma.pdf", status_sex_dir), width = 8, height = 6)
  max_x = max(degs$baseMean)
  plot(ggplot(degs) + 
    geom_point(aes(x = baseMean, y = log2FoldChange, col = sign)) +
    scale_x_log10() +
    scale_color_manual(values = c("gray", "blue", "red")) + 
    ggtitle("Naive vs MMF"))
  dev.off()
  
  write.csv(
    degs[order(degs$log2FoldChange, decreasing = TRUE), ],
    sprintf("%s/degs.csv", status_sex_dir)
  )
  
  sign.only <- degs[(degs$padj < 0.05), ]
  write.csv(
    sign.only[order(sign.only$log2FoldChange, decreasing = TRUE), ],
    sprintf("%s/degs-sign.csv", status_sex_dir)
  )
  
  print("done")
}

cell_types <- list.files(BASE)
for (cell_type in cell_types) {
  process_cell_type(cell_type)
}

