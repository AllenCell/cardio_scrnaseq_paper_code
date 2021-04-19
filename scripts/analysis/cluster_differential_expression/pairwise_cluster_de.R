#!/usr/bin/env Rscript

library(Seurat)
library(edgeR)

# load seurat object
load("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/D0_D12_D24_D90.RDadta")

setwd("/allen/aics/gene-editing/Manuscripts/Transcriptomics_2019/scrnaseq_supplement/paper_plots/small_molecule_cardio_timepoints_indepth/de_analysis_calcnorm/")


get_gene_name <- function(name){
  gene_name <- unlist(strsplit(name, '_'))
  name_length <- length(gene_name)
  return(paste((gene_name[-name_length]), collapse = "_"))
}


save_dir <- "cluster_de_genes/"

clusters <- unique(cardio@meta.data$res.0.5)
cluster_count <- length(clusters)

# pairwise cluster comparisons; resolution=0.5
cluster_index <- 2
for (i in 1:14) {
  for (j in cluster_index:14) {
    if (j > cluster_count) {
      break
    }
    print(c(i, j))
    current.cells <- rownames(cardio@meta.data)[cardio@meta.data$res.0.5 %in% c(clusters[i], clusters[j])]
    current.meta.data <- cardio@meta.data[current.cells, ]
    raw_counts = cardio@raw.data[, current.cells]

    keep.genes = raw_counts > 0                                                                                                                                                
    keep.genes = rowSums(keep.genes) >= 10

    raw_counts = raw_counts[keep.genes, current.cells]

    dge <- DGEList(
      counts = raw_counts,
      group = factor(current.meta.data$res.0.5)
    )
    dge <- calcNormFactors(dge)
    group_cluster <- factor(current.meta.data$res.0.5)
    design <- model.matrix(~0+group_cluster)
    
    cluster_pair <- levels(group_cluster)
    dir.create(paste0(save_dir, "/cluster", cluster_pair[1], "v", cluster_pair[2]), recursive = TRUE)
#     print(design)
    
    
    dge <- estimateDisp(dge, design = design)
    
    fit <- glmFit(dge, design)
    
    
    save(fit, file = paste0(save_dir, "fit_cluster", cluster_pair[1], "v", cluster_pair[2], ".Robj"))
    
    res <- glmLRT(fit, contrast = c(1, -1))
    
    print(res$comparison)
    
    results_table <- topTags(res, n = 5000L, sort.by = "logFC", p.value = 0.05)
    
    results_table <- results_table$table[abs(results_table$table[["logFC"]]) >= 1, ]
    
    pos_cat <- cluster_pair[1]
    neg_cat <- cluster_pair[2]
    up.cluster.cells <- rownames(current.meta.data)[current.meta.data[["res.0.5"]] == pos_cat]
    down.cluster.cells <- rownames(current.meta.data)[current.meta.data[["res.0.5"]] == neg_cat]
    
    count_mat <- cardio@raw.data[, current.cells]
    
    genes.expFraction.up <- sapply(X = rownames(results_table), 
                                   FUN = function(x) {
                                     gene = x
                                     count_mat = count_mat
                                     cell_metadata = current.meta.data
                                     gene_lfc = results_table["logFC"][x,]
                                     
                                     cluster.cells <- c()
                                     if (gene_lfc > 0) {
                                       cluster.cells <- up.cluster.cells
                                     } else if (gene_lfc < 0) {
                                       cluster.cells <- down.cluster.cells
                                     }
                                     gene_counts <- count_mat[gene, cluster.cells]
                                     detected_cells <- gene_counts > 0
                                     return(sum(detected_cells) / length(detected_cells))
                                   },
                                   
                                   simplify = TRUE)
    
    genes.expFraction.down <- sapply(X = rownames(results_table), 
                                     FUN = function(x) {
                                       gene = x
                                       count_mat = count_mat
                                       cell_metadata = current.meta.data
                                       gene_lfc = results_table["logFC"][x,]
                                       
                                       cluster.cells <- c()
                                       if (gene_lfc < 0) {
                                         cluster.cells <- up.cluster.cells
                                       } else if (gene_lfc > 0) {
                                         cluster.cells <- down.cluster.cells
                                       }
                                       gene_counts <- count_mat[gene, cluster.cells]
                                       detected_cells <- gene_counts > 0
                                       return(sum(detected_cells) / length(detected_cells))
                                     },
                                     
                                     simplify = TRUE)
    
    temp_df <- data.frame(up = genes.expFraction.up, down = genes.expFraction.down)
    
    all.equal(rownames(temp_df), rownames(results_table))
    
    results_table <- cbind(results_table, temp_df)
    
    results_table.final <- results_table[results_table$up > 0.30, ]
    
    gene_df <- data.frame(gene = sapply(rownames(results_table.final), get_gene_name, USE.NAMES = TRUE, simplify = TRUE))
    
    write.csv(results_table.final, file = paste0(save_dir, "/cluster", pos_cat, "v", neg_cat, "/cluster", pos_cat, "v", neg_cat, ".csv"))
  }
  cluster_index <- cluster_index + 1 
}
