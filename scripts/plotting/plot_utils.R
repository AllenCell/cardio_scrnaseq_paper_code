# Functions used to make plots


# get sample # for each cell
get_sample_num = function(barcode) {
    sample_id = as.integer(unlist(strsplit(barcode, '_'))[[3]])
    return(paste0("E", sample_id + 1))
}


# umap plots color coded by a cell variable (ex. day, protocol)
make_umap_subset = function(embeddings,
                            color,
                            values,
                            labels,
                            rasterize=FALSE,
                            legend="",
                            point_size=0.05,
                            font_size=13,
                            guide_size=4
                            ) {
    p.umap = ggplot(embeddings, aes_string("UMAP1", "UMAP2", color = color))

    if (rasterize) {
        p.umap = p.umap + geom_point_rast(alpha = 1, size = point_size)
    } else {

        p.umap = p.umap + geom_point(alpha = 1, size = point_size)
    }

    p.umap = p.umap +
        scale_color_manual(name=legend, values=values) +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        guides(color = guide_legend(nrow = 1, override.aes = list(size=guide_size))) +
        theme(legend.text = element_text(size = font_size)) +
        theme(legend.position = "top", legend.justification = "center") +
        theme(text = element_text(size = font_size), axis.text = element_text(size = font_size)) 

    return(p.umap)
}


# umap plots color coded by gene normalized transcript abundance
make_umap_gene = function(embeddings,
                          gene,
                          label,
                          rasterize=FALSE,
                          point_size=0.5,
                          font_size=13,
                          guide_size=4,
                          legend_pos=c(0.8, 0.95),
                          legend_dir="horizontal"
                          ) {
    p.gene = ggplot(embeddings, aes_string("UMAP1", "UMAP2", color = gene))

    if (rasterize) {
        p.gene = p.gene + geom_point_rast(alpha = 1, size = point_size)
    } else {
        p.gene = p.gene + geom_point(alpha = 1, size = point_size)
    }
        
    p.gene = p.gene +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        scale_color_gradient(name = label, low = "gray", high = "red") +
        theme(legend.text = element_text(size = font_size)) +
        theme(legend.title = element_text(size = font_size)) + 
        theme(text = element_text(size = font_size), axis.text = element_text(size = font_size)) +
        theme(legend.position = legend_pos) +
        theme(legend.direction=legend_dir)

    return(p.gene)
}


# group violin plot by cluster
cluster_group_violin = function(seurat_obj,
                                cluster_color_df,
                                group_order = c("0", "1", "2", "3", "4", "5", "6", "7"),
                                resolution = "res.0.4",
                                genes.plot = c("MKI67", "TNNT2", "MYH6", "MYH7"),
                                font_size = 11
                                ) {
    cells.plot = rownames(seurat_obj@meta.data[seurat_obj@meta.data[, resolution] %in% group_order, ])

    gene_data = data.frame("sample_name" = cells.plot)
    gene_df = sapply(X = genes.plot, FUN = function(x) { return(seurat_obj@data[x, cells.plot]) })
    gene_data = cbind(gene_data, gene_df)


    tail(genes.plot, n = 1)
    length(genes.plot)

    anno = data.frame("sample_name" = cells.plot)
    anno = cbind(anno, seurat_obj@meta.data[cells.plot, ])

    anno$cluster_id = anno[, resolution]
    anno$cluster_label = sapply(X = anno[, resolution], FUN = function(x) { return(paste0("C", x))})

    anno$cluster = anno[, resolution]
    anno = merge(anno, cluster_color_df, by = "cluster")

    cluster_cell_types = group_violin_plot(data = gene_data,
                                           anno = anno,
                                           genes = genes.plot,
                                           grouping = "cluster",
                                           group_order = factor(group_order, levels = group_order),
                                           log_scale = FALSE,
                                           font_size = font_size,
                                           label_height = 6,
                                           show_counts = FALSE,
                                           rotate_counts = FALSE,
                                           max_width = 15
                                           ) 
    return(cluster_cell_types)
}


# group violin plot by differentiation experiment
experiment_group_violin = function(seurat_obj,
                                   group_order = c("7_10_2017", "7_13_2017", "7_20_2017", "7_24_2017", "7_27_2017"),
                                   resolution = "diff_exp",
                                   genes.plot = c("MKI67", "TNNT2", "MYH6", "MYH7")) {

    values = c("7_10_2017" = "Exp1", "7_13_2017" = "Exp2", "7_20_2017" = "Exp3", "7_24_2017" = "Exp4", "7_27_2017" = "Exp5")

    cells.plot = rownames(get(seurat_obj)@meta.data[get(seurat_obj)@meta.data$diff_exp %in% group_order, ])

    gene_data = data.frame("sample_name" = cells.plot)
    gene_df = sapply(X = genes.plot,
                     FUN = function(x) {
                         return(get(seurat_obj)@data[x, cells.plot])
                     })
    gene_data = cbind(gene_data, gene_df)

    # tail(genes.plot, n = 1)
    # length(genes.plot)

    anno = data.frame("sample_name" = cells.plot)
    anno = cbind(anno, get(seurat_obj)@meta.data[cells.plot, ])

    anno$diff_exp_id = anno$diff_exp
    anno$diff_exp_label = anno$diff_exp

    anno$cluster = anno$diff_exp
    anno = merge(anno, diff_exp_color_df, by = "diff_exp")

    cluster_cell_types = group_violin_plot(data = gene_data,
                                           anno = anno,
                                           genes = genes.plot,
                                           grouping = "diff_exp",
                                           group_order = factor(group_order, levels = group_order),
                                           log_scale = FALSE,
                                           font_size = 11,
                                           label_height = 6,
                                           show_counts = FALSE,
                                           rotate_counts = FALSE,
                                           max_width = 15)

    return(cluster_cell_types)
}


# umap highlighting one experiment at a time
plot_single_exp = function(embedding_df,
                           current_exp="7_10_2017",
                           point_size=0.5,
                           font_size=2
                          ) {

    p = ggplot(embedding_df %>% filter(diff_exp!=current_exp)) +
        geom_point(aes(x=UMAP1, y=UMAP2), size = point_size, color = "gray") +
        geom_point(data = embedding_df %>% filter(diff_exp==current_exp), aes(x=UMAP1, y=UMAP2), size=point_size, color="red") +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        guides(color = guide_legend(nrow = 1, override.aes = list(size=2))) +
        theme(legend.text = element_text(size = font_size)) + 
        theme(legend.position = "top", legend.justification = "center") +
        theme(text = element_text(size = font_size), axis.text = element_text(size = font_size))

    return(p)
}


# umap highlighting one cell line at a time
plot_single_cell_line = function(embedding_df,
                                 current_cell_line="AICS11",
                                 point_size=0.5,
                                 font_size=2
                                 ) {

    p = ggplot(embedding_df %>% filter(cell_line!=current_cell_line)) +
        geom_point(aes(x=UMAP1, y=UMAP2), size = point_size, color = "gray") +
        geom_point(data = embedding_df %>% filter(cell_line==current_cell_line), aes(x=UMAP1, y=UMAP2), size=point_size, color="red") +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        guides(color = guide_legend(nrow = 1, override.aes = list(size=2))) +
        theme(legend.text = element_text(size = font_size)) + 
        theme(legend.position = "top", legend.justification = "center") +
        theme(text = element_text(size = font_size), axis.text = element_text(size = font_size))

    return(p)
}


# umap highlighting one sample at a time
plot_single_sample = function(embedding_df,
                              sample="1",
                              point_size=0.05,
                              font_size=2
                              ) {

    p = ggplot(embedding_df %>% filter(sample_num!=sample)) +
        geom_point(aes(x=UMAP1, y=UMAP2), size = point_size, color = "gray") +
        geom_point(data = embedding_df %>% filter(sample_num==sample), aes(x=UMAP1, y=UMAP2), size=point_size, color="red") +
        theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        guides(color = guide_legend(nrow = 1, override.aes = list(size=2))) +
        theme(legend.text = element_text(size = font_size)) + 
        theme(legend.position = "top", legend.justification = "center") +
        theme(text = element_text(size = font_size), axis.text = element_text(size = font_size))

    return(p)
}
