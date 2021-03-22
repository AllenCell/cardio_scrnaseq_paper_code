# get sample # for each cell
get_sample_num = function(barcode) {
        sample_id = as.integer(unlist(strsplit(barcode, '_'))[[3]])
    return(paste0("E", sample_id + 1))
        }

# umap plots D0, D12, D24, D90
make_umap_subset = function(embeddings,
                            color,
                            values,
                            labels,
                            rasterize=FALSE,
                            legend="",
                            point_size=0.05,
                            font_size=13,
                            guide_size=4
                            )
{
    p.umap = ggplot(embeddings, aes_string("UMAP1", "UMAP2", color = color))

    if (rasterize) {
        p.umap = p.umap + geom_point(alpha = 1, size = point_size)
    } else {

        p.umap = p.umap + geom_point_rast(alpha = 1, size = point_size)
    }

    p.umap = p.umap + scale_color_manual(name=legend,
                                 values=values
                                 ) +
                theme_bw() +
                theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
                guides(color = guide_legend(nrow = 1, override.aes = list(size=guide_size))) +
                theme(legend.text = element_text(size = font_size)) +
                theme(legend.position = "top", legend.justification = "center") +
                theme(text = element_text(size = font_size), axis.text = element_text(size = font_size)) 

    return(p.umap)
}


make_umap_gene = function(embeddings,
                     gene,
                     label,
                     rasterize=FALSE,
                     point_size=0.5,
                     font_size=13,
                     guide_size=4,
                     legend_pos=c(0.8, 0.95),
                     legend_dir="horizontal"
                     )
{
    p.gene = ggplot(embeddings, aes_string("UMAP1", "UMAP2", color = gene))

    if (rasterize) {
        p.gene = p.gene + geom_point_rast(alpha = 1, size = point_size)
    } else {
        p.gene = p.gene + geom_point(alpha = 1, size = point_size)
        }
        
    p.gene = p.gene + theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        scale_color_gradient(name = label, low = "gray", high = "red") +
        theme(legend.text = element_text(size = font_size)) +
        theme(legend.title = element_text(size = font_size)) + 
        theme(text = element_text(size = font_size), axis.text = element_text(size = font_size)) +
        theme(legend.position = legend_pos) +
        theme(legend.direction=legend_dir)

    return(p.gene)
}


cluster_group_violin = function(seurat_obj,
                                cluster_color_df,
                                group_order = c("0", "1", "2", "3", "4", "5", "6", "7"),
                                resolution = "res.0.4",
                                genes.plot = c("MKI67", "TNNT2", "MYH6", "MYH7")
                                )
{


    cells.plot = rownames(seurat_obj@meta.data[seurat_obj@meta.data[, resolution] %in% group_order, ])

    gene_data = data.frame("sample_name" = cells.plot)
    gene_df = sapply(X = genes.plot, FUN = function(x) { return(seurat_obj@data[x, cells.plot]) })
    gene_data = cbind(gene_data, gene_df)


    tail(genes.plot, n = 1)
    length(genes.plot)

    anno = data.frame("sample_name" = cells.plot)
    anno = cbind(anno, seurat_obj@meta.data[cells.plot, ])

    anno$cluster_id = anno[, resolution]
    anno$cluster_label = sapply(X = anno[, resolution], FUN = function(x) { return(paste0("C", x))
                                                                                                                               }
                                                                                                                               )

    anno$cluster = anno[, resolution]
    anno = merge(anno, cluster_color_df, by = "cluster")

    cluster_cell_types = group_violin_plot(data = gene_data, anno = anno, genes = genes.plot, grouping = "cluster", group_order = factor(group_order, levels = group_order), log_scale = FALSE, font_size = 11, label_height = 6, show_counts = FALSE, rotate_counts = FALSE, max_width = 15) 

    return(cluster_cell_types)
}

make_heatmap = function(seurat_obj,
                        annotation_colors,
                        gene_list=c(),
                        cluster_res="res.0.4",
                        cluster_ids=c(),
                        cluster_heatmap_rows=FALSE,
                        cluster_heatmap_cols=FALSE,
                        add_annotations=c(),
                        max_exp=4
                        )
{
    keep.cells = rownames(cardio@meta.data[seurat_obj@meta.data[, cluster_res] %in% cluster_ids, ])
    top.marker.mat = seurat_obj@scale.data[unique(gene_list), keep.cells]

    marker.gene.names = rownames(top.marker.mat)
    top.marker.metdata = seurat_obj@meta.data[keep.cells, c(grouping, add_annotations)]

    top.marker.metdata$cell = rownames(top.marker.metdata)
    
    top.marker.metdata$Cluster = factor(top.marker.metdata[, cluster_res], levels = cluster_ids)

    # order by cluster
    top.marker.metdata = top.marker.metdata %>% arrange(Cluster)

    cells <- top.marker.metdata$cell
    top.marker.metdata <- top.marker.metdata[, c("Cluster", add_annotation)]
    rownames(top.marker.metdata) <- cells

    top_marker.mat.ordered <- top.marker.mat[, cells]
    top_marker.mat.ordered[top_marker.mat.ordered > max_exp] <- max_exp 

    palette_length <- 1000
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = palette_length)

    my_breaks <- c(seq(min(top_marker.mat.ordered), 0, length.out=ceiling(palette_length/2) + 1), seq(max(top_marker.mat.ordered)/palette_length, max(top_marker.mat.ordered), length.out=floor(palette_length/2))) 

    p.heatmap = pheatmap(
            top_marker.mat.ordered,
            scale = "none",
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            annotation_col = top.marker.metdata[, c("Cluster", add_annotation)],
            labels_col = top.marker.metdata$Cluster,
            cellheight = 8,
            cellwidth = 0.02,
            fontsize_row = 8,
            color = my_palette,
            breaks = my_breaks,
            annotation_names_col = FALSE,
            show_colnames = FALSE,
            annotation_colors = ann_colors,
            border_color = NA,
            fontsize = 16,
            treeheight_row = 20
            )

    return(p.heatmap)
}
