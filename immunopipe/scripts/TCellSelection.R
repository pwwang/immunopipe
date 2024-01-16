source("{{biopipen_dir}}/utils/misc.R")
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)
library(ggplot2)
library(ggprism)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(factoextra)
library(logger)
library(glue)

srtfile = {{in.srtobj | r}}
immfile = {{in.immdata | r}}
outdir = {{out.outdir | r}}
joboutdir = {{job.outdir | r}}
rdsfile = {{out.rdsfile | r}}
ignore_tcr = {{envs.ignore_tcr | r}}
indicator_genes = {{envs.indicator_genes | r}}
tcell_selector = {{envs.tcell_selector | r}}
kmeans_args = {{envs.kmeans | r: todot = "-"}}

# Check if we have enough genes for k-means clustering
if (
    (is.null(tcell_selector) || length(tcell_selector) == 0) &&
    length(indicator_genes) < 2 &&
    isTRUE(ignore_tcr)
) {
    stop("You need at least 2 markers to perform k-means clustering with TCR being ignored.")
}

log_info("Reading data from input files ...")
sobj = readRDS(srtfile)
if (ignore_tcr || is.null(ignore_tcr) || is.null(immfile)) {
    log_info("Ignoring TCR information ...")
    immdata <- NULL
} else {
    immdata <- read.table(immfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
}

# Cluster   CD3E   ...
# 0         1.1    ...
# 1         0      ...
# ...
log_info("Fetching indicator gene expression ...")
assay <- DefaultAssay(sobj)
indicators = AverageExpression(sobj, features = indicator_genes, assays = assay)[[assay]] %>%
    # It's a sparse matrix, convert it to a dense matrix
    as.matrix() %>%
    t() %>%
    # scale() %>%  # scale on each gene
    as.data.frame()

# colnames(indicators) = indicator_genes
indicators = indicators %>% rownames_to_column("Cluster") %>%
    mutate(Cluster = sub("^RNA\\.", "", Cluster))

cluster_sizes = table(Idents(sobj))[indicators$Cluster]
indicators$Cluster_Size = as.numeric(cluster_sizes)

if (!is.null(immdata)) {
    log_info("Fetching clonotype percentage for each cluster ...")
    tcells = rownames(immdata)

    clonotype_pct_ident = function(ident) {
        ident_cells = WhichCells(sobj, idents = ident)
        has_clonotypes = table(ident_cells %in% tcells)
        unname(has_clonotypes["TRUE"] / sum(has_clonotypes))
    }

    clonotype.pct = list()
    for (ident in unique(Idents(sobj))) {
        clonotype.pct[[ident]] = clonotype_pct_ident(ident)
        if (is.na(clonotype.pct[[ident]])) {
            clonotype.pct[[ident]] = 0
        }
    }

    if (sum(unlist(clonotype.pct)) == 0) {
        msg = glue("No clonotype information found in the Seurat object

    The barcode information from the Seurat object does not match the barcode information from the scTCR-seq data
    The barcode from scRNA-seq data is like: {rownames(sobj@meta.data)[1]}
    The barcode from scTCR-seq data is like: {tcells[1]}

    Have you specified correct `envs.prefix` for `ImmunarchLoading` process?")
        stop(msg)
    }

    # Cluster   CD3E   ...   Clonotype_Pct
    # 0         1.1    ...   0.8
    # 1         0      ...   0.1
    # ...
    indicators = add_column(
        indicators,
        Clonotype_Pct = unlist(clonotype.pct[indicators$Cluster])
    ) %>% replace_na(list(Clonotype_Pct = 0))
}

km_plot = NULL
if (!is.null(tcell_selector) || isFALSE(tcell_selector)) {
    # mutate_code = paste0("mutate(indicators, is_TCell = ", tcell_selector, ")")
    # indicators = eval(parse(text = mutate_code))
    log_info("Using user-defined T cell selector ...")
    indicators = indicators %>% mutate(is_TCell = !!parse_expr(tcell_selector))
} else {
    # Use k-means to determine T cell clusters
    # based on the indicators and clonotype percentage
    log_info("Using k-means to determine T cell clusters ...")

    km_df = select(indicators, -c("Cluster", "Cluster_Size"))
    kmeans_args$x = km_df
    kmeans_args$centers = 2
    km = do_call(kmeans, kmeans_args)
    indicators = indicators %>% mutate(km_cluster = km$cluster)
    tcell_kmcluster = indicators %>%
        group_by(km_cluster)

    if (!is.null(indicators$Clonotype_Pct)) {
        # Summarise the clonotype percentage for each km cluster
        tcell_kmcluster = tcell_kmcluster %>%
            summarise(km_cluster_clono_pct = mean(Clonotype_Pct)) %>%
            arrange(desc(km_cluster_clono_pct)) %>%
            # The clusters with the highest clonotype percentage
            # are likely to be T cells
            head(1) %>%
            pull(km_cluster) %>%
            unname()
    } else {
        tcell_kmcluster = tcell_kmcluster %>%
            summarise(km_pos_gene_expr = mean(!!sym(indicator_genes[1]))) %>%
            arrange(desc(km_pos_gene_expr)) %>%
            # The clusters with the highest expression of the first
            # indicator gene are likely to be T cells
            head(1) %>%
            pull(km_cluster) %>%
            unname()
    }
    indicators = indicators %>% mutate(is_TCell = km_cluster == tcell_kmcluster)
    tcell_kmcluster_id = indicators %>%
        filter(is_TCell) %>% slice_head(n = 1) %>% pull(km_cluster)
    if (tcell_kmcluster_id == 1) {
        km$cluster = c("T cell", "Other")[km$cluster]
    } else {
        km$cluster = c("Other", "T cell")[km$cluster]
    }

    # Plot the k-means clusters
    km_plot = fviz_cluster(km, data = km_df, ellipse.type = "convex") +
        labs(title = "K-means clustering of T cell indicators") +
        theme_prism(base_size = 16) +
        theme(legend.title = element_text(size=12)) +
        scale_color_manual(values=pal_biopipen()(2), name="cluster")
    png(
        file.path(outdir, "kmeans.png"),
        res=70,
        height=600,
        width=800
    )
    print(km_plot)
    dev.off()

    add_report(
        list(
            kind = "descr",
            content = ifelse(
                is.null(indicators$Clonotype_Pct),
                paste0(
                    "When no `tcell_selector` specified, a k-means clustering is done with `n=2` ",
                    "using the expressions of indicator genes for each cluster. ",
                    "Then the mean expression of ",
                    indicator_genes[1],
                    " is calculated for each kmeans cluster and ",
                    "the cluster with the highest mean expression is selected, ",
                    "in which the Seurat clusters are assigned as T cells."
                ),
                paste0(
                    "When no `tcell_selector` specified, a k-means clustering is done with `n=2` ",
                    "using the expressions of indicator genes and clonotype percentage for each cluster. ",
                    "Then the mean clonotype percentage is calculated for each kmeans cluster and ",
                    "the cluster with the highest mean clonotype percentage is selected, ",
                    "in which the Seurat clusters are assigned as T cells."
                )
            )
        ),
        list(
            kind = "image",
            src = file.path(outdir, "kmeans.png")
        ),
        h1 = "K-means clustering"
    )
}

# Write the indicator gene expression and
# clonotype percentage to a file
log_info("Writing indicator gene expression (and clonotype percentage) to file ...")
write.table(
    indicators,
    file.path(outdir, 'data.txt'),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

add_report(
    list(
        kind = "descr",
        content = paste0(
            "The data that is used to determine the T-cell clusters, ",
            "including the gene expression values for indicator genes",
            ifelse(is.null(indicators$Clonotype_Pct), ". ", " and the clonotype percentage. "),
            "The cluster size and the final determined T-cell flags. ",
            "The expression values are normalized with mean 0 and standard deviation 1."
        )
    ),
    list(
        kind = "table",
        src = file.path(outdir, "data.txt")
    ),
    h1 = "Table of indicator gene expression and clonotype percentage"
)

# Plot the indicator gene expression and
# clonotype percentage and mark the T cell clusters
plot_indicator_gene = function(gene, x) {
    log_info("- Gene: {gene}")
    if (gene == x) { return() }
    x_lab = ifelse(x == "Clonotype_Pct", "clonotype percentage", paste0(x, " expression"))
    p = ggplot(indicators, aes(x=!!sym(x), y=!!sym(gene))) +
        geom_point(aes(color=is_TCell, size=Cluster_Size), shape=19) +
        geom_label_repel(
            aes(label=Cluster),
            box.padding   = 0.35,
            point.padding = 0.5,
            segment.color = 'grey50'
        ) +
        labs(x=x_lab, y=paste0(gene, " mean expression")) +
        theme_prism(base_size = 16) +
        theme(legend.title = element_text(size=12)) +
        scale_color_manual(values=pal_biopipen()(2), name="is_TCell")

    plotfile = file.path(outdir, paste0(slugify(gene), "-vs-", x, ".png"))
    png(plotfile, res=100, height=600, width=1200)
    print(p)
    dev.off()

    add_report(
        list(src = plotfile, name = gene),
        h1 = paste0("Indicator gene expression vs ", x_lab),
        ui = "table_of_images"
    )
}

log_info("Plotting indicator gene expression vs clonotype percentage/first gene ...")


if (length(indicator_genes) > 0) {
    add_report(
        list(
            kind = "descr",
            content = ifelse(
                is.null(indicators$Clonotype_Pct),
                paste0(
                    "Scatter plots showing the average gene expression against expression of ",
                    indicator_genes[1],
                    " for each cluster. ",
                    "The x-axis is the expression of ",
                    indicator_genes[1],
                    " The y-axis is the average gene expression of the given gene. ",
                    "The expression values are normalized with mean 0 and standard deviation 1."
                ),
                paste0(
                    "Scatter plots showing the average gene expression against the clonotype percentage for each cluster. ",
                    "The x-axis is the clonotype percentage, which is the percentage of the cells with clonotype information detected by scTCR-seq data. ",
                    "The y-axis is the average gene expression of the given gene. ",
                    "The expression values are normalized with mean 0 and standard deviation 1."
                )
            )
        ),
        h1 = paste0("Indicator gene expression vs ", ifelse(
            is.null(indicators$Clonotype_Pct),
            indicator_genes[1],
            "clonotype percentage"
        )),
        ui = "flat"
    )

    sapply(
        indicator_genes,
        plot_indicator_gene,
        x = ifelse(!is.null(indicators$Clonotype_Pct), "Clonotype_Pct", indicator_genes[1])
    )
}

# Plot the percentage of T cells in each sample
log_info("Plotting (selected) T cell composition per sample ...")
{
    plot_data <- sobj@meta.data %>%
        left_join(indicators, by=c(seurat_clusters = "Cluster")) %>%
        group_by(Sample, is_TCell) %>%
        summarise(nCells = n(), .groups = "drop") %>%
        mutate(is_TCell = factor(ifelse(is_TCell, "T cell", "Other"), levels=c("T cell", "Other")))

    p <- ggplot(
            plot_data %>% mutate(is_TCell = factor(is_TCell, levels = c("Other", "T cell"))),
            aes(x=Sample, y=nCells, fill=is_TCell)
        ) +
        geom_bar(stat="identity", position="stack") +
        labs(x="Sample", y="Number of cells") +
        theme_prism() +
        # rotate x-axis labels
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(values=rev(pal_biopipen()(2)), name="is_TCell")
    p_perc <- ggplot(
            plot_data %>%
                group_by(Sample) %>%
                mutate(percCells = nCells / sum(nCells) * 100) %>%
                ungroup() %>%
                mutate(is_TCell = factor(is_TCell, levels = c("Other", "T cell"))),
            aes(x=Sample, y=percCells, fill=is_TCell)
        ) +
        geom_bar(stat="identity", position="stack") +
        labs(x="Sample", y="Percentage of cells") +
        theme_prism() +
        # rotate x-axis labels
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        scale_fill_manual(values=rev(pal_biopipen()(2)), name="is_TCell")

    png(
        file.path(outdir, "tcell_per_sample.png"),
        res=100,
        height=600,
        width=2 * (length(unique(sobj@meta.data$Sample)) * 50 + 100)
    )
    print(p | p_perc)
    dev.off()
}

# Plot a pie chart of the T cell composition
log_info("Plotting T cell composition ...")
{
    p_pie = ggplot(
            plot_data %>%
                group_by(is_TCell) %>%
                summarise(nCells = sum(nCells), .groups = "drop") %>%
                mutate(
                    percCells = nCells / sum(nCells) * 100,
                    percCells = paste0(
                        format(nCells, big.mark=","),
                        " (", round(percCells, 1), "%)"
                    )
                ),
            aes(x="", y=nCells, fill=is_TCell, label=percCells)
        ) +
        geom_bar(width = 1, stat = "identity", position = position_stack(reverse = TRUE)) +
        coord_polar("y", start=0) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_label_repel(
            position = position_stack(vjust = 0.5),
            color = "#333333",
            fill = "#eeeeee",
            size = 4
        ) +
        scale_fill_manual(values=pal_biopipen()(2), name = "")
    png(
        file.path(outdir, "tcell_pie.png"),
        res=100,
        height=600,
        width=800
    )
    print(p_pie)
    dev.off()
}

add_report(
    list(
        src = file.path(outdir, "tcell_pie.png"),
        name = "T cells in pie charts",
        desc = "The number and percentage of T cells and other cells in pie charts"
    ),
    list(
        src = file.path(outdir, "tcell_per_sample.png"),
        name = "T cells in each sample",
        desc = "Absolute number of T cells and other cells in each sample (left) and percentage of T cells and other cells in each sample (right)."
    ),
    h1 = "T cell composition",
    ui = "table_of_images"
)

# Save first. Do not save the other features for visualization to the Seurat object
log_info("Saving selected T cells as Seurat object ...")
outobj = subset(sobj, idents = indicators$Cluster[indicators$is_TCell])

saveRDS(outobj, rdsfile)

# feature plots
log_info("Plotting feature plots ...")
add_report(
    list(
        kind = "descr",
        content = "Showing feature plots of indicator genes and selected T cells. "
    ),
    h1 = "Feature plots"
)
# adding TCR information
features = indicator_genes
if (!is.null(immdata)) {
    features = c(features, "Clonotype_Pct", "TCR")
    sobj@meta.data$TCR = rownames(sobj@meta.data) %in% rownames(immdata)
    sobj@meta.data$Clonotype_Pct = sobj@meta.data %>% left_join(
        indicators,
        by = c(seurat_clusters = "Cluster")
    ) %>% pull(Clonotype_Pct)
}
features = c(features, "Selected_TCells")
sobj@meta.data$Selected_TCells = sobj@meta.data %>% left_join(
    indicators,
    by = c(seurat_clusters = "Cluster")
) %>% pull(is_TCell)
for (feature in features) {
    log_info("- Feature: {feature}")
    plotfile = file.path(outdir, paste0("feature_", slugify(feature), ".png"))
    png(plotfile, res=100, height=600, width=800)
    p <- FeaturePlot(sobj, features = feature)
    print(p)
    dev.off()

    add_report(
        list(src = plotfile, name = feature),
        h1 = "Feature plots",
        ui = "table_of_images"
    )
}

save_report(joboutdir)
