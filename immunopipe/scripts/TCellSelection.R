library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(plotthis)
library(scplotter)
library(biopipen.utils)

srtfile = {{in.srtobj | r}}
immfile = {{in.immdata | r}}
outdir = {{out.outdir | r}}
joboutdir = {{job.outdir | r}}
outfile = {{out.outfile | r}}
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

log <- get_logger()
reporter <- get_reporter()

log$info("Reading seurat object for scRNA-seq data ...")
sobj = read_obj(srtfile)
if (ignore_tcr || is.null(ignore_tcr) || is.null(immfile)) {
    log$info("Ignoring TCR information ...")
    immdata <- NULL
} else {
    log$info("Reading scTCR-seq data ...")
    immdata <- read_obj(immfile)
    immdata <- do_call(rbind, immdata)
}

# Cluster   CD3E   ...
# 0         1.1    ...
# 1         0      ...
# ...
log$info("Fetching indicator gene expression ...")
assay <- DefaultAssay(sobj)
indicators = AverageExpression(sobj, features = indicator_genes, assays = assay)[[assay]] %>%
    # It's a sparse matrix, convert it to a dense matrix
    as.matrix() %>%
    t() %>%
    # scale() %>%  # scale on each gene
    as.data.frame()

if (ncol(indicators) == 1) {
    colnames(indicators) = indicator_genes[1]
}

# colnames(indicators) = indicator_genes
indicators = indicators %>% rownames_to_column("Cluster") %>%
    mutate(Cluster = sub("^RNA\\.", "", Cluster))

cluster_sizes = table(Idents(sobj))[indicators$Cluster]
if (all(is.na(cluster_sizes)) && startsWith(indicators$Cluster[1], "g")) {
    log$warn("- Cluster Idents mismatch, try without prefix 'g'")
    indicators$Cluster <- sub("^g", "", indicators$Cluster)
    cluster_sizes = table(Idents(sobj))[indicators$Cluster]
}
indicators$Cluster_Size = as.numeric(cluster_sizes)

if (!is.null(immdata)) {
    log$info("Fetching clonotype percentage for each cluster ...")
    tcells = immdata$barcode %||% immdata[, 1, drop = TRUE]

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
    The barcode from scTCR-seq data is like: {tcells[1]}")
        stop(msg)
    }

    # Cluster   CD3E   ...   Clonotype_Pct
    # 0         1.1    ...   0.8
    # 1         0      ...   0.1
    # ...
    indicators$Clonotype_Pct <- unlist(clonotype.pct[indicators$Cluster])
    indicators$Clonotype_Pct[is.na(indicators$Clonotype_Pct)] <- 0
}

km_plot = NULL
if (!is.null(tcell_selector) || isFALSE(tcell_selector)) {
    # mutate_code = paste0("mutate(indicators, is_TCell = ", tcell_selector, ")")
    # indicators = eval(parse(text = mutate_code))
    log$info("Using user-defined T cell selector ...")
    indicators = indicators %>% mutate(is_TCell = !!parse_expr(tcell_selector))
} else {
    # Use k-means to determine T cell clusters
    # based on the indicators and clonotype percentage
    log$info("Using k-means to determine T cell clusters ...")

    km_df = select(indicators, -c("Cluster", "Cluster_Size"))
    kmeans_args$x = km_df
    kmeans_args$centers = 2
    km = do_call(kmeans, kmeans_args)
    indicators = indicators %>% mutate(km_cluster = km$cluster)
    tcell_kmcluster = indicators %>% group_by(km_cluster)

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

    pca <- prcomp(km_df, scale. = TRUE)
    pca <- as.data.frame(pca$x[, 1:2])
    pca$group <- factor(km$cluster, levels = c("T cell", "Other"))

    km_plot <- plotthis::DimPlot(
        pca, group_by = "group",
        add_mark = TRUE, mark_type = "ellipse", palette = "Set1"
    )

    png(
        file.path(outdir, "kmeans.png"),
        res = 100,
        height = attr(km_plot, "height") * 100,
        width = attr(km_plot, "width") * 100
    )
    print(km_plot)
    dev.off()

    reporter$add(
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
log$info("Writing indicator gene expression (and clonotype percentage) to file ...")
write.table(
    indicators,
    file.path(outdir, 'data.txt'),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

reporter$add(
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
    log$info("- Gene: {gene}")
    if (gene == x) { return() }
    x_lab = ifelse(x == "Clonotype_Pct", "clonotype percentage", paste0(x, " expression"))

    p <- plotthis::ScatterPlot(
        indicators,
        x = x,
        y = gene,
        color_by = "is_TCell",
        size_by = "Cluster_Size",
        palette = "Set1"
    )

    plotfile = file.path(outdir, paste0(slugify(gene), "-vs-", x, ".png"))
    png(plotfile, res=100, height=attr(p, "height") * 100, width=attr(p, "width") * 100)
    print(p)
    dev.off()

    reporter$add(
        list(src = plotfile, name = gene),
        h1 = paste0("Indicator gene expression vs ", x_lab),
        ui = "table_of_images"
    )
}

log$info("Plotting indicator gene expression vs clonotype percentage/first gene ...")


if (length(indicator_genes) > 0) {
    reporter$add(
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
    indicators$is_TCell <- factor(
        indicators$is_TCell,
        levels = c("TRUE", "FALSE"),
        labels = c("T cell", "Other")
    )
    sapply(
        indicator_genes,
        plot_indicator_gene,
        x = ifelse(!is.null(indicators$Clonotype_Pct), "Clonotype_Pct", indicator_genes[1])
    )
}

# Plot the percentage of T cells in each sample
log$info("Plotting (selected) T cell composition per sample ...")
{
    plot_data <- sobj@meta.data %>%
        left_join(indicators, by=c(seurat_clusters = "Cluster")) %>%
        group_by(Sample, is_TCell) %>%
        summarise(nCells = n(), .groups = "drop")

    p_num <- plotthis::BarPlot(
        plot_data,
        x = "Sample",
        y = "nCells",
        group_by = "is_TCell",
        palette = "Set1",
        position = "stack",
        x_text_angle = 90
    )
    p_frac <- plotthis::BarPlot(
        plot_data %>%
            group_by(Sample) %>%
            mutate(percCells = nCells / sum(nCells) * 100) %>%
            ungroup(),
        x = "Sample",
        y = "percCells",
        group_by = "is_TCell",
        palette = "Set1",
        position = "stack",
        x_text_angle = 90
    )
    p <- plotthis:::combine_plots(list(p_num, p_frac), ncol = 2)

    png(
        file.path(outdir, "tcell_per_sample.png"),
        res = 100,
        height = attr(p, "height") * 100,
        width = attr(p, "width") * 100
    )
    print(p)
    dev.off()
}

# Plot a pie chart of the T cell composition
log$info("Plotting T cell composition ...")
{
    p_pie <- plotthis::PieChart(
        plot_data %>%
            group_by(is_TCell) %>%
            summarise(nCells = sum(nCells), .groups = "drop"),
        x = "is_TCell",
        y = "nCells",
        palette = "Set1"
    )
    png(
        file.path(outdir, "tcell_pie.png"),
        res = 100,
        height = attr(p_pie, "height") * 100,
        width = attr(p_pie, "width") * 100
    )
    print(p_pie)
    dev.off()
}

reporter$add(
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
log$info("Saving selected T cells as Seurat object ...")
outobj = subset(sobj, idents = indicators$Cluster[as.character(indicators$is_TCell) == "T cell"])

save_obj(outobj, outfile)

# feature plots
log$info("Plotting feature plots ...")
reporter$add(
    list(
        kind = "descr",
        content = "Showing feature plots of indicator genes and selected T cells. "
    ),
    h1 = "Feature plots"
)
p <- scplotter::FeatureStatPlot(
    sobj,
    features = indicator_genes,
    plot_type = "dim"
)
png(
    file.path(outdir, "feature-plots.png"),
    res = 100,
    height = attr(p, "height") * 100,
    width = attr(p, "width") * 100
)
print(p)
dev.off()
reporter$add(
    list(kind = "image", src = file.path(outdir, "feature-plots.png"), name = "Feature plots"),
    h1 = "Feature plots"
)

# Plot dim plots for TCR and selected T cells
log$info("Plotting dim plots for cells with TCR data and final selected T cells ...")

tcr_plots <- list()
if (!is.null(immdata)) {
    sobj@meta.data$TCR = factor(
        rownames(sobj@meta.data) %in% tcells,
        levels = c(TRUE, FALSE),
        labels = c("TCR", "No TCR")
    )
    sobj@meta.data$Clonotype_Pct = sobj@meta.data %>% left_join(
        indicators,
        by = c(seurat_clusters = "Cluster")
    ) %>% pull(Clonotype_Pct)

    p1 <- scplotter::FeatureStatPlot(
        sobj,
        features = "Clonotype_Pct",
        plot_type = "dim",
        subtitle = "Clonotype Percentage"
    )

    p2 <- scplotter::CellDimPlot(
        sobj,
        group_by = "TCR",
        palette = "Set1"
    )

    tcr_plots <- list(p1, p2)
}

sobj@meta.data$Selected_TCells = sobj@meta.data %>% left_join(
    indicators,
    by = c(seurat_clusters = "Cluster")
) %>% pull(is_TCell)

p3 <- scplotter::CellDimPlot(
    sobj,
    group_by = "Selected_TCells",
    palette = "Set1"
)

tcr_plots <- c(tcr_plots, list(p3))
p <- plotthis:::combine_plots(tcr_plots, ncol = 2)
png(
    file.path(outdir, "tcr-selected-tcells.png"),
    res = 100,
    height = attr(p, "height") * 100,
    width = attr(p, "width") * 100
)
print(p)
dev.off()
reporter$add(
    list(kind = "image", src = file.path(outdir, "tcr-selected-tcells.png"),
        name = "TCR and selected T cells"),
    h1 = "TCR and selected T cells"
)

reporter$save(joboutdir)
