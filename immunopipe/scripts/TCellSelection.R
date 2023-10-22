library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)
library(ggplot2)
library(ggprism)
library(ggrepel)
library(ggpubr)
library(factoextra)

srtfile = {{in.srtobj | quote}}
immfile = {{in.immdata | quote}}
outdir = {{out.outdir | quote}}
rdsfile = {{out.rdsfile | quote}}
indicator_genes = {{envs.indicator_genes | r}}
tcell_indicator = {{envs.tcell_indicator | r}}

sobj = readRDS(srtfile)
immdata = readRDS(immfile)

# Cluster   CD3E   ...
# 0         1.1    ...
# 1         0      ...
# ...
indicators = AverageExpression(sobj, features = indicator_genes, assays = "RNA") %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()

colnames(indicators) = indicator_genes
indicators = indicators %>% rownames_to_column("Cluster") %>%
    mutate(Cluster = sub("^RNA\\.", "", Cluster))

cluster_sizes = table(Idents(sobj))[indicators$Cluster]
indicators$Cluster_Size = as.numeric(cluster_sizes)

tcells = c()
for (sample in names(immdata$data)) {
    cells = immdata$data[[sample]]$Barcode %>% strsplit(";") %>% unlist()
    cells = paste(sample, cells, sep="_")
    tcells = c(tcells, cells)
}

clonotype_pct_ident = function(ident) {
    ident_cells = WhichCells(sobj, idents = ident)
    has_clonotypes = table(ident_cells %in% tcells)
    unname(has_clonotypes["TRUE"] / sum(has_clonotypes))
}

clonotype.pct = list()
for (ident in unique(Idents(sobj))) {
    clonotype.pct[[ident]] = clonotype_pct_ident(ident)
}

# Cluster   CD3E   ...   Clonotype_Pct
# 0         1.1    ...   0.8
# 1         0      ...   0.1
# ...
indicators = add_column(
    indicators,
    Clonotype_Pct = unlist(clonotype.pct[indicators$Cluster])
) %>% replace_na(list(Clonotype_Pct = 0))

if (!is.null(tcell_indicator) || isFALSE(tcell_indicator)) {
    mutate_code = paste0("mutate(indicators, is_TCell = ", tcell_indicator, ")")
    indicators = eval(parse(text = mutate_code))
} else {
    # Use k-means to determine T cell clusters
    # based on the indicators and clonotype percentage
    km_df = select(indicators, -c("Cluster", "Cluster_Size"))
    km = kmeans(km_df, 2)
    indicators = indicators %>% mutate(km_cluster = km$cluster)
    tcell_kmcluster = indicators %>%
        group_by(km_cluster) %>%
        # Summarise the clonotype percentage for each km cluster
        summarise(km_cluster_clono_pct = mean(Clonotype_Pct)) %>%
        arrange(desc(km_cluster_clono_pct)) %>%
        # The clusters with the highest clonotype percentage
        # are likely to be T cells
        head(1) %>%
        pull(km_cluster) %>%
        unname()
    indicators = indicators %>% mutate(is_TCell = km_cluster == tcell_kmcluster)

    # Plot the k-means clusters
    km_plot = fviz_cluster(km, data = km_df, ellipse.type = "convex") +
        labs(title = "K-means clustering of T cell indicators") +
        theme_prism(base_size = 16) +
        theme(legend.title = element_text(size=12))
    png(
        file.path(outdir, "kmeans.png"),
        res=70,
        height=600,
        width=800
    )
    print(km_plot)
    dev.off()
}

write.table(
    indicators,
    file.path(outdir, 'data.txt'),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

# Plot the indicator gene expression and
# clonotype percentage and mark the T cell clusters
plot_indicator_gene = function(gene) {
    p = ggplot(indicators, aes(x=Clonotype_Pct, y=!!sym(gene))) +
        geom_point(aes(color=is_TCell, size=Cluster_Size), shape=19) +
        geom_label_repel(
            aes(label=Cluster),
            box.padding   = 0.35,
            point.padding = 0.5,
            segment.color = 'grey50'
        ) +
        labs(x="Clonotype Percentage", y=paste0(gene, " mean expression")) +
        theme_prism(base_size = 16) +
        theme(legend.title = element_text(size=12))

    png(
        file.path(outdir, paste0(gene, "-vs-clonopct.png")),
        res=100,
        height=600,
        width=1200
    )
    print(p)
    dev.off()
}

sapply(indicator_genes, plot_indicator_gene)

out = subset(sobj, idents = indicators$Cluster[indicators$is_TCell])
saveRDS(out, rdsfile)
