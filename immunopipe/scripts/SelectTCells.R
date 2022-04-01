library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)
library(ggrepel)

srtfile = {{in.srtobj | quote}}
immfile = {{in.immdata | quote}}
outdir = {{out.outdir | quote}}
rdsfile = {{out.rdsfile | quote}}
indicator_gene = {{envs.indicator_gene | quote}}

sobj = readRDS(srtfile)
immdata = readRDS(immfile)

cd3e.expr = as.matrix(sobj@assays$RNA[indicator_gene,])
cd3e.means = lapply(split(cd3e.expr, Idents(sobj)), mean)

tcells = c()
for (sample in names(immdata$data)) {
    cells = immdata$data[[sample]]$Barcode %>% strsplit(";") %>% unlist()
    cells = paste(sample, cells, sep="_")
    tcells = c(tcells, cells)
}

clonotype_pct_ident = function(ident) {
    ident_cells = WhichCells(sobj, idents = ident)
    has_clonotypes = table(ident_cells %in% tcells)
    has_clonotypes["TRUE"] / sum(has_clonotypes)
}

clonotype.pct = list()
for (ident in unique(Idents(sobj))) {
    clonotype.pct[[ident]] = clonotype_pct_ident(ident)
}

df = bind_rows(
    as.data.frame(cd3e.means),
    as.data.frame(clonotype.pct)
)

rownames(df) = c(paste0(indicator_gene, "_means"), "Clonotype_pct")
df = t(df) %>% as.data.frame() %>%
    rownames_to_column("ClusterID") %>%
    mutate(
        ClusterID = sub("X", "", ClusterID),
        Cluster = paste0("Cluster", ClusterID)
    )
tcell_clusters = df %>% filter({{envs.tcell_filter}}) %>% pull(ClusterID)
df = df %>% mutate(
    Group = if_else(ClusterID %in% tcell_clusters, "T Cell", "Non-T Cell")
)
write.table(
    df,
    file.path(outdir, 'clusters.txt'),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)
png(
    file.path(outdir, "subdivision-rationale.png"),
    res=100,
    height=600,
    width=1200
)

ggplot(df, aes_string(x="Clonotype_pct", y=paste0(indicator_gene, "_means"))) +
    geom_point(aes(color=Group)) +
    geom_label_repel(
        aes(label=Cluster),
        box.padding   = 0.35,
        point.padding = 0.5,
        segment.color = 'grey50'
    ) +
    theme_prism(base_size = 16)

dev.off()

out = subset(sobj, idents = tcell_clusters)
saveRDS(out, rdsfile)
