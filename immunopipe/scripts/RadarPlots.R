library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggradar)

srtfile = {{in.srtobj | quote}}
groupfile = {{in.groupfile | r}}
direction = {{in.direction | r}}
breaks = as.integer({{in.breaks | split: "," | r}})
outfile = {{out.outfile | quote}}

sobj = readRDS(srtfile)
allgroups = read.table(
    groupfile,
    row.names=1,
    header=T,
    sep="\t",
    check.names = F
) |> select(last_col())

group = colnames(allgroups)
cells = intersect(rownames(sobj@meta.data), rownames(allgroups))
data = cbind(sobj@meta.data[cells, ], allgroups[cells,, drop=FALSE])

data = data |>
    mutate(
        Cluster = if_else(
            !is.na(as.integer(seurat_clusters)),
            factor(paste0("Cluster", seurat_clusters)),
            seurat_clusters
        )
    ) |>
    group_by(Cluster, !!sym(group)) |>
    count() |>
    pivot_wider(id_cols = group, names_from = "Cluster", values_from = "n")

counts = data[2:ncol(data)]
if (direction == "inter-cluster") {
    data[2:ncol(data)] = t(t(counts) / rowSums(t(counts)))
} else {
    data[2:ncol(data)] = counts / rowSums(counts)
}
p = ggradar(
    data,
    values.radar = paste0(breaks, "%"),
    grid.min = breaks[1] / 100,
    grid.mid = breaks[2] / 100,
    grid.max = breaks[3] / 100
)
png(outfile, res=100, width=1500, height=1000)
print(p)
dev.off()
