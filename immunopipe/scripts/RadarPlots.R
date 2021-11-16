library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggradar)

srtfile = {{in.srtobj | quote}}
groupfiles = {{in.groupfiles | r}}
direction = {{in.direction | r}}
breaks = as.integer({{in.breaks | split: "," | r}})
outfile = {{out.outfile | quote}}

sobj = readRDS(srtfile)
allgroups = NULL
for (groupfile in groupfiles) {
    groups = read.table(groupfile, row.names=NULL, header=T, sep="\t", check.names = F)
    n_samples = ncol(groups) - 1
    df = groups %>% rowwise() %>%
        mutate(
            across(2:(n_samples+1),
            ~ sapply(
                strsplit(.x, ";", fixed=TRUE),
                function(x) paste(cur_column(), x, sep="_", collapse=";")
            )
        )) %>%
        unite("ALL", 2:(n_samples+1)) %>%
        rename(group = 1)
    if (is.null(allgroups)) {
        allgroups = df
    } else {
        allgroups = bind_rows(allgroups, df)
    }
}

for (ident in unique(Idents(sobj))) {
    cells = WhichCells(sobj, ident = ident)
    if (!is.na(as.integer(ident))) {
        cluster = paste0("Cluster", ident)
    } else {
        cluster = ident
    }
    allgroups = allgroups %>% rowwise() %>% mutate(
        .x = length(intersect(unlist(strsplit(ALL, ";")), cells))
    )
    colnames(allgroups)[ncol(allgroups)] = cluster
}

allgroups = allgroups %>% select(-"ALL") %>% arrange(group)
counts = allgroups[2:ncol(allgroups)]
if (direction == "inter-cluster") {
    allgroups[2:ncol(allgroups)] = t(t(counts) / rowSums(t(counts)))
} else {
    allgroups[2:ncol(allgroups)] = counts / rowSums(counts)
}
p = ggradar(
    allgroups,
    values.radar = paste0(breaks, "%"),
    grid.min = breaks[1] / 100,
    grid.mid = breaks[2] / 100,
    grid.max = breaks[3] / 100
)
png(outfile, res=100, width=1500, height=1000)
print(p)
dev.off()


