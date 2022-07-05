source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggradar)

srtfile = {{in.srtobj | quote}}
outdir = {{out.outdir | quote}}
cases = {{envs.cases | config: "toml" | r}}

sobj = readRDS(srtfile)

mutate_meta = function(meta, mutaters) {
    if (is.null(mutaters)) {
        return(meta)
    }
    expr = list()
    for (key in names(mutaters)) {
        expr[[key]] = parse_expr(mutaters[[key]])
    }

    meta |> mutate(!!!expr)
}

do_case = function(case) {
    casepms = cases[[case]]
    outfile = file.path(outdir, paste0(slugify(case), ".png"))

    meta = mutate_meta(sobj@meta.data, casepms$mutaters)
    if (!is.null(casepms$filter)) {
        meta = meta |> filter(eval(parse(text=casepms$filter)))
    }
    meta = meta |>
        mutate(
            Cluster = if_else(
                !is.na(as.integer(seurat_clusters)),
                factor(paste0("Cluster", seurat_clusters)),
                seurat_clusters
            )
        ) |>
        group_by(Cluster, !!sym(casepms$by)) |>
        count() |>
        pivot_wider(id_cols = casepms$by, names_from = "Cluster", values_from = "n")

    meta[is.na(meta)] = 0
    if (!is.null(casepms$order)) {
        meta[[casepms$by]] = factor(meta[[casepms$by]], levels = casepms$order)
    }
    perc_with_na = casepms$perc_with_na
    if (is.null(perc_with_na)) {
        perc_with_na = FALSE
    }

    counts = if (perc_with_na) meta[, 2:ncol(meta)] else meta[!is.na(meta[[casepms$by]]), 2:ncol(meta)]
    if (!is.null(casepms$direction) && casepms$direction == "inter-cluster") {
        meta[2:ncol(meta)] = t(t(counts) / rowSums(t(counts)))
    } else {
        meta[2:ncol(meta)] = counts / rowSums(counts)
    }

    meta = meta  |> filter(!is.na(!!sym(casepms$by)))

    p = ggradar(
        meta,
        values.radar = paste0(casepms$breaks, "%"),
        grid.min = casepms$breaks[1] / 100,
        grid.mid = casepms$breaks[2] / 100,
        grid.max = casepms$breaks[3] / 100,
        plot.title = case
    )
    devpars = casepms$devpars
    if (is.null(devpars)) {
        devpars = list(res=100, width=1500, height=1000)
    }
    devpars$filename = outfile
    do.call(png, devpars)
    print(p)
    dev.off()
}

sapply(names(cases), do_case)
