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
    print("Processing case: " %>% paste(case))
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

    meta = meta  |> filter(!is.na(!!sym(casepms$by)))
    counts = if (perc_with_na) meta[, 2:ncol(meta)] else meta[!is.na(meta[[casepms$by]]), 2:ncol(meta)]
    if (!is.null(casepms$direction) && casepms$direction == "inter-cluster") {
        meta[2:ncol(meta)] = t(t(counts) / rowSums(t(counts)))
    } else {
        meta[2:ncol(meta)] = counts / rowSums(counts)
    }

    # group            mpg   cyl  disp    hp  drat    wt   qsec    vs    am
    # <chr>          <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>
    # Ford Pantera L 0.230   1   0.698 0.749 0.673 0.424 0          0     1
    # Ferrari Dino   0.396   0.5 0.184 0.435 0.396 0.321 0.119      0     1
    # Maserati Bora  0.196   1   0.573 1     0.359 0.526 0.0119     0     1
    # Volvo 142E     0.468   0   0.124 0.201 0.622 0.324 0.488      1     1
    if (casepms$breaks == "auto" || is.null(casepms$breaks)) {
        maxval = max(meta[, 2:ncol(meta)])
        if (maxval <= 0.1) {  # 10%
            casepms$breaks = c(0, 5, 10)
        } else if (maxval <= 0.2) {
            casepms$breaks = c(0, 10, 20)
        } else if (maxval <= 0.3) {
            casepms$breaks = c(0, 15, 30)
        } else if (maxval <= 0.4) {
            casepms$breaks = c(0, 20, 40)
        } else if (maxval <= 0.5) {
            casepms$breaks = c(0, 25, 50)
        } else if (maxval <= 0.6) {
            casepms$breaks = c(0, 30, 60)
        } else if (maxval <= 0.7) {
            casepms$breaks = c(0, 35, 70)
        } else if (maxval <= 0.8) {
            casepms$breaks = c(0, 40, 80)
        } else if (maxval <= 0.9) {
            casepms$breaks = c(0, 45, 90)
        } else {
            casepms$breaks = c(0, 50, 100)
        }
    }
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
