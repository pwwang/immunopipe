source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(rlang)
library(dplyr)
library(tibble)
library(parallel)
library(Seurat)
library(ggprism)
library(tidyseurat)

srtfile = {{in.srtobj | r}}
cases = {{envs.cases | r}}
outdir = {{out.outdir | r}}
ncores = {{envs.ncores | r}}


print("- Reading seurat object ...")
srtobj = readRDS(srtfile)


mutate_meta = function(obj, mutaters) {
    meta = obj@meta.data
    if (!is.null(mutaters)) {
        expr = list()
        for (key in names(mutaters)) {
            expr[[key]] = parse_expr(mutaters[[key]])
        }
        obj@meta.data = meta |> mutate(!!!expr)
    }
    return(obj)
}


do_case = function(case) {
    casepms = cases[[case]]

    if (!is.null(casepms$mutaters)) {
        sobj = mutate_meta(srtobj, casepms$mutaters)
    } else {
        sobj = srtobj
    }

    if (!is.null(casepms$filter)) {
        sobj = sobj |> filter(eval(parse(text=casepms$filter)))
    }

    odir = file.path(outdir, slugify(case))
    dir.create(odir, showWarnings = FALSE)
    cat(case, file = file.path(odir, "case.txt"))

    # rows: genes
    # cols: cells
    exprs = as.data.frame(
        GetAssayData(sobj, slot = "data", assay = "RNA")
    )
    genes = rownames(exprs)

    get_gene_exprs = function(gene) {
        # Cell  Expr  Group
        gdata = exprs[gene, , drop=FALSE] |>
            rownames_to_column("Gene") |>
            pivot_longer(!Gene, names_to = "Cell", values_to = "Expr") |>
            select(-"Gene")
        gdata = cbind(gdata, sobj@meta.data[gdata$Cell,,drop=FALSE])
        colnames(gdata)[ncol(gdata)] = "Group"
        gdata
    }

    do_gene = function(gene) {
        gdata = get_gene_exprs(gene)
        tryCatch({
            aov_ret = aov(Expr ~ Group, data=gdata)
            pval = summary(aov_ret)[[1]][["Pr(>F)"]][1]
            # if (!is.na(pval) && pval < 0.05) {
            #     write.table(
            #         gexprs,
            #         file.path(odir, paste(gene, "plotdata", "txt", sep=".")),
            #         row.names=FALSE,
            #         col.names=TRUE,
            #         sep="\t",
            #         quote=FALSE
            #     )
            # }
            pval
        }, error = function(e) {
            NA
        })
    }

    # Do ANOVA for each gene
    pvals = mclapply(genes, do_gene, mc.cores = ncores) %>% unlist()
    adjpvals = p.adjust(pvals, method = "BH")

    metamarkers = data.frame(Gene=genes, p_val=pvals, p_val_adj=adjpvals)
    metamarkers = metamarkers[complete.cases(metamarkers), , drop=FALSE] |>
        # filter(p_val_adj < 0.05) |>
        arrange(p_val_adj)

    # Save the results
    write.table(
        metamarkers,
        file.path(odir, "meta_makers.txt"),
        col.names=T,
        row.names=F,
        sep="\t",
        quote=FALSE
    )

    # Plot the top 10 genes in each group with violin plots
    for (gene in metamarkers |> slice_head(n=10) |> pull(Gene)) {
        plotdata = get_gene_exprs(gene)
        plotViolin(
            plotdata,
            args = list(aes(x=Group, y=Expr)),
            ggs = c("theme_prism()"),
            outfile = file.path(odir, paste(gene, "png", sep="."))
        )
    }
}


sapply(names(cases), do_case)
