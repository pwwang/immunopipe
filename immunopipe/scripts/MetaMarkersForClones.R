library(dplyr)
library(tibble)
library(parallel)
library(Seurat)
library(ggprism)
source("{{biopipen_dir}}/utils/plot.R")

srtobj = {{in.srtobj | r}}
groupfile = {{in.groupfile | r}}
markersdir = {{in.markersdir | r}}
config = {{in.configfile | read | toml_loads | r}}
outdir = {{out.outdir | r}}
ncores = {{envs.ncores | r}}
min_cells = {{envs.min_cells | r}}
venn_devpars = {{envs.venn_devpars | r}}
upset_devpars = {{envs.upset_devpars | r}}

# markers overlapping
markers_ov_dir = file.path(outdir, "markers_overlapping")
dir.create(markers_ov_dir, showWarnings = FALSE)
if (length(config$overlap) > 0) {
    print("- Overlapping markers from different comparisons")
    for (ovname in names(config$overlap)) {
        print(paste0("  ", ovname))
        ovmarkers = list()
        ovgenes = c()

        modir = file.path(markers_ov_dir, ovname)
        dir.create(modir, showWarnings = FALSE)

        for (ovgroup in config$overlap[[ovname]]) {
            markerfile = file.path(markersdir, ovgroup, "markers.txt")
            if (!file.exists(markerfile)) {
                ovmarkers[[ovgroup]] = c()
                next
            }
            con = file(markerfile)
            first_line = readLines(con, n=1)
            close(con)
            if (grepl("Error", first_line)) {
                ovmarkers[[ovgroup]] = c()
                next
            }
            markers = read.table(
                markerfile,
                header = TRUE,
                row.names = NULL,
                sep = "\t",
                check.names = FALSE
            ) %>% filter(p_val_adj < 0.05)
            ovmarkers[[ovgroup]] = markers$gene
            ovgenes = unique(c(ovgenes, ovmarkers[[ovgroup]]))
        }
        if (length(ovmarkers) <= 4) {
            plotVenn(
                ovmarkers,
                devpars = venn_devpars,
                outfile = file.path(modir, "overlapping_markers.png")
            )
        } else {
            plotUpset(
                ovmarkers,
                devpars = upset_devpars,
                ggs = c("geom_bar(aes(x=V1))", "xlab(NULL)"),
                outfile = file.path(modir, "overlapping_markers.png")
            )
        }

        ovmarkers_table = data.frame(Gene=ovgenes)
        for (ovname in names(ovmarkers)) {
            ovmarkers_table[[ovname]] = FALSE
            ovmarkers_table[ovmarkers_table$Gene %in% ovmarkers[[ovname]], ovname] = TRUE
        }

        write.table(
            ovmarkers_table,
            file.path(modir, "overlapping_markers.txt"),
            row.names=F,
            col.names=T,
            sep="\t",
            quote=FALSE)
    }

}

# meta markers
meta_markers_dir = file.path(outdir, "meta_markers")
dir.create(meta_markers_dir, showWarnings = FALSE)
if (length(config$meta) > 0) {
    print("- Meta markers for different groups")
    print("  Reading seurat object ...")
    seurat_obj = readRDS(srtobj)
    # Cell_Barcode  Group
    groups = read.table(groupfile, header=T, row.names=1, sep="\t", check.names=F) |>
        select(last_col())

    DefaultAssay(seurat_obj) <- "RNA"
    seurat_obj = NormalizeData(seurat_obj)
    allcells = intersect(rownames(groups), colnames(seurat_obj))
    # rows: genes
    # cols: cells
    exprs = as.data.frame(
        GetAssayData(seurat_obj, slot = "data", assay = "RNA")
    )[, allcells, drop=FALSE]
    genes = rownames(exprs)

    do_case = function(mmname) {
        print(paste("  Handling case", mmname))
        mmdir = file.path(meta_markers_dir, mmname)
        dir.create(mmdir, showWarnings = FALSE)

        # Pull the cells of this group
        subsets = config$meta[[mmname]]
        mmgroups = groups[groups[, 1] %in% subsets, , drop=FALSE]
        cells = intersect(rownames(mmgroups), allcells)
        exprmat = exprs[, cells, drop=FALSE]

        get_gene_exprs = function(gene) {
            gdata = exprmat[gene, , drop=FALSE] |>
                rownames_to_column("Gene") |>
                pivot_longer(!Gene, names_to = "Cell", values_to = "Expr") |>
                select(-"Gene")
            gdata = cbind(gdata, mmgroups[gdata$Cell,,drop=FALSE])
            colnames(gdata)[ncol(gdata)] = "Group"
            gdata
        }

        do_gene = function(gene) {
            gdata = get_gene_exprs(gene)
            tryCatch({
                aov_ret = aov(Expr ~ Group, data=gdata)
                pval = summary(aov_ret)[[1]][["Pr(>F)"]][1]
                if (!is.na(pval) && pval < 0.05) {
                    write.table(
                        gexprs,
                        file.path(mmdir, paste(gene, "plotdata", "txt", sep=".")),
                        row.names=FALSE,
                        col.names=TRUE,
                        sep="\t",
                        quote=FALSE
                    )
                }
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
            file.path(mmdir, "meta_makers.txt"),
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
                outfile = file.path(mmdir, paste(gene, "png", sep="."))
            )
        }
    }

    for (mmname in names(config$meta)) {
        do_case(mmname)
    }

}
