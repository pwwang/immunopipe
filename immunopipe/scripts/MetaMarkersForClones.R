library(dplyr)
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
            ovmarkers[[ovgroup]] = markers$Gene
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
    groups = read.table(groupfile, header=T, row.names=1, sep="\t", check.names=F)

    do_one_metamarkers_case = function(mmname) {
        print(paste("  Handling case", mmname))
        mmdir = file.path(meta_markers_dir, mmname)
        dir.create(mmdir, showWarnings = FALSE)

        metagroups = config$meta[[mmname]]

        cells = lapply(metagroups, function(mgroup) {
            groups[mgroup, 1] %>% strsplit(";", fixed=T) %>% unlist()
        })
        names(cells) = metagroups

        exprs = lapply(metagroups, function(mgroup) {
            sobj = subset(seurat_obj, cells = cells[[mgroup]])
            DefaultAssay(sobj) <- "RNA"
            sobj = NormalizeData(sobj)
            as.data.frame(
                GetAssayData(sobj, slot = "data", assay = "RNA")
            )
        })
        names(exprs) = metagroups

        do_gene = function(gene) {
            gexprs = lapply(metagroups, function(mgroup) {
                expr = unlist(exprs[[mgroup]][gene,,drop=T])
                if (length(expr) < min_cells) {
                    return(NULL)
                } else {
                    return(data.frame(Expr=expr, Group=mgroup))
                }
            })
            gexprs = do.call(rbind, gexprs)
            tryCatch({
                aov_ret = aov(Expr ~ Group, data=gexprs)
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
        genes = rownames(exprs[[metagroups[1]]])
        pvals = mclapply(genes, do_gene, mc.cores = ncores) %>% unlist()
        adjpvals = p.adjust(pvals, method = "BH")

        metamarkers = data.frame(Gene=genes, p_val=pvals, p_val_adj=adjpvals)
        metamarkers = metamarkers[complete.cases(metamarkers), , drop=FALSE] %>%
            filter(p_val_adj < 0.05) %>%
            arrange(p_val_adj)

        write.table(
            metamarkers,
            file.path(mmdir, "meta_makers.txt"),
            col.names=T,
            row.names=F,
            sep="\t",
            quote=FALSE
        )

        for (gene in metamarkers %>% slice_head(n=10) %>% pull(Gene)) {
            plotdata = read.table(
                file.path(mmdir, paste(gene, "plotdata", "txt", sep=".")),
                row.names=NULL,
                header=TRUE,
                sep="\t",
                check.names = FALSE
            )
            plotViolin(
                plotdata,
                args = list(aes(x=Group, y=Expr)),
                ggs = c("theme_prism()"),
                outfile = file.path(mmdir, paste(gene, "png", sep="."))
            )
        }
    }

    for (mmname in names(config$meta)) {
        do_one_metamarkers_case(mmname)
    }
}
