source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")

library(dplyr)

mfdir = {{in.mfdir | r}}
outdir = {{out.outdir | r}}
cases = {{envs.cases | r}}

do_case = function(case) {
    overlaps = cases[[case]]$overlaps
    devpars = list_update(
        list(res = 100, width = 1200, height = 1000),
        cases[[case]]$devpars
    )
    odir = file.path(outdir, case)
    dir.create(odir, showWarnings = FALSE)

    ovmarkers = list()
    ovgenes = c()
    for (ovgroup in overlaps) {
        markerfile = file.path(mfdir, ovgroup, "markers.txt")
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
        ) |> filter(p_val_adj < 0.05)
        ovmarkers[[ovgroup]] = markers$gene
        ovgenes = unique(c(ovgenes, ovmarkers[[ovgroup]]))
    }
    if (length(ovmarkers) <= 4) {
        plotVenn(
            ovmarkers,
            devpars = devpars,
            outfile = file.path(odir, "overlapping_markers.png")
        )
    } else {
        plotUpset(
            ovmarkers,
            devpars = devpars,
            ggs = c("geom_bar(aes(x=V1))", "xlab(NULL)"),
            outfile = file.path(odir, "overlapping_markers.png")
        )
    }

    ovmarkers_table = data.frame(Gene=ovgenes)
    for (ovname in names(ovmarkers)) {
        ovmarkers_table[[ovname]] = FALSE
        ovmarkers_table[ovmarkers_table$Gene %in% ovmarkers[[ovname]], ovname] = TRUE
    }

    write.table(
        ovmarkers_table,
        file.path(odir, "overlapping_markers.txt"),
        row.names=F,
        col.names=T,
        sep="\t",
        quote=FALSE)
}

sapply(names(cases), do_case)
