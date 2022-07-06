library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggprism)

theme_set(theme_prism())

sobjfile = {{ in.sobjfile | r }}
outdir = {{ out.outdir | r }}
cases = {{ envs.cases | r }}

set.seed(8525)
sobj = readRDS(sobjfile)

clusters = unique(Idents(sobj))

fisher_pvals = NULL
chisq_pvals = NULL

do_one = function(case, cluster) {
    print(paste("Doing", case, "and cluster", cluster, "..."))
    cluster_name = if (is.na(as.integer(cluster))) cluster else paste0("Cluster", cluster)
    casevals = cases[[case]]

    if (is.character(casevals$cut) && casevals$cut == "quantile") {
        probs = seq(0, 1, 0.25)
    } else if (is.character(casevals$cut) && casevals$cut == "median") {
        probs = seq(0, 1, 0.5)
    } else {
        probs = casevals$cut
    }

    breaks = quantile(
        sobj@meta.data$Proportion,
        probs = probs,
        na.rm = TRUE,
    )
    labels = paste0("q", 1:(length(breaks)-1))

    metadata = sobj@meta.data |>
        filter(seurat_clusters == cluster) |>
        mutate(
            Quantile = cut(
                Proportion,
                breaks=breaks,
                labels=labels,
                include.lowest=TRUE
            )
        ) |>
        filter(!is.na(Quantile))

    designs = casevals$design
    subsettings = casevals$subsetting
    for (dname in names(designs)) {
        print(paste("  - Design", dname, "..."))
        output_dir = file.path(outdir, case, cluster_name, dname)
        dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
        contfile = file.path(output_dir, "contingency.txt")

        groups = designs[[dname]]
        counts = lapply(groups, function(group) {
            expect = tibble(Quantile=labels, n=0)
            out = metadata |>
                filter(eval(parse(text=subsettings[[group]]))) |>
                count(Quantile)

            expect |>
                left_join(out, by="Quantile") |>
                replace_na(list(n.y = 0)) |>
                rowwise() |>
                mutate(n = sum(c_across(n.x:n.y))) |>
                select(n) |>
                ungroup()
        })

        conttable = do.call(cbind, counts)
        rownames(conttable) = labels
        colnames(conttable) = groups
        conttable = conttable |>
            rownames_to_column("Quantile") |>
            relocate(Quantile, .before=1)

        write.table(
            conttable,
            contfile,
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t"
        )

        # box plot showing proportions
        count_fig = file.path(output_dir, "counts.boxplot.png")
        plotdata = conttable |>
            pivot_longer(-"Quantile", names_to="Group", values_to="Count")
        p = ggplot(plotdata) +
            geom_col(
                aes(x=Group, y=Count, group=Quantile, fill=Quantile),
                position="dodge2"
            )
        png(count_fig, res=100, height=1000, width=1000)
        print(p)
        dev.off()

        chisq_file = file.path(output_dir, "chisq.txt")
        conttable = conttable |> select(-"Quantile")
        conttable = conttable[rowSums(conttable) > 0, , drop=FALSE]
        chisq.out = chisq.test(conttable)
        writeLines(
            capture.output(chisq.out),
            chisq_file
        )
        chisq_pvals <<- bind_rows(
            chisq_pvals,
            list(
                Cluster = cluster_name,
                Case = dname,
                p_value = chisq.out$p.value
            )
        )

        fisher_file = file.path(output_dir, "fisher.txt")
        if (nrow(conttable) < 2) {
            fisher.out = list(p.value = NA)
        } else {
            fisher.out = tryCatch({
                fisher.test(conttable)
            }, error=function(e) {
                fisher.test(conttable, simulate.p.value=TRUE, B=1e4)
            })
        }
        writeLines(
            capture.output(fisher.out),
            fisher_file
        )

        fisher_pvals <<- bind_rows(
            fisher_pvals,
            list(
                Cluster = cluster_name,
                Case = dname,
                p_value = fisher.out$p.value
            )
        )
    }
}

for (case in names(cases)) {
    for (cluster in clusters) {
        do_one(case, cluster)
    }

    write.table(
        fisher_pvals,
        file.path(outdir, case, "fisher.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    write.table(
        chisq_pvals,
        file.path(outdir, case, "chisq.txt"),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE
    )

    fisher_pvals = NULL
    chisq_pvals = NULL
}
