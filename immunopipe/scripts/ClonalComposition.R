
library(RColorBrewer)
library(RcppTOML)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

pccdir = "{{ in.pccdir }}"
tcrdir = "{{ in.tcrdir }}"
outdir = "{{ out.outdir }}"
tclusters = '{{ args.tclusters }}'
multipt_samples = '{{ args.multipt_samples }}'

tclusters = parseTOML(tclusters, fromFile=FALSE)
multipt_samples = parseTOML(multipt_samples, fromFile=FALSE)

idents.tcell.pal = unlist(tclusters$colors)
names(idents.tcell.pal) = unlist(tclusters$names)

dir.create(outdir, showWarnings = FALSE)

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(
        toupper(substring(s, 1, 1)),
        {s <- substring(s, 2); if(strict) tolower(s) else s},
        sep = "",
        collapse = " "
    )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

load(file.path(pccdir, "global.clonotype.ident.RData"))
load(file.path(pccdir, "global.counts.RData"))

biggest.tumor.clones = names(head(sort(clonotypes.uniq.tcounts[rownames(clonotype.tumor.idents)], decr=T), 20))
biggest.normal.clones = names(head(sort(clonotypes.uniq.ncounts[rownames(clonotype.normal.idents)], decr=T), 20))

# layout.mat <- matrix(
#     c(1,2,3,4,
#       5,6,7,8,
#       9,10,11,12,
#       13,14,15,16,
#       17,18,0,0,
#       19,20,0,0),
#     6, 4,
#     byrow=T
# )
# layout.widths <- c(1, 1, 1, 1)
# layout.heights <- c(1, 1, 1, 1, 1, 1)


# layout(layout.mat, layout.widths, layout.heights)
# par(mar=c(0.3, 0, 1.0, 0))
#par(lheight=0.8)

plot.clonal.composition = function(biggest.clones, counts, idents, plotfile) {
    plots = list()
    i = 1
    for (cl in biggest.clones) {
        size <- counts[cl]

        clone <- paste(capwords(strsplit(cl,"[.]")[[1]][1]), strsplit(cl,"[.]")[[1]][2], sep=".")
        df = idents[cl,,drop=F] %>%
            as.data.frame() %>%
            rownames_to_column('Clone') %>%
            pivot_longer(2:ncol(.), names_to="Cluster", values_to="Count")

        p = ggplot(df, aes(x="", y=Count, fill=Cluster)) +
            geom_bar(stat="identity", width=0.038*sqrt(size)) +
            coord_polar("y", start=0) +
            scale_fill_manual(values=idents.tcell.pal) +
            ggtitle(bquote(.(clone) ~~ italic(n) == .(size))) +
            theme_void() +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

        plots[[i]] = p
        i = i + 1
    }
    png(plotfile, res=100, width=700, height=2000)
    plots$ncol = 2
    plots$scale = sqrt(counts[biggest.clones])
    # scale to .3 ~ 1
    plots$scale = .3 + .7 * (plots$scale - min(plots$scale)) / (max(plots$scale) - min(plots$scale))
    print(do.call(plot_grid, plots))
    dev.off()
}

plot.clonal.composition(
    biggest.tumor.clones,
    clonotypes.uniq.tcounts,
    clonotype.tumor.idents,
    file.path(outdir, "clonal-pies.tumor.png")
)

plot.clonal.composition(
    biggest.normal.clones,
    clonotypes.uniq.ncounts,
    clonotype.normal.idents,
    file.path(outdir, "clonal-pies.normal.png")
)


plot.ms.clonal.composition = function(sc_df, all_idents, by, plotfile, ncol=2, top=5) {
    plots = list()
    filter_idx = T
    for (i in seq_len(ncol)) {
        filter_idx = filter_idx & (sc_df[, i] %in% rownames(all_idents))
    }
    sc_df = sc_df %>%
        filter(filter_idx) %>%
        arrange(desc(.[[by]])) %>%
        slice_head(n=top)
    i = 1
    sizes = c()
    for (r in seq_len(nrow(sc_df))) {
        cls = as.character(unlist(sc_df[r, 1:ncol, drop=TRUE]))
        counts = unname(unlist(sc_df[r, (ncol+1):(2*ncol), drop=TRUE]))
        names(counts) = cls
        sizes = c(sizes, counts[cls])
        for (cl in cls) {
            size <- counts[cl]
            clone <- paste(capwords(strsplit(cl,"[.]")[[1]][1]), strsplit(cl,"[.]")[[1]][2], sep=".")

            df = all_idents[cl,,drop=F] %>%
                as.data.frame() %>%
                rownames_to_column('Clone') %>%
                pivot_longer(2:ncol(.), names_to="Cluster", values_to="Count")

            p = ggplot(df, aes(x="", y=Count, fill=Cluster)) +
                geom_bar(stat="identity", width=0.038*sqrt(size)) +
                coord_polar("y", start=0) +
                scale_fill_manual(values=idents.tcell.pal) +
                ggtitle(bquote(.(clone) ~~ italic(n) == .(size))) +
                theme_void() +
                theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

            plots[[i]] = p
            i = i + 1
        }
    }
    png(plotfile, res=100, width=ncol*350, height=1200)
    plots$ncol = ncol
    plots$scale = sqrt(sizes)
    plots$scale = .3 + .7 * (plots$scale - min(plots$scale)) / (max(plots$scale) - min(plots$scale))
    print(do.call(plot_grid, plots))
    dev.off()
}

multiptsamples_dir = file.path(outdir, "MultiPatient-Samples")
dir.create(multiptsamples_dir, showWarnings=FALSE)

all_shared_clones = list()
x = 1
ms_patients = c()
for (mspatient in names(multipt_samples)) {
    if (mspatient == "COMPARING") next
    ms_patients = c(ms_patients, mspatient)

    msp_dir = file.path(multiptsamples_dir, mspatient)
    dir.create(msp_dir, showWarnings=FALSE)

    # tcounts_idx = TRUE
    # ncounts_idx = TRUE
    masters = list()
    i = 1
    for (ms in multipt_samples[[mspatient]]) {
        # tcounts_idx = tcounts_idx | startsWith(names(clonotypes.uniq.tcounts), ms)
        # ncounts_idx = ncounts_idx | startsWith(names(clonotypes.uniq.tcounts), ms)

        masters[[i]] = read.table(
            file.path(tcrdir, paste0(ms, ".master")),
            header=FALSE,
            row.names=NULL,
            sep="\t"
        )
        i = i + 1
    }
    masters$by = c("V2", "V4")
    shared_clones = do.call(inner_join, masters) %>% select(
            Clone1=V1.x,
            Clone2=V1.y,
            alpha_consensus=V2,
            beta_consensus=V4
        ) %>% mutate(
            Count1_Tumor=clonotypes.uniq.tcounts[Clone1],
            Count1_Normal=clonotypes.uniq.ncounts[Clone1],
            Count2_Tumor=clonotypes.uniq.tcounts[Clone2],
            Count2_Normal=clonotypes.uniq.tcounts[Clone2],
            Count1=Count1_Tumor + Count1_Normal,
            Count2=Count2_Tumor + Count2_Normal,
            Total=Count1 + Count2,
            .before=alpha_consensus
        ) %>% filter(
            Clone1 %in% rownames(clonotype.tumor.idents) |
            Clone1 %in% rownames(clonotype.tumor.idents) |
            Clone2 %in% rownames(clonotype.normal.idents) |
            Clone2 %in% rownames(clonotype.normal.idents)
        ) %>% arrange(desc(Total))
    all_shared_clones[[x]] = shared_clones
    x = x + 1

    write.table(
        shared_clones,
        file.path(msp_dir, "Shared_Clones.txt"),
        row.names=FALSE,
        col.names=TRUE,
        sep="\t",
        quote=FALSE
    )

    # tumor and normal merged clones
    plot.ms.clonal.composition(
        shared_clones %>% select(Clone1, Clone2, Count1, Count2),
        all_idents=clonotype.idents,
        by="Count1",
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Merged_OrderBy.", multipt_samples[[mspatient]][1] ,".png")
        )
    )

    plot.ms.clonal.composition(
        shared_clones %>% select(Clone1, Clone2, Count1, Count2),
        all_idents=clonotype.idents,
        by="Count2",
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Merged_OrderBy.", multipt_samples[[mspatient]][2] ,".png")
        )
    )

    # shared between tumor/normal/pre/post samples
    plot.ms.clonal.composition(
        shared_clones %>% mutate(
            Clone1_Tumor=paste0(Clone1, ".Tumor"),
            Clone1_Normal=paste0(Clone1, ".Normal"),
            Clone2_Tumor=paste0(Clone2, ".Tumor"),
            Clone2_Normal=paste0(Clone2, ".Normal"),
            .before=1
        ) %>% select(-Clone1, -Clone2),
        all_idents=bind_rows(
            clonotype.tumor.idents %>%
                as.data.frame() %>%
                rownames_to_column('Clone') %>%
                mutate(Clone=paste0(Clone, ".Tumor")) %>%
                column_to_rownames('Clone'),
            clonotype.normal.idents %>%
                as.data.frame() %>%
                rownames_to_column('Clone') %>%
                mutate(Clone=paste0(Clone, ".Normal")) %>%
                column_to_rownames('Clone')
        ),
        by="Total",
        ncol=4,
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Across_OrderBy.Total.png")
        )
    )

    # tumor only
    plot.ms.clonal.composition(
        shared_clones %>% select(Clone1, Clone2, Count1_Tumor, Count2_Tumor),
        all_idents=clonotype.tumor.idents,
        by="Count1_Tumor",
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Tumor_OrderBy.", multipt_samples[[mspatient]][1] ,".png")
        )
    )

    plot.ms.clonal.composition(
        shared_clones %>% select(Clone1, Clone2, Count1_Tumor, Count2_Tumor),
        all_idents=clonotype.tumor.idents,
        ncol=2,
        by="Count2_Tumor",
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Tumor_OrderBy.", multipt_samples[[mspatient]][2] ,".png")
        )
    )

    # normal only
    plot.ms.clonal.composition(
        shared_clones %>% select(Clone1, Clone2, Count1_Normal, Count2_Normal),
        all_idents=clonotype.normal.idents,
        by="Count1_Normal",
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Normal_OrderBy.", multipt_samples[[mspatient]][1] ,".png")
        )
    )

    plot.ms.clonal.composition(
        shared_clones %>% select(Clone1, Clone2, Count1_Normal, Count2_Normal),
        all_idents=clonotype.normal.idents,
        by="Count2_Normal",
        plotfile=file.path(
            msp_dir,
            paste0("Clonal_Composition_Normal_OrderBy.", multipt_samples[[mspatient]][2] ,".png")
        )
    )

}

# # All multi-sample patients
# msp_dir = file.path(multiptsamples_dir, "All_MultiSample_Patients")
# dir.create(msp_dir, showWarnings=FALSE)

# all_shared_clones$by = c('alpha_consensus', 'beta_consensus')
# all_shared_clones$suffix = ms_patients

# all_sc_df = do.call(inner_join, all_shared_clones) %>%
#     select(
#         startsWith("Clone1.") |
#         startsWith("Clone2.") |
#         startsWith("Count1_") |
#         startsWith("Count2_")
#     ) %>% rowwise() %>% mutate(
#         Total=sum(c_across(startsWith("Count1_") | startsWith("Count2_"))),
#         Total_Normal=sum(c_across(startsWith("Count1_Tumor") | startsWith("Count2_Tumor"))),
#         Total_Tumor=sum(c_across(startsWith("Count1_Normal") | startsWith("Count2_Normal"))),
#     )

# plot.ms.clonal.composition(
#     all_sc_df,
#     all_idents=all_idents,
#     by="Total",
#     ncol=len(all_ms_samples),
#     plotfile=file.path(
#         msp_dir,
#         paste0("Clonal_Composition_Merged_OrderBy.Total.png")
#     )
# )

# plot.ms.clonal.composition(
#     all_sc_df,
#     all_idents=clonotype.tumor.idents,
#     by="Total_Tumor",
#     ncol=len(all_ms_samples),
#     plotfile=file.path(
#         msp_dir,
#         paste0("Clonal_Composition_Tumor_OrderBy.Total.png")
#     )
# )

# plot.ms.clonal.composition(
#     all_sc_df,
#     all_idents=clonotype.normal.idents,
#     by="Total_Normal",
#     ncol=len(all_ms_samples),
#     plotfile=file.path(
#         msp_dir,
#         paste0("Clonal_Composition_Normal_OrderBy.Total.png")
#     )
# )
