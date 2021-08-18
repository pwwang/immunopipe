library(dplyr)
library(tibble)
library(ggplot2)
library(ggprism)
library(ggrepel)
library(RcppTOML)
library(doParallel)
library(Seurat)
library(enrichR)

itgdir = "{{ in.itgdir }}"
outdir = "{{ out.outdir }}"
ncores = {{ args.ncores }}
commoncfg = '{{ args.commoncfg }}'

setEnrichrSite("Enrichr")
registerDoParallel(min(2, ncores))
gsea_dbs = parseTOML(commoncfg, fromFile=FALSE)$GSEA_DBs

dir.create(outdir, showWarnings = FALSE)

# plot CD3E/clonotype proportion to separate T and non-T cells
load(file.path(itgdir, "global.0.8.RData"))
cd3e.expr <- as.matrix(global.obj@assays$integrated@data)["CD3E",]
p2 <- read.csv(file.path(itgdir, "global.0.8.umap"), row.names=1, sep="\t")

oo <- sample(nrow(p2))
print(table(p2$ident)) # print to stdout

clonotype.pct <- list()
for (ident in names(table(p2$ident))) {
    x <- table(is.na(p2$clonotype[which(p2$ident == ident)]))
    clonotype.pct[[ident]] <- x[1]/(x[1]+x[2])
}
clonotype.pct <- unlist(clonotype.pct)
names(clonotype.pct) <- names(table(p2$ident))

cd3e.means <- lapply(split(cd3e.expr, p2$ident), mean)
names(cd3e.means) <- names(table(p2$ident))

# plot

tcell.ident <- names((which(clonotype.pct > .2 & cd3e.means > 0)))
nont.ident <- setdiff(names(cd3e.means), tcell.ident)

col <- rep(NA,length(cd3e.means))
names(col) <- names(table(p2$ident))
col[tcell.ident] <- "T cells"
col[nont.ident] <- "Non-T cells"

df = data.frame(
    clonotype.pct=unlist(clonotype.pct),
    cd3e.means=unlist(cd3e.means),
    cluster=names(col),
    group=col
)
write.table(
    df,
    file.path(outdir, 'clusters.txt'),
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

png(
    file.path(outdir, "subdivision-rationale.png"),
    res=100,
    height=600,
    width=1200
)

ggplot(df, aes(x=clonotype.pct, y=cd3e.means)) +
    geom_point(aes(color=group)) +
    geom_label_repel(
        aes(label=cluster),
        box.padding   = 0.35,
        point.padding = 0.5,
        segment.color = 'grey50'
    ) +
    theme_prism(base_size = 16)

dev.off()


write.table(p2[which(p2$ident %in% tcell.ident),],
            file=file.path(outdir, "tcell.main.xls"),
            quote=F, sep="\t", row.names=T, col.names=NA)

write.table(p2[which(p2$ident %in% nont.ident),],
            file=file.path(outdir, "nont.main.xls"),
            quote=F, sep="\t", row.names=T, col.names=NA)


# Markers
markers_dir = file.path(outdir, "Markers")
dir.create(markers_dir, showWarnings = FALSE)
markers_gsea_dir = file.path(outdir, "Markers-GSEA")
dir.create(markers_gsea_dir, showWarnings = FALSE)

foreach (markfile = Sys.glob(file.path(itgdir, "Markers", "*.markers.RData"))) %dopar% {
    print(paste("* For ident", basename(markfile), "..."))
    load(markfile)
    write.table(
        markers,
        file.path(
            markers_dir,
            sub(".RData", ".txt", basename(markfile), fixed=T)
        ),
        col.names = TRUE,
        row.names = TRUE,
        sep="\t",
        quote = FALSE
    )

    # GSEA
    ident = unlist(strsplit(basename(markfile), ".", fixed=TRUE))[2]
    ident_gsea_dir = file.path(markers_gsea_dir, ident)
    dir.create(ident_gsea_dir, showWarnings = FALSE)
    genes = markers %>% filter(p_val_adj < 0.05) %>% rownames()
    enriched = enrichr(genes, gsea_dbs)

    for (db in gsea_dbs) {
        outtable = file.path(ident_gsea_dir, paste0('enrichr_', db, '.txt'))
        outfig = file.path(ident_gsea_dir, paste0('enrichr_', db, '.png'))

        write.table(enriched[[db]], outtable, col.names=T, row.names=F, sep="\t", quote=F)

        png(outfig, width=1000, height=1000, res=100)
        print(plotEnrich(enriched[[db]], title=db))
        dev.off()
    }
}
