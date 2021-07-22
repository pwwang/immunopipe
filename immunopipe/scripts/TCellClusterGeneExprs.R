library(Seurat)
library(dplyr)
library(RcppTOML)
library(RColorBrewer)

cldir = "{{ in.cldir }}"
outdir = "{{ out.outdir }}"
genes = '{{ args.genes }}'
tclusters = '{{ args.tclusters }}'

dir.create(outdir, showWarnings = FALSE)
genes = parseTOML(genes, fromFile=FALSE)
tclusters = parseTOML(tclusters, fromFile=FALSE)

idents.tcell.ordered = names(tclusters)
idents.tcell.labels = unname(unlist(tclusters))
# colors configurable?
idents.tcell.pal <- c(brewer.pal(12,"Paired"), brewer.pal(8, "Dark2"))[
    1:length(idents.tcell.ordered)
]

load(file.path(cldir, "Seurat-1", 'global.tcell.RData'))
# Scale the RNA assay so we can have all genes
integrated.features = rownames(global.obj@assays$integrated@scale.data)
integrated.features = unique(c(integrated.features, unlist(genes)))
scale.data = ScaleData(global.obj, assay="RNA", features=integrated.features)
# scale.data <- global.obj@assays$integrated@scale.data
scale.data <- scale.data@assays$RNA@scale.data

idents <- factor(Idents(global.obj), levels=idents.tcell.ordered, labels=idents.tcell.labels)
nidents = length(levels(idents))

large.space <- 1.4
small.space <- ((nidents-1) - 2*large.space)/(nidents-3)
seq1 <- seq(from=1, by=small.space, length.out=6)
seq2 <- seq(from=seq1[length(seq1)] + large.space, by=small.space, length.out=7)
# space <- c(seq1,seq2,seq2[length(seq2)] + large.space)

space = seq(from=1, by=small.space, length.out=nidents)

plot_one_gene = function(alias, gene) {
    split.counts <- split(scale.data[gene,], idents)
    # boxplot(
    tryCatch(
        vioplot(
            split.counts,
            names=rep("",nidents),
            col=idents.tcell.pal,
            yaxt="n",
            lwd=0.5,
            outline=F,
            ylim=c(-1,2),
            range=1,
            notch=T,
            at=space
        ),
        error=function(e)
        boxplot(
            split.counts,
            names=rep("",nidents),
            col=idents.tcell.pal,
            yaxt="n",
            lwd=0.5,
            outline=F,
            ylim=c(-1,2),
            range=1,
            notch=T,
            at=space
        )
    )
    label = ifelse (gene == alias, gene, paste0(gene, " (",alias,")"))
    mtext(label, side=4, line=0.3, cex=0.5, las=2)
    1 # don't do mtext twice
}
png(file.path(outdir, "biomarker-boxplots.png"), res=150, width=8.1, height=14.25, units='in')
par(mfcol=c(floor(length(genes)/2)+2,2))
par(mgp=c(3.4, 0.3, 0))
par(mar=c(0, 2.5, 0, 4.1))  # Needs to be large enough for (Lymphotactin)
par(tck=0)

for (alias in names(genes)) {
    print(paste("Plotting", alias, genes[[alias]]))
    plot_one_gene(alias, genes[[alias]])
}

for (i in 1:nidents) {
    mtext(idents.tcell.labels[i], side=1, line=0.4, at=space[i], las=2, cex=0.6)
}

ngenes = length(genes)
mtext("Expression level (log scaled value)",
    side=2, line=-1.8,
    at=(ngenes-(ngenes-2)/2)/ngenes,
    outer=T, cex=0.6)

dev.off()
