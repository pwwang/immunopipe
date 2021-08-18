library(colorspace)
library(RcppTOML)

tcrdir = "{{ in.tcrdir }}"
tcdir = "{{ in.tcdir }}"
itdir = "{{ in.itdir }}"
outdir = "{{ out.outdir }}"
tclusters = '{{ args.tclusters }}'

dir.create(outdir, showWarnings = FALSE)
tclusters = parseTOML(tclusters, fromFile=FALSE)$names

idents.tcell.ordered = names(tclusters)
idents.tcell.labels = unname(unlist(tclusters))

C <- matrix(0,8,8)
for (i in 1:8) {
    for (j in 1:8) {
        lum <- 85 - 85*((i-1)/7)^.5*((j-1)/7)^.5
        a <- 100*(i-1)/7
        b <- -100*(j-1)/7
        C[i,j] <- hex(LAB(lum,a,b),fixup=TRUE)
    }
}

arg.max <- function (x) {
    l <- which(x == max(x))
    if (length(l) == 1) {
        l
    } else {
        NA
    }
}

is.normal = function(src, controls = c("wbc", "blood", "normal", "control")) {
    tolower(src) %in% controls
}

compute.clonotype.ident <- function (
    clonotype.tcounts.presep,
    clonotype.ncounts.presep,
    outfile
) {

    p2.tcell <- read.csv(file.path(tcdir, "Seurat-1", "global.tcell.umap"), row.names=1, sep="\t")
    idents <- as.character(p2.tcell$ident)
    idents <- factor(idents, levels=idents.tcell.ordered, labels=idents.tcell.labels)
    clonotypes <- as.character(p2.tcell$clonotype)

    nidents <- length(levels(idents))

    clonotypes.uniq <- unique(clonotypes)
    clonotypes.uniq <- clonotypes.uniq[!is.na(clonotypes.uniq)]  # 55260 or 56974

    clonotype.tcounts.presep = clonotype.tcounts.presep[
        intersect(clonotypes.uniq, names(clonotype.tcounts.presep))
    ]
    clonotype.ncounts.presep = clonotype.ncounts.presep[
        intersect(clonotypes.uniq, names(clonotype.ncounts.presep))
    ]

    t_idx = sort(match(clonotypes.uniq, names(clonotype.tcounts.presep)))
    n_idx = sort(match(clonotypes.uniq, names(clonotype.ncounts.presep)))

    clonotypes.tumor.uniq <- names(clonotype.tcounts.presep)[t_idx]
    nclonotypes.tumor <- length(clonotypes.tumor.uniq)
    clonotypes.normal.uniq <- names(clonotype.ncounts.presep)[n_idx]
    nclonotypes.normal <- length(clonotypes.normal.uniq)

    clonotypes.tumor.sample <- unlist(lapply(
        strsplit(clonotypes.tumor.uniq,"[.]"),
        function (x) {paste(x[1],x[2],sep=".")}
    ))
    clonotypes.normal.sample <- unlist(lapply(
        strsplit(clonotypes.normal.uniq,"[.]"),
        function (x) {paste(x[1],x[2],sep=".")}
    ))

    nclonotypes = length(clonotypes.uniq)
    clonotype.idents <- matrix(NA, nclonotypes, nidents)
    rownames(clonotype.idents) <- clonotypes.uniq
    colnames(clonotype.idents) <- levels(idents)
    tab <- table(clonotypes, idents)  # 55260 x 16 or 56974 x 17
    clonotype.idents[rownames(tab),] <- tab

    clonotype.tumor.idents <- matrix(NA, nclonotypes.tumor, nidents)
    clonotype.normal.idents <- matrix(NA, nclonotypes.normal, nidents)
    rownames(clonotype.tumor.idents) <- clonotypes.tumor.uniq
    colnames(clonotype.tumor.idents) <- levels(idents)
    rownames(clonotype.normal.idents) <- clonotypes.normal.uniq
    colnames(clonotype.normal.idents) <- levels(idents)

    normal_ind = is.normal(p2.tcell$source)
    tab.tumor = table(clonotypes[!normal_ind], idents[!normal_ind])
    clonotype.tumor.idents[rownames(tab.tumor),] <- tab.tumor

    tab.normal = table(clonotypes[normal_ind], idents[normal_ind])
    clonotype.normal.idents[rownames(tab.normal),] <- tab.normal

    clonotype.tumor.ident <- factor(
        levels(idents)[apply(clonotype.tumor.idents, 1, arg.max)],
        levels=idents.tcell.labels
    )
    clonotype.normal.ident <- factor(
        levels(idents)[apply(clonotype.normal.idents, 1, arg.max)],
        levels=idents.tcell.labels
    )

    names(clonotype.tumor.ident) <- clonotypes.tumor.uniq
    names(clonotype.normal.ident) <- clonotypes.normal.uniq

    clonotype.tumor.idents = clonotype.tumor.idents[complete.cases(clonotype.tumor.idents), ]
    clonotype.normal.idents = clonotype.normal.idents[complete.cases(clonotype.normal.idents), ]
    # # Each 55260 x 16 or 56974 x 17
    # clonotype.tumor.idents <- matrix(0, nclonotypes, nidents)
    # clonotype.normal.idents <- matrix(0, nclonotypes, nidents)

    # # Not every sample had blood, so set NA accordingly
    # # clonotype.blood.idents <- matrix(NA, nclonotypes, nidents)
    # # clonotype.blood.idents[clonotypes.sample %in% tnb.samples,] <- 0   # For samples with blood

    # rownames(clonotype.tumor.idents) <- clonotypes.uniq
    # rownames(clonotype.normal.idents) <- clonotypes.uniq
    # # rownames(clonotype.blood.idents) <- clonotypes.uniq
    # colnames(clonotype.tumor.idents) <- levels(idents)
    # colnames(clonotype.normal.idents) <- levels(idents)
    # # colnames(clonotype.blood.idents) <- levels(idents)

    # p2.source <- as.character(p2.tcell$source)

    # tab <- table(clonotypes[!is.normal(p2.source)], idents[!is.normal(p2.source)])
    # clonotype.tumor.idents[rownames(tab),] <- tab

    # tab <- table(clonotypes[is.normal(p2.source)], idents[is.normal(p2.source)])
    # clonotype.normal.idents[rownames(tab),] <- tab

    # tab <- table(clonotypes[p2.source == "Blood"], idents[p2.source == "Blood"])
    # clonotype.blood.idents[rownames(tab),] <- tab

    save(
        clonotypes.tumor.uniq,
        clonotypes.normal.uniq,
        clonotypes.tumor.sample,
        clonotypes.normal.sample,
        clonotype.idents,
        clonotype.tumor.ident,
        clonotype.normal.ident,
        clonotype.tumor.idents,
        clonotype.normal.idents, #clonotype.blood.idents,
        file=outfile
    )

    clonotypes.uniq.tcounts <- rowSums(clonotype.tumor.idents)
    clonotypes.uniq.ncounts <- rowSums(clonotype.normal.idents)
    # clonotypes.uniq.bcounts <- rowSums(clonotype.blood.idents)

    list(clonotypes.uniq=clonotypes.uniq,
         tcounts=clonotypes.uniq.tcounts,
         ncounts=clonotypes.uniq.ncounts)
}


compute.residency <- function (
    clonotypes.uniq.tcounts,
    clonotypes.uniq.ncounts,
    # clonotypes.uniq.bcounts,
    # clonotypes.uniq,
    countfile,
    resdfile
) {
    # nclonotypes <- length(clonotypes.uniq)
    common_clonotypes = intersect(
        names(clonotypes.uniq.tcounts),
        names(clonotypes.uniq.ncounts)
    )
    nclonotypes <- length(common_clonotypes)

    cl.colors <- rep(NA,nclonotypes)
    clonotype.residency <- rep(NA,nclonotypes)
    names(cl.colors) <- as.character(common_clonotypes)
    names(clonotype.residency) <- as.character(common_clonotypes)

    for (cl in common_clonotypes) {
        tumor.count <- clonotypes.uniq.tcounts[cl]
        normal.count <- clonotypes.uniq.ncounts[cl]
        # blood.count <- clonotypes.uniq.bcounts[cl]
        # blood.count <- 1

        if (tumor.count == 0 && normal.count == 0) {
            stop("Both tumor and normal counts are zero")
            # if (is.na(blood.count) || blood.count == 0) {
            #     stop("tumor, normal, and blood counts are all zero")
            # } else if (blood.count == 1) {
            #     cl.colors[cl] <- residency.pal[6]
            #     clonotype.residency[cl] <- "b"
            # } else {
            #     cl.colors[cl] <- residency.pal[7]
            #     clonotype.residency[cl] <- "B"
            # }
        } else if (tumor.count == 1 && normal.count == 0) {
            cl.colors[cl] <- orange
            clonotype.residency[cl] <- "t"
        } else if (normal.count == 1 && tumor.count == 0) {
            cl.colors[cl] <- yellow
            clonotype.residency[cl] <- "n"
        } else {
            if (tumor.count == 0) {
                clonotype.residency[cl] <- "N"
            } else if (normal.count == 0) {
                clonotype.residency[cl] <- "T"
            } else {
                clonotype.residency[cl] <- "D"
            }

            if (tumor.count <= 1) {
                tumor.level <- 1
            } else {
                tumor.level <- ceiling(9/5 * log(tumor.count))
                if (tumor.level > 8) { tumor.level = 8 }
            }
            if (normal.count <= 1) {
                normal.level <- 1
            } else {
                normal.level <- ceiling(9/5 * log(normal.count))
                if (normal.level > 8) { normal.level = 8 }
            }
            cl.colors[cl] <- C[tumor.level,normal.level]
        }
    }

    # Compute and all and restrict to T cells
    p2.all <- read.csv(
        file.path(itdir, "global.0.8.umap"),
        row.names=1,
        sep="\t"
    )
    p2.tcell <- read.csv(
        file.path(tcdir, "Seurat-1", "global.tcell.umap"),
        row.names=1,
        sep="\t"
    )

    tcell.colors <- cl.colors[
        intersect(
            names(cl.colors),
            as.character(p2.tcell$clonotype)
        )
    ]  # 141623
    # tcell.colors <- ifelse(is.na(tcell.colors), "white", tcell.colors)   # leave as NA
    #cell.colors.tn <- cell.colors
    #cell.colors.tn[which(p2$source == "Blood")] <- NA

    clonotype.residency.f <- factor(clonotype.residency, levels=c("n","N","D","T","t"))  # 56974
    names(clonotype.residency.f) <- common_clonotypes

    all.residency <- clonotype.residency.f[
        intersect(
            names(clonotype.residency.f),
            as.character(p2.all$clonotype)
        )
    ]
    all.residency.f <- factor(all.residency, levels=c("n","N","D","T","t"))  # 200626
    tcell.residency <- clonotype.residency.f[
        intersect(
            names(clonotype.residency.f),
            as.character(p2.tcell$clonotype)
        )
    ]
    tcell.residency.f <- factor(tcell.residency, levels=c("n","N","D","T","t"))  # 141623

    # clonotype.bloodexp <- rep(NA, length(clonotypes.uniq))
    # clonotype.bloodexp[which(clonotypes.uniq.bcounts == 0)] <- "Ind"
    # clonotype.bloodexp[which(clonotypes.uniq.bcounts == 1)] <- "Non"
    # clonotype.bloodexp[which(clonotypes.uniq.bcounts > 1)] <- "Exp"
    # clonotype.bloodexp.f <- factor(clonotype.bloodexp, levels=c("Ind","Non","Exp"))
    # names(clonotype.bloodexp.f) <- clonotypes.uniq

    # all.bloodexp <- clonotype.bloodexp.f[as.character(p2.all$clonotype)]
    # all.bloodexp.f <- factor(all.bloodexp, levels=c("Ind","Non","Exp"))
    # tcell.bloodexp <- clonotype.bloodexp.f[as.character(p2.tcell$clonotype)]
    # tcell.bloodexp.f <- factor(tcell.bloodexp, levels=c("Ind","Non","Exp"))

    save(tcell.colors, cl.colors,
         clonotype.residency.f, all.residency.f, tcell.residency.f,
        #  clonotype.bloodexp.f, all.bloodexp.f, tcell.bloodexp.f,
         file=resdfile)

    all.tcounts <- clonotypes.uniq.tcounts[
        intersect(
            names(clonotypes.uniq.tcounts),
            as.character(p2.all$clonotype)
        )
    ]
    all.ncounts <- clonotypes.uniq.ncounts[
        intersect(
            names(clonotypes.uniq.ncounts),
            as.character(p2.all$clonotype)
        )
    ]
    # all.bcounts <- clonotypes.uniq.bcounts[as.character(p2.all$clonotype)]

    tcell.tcounts <- clonotypes.uniq.tcounts[
        intersect(
            names(clonotypes.uniq.tcounts),
            as.character(p2.tcell$clonotype)
        )
    ]
    tcell.ncounts <- clonotypes.uniq.ncounts[
        intersect(
            names(clonotypes.uniq.ncounts),
            as.character(p2.tcell$clonotype)
        )
    ]
    # tcell.bcounts <- clonotypes.uniq.bcounts[as.character(p2.tcell$clonotype)]

    save(clonotypes.uniq.tcounts, clonotypes.uniq.ncounts, #clonotypes.uniq.bcounts,
         all.tcounts, all.ncounts, #all.bcounts,
         tcell.tcounts, tcell.ncounts, #tcell.bcounts,
         file=countfile)
}


clonotype.tcounts.presep <- c()
clonotype.ncounts.presep <- c()
for (cdata in Sys.glob(file.path(tcrdir, "*.clonotype.counts.RData"))) {
    load(cdata)
    tcounts = if (is.normal(sources[1])) clonotype.counts2 else clonotype.counts1
    ncounts = if (is.normal(sources[1])) clonotype.counts1 else clonotype.counts2
    clonotype.tcounts.presep <- c(clonotype.tcounts.presep, tcounts)
    clonotype.ncounts.presep <- c(clonotype.ncounts.presep, ncounts)
}

counts <- compute.clonotype.ident(
    clonotype.tcounts.presep,
    clonotype.ncounts.presep,
    outfile = file.path(outdir, "global.clonotype.ident.RData")
)

compute.residency(
    counts$tcounts,
    counts$ncounts,
    # counts$clonotypes.uniq,
    countfile = file.path(outdir, "global.counts.RData"),
    resdfile = file.path(outdir, "global.residency.RData")
)
