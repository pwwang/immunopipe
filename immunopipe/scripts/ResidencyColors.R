library(dplyr)
library(colorspace)

samples = "{{ in.samples }}"
tcrcount_dir = "{{ in.tcr_counts }}"
outdir = "{{ out.outdir }}"

colpat = "{{ args.colpat }}"

source(colpat)

dir.create(outdir, showWarnings = FALSE)

# Generate gradient colors
# at most 64?
# more colors and more efficient?
C <- matrix(0,8,8)
for (i in 1:8) {
    for (j in 1:8) {
        lum <- 85 - 85*((i-1)/7)^.5*((j-1)/7)^.5
        a <- 100*(i-1)/7
        b <- -100*(j-1)/7
        C[i,j] <- hex(LAB(lum,a,b),fixup=TRUE)
    }
}

residency.colors <- function (sample) {
    load(file.path(tcrcount_dir, paste0(sample, ".clonotype.counts.RData")))
    # save(clonotypes, clonotype.counts1, clonotype.counts2, sources,
    #      clonotype.residency, cell.residency, cell.source,
    #      file=file.path(outdir, paste0(sample,".clonotype.counts.RData")))

    nclonotypes <- length(clonotype.counts1)
    cl.colors <- rep(NA,nclonotypes)
    clonotype.residency <- rep(NA,nclonotypes)
    for (i in 1:nclonotypes) {
        count1 <- clonotype.counts1[i]
        count2 <- clonotype.counts2[i]
        # if (count1 == 0 && count2 == 0) { # no blood
        #     blood.count <- clonotype.bcounts[i]
        #     if (blood.count == 1) {
        #         cl.colors[i] <- residency.pal[6]
        #         clonotype.residency[i] <- "b"
        #     } else {
        #         cl.colors[i] <- residency.pal[7]
        #         clonotype.residency[i] <- "B"
        #     }
        if (count1 == 1 && count2 == 0) {
            cl.colors[i] <- orange
            clonotype.residency[i] <- tolower(sources[1])
        } else if (count2 == 1 && count1 == 0) {
            cl.colors[i] <- yellow
            clonotype.residency[i] <- tolower(sources[2])
        } else {
            if (count1 == 0) {
                clonotype.residency[i] <- toupper(sources[2])
            } else if (count2 == 0) {
                clonotype.residency[i] <- toupper(sources[1])
            } else {
                clonotype.residency[i] <- "Dual"
            }

            if (count1 <= 1) {
                level1 <- 1
            } else {
                level1 <- ceiling(9/5 * log(count1))
                if (level1 > 8) { level1 = 8 }
            }
            if (count2 <= 1) {
                level2 <- 1
            } else {
                level2 <- ceiling(9/5 * log(count2))
                if (level2 > 8) { level2 = 8 }
            }
            cl.colors[i] <- C[level1,level2]
        }
    }

    names(cl.colors) <- names(clonotype.counts1)
    cell.colors <- cl.colors[as.character(clonotypes)]
    cell.colors <- ifelse(is.na(cell.colors), gray, cell.colors)

    clonotype.residency.f <- factor(clonotype.residency, levels=c(
        tolower(sources[2]),
        tolower(sources[1]),
        toupper(sources[2]),
        toupper(sources[1]),
        "Dual"))  # 56974
    names(clonotype.residency.f) <- names(clonotype.counts1)

    cell.residency <- clonotype.residency.f[as.character(clonotypes)]
    cell.residency.f <- factor(cell.residency, levels=c(
        tolower(sources[2]),
        tolower(sources[1]),
        toupper(sources[2]),
        toupper(sources[1]),
        "Dual"))


    save(cell.colors, cl.colors, clonotype.residency.f, cell.residency.f,
         file=file.path(outdir, paste0(sample,".residency.colors.RData")))

}


load(samples)

# Only look at those who have both BM and WBC TCR data
# MM003BM-earlier	scTCR	MM003-earlier	BM	TCR_MM003BMCD138neg042516/
# MM003WBC-earlier	scTCR	MM003-earlier	WBC	TCR_MM003WBC042516CD15neg/
patients = samples %>%
    filter(Type == "scTCR") %>%
    group_by(Patient) %>%
    filter(n() == 2) %>%
    slice_head(n=1) %>%
    pull(Patient)

for (patient in patients) {
    residency.colors(patient)
}
