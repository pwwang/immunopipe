library(dplyr)
library(ggplot2)
library(cowplot)
library(tibble)
library(ggVennDiagram)
library(ggprism)

samples = "{{ in.samples }}"
tcrcount_dir = "{{ in.tcr_counts }}"
rc_colors = "{{ in.rc_colors }}"
outdir = "{{ out.outdir }}"

dir.create(outdir, showWarnings = FALSE)


exponent <- function (x) {
    floor(log10(abs(x)))
}

mantissa <- function (x) {
    mant <- log10(abs(x))
    10^(mant - floor(mant))
}


residency_one_patient = function(patient) {
    print(paste('Handling', patient, '...'))
    countfile = file.path(tcrcount_dir, paste0(patient, '.clonotype.counts.RData'))
    colorfile = file.path(rc_colors, paste0(patient, '.residency.colors.RData'))

    load(countfile)
    load(colorfile)

    tibble(
        Subject = patient,
        Source1 = sources[1],
        Source2 = sources[2],
        Source1_Unique = sum( tolower(clonotype.residency) == tolower(sources[1]) ),
        Source2_Unique = sum( tolower(clonotype.residency) == tolower(sources[2]) ),
        Shared = sum( tolower(clonotype.residency) == tolower('Dual') ),
        Total = length(clonotype.residency),
        Source1_Unique_Perc = Source1_Unique / Total,
        Source2_Unique_Perc = Source2_Unique / Total,
        Shared_Perc = Shared / Total,
    )
}


load(samples)

# Generate residency table
# Only look at those who have both BM and WBC TCR data
# MM003BM-earlier	scTCR	MM003-earlier	BM	TCR_MM003BMCD138neg042516/
# MM003WBC-earlier	scTCR	MM003-earlier	WBC	TCR_MM003WBC042516CD15neg/
patients = samples %>%
    filter(Type == "scTCR") %>%
    group_by(Patient) %>%
    filter(n() == 2) %>%
    slice_head(n=1) %>%
    pull(Patient)

residency_table = c()
for (patient in patients) {
    residency_table = bind_rows(residency_table, residency_one_patient(patient))
}

# Subject	Source1	Source2	Source1_Unique	Source2_Unique	Shared	Total	Source1_Unique_Perc	Source2_Unique_Perc	Shared_Perc
# MM003-earlier	BM	WBC	31	536	3	570	0.0543859649122807	0.940350877192982	0.00526315789473684

write.table(
    residency_table,
    file.path(outdir, 'residency_table.txt'),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)

# Venn diagram

plotlist = list()
for (patient in patients) {
    sources = residency_table %>%
        filter(Subject == patient) %>%
        select(matches('Source\\d+$')) %>%
        unlist()

    ndata = residency_table %>%
        filter(Subject == patient) %>%
        select(4:6) %>%
        unlist()

    venndata = list()
    if (ndata[3] == 0) {
        venndata[[sources[1]]] = c(
            paste0(sources[1], seq_len(ndata[1]))
        )
        venndata[[sources[2]]] = c(
            paste0(sources[2], seq_len(ndata[2]))
        )
    } else {
        venndata[[sources[1]]] = c(
            paste0('Shared', seq_len(ndata[3])),
            paste0(sources[1], seq_len(ndata[1]))
        )
        venndata[[sources[2]]] = c(
            paste0('Shared', seq_len(ndata[3])),
            paste0(sources[2], seq_len(ndata[2]))
        )
    }
    g = ggVennDiagram(venndata, label_alpha=0) +
        ggtitle(patient) +
        theme(plot.title = element_text(hjust = 0.5))
    plotlist[[length(plotlist) + 1]] = g
}

venn_fig = file.path(outdir, 'Venn.png')
height = ceiling(length(patients) / 4) * 500
png(venn_fig, width=2400, height=height, res=100)
g = plot_grid(plotlist=plotlist, ncol=4, align = "hv")
print(g)
dev.off()

# Pie charts
plotlist = list()
for (patient in patients) {
    sources = residency_table %>%
        filter(Subject == patient) %>%
        select(matches('Source\\\\d+$')) %>%
        unlist()

    sdata = residency_table %>%
        filter(Subject == patient) %>%
        select(8:10) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column('Source') %>%
        mutate(
            Source=c(
                paste(sources[1], "Unique %"),
                paste(sources[2], "Unique %"),
                "Shared %"
            )
            # ypos=cumsum(V1) - .5*V1
        )

    g = ggplot(sdata, aes(x="", y=V1, fill=Source)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + ggtitle(patient) +
      theme(plot.title = element_text(hjust = 0.5))
    plotlist[[length(plotlist) + 1]] = g
}

pie_fig = file.path(outdir, 'Pie.png')
height = ceiling(length(patients) / 4) * 500
png(pie_fig, width=2400, height=height, res=100)
g = plot_grid(plotlist=plotlist, ncol=4, align = "hv")
print(g)
dev.off()

# Scatter plots
plot_scatter = function(patient) {
    countfile = file.path(tcrcount_dir, paste0(patient, '.clonotype.counts.RData'))
    colorfile = file.path(rc_colors, paste0(patient, '.residency.colors.RData'))

    load(countfile)
    load(colorfile)

    dual = which(clonotype.counts1 > 0 & clonotype.counts2 > 0)
    if (length(dual) <= 2) {
        test = list(estimate=NA, p.value=NA)
    } else {
        test <- cor.test(log(clonotype.counts1[dual]),log(clonotype.counts2[dual]))
    }
    sum.counts1 <- sum(clonotype.counts1)
    sum.counts2 <- sum(clonotype.counts2)

    counts1.norm <- jitter(1+clonotype.counts1, amount=0.25)/sum.counts1
    counts2.norm <- jitter(1+clonotype.counts2, amount=0.25)/sum.counts2

    oo <- sample(length(counts1.norm))
    plotdata = data.frame(x=counts1.norm[oo], y=counts2.norm[oo])
    names(plotdata) = sources
    plotdata$color = cl.colors[oo]

    if (tolower(sources[1]) %in% c('wbc', 'blood', 'normal', 'control')) {
        x = sources[1]
        y = sources[2]
        sum_counts1 = sum.counts1
        sum_counts2 = sum.counts2

    } else {
        x = sources[2]
        y = sources[1]
        sum_counts2 = sum.counts1
        sum_counts1 = sum.counts2
    }
    xbreaks = c(1/sum_counts1, 0.001+1/sum_counts1, 0.01+1/sum_counts1, 0.1+1/sum_counts1)
    ybreaks = c(1/sum_counts2, 0.001+1/sum_counts2, 0.01+1/sum_counts2, 0.1+1/sum_counts2)

    minx = min(plotdata[[x]])
    miny = min(plotdata[[y]])
    maxx = max(plotdata[[x]])
    maxy = max(plotdata[[y]])
    color = plotdata$color
    names(color) = color
    patient = as.character(patient)
    n.formatted <- formatC(length(oo), format="f", big.mark=",", digits=0)
    r.formatted <- format(test$estimate,digits=2,scientific=F)
    if (is.na(test$p.value)) {
        subtitle = bquote(italic(n)[D] == .(length(dual)) ~~ italic(r) == .(r.formatted) ~~ italic(P) == "NA")
    } else if (test$p.value < 1e-4) {
        P.mant <- format(mantissa(test$p.value),digits=2)
        P.exp <- exponent(test$p.value)
        subtitle = bquote(italic(n)[D] == .(length(dual)) ~~ italic(r) == .(r.formatted) ~~ italic(P) == .(P.mant) %*% 10^.(P.exp))
    } else {
        P.formatted <- format(test$p.value,digits=2)
        subtitle = bquote(italic(n)[D] == .(length(dual)) ~~ italic(r) == .(r.formatted) ~~ italic(P) == .(P.formatted))
    }
    ggplot(plotdata) +
        geom_point(aes_string(x=x, y=y, color='color'), shape=1) +
        scale_color_manual(values=color) +
        scale_x_continuous(
            trans="log2",
            limits=c(minx, maxx),
            breaks=xbreaks,
            labels=c("0","0.001","0.01","0.1")
        ) +
        scale_y_continuous(
            trans="log2",
            limits=c(miny, maxy),
            breaks=ybreaks,
            labels=c("0","0.001","0.01","0.1")
        ) +
        theme_prism(base_size = 16) +
        theme(legend.position = "none") +
        labs(
            title=bquote(.(patient)~(italic(n) == .(n.formatted))),
            subtitle=subtitle
        ) +
        geom_segment(
            data=data.frame(
                x=c(
                    1.5/sum_counts1,
                    minx,
                    1.5/sum_counts1,
                    minx,
                    2.5/sum_counts1
                ),
                xend=c(
                    maxx, # diagnal
                    maxx,  # horizontal
                    1.5/sum_counts1,   # vertical
                    1.5/sum_counts1,  # horizontal short
                    2.5/sum_counts1  # vertical short

                ),
                y=c(
                    1.5/sum_counts2,
                    1.5/sum_counts2,
                    miny,
                    2.5/sum_counts2,
                    miny
                ),
                yend=c(
                    maxy,
                    1.5/sum_counts2,
                    maxy,
                    2.5/sum_counts2,
                    1.5/sum_counts2
                )
            ),
            aes(x=x, y=y, xend=xend, yend=yend), color='gray'
        )

}

plotlist = list()
for (patient in patients) {
    plotlist[[length(plotlist) + 1]] = plot_scatter(patient)
}

scatter_fig = file.path(outdir, 'Scatter.png')
height = ceiling(length(patients) / 4) * 500
png(scatter_fig, width=2000, height=height, res=100)
g = plot_grid(plotlist=plotlist, ncol=4, align = "hv")
print(g)
dev.off()
