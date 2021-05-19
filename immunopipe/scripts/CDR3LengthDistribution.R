library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)

immdata = "{{ in.immdata }}"
outdir = "{{ out.outdir }}"

dir.create(outdir, showWarnings = FALSE)

load(immdata)

data = c()
for (sample in names(immdata$data)) {
    sample_data = immdata$data[[sample]] %>%
        select(starts_with('CDR3.')) %>%
        separate_rows(1, 2, sep=";") %>%
        transmute(across(.fns=nchar, .names='{.col}.length'),
                  Sample=sample)

    dir.create(file.path(outdir, sample), showWarnings = FALSE)
    # nt
    png(file.path(outdir, sample, 'nt.png'), width=1200, height=1000, res=100)
    g = sample_data %>%
        group_by(Sample) %>%
        summarise(
            Count=as.numeric(table(CDR3.nt.length)),
            CDR3.nt.length=as.numeric(names(table(CDR3.nt.length)))
        ) %>%
        ggplot(aes(x=CDR3.nt.length, y=Count)) +
        geom_col(aes(fill=Sample), position=position_dodge()) +
        theme_prism(base_size = 16)
    print(g)
    dev.off()
    # aa
    png(file.path(outdir, sample, 'aa.png'), width=1200, height=1000, res=100)
    g = sample_data %>%
        group_by(Sample) %>%
        summarise(
            Count=as.numeric(table(CDR3.aa.length)),
            CDR3.aa.length=as.numeric(names(table(CDR3.aa.length)))
        ) %>%
        ggplot(aes(x=CDR3.aa.length, y=Count)) +
        geom_col(aes(fill=Sample), position=position_dodge()) +
        theme_prism(base_size = 16)
    print(g)
    dev.off()

    data = data %>% bind_rows(sample_data)
}

# all samples
dir.create(file.path(outdir, 'All_samples'), showWarnings = FALSE)
# nt
png(file.path(outdir, 'All_samples', 'all_samples.nt.png'),
    width=1200, height=1000, res=100)
data %>%
    group_by(Sample) %>%
    summarise(
        Count=as.numeric(table(CDR3.nt.length)),
        CDR3.nt.length=as.numeric(names(table(CDR3.nt.length)))
    ) %>%
    ggplot(aes(x=CDR3.nt.length, y=Count)) +
    geom_col(aes(fill=Sample), position=position_dodge()) +
    theme_prism(base_size = 16)
dev.off()
# aa
png(file.path(outdir, 'All_samples', 'all_samples.aa.png'),
    width=1200, height=1000, res=100)
data %>%
    group_by(Sample) %>%
    summarise(
        Count=as.numeric(table(CDR3.aa.length)),
        CDR3.aa.length=as.numeric(names(table(CDR3.aa.length)))
    ) %>%
    ggplot(aes(x=CDR3.aa.length, y=Count)) +
    geom_col(aes(fill=Sample), position=position_dodge()) +
    theme_prism(base_size = 16)
dev.off()
