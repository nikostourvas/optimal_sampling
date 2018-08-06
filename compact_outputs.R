# Plots of common sample sizes ####  

pdf(paste(id, "100_repl_compact.pdf", sep = "_"), 
    width = 24, height = 13.5, compress = FALSE)

my_palette <- brewer.pal(12, "Set3") # create a new palette
my_palette <- colorRampPalette(my_palette)(19) # how many colors this palette will have

Hobs_means_tidy_compact <- subset(Hobs_means_tidy, 
                            samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)

Hobs_means_tidy_compact$marker_num <- 
  factor(Hobs_means_tidy_compact$marker_num, levels = unique(
    as.character(Hobs_means_tidy_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_Ho <- expression(paste(
    "Mean Observed Heterozygosity (", H[o], ") by sample size & marker number"))}


y_axis_Ho <- seq(0.05, 0.95, 0.05)

p_Ho_tidy_compact <- ggplot(Hobs_means_tidy_compact, aes(x = samp_size, y = value)) + 
  geom_boxplot(aes(fill = samp_size)) + 
  facet_wrap(~ marker_num, nrow = 2)

p_Ho_tidy_compact + ggtitle(title_Ho) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = expression(paste("Mean Observed Heterozygosity (", H[o], ")" )), breaks = y_axis_Ho) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")





Hexp_means_tidy_compact <- subset(Hexp_means_tidy, 
                            samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)

Hexp_means_tidy_compact$marker_num <- 
  factor(Hexp_means_tidy_compact$marker_num, levels = unique(
    as.character(Hexp_means_tidy_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_He <- expression(paste(
    "Mean Expected Heterozygosity (", H[e], ") by sample size & marker number"))}

y_axis_He <- seq(0.05, 0.95, 0.05)

p_He_tidy_compact <- ggplot(Hexp_means_tidy_compact, aes(x = samp_size, y = value)) + 
  geom_boxplot(aes(fill = samp_size)) + 
  facet_wrap(~ marker_num, nrow = 2)

p_He_tidy_compact + ggtitle(title_He) + xlab("Sample Size") + 
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = expression(paste("Mean Expected Heterozygosity (", H[e], ")" )), breaks = y_axis_He) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")





ar_means_tidy_compact <- subset(ar_means_tidy, 
                                samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)

ar_means_tidy_compact$marker_num <- 
  factor(ar_means_tidy_compact$marker_num, levels = unique(
    as.character(ar_means_tidy_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
  
}else{
  title_ar <- expression(paste(
    "Mean Allelic richness (Ar) by sample size & marker number"))}
y_axis_ar <- 1:30

p_ar_tidy_compact <- ggplot(ar_means_tidy_compact, aes(x = samp_size, y = value)) + 
  geom_boxplot(aes(fill = samp_size)) + 
  facet_wrap(~ marker_num, nrow = 2)

p_ar_tidy_compact + ggtitle(title_ar) + xlab("Sample Size") + 
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Mean Allelic richness (Ar)", breaks = y_axis_ar) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")




perc_repl_detect_compact <- subset(perc_repl_detect, 
                              samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)

perc_repl_detect_compact$marker_num <- 
  factor(perc_repl_detect_compact$marker_num, levels = unique(
    as.character(perc_repl_detect_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_perc <- expression(paste(
    "Allele detection by sample size & marker number"))}


y_axis_perc <- seq(0, 100, 5)

p_perc_compact <- ggplot(perc_repl_detect_compact, aes(x = samp_size, group = 1)) + 
  geom_point(aes(y = percent_f1, colour = "percent_f1")) + 
  geom_point(aes(y = percent_f2, colour = "percent_f2")) +
  geom_line(aes(y = percent_f1, colour = "percent_f1", linetype = "percent_f2")) + 
  geom_line(aes(y = percent_f2, colour = "percent_f2", linetype = "percent_f1")) +
  geom_hline(yintercept = 95, linetype = "dashed") + 
  facet_wrap(~ marker_num, nrow =2)

p_perc_compact + ggtitle(title_perc) + xlab("Sample Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(
    name = "Replicates with all alleles > 0.05 (cyan solid line) & > 0.01 (red dashed line) detected (%)",
    breaks = y_axis_perc) +
  theme(legend.position = "none")


dev.off()
