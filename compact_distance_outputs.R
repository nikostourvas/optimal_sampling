# Plots of common sample sizes ####  

pdf(paste(id, "distances_100_repl_compact.pdf", sep = "_"), 
    width = 24, height = 13.5, compress = FALSE)

my_palette <- brewer.pal(12, "Set3") # create a new palette
my_palette <- colorRampPalette(my_palette)(19) # how many colors this palette will have


fst_tidy_compact <- subset(fst_tidy,
                           samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)
fst_tidy_compact$marker_num <- 
  factor(fst_tidy_compact$marker_num, levels = unique(
    as.character(fst_tidy_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_Fst <- expression(paste(
    "Mean pairwise Fst by sample size & marker number"))}

y_axis_fst <- seq(-0.200, 0.950, 0.005)

p_fst_tidy <- ggplot(fst_tidy_compact, aes(x = samp_size, y = original_values)) +
  geom_boxplot(aes(fill = samp_size)) +
  facet_wrap(~ marker_num, nrow = 2)

p_fst_tidy + ggtitle(title_Fst) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Mean pairwise Fst", breaks = y_axis_fst) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")







dch_tidy_compact <- subset(dch_tidy,
                           samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)

dch_tidy_compact$marker_num <- 
  factor(dch_tidy_compact$marker_num, levels = unique(
    as.character(dch_tidy_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_dch <- expression(paste(
    "Cavalli-Sforza and Edwards Chord distance by sample size & marker number"))}

y_axis_dch <- seq(-0.10, 0.95, 0.02)

p_dch_tidy <- ggplot(dch_tidy_compact, aes(x = samp_size, y = Dch)) +
  geom_boxplot(aes(fill = samp_size)) +
  facet_wrap(~ marker_num, nrow = 2)

p_dch_tidy + ggtitle(title_dch) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Cavalli-Sforza and Edwards Chord distance", breaks = y_axis_dch) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")



jost_tidy_compact <- subset(jost_tidy,
                            samp_size == 25 | samp_size == 30 | samp_size == 50 | samp_size == 75)
jost_tidy_compact$marker_num <- 
  factor(jost_tidy_compact$marker_num, levels = unique(
    as.character(jost_tidy_compact$marker_num)))

if(id == "Abies_DE_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek adult population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian adult population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek regeneration population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian regeneration population of ", 
    italic("A. alba")))
}else if (id == "Abies_DE_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german seed population of ", 
    italic("A. alba")))
}else if (id == "Abies_GR_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek seed population of ", 
    italic("A. borisii-regis")))
}else if (id == "Abies_SL_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian seed population of ", 
    italic("A. alba")))
  
}else if (id == "Fagus_DE_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Adult"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian adult population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Regen"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian regeneration population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_DE_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for german seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_GR_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for greek seed population of ", 
    italic("F. sylvatica")))
}else if (id == "Fagus_SL_Seed"){
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number for slovenian seed population of ", 
    italic("F. sylvatica")))
  
}else{
  title_jost <- expression(paste(
    "Jost's D by sample size & marker number"))}

y_axis_jost <- seq(-0.20, 0.95, 0.02)

p_jost_tidy <- ggplot(jost_tidy_compact, aes(x = samp_size, y = original_values)) +
  geom_boxplot(aes(fill = samp_size)) +
  facet_wrap(~ marker_num, nrow = 2)

p_jost_tidy + ggtitle(title_jost) + xlab("Sample Size") +
  scale_fill_manual(values = my_palette) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = "Jost's D", breaks = y_axis_jost) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")

dev.off()