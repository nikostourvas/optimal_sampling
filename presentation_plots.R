# First run the plots section from optimal_sampling.R

# Filter results for 3 loci classes

library(dplyr)
Hobs_means_tidy3 <- Hobs_means_tidy %>% 
  filter(marker_num %in% c('5 markers', '8 markers', '11 markers'))


my_palette <- brewer.pal(12, "Set3") # create a new palette
my_palette <- colorRampPalette(my_palette)(19) # how many colors this palette will have

y_axis_Ho <- seq(0.05, 0.95, 0.05)

p_Ho_tidy <- ggplot(Hobs_means_tidy3, aes(x = samp_size, y = value)) + 
  geom_boxplot(aes(fill = samp_size)) + 
  facet_wrap(~ marker_num, nrow = 1)

p_Ho <- p_Ho_tidy + labs(x = NULL) +
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = expression(paste(" ", H[o], "" )), 
                     breaks = y_axis_Ho) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")

p_Ho
  
# He ----------------------------------------------------------------------

Hexp_means_tidy3 <- Hexp_means_tidy %>% 
  filter(marker_num %in% c('5 markers', '8 markers', '11 markers'))

y_axis_He <- seq(0.05, 0.95, 0.05)

p_He_tidy <- ggplot(Hexp_means_tidy3, aes(x = samp_size, y = value)) + 
  geom_boxplot(aes(fill = samp_size)) + 
  facet_wrap(~ marker_num, nrow = 1)

p_He <- p_He_tidy + xlab("Sample Size") + 
  scale_fill_manual(values = my_palette) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) + 
  theme(text = element_text(size = 18)) +
  theme(title = element_text(size = 18)) +
  scale_y_continuous(name = expression(paste("",H[e],"") ), 
                     breaks = y_axis_He) +
  stat_summary(fun.y=mean, geom="point", shape=4, size=2, color="black", fill="black")



# combine -----------------------------------------------------------------

library(gridExtra)

png("DE_A_Ho_He.png",
    width = 720,
    height = 1280,
    units = "px")

grid.arrange(p_Ho, p_He,                
             ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,1), c(2,2)))

dev.off()


# AR ----------------------------------------------------------------------


