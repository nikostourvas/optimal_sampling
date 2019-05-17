# CV function for all parameters

cv_opt <- function(input){
  sd(input, na.rm = TRUE) /
    mean(input, na.rm = TRUE) * 100
}


library(dplyr)

test0 <- ar_means_tidy %>% 
  group_by(samp_size, marker_num) %>%
  summarise(cv_opt(value))

library(ggplot2)
p_test_tidy <- ggplot(test0, 
                      aes(x = samp_size, 
                          y = `cv_opt(value)`)) + 
  geom_point() + 
  facet_wrap(~ marker_num, nrow = 2)

p_test_tidy