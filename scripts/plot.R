library(ggplot2)
library(tidyverse)

mytheme <- theme(axis.line = element_line(), legend.key=element_rect(fill = NA),
                 text = element_text(size=22),# family = 'PT Sans'),
                 legend.key.width = unit(2, 'cm'),
                 legend.key.size = unit(1.5, 'lines'),
                 panel.background = element_rect(fill = "white"))

plot_prevalences <- function(dynamics_csv = "test_dynamics.csv", 
                             write_path = "figures/species_dynamics.pdf") {
  
  read_csv(dynamics_csv) %>%
    pivot_longer(!step, names_to = "variable", values_to = "value") %>%
    filter((variable != "total_grass") & 
           (variable != "total_fenced_area") &
           (variable != "landholder")) %>%
    ggplot(aes(x=step, y=value, color=variable)) + geom_line() + mytheme
  
  ggsave(write_path, width = 7, height = 4.5)
  
}


grass_tile <- function(grass_df_or_csv = "test_grass.csv", 
                       write_path = "figures/grass_test.pdf") {
  
  if (is_character(grass_df_or_csv)) {
    grass_df <- read_csv(grass_df_or_csv)
  } else {
    grass_df <- grass_df_or_csv
  }
  
  ggplot(grass_df, aes(x = x, y = y, fill = grass)) + 
    geom_tile() + 
    scale_fill_gradient2(low = "#000000", mid = "#333300", high = "#63AF03") +
    mytheme
  
  ggsave(write_path, width = 7, height = 5.5)
}
