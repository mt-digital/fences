library(tidyr)

mytheme <- theme(axis.line = element_line(), legend.key=element_rect(fill = NA),
                 text = element_text(size=22),# family = 'PT Sans'),
                 legend.key.width = unit(2, 'cm'),
                 legend.key.size = unit(1.5, 'lines'),
                 panel.background = element_rect(fill = "white"))

plot_prevalences <- function(dynamics_csv) {
  read_csv(dynamics_csv) %>%
    pivot_longer(!step, names_to = "variable", values_to = "value") %>%
    filter((variable != "total_grass") & (variable != "total_fenced_area") & (variable != "landholder")) %>%
    ggplot(aes(x=step, y=value, color=variable)) + geom_line() + mytheme
}