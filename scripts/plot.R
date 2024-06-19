library(ggplot2)
library(tidyverse)


mytheme <- theme(axis.line = element_line(), legend.key=element_rect(fill = NA),
                 text = element_text(size=22),# family = 'PT Sans'),
                 legend.key.width = unit(2, 'cm'),
                 legend.key.size = unit(1.5, 'lines'),
                 panel.background = element_rect(fill = "white"))


plot_prevalences <- function(dynamics_csv = "prevalence.csv", 
                             write_path = "figures/prevalence_dynamics.pdf") {
  
  read_csv(dynamics_csv) %>%
    pivot_longer(!step, names_to = "variable", values_to = "value") %>%
    
    # Remove unwanted columns that are not species prevalence.
    filter((variable != "total_grass") & 
           (variable != "total_fenced_area") &
           (variable != "landholder")) %>%
    
    ggplot(aes(x=step, y=value, color=variable)) + geom_line() + mytheme
  
  ggsave(write_path, width = 7, height = 4.5)
  
}


plot_map <- function(grass_df_or_csv, fence_csv, animal_csv,
                     write_path = "figures/map.pdf") {
  
  if (is_character(grass_df_or_csv)) {
    grass_df <- read_csv(grass_df_or_csv)
  } else {
    grass_df <- grass_df_or_csv
  }

  if (file.exists(fence_csv)) {
    fence_df = read_csv(fence_csv) |>
      mutate(landholder_id = as.factor(landholder_id)) |>
      group_by(landholder_id) |>
      # Calculate convex hulls, which may have NA if hull cardinality is less.
      mutate(chull_x = chull_w_backfill(x, y, "x"),
             chull_y = chull_w_backfill(x, y, "y")) |>
      # Remove chull rows with NA, drop old x & y cols, replace with chull.
      drop_na() |> mutate(x = chull_x, y = chull_y) |> select(!c(chull_x, chull_y))

    animal_df = read_csv(animal_csv)
    
    ggplot(grass_df, aes(x = x, y = y, fill = grass_layer)) + 
      geom_tile() + 
      scale_fill_gradient2(low = "#000000", mid = "#333300", high = "#63AF03") +
      geom_polygon(aes(x = x, y = y, group = landholder_id, color = "brown"), 
                   fence_df, 
                   fill = NA, size=2, show.legend = FALSE) +
      mytheme
  } else {
    ggplot(grass_df, aes(x = x, y = y, fill = grass_layer)) + 
      geom_tile() + 
      scale_fill_gradient2(low = "#000000", mid = "#333300", high = "#63AF03") +
      mytheme
  }
    
  ggsave(write_path, width = 7, height = 5.5)
}


  ##
  # Calculate convex hull and backfill with NULL to do efficient transformations
  # in call to mutate in plot_map above.
  #
  chull_w_backfill <- function(xvec, yvec, x_or_y) {
    
    stopifnot(length(xvec) == length(yvec))
    stopifnot(x_or_y %in% c("x", "y"))

    if (x_or_y == "x")
      vec <- xvec
    else if (x_or_y == "y")
      vec <- yvec

    retvec = rep(NA, length(xvec))
    chull_idxs = chull(xvec, yvec)
    retvec[1:length(chull_idxs)] = vec[chull_idxs]

    return (retvec)
  }
