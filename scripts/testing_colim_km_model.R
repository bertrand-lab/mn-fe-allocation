library(ggplot2)
library(magrittr)

mu_sub <- 0.1

mu_no_sub <- 1 - mu_sub

km1 <- 10
km2 <- 10
km3 <- 10
km4 <- 10

colim_comb <- function(mu_sub, s1, s2, 
                       km1 = 10, km2 = 10, 
                       km3 = 10, km4 = 10){
  mu_no_sub <- 1 - mu_sub
  part1 <- mu_no_sub*(s1/(s1 + km1)*s2/(s2 + km2))
  # part1 <- mu_no_sub*min(c(s1/(s1 + km1), s2/(s2 + km2)))
  part2 <- mu_sub*((s1/(km3 + s1 + (s2*km3/km4))) + (s2/(km4 + s2 + (s1*km4/km3))))
  
  return(part1 + part2)
}

range_colim <- function(mu_sub, s_range, 
                        km1 = 10, km2 = 10, km3 = 10, km4 = 10){
  
  df_out <- data.frame(growth_rate = numeric(),
                       s1 = numeric(),
                       s2 = numeric(),
                       mu_sub_val = numeric())
  
  for(mu_sub_val_i in 1:length(mu_sub)){
    
    for(s_range_i in 1:length(s_range)){
      
      growth_val_out <- colim_comb(mu_sub = mu_sub[mu_sub_val_i],
                                   s1 = s_range,
                                   s2 = s_range[s_range_i],
                                   km1 = km1, 
                                   km2 = km2, 
                                   km3 = km3,
                                   km4 = km4)
      temp_df <- data.frame(growth_rate = growth_val_out,
                            s1 = s_range,
                            s2 = rep(s_range[s_range_i], 
                                     length(s_range)),
                            mu_sub_val = rep(mu_sub[mu_sub_val_i], 
                                             length(s_range)))
      
      df_out <- rbind(temp_df, 
                      df_out)
    }
  }
  
  
  return(df_out)
}




diagnostic_plots_range_colim <- function(mu_sub,
                                         s_range,
                                         km1 = 10,
                                         km2 = 10, 
                                         km3 = 10,
                                         km4 = 10){
  
  out01 <- range_colim(mu_sub, s_range = s_range,
                       km1= km1, km2 = km2, km3 = km3, km4 = km4)
  
  contour_plots <- out01 %>% 
    ggplot(aes(x = s1, y = s2, 
               fill = growth_rate)) +
    # geom_point(size = 4) +
    geom_raster(interpolate = TRUE) +
    facet_grid(~factor(mu_sub_val)) +
    geom_contour(aes(z = growth_rate),
    colour = 'black', lwd = 0.4,
    alpha = 0.7) +
    theme_bw()
  
  
  colim_experiments <- range_colim(mu_sub = mu_sub, 
                                   s_range = c(km1, max(s_range)),km1 = km1, 
                                   km2 = km2, km3 = km3, km4 = km4)
  
  bar_plots <- colim_experiments %>% 
    ggplot(aes(x = interaction(s1, s2), y = growth_rate)) +
    geom_col() +
    theme_bw() +
    facet_grid(~factor(mu_sub_val))
  
  
  mean_colim_exp <- colim_experiments %>%
    filter(s1 != max(s_range) | s2 != max(s_range),
           s1 != km1 | s2 != km1) %>% 
    group_by(mu_sub_val) %>% 
    summarize(u_min = min(growth_rate),
              u_mean = mean(growth_rate))
  
  low_colim_exp <- colim_experiments %>% 
    filter(s1 == km1 & s2 == km1)
  
  inter_scores <- low_colim_exp %>% 
    inner_join(mean_colim_exp, 
               by = c('mu_sub_val')) %>% 
    mutate(int_index_min = u_min - growth_rate,
           int_index_mean = u_mean - growth_rate)
  
  min_index_score_plot <- inter_scores %>% 
    ggplot(aes(x = s1, y = int_index_min)) +
    geom_point() +
    theme_bw() +
    facet_grid(~mu_sub_val)
  
  mean_index_score_plot <- inter_scores %>% 
    ggplot(aes(x = s1, y = int_index_mean)) +
    geom_point() +
    theme_bw() +
    facet_grid(~mu_sub_val)
  
  ggarrange(contour_plots,
            bar_plots,
            ggarrange(min_index_score_plot, mean_index_score_plot), 
            nrow = 3, heights = c(2, 2, 1))
  
  # inter_scores$id_val <- rep(input_val, nrow(inter_scores))
  
  # return(inter_scores)
}


km_dif_plot <- diagnostic_plots_range_colim(mu_sub = c(0, 0.5, 0.9),
                             s_range = seq(0, 50, 0.5),
                             km2 = 1, km4 = 10)
km_dif_plot_both <- diagnostic_plots_range_colim(mu_sub = c(0, 0.5, 0.9),
                                            s_range = seq(0, 50, 0.5),
                                            km2 = 1, km1 = 1)
km_same_plot <- diagnostic_plots_range_colim(mu_sub = c(0, 0.5, 0.9),
                             s_range = seq(0, 50, 0.5),
                             km1 = 10, km3 = 10)

ggarrange(km_same_plot, 
          km_dif_plot, 
          km_dif_plot_both, nrow = 1)

out01 <- range_colim(c(0, 0.5, 0.9), 
                     s_range = seq(0, 50, 0.5), km1 = 1, km2 = 10)
out02 <- range_colim(c(0, 0.5, 0.9), 
                     s_range = seq(0, 50, 0.5), km1 = 1, km2 = 1)
out03 <- range_colim(c(0, 0.5, 0.9), 
                     s_range = seq(0, 50, 0.5), km1 = 10, km2 = 10)

out01$km_val <- rep('Unbalanced', nrow(out01))
out02$km_val <- rep('Balanced but Equally Decreased', nrow(out02))
out03$km_val <- rep('Balanced, Base Model', nrow(out03))

out_all <- rbind(out01, out02, out03)

mean_colim_exp <- out_all %>%
  mutate(s3 = s1 + s2) %>% 
  filter(s3 == 60, s1 == 10 | s2 == 10) %>%
  group_by(mu_sub_val,
           km_val) %>%
  summarize(u_min = min(growth_rate),
            u_mean = mean(growth_rate))

low_colim_exp <- out_all %>%
  filter(s1 == 10 & s2 == 10)

inter_scores <- low_colim_exp %>%
  inner_join(mean_colim_exp,
             by = c('mu_sub_val',
                    'km_val')) %>%
  mutate(int_index_min = u_min - growth_rate,
         int_index_mean = u_mean - growth_rate)

min_comparison_inter <- inter_scores %>% 
  ggplot(aes(x = factor(mu_sub_val), 
             y = int_index_min)) +
  geom_point(aes(shape = km_val), 
             size = 7, 
             alpha = 0.5) +
  theme_bw() +
  ylab('Interdependence Index (minimum)') +
  xlab('Degree of Colimitation/Interdependence')

mean_comparison_inter <- inter_scores %>% 
  ggplot(aes(x = factor(mu_sub_val), 
             y = int_index_mean)) +
  geom_point(aes(shape = km_val), 
             size = 7, 
             alpha = 0.5) +
  theme_bw() +
  ylab('Interdependence Index (mean)') +
  xlab('Degree of Colimitation/Interdependence')

metric_comparison <- ggarrange(min_comparison_inter, mean_comparison_inter, 
          common.legend = TRUE)
metric_comparison

contour_comparison_out <- out_all %>% 
  ggplot(aes(x = s1, y = s2, 
             fill = growth_rate)) +
  # geom_point(size = 4) +
  geom_raster(interpolate = TRUE) +
  facet_grid(km_val~factor(mu_sub_val)) +
  geom_contour(aes(z = growth_rate),
               colour = 'black', lwd = 0.4,
               alpha = 0.7) +
  theme_bw()


# out01 <- range_colim(c(0, 0.5, 0.9), s_range = seq(0, 50, 0.5))
