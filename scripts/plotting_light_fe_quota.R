# light relationship with TM quota

library(ggplot2)
library(magrittr)

# must first run model_output_contour_plots.R

light_values <- read_in_data_loop(data_name_string = 'mnfe4_model4_aug19_post_inference_base_light', 
                  dir_to_look = 'data/model_output/', 
                  string_to_not_include = 'baseline')


light_out_mean <- light_values %>% 
  group_by(Fex, Mnx, I) %>% 
  summarize_all(mean)# %>% filter(Mnx < 1001, Fex < 1001)

light_out_plot <- light_out_mean %>%
  ggplot(aes(x = I, y = total_fe_amol)) +
  geom_point(aes(shape = factor(Mnx)),
             size = 4, alpha = 0.7) +
  facet_grid(~Fex) +
  theme_bw() +
  scale_shape_discrete('dMn (pm)') +
  theme(legend.position = c(0.8, 0.8),
        strip.background = element_rect(fill = 'white')) +
  # xlab('umol Einsteins * m^-2 *s^-1') +
  ylab('Fe Quota (amol Fe per cell)') +
  xlab(label = expression(paste('uEin m'^-2, ' s'^-1)))

ggsave(light_out_plot, 
       filename = 'figures/light_fe_quota.png',
       width = 7.67, 
       height = 5.64)
