# plotting model introduction

library(ggplot2)
library(dplyr)

model_out_mean <- model_out_mnfe4_model1 %>% 
  group_by(Fex, Mnx) %>% 
  summarize_all(mean)# %>% filter(Fex > 49, Mnx > 49)

# empty plot for model schematic ------------------------------------------

empty_plot <- model_out_mean %>% ggplot() + 
  theme(plot.background = element_rect(fill = 'white', 
                                       colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'white'));empty_plot


# plotting parameter posteriors -------------------------------------------

# cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_h_10_summarized.csv')
# eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_h_10_summarized.csv')
# avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_h_10_summarized.csv')
# cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_h_11_10000_summarized.csv')
# eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_h_11_10000_summarized.csv')
# avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_h_11_10000_summarized.csv')
# cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_h_12_20000_summarized.csv')
# eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_h_12_20000_summarized.csv')
# avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_h_12_20000_summarized.csv')
# cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_h_115_30000_summarized.csv')
# eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_h_115_30000_summarized.csv')
# # avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_h_115_30000_summarized.csv')
# cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_h_115_10000_dyntfe_summarized.csv')
# eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_h_115_10000_dyntfe_summarized.csv')
# avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_h_115_10000_dyntfe_summarized.csv')
# cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_combined_h40_summarized.csv')
# eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_combined_h40_summarized.csv')
# avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_combined_h40_summarized.csv')

cost_par_meta <- read.csv('data/abc_intermediate/cost_par_meta_combined_meta_cohen_nunn_h2_summarized.csv')
avail_meta <- read.csv('data/abc_intermediate/avail_space_meta_combined_meta_cohen_nunn_h2_summarized.csv')
eps_meta <- read.csv('data/abc_intermediate/epsilon_a_meta_combined_meta_cohen_nunn_h2_summarized.csv')

cost_meta_p1 <- cost_par_meta %>% 
  ggplot(aes(x = cost_par, 
             y = probability)) +
  geom_col() +
  # geom_point() +
  theme_bw() +
  ylab('Posterior\nProbability') +
  xlab('Fe / Mn Internal Cost\nParameter');cost_meta_p1

avail_meta_p2 <- avail_meta %>% 
  ggplot(aes(x = avail_space, y = probability)) +
  geom_col() +
  # geom_smooth() +
  theme_bw() +
  ylab('Posterior\nProbability') +
  xlab('Available Membrane\nSpace Parameter');avail_meta_p2

eps_meta_p2 <- eps_meta %>% 
  ggplot(aes(x = epsilon_a, y = probability)) +
  geom_col() +
  theme_bw() +
  ylab('') +
  xlab('MnSOD Efficacy Parameter');eps_meta_p2

# plotting fe transporters ------------------------------------------------

fe_transporters_fex <- model_out_mean %>% 
  filter(Mnx == 500) %>%
  ggplot(aes(x  = Fex, y = Tfe)) +
  geom_point(size = 3, alpha = 0.9, aes(fill = u_trans),
             pch = 21) +
  ylab('\nFe Transporters\n(per cell)') +
  labs(x = '[dFe] (pM)') +
  theme_bw() +
  theme(legend.position = c(0.8, 0.6),
        legend.box.background = element_rect(colour = 'black', fill = alpha('white', 0)),
        legend.background = element_rect(fill=alpha('white', 0.7)),
        legend.text = element_text(angle = 45));fe_transporters_fex


# growth rate plot --------------------------------------------------------

### this section required model_output_contour_plots.R to be run first

# growth_rate_plot <- plots10[[3]] +
#   ylim(0, 501) +
#   scale_fill_distiller(name = 'Growth Rate\n (per day)') +
#   theme(panel.grid = element_blank(),
#         legend.key.size = unit(0.3, "cm"),
#         legend.key.width = unit(0.3,"cm"),
#         legend.position = c(0.65, 0.6),
#         legend.text = element_text(angle = 45),
#         legend.box.background = element_rect(colour = 'black', fill = alpha('white', 0)),
#         legend.background = element_rect(fill=alpha('white', 
#                                                     0.8))) +
#   labs(title = NULL);growth_rate_plot

growth_rate_plot_2 <- plots10[[3]] +
  ylim(0, 501) +
  scale_fill_distiller(name = 'Growth Rate\n (per day)') +
  theme(panel.grid = element_blank(),
        legend.position = 'top',
        legend.text = element_text(angle = 45)) +
  labs(title = NULL);growth_rate_plot_2
# posterior predictive check metaproteome ---------------------------------

### reading in posteriors for protein pool relative change
a_post_meta <- read.csv('data/abc_intermediate/A_meta_combined_meta_cohen_nunn_h2_summarized.csv')
p_post_meta <- read.csv('data/abc_intermediate/P_meta_combined_meta_cohen_nunn_h2_summarized.csv')
tn_post_meta <- read.csv('data/abc_intermediate/Tn_meta_combined_meta_cohen_nunn_h2_summarized.csv')
tfe_post_meta <- read.csv('data/abc_intermediate/Tfe_meta_combined_meta_cohen_nunn_h2_summarized.csv')
r_post_meta <- read.csv('data/abc_intermediate/R_meta_combined_meta_cohen_nunn_h2_summarized.csv')

a_post_meta$coarse_grain <- rep('Antioxidants', nrow(a_post_meta))
p_post_meta$coarse_grain <- rep('Photosystem Units', nrow(p_post_meta))
tn_post_meta$coarse_grain <- rep('Nitrogen Uptake and \nBiosynthesis', nrow(tn_post_meta))
r_post_meta$coarse_grain <- rep('Ribosomes', nrow(tn_post_meta))
tfe_post_meta$coarse_grain <- rep('Iron Transporters and \nInternal Cost', nrow(tfe_post_meta))

names(tfe_post_meta)[3] <- 'fold_change'
names(tn_post_meta)[3] <- 'fold_change'
names(p_post_meta)[3] <- 'fold_change'
names(r_post_meta)[3] <- 'fold_change'
names(a_post_meta)[3] <- 'fold_change'

post_meta_prots <- rbind(a_post_meta, p_post_meta, 
                         tn_post_meta, tfe_post_meta, 
                         r_post_meta)


# coarse_diatoms_specific_norm_relative is a dataframe that needs to be read in using meta_abc_relative.R

post_pred_figure_df_2 <- coarse_diatoms_specific_norm_relative %>% 
  filter(day == 3) %>% 
  filter(coarse_grains != 'U')

post_meta_prots$model <- rep('Model', nrow(post_meta_prots))

post_pred_figure_df_2 <- coarse_diatoms_specific_norm_relative %>% 
  filter(day == 3) %>% 
  filter(coarse_grains != 'U')
post_pred_figure_df_2$model <- rep('Observations', nrow(post_pred_figure_df_2))

levels(post_pred_figure_df_2$coarse_grains) <- c('Antioxidants', 
                                                 'Photosystem Units', 
                                                 'Ribosomes',
                                                 'Iron Transporters and \nInternal Cost', 
                                                 'Nitrogen Uptake and \nBiosynthesis',
                                                 'U')


post_pred_figure_df_2$model <- rep('Observations', nrow(post_pred_figure_df_2))
post_pred_figure_df_2$probability <- rep(0.5, nrow(post_pred_figure_df_2))

post_predi_value <- post_meta_prots %>% 
  dplyr::filter(fold_change < 2) %>%
  ggplot(aes(x = fold_change, 
             y = probability, 
             fill = model)) +
  geom_col(fill = 'seagreen2', colour = 'grey30') +
  facet_grid(~coarse_grain) +
  # xlim(0.5, 1.5) +
  # coord_flip() +
  ylab('Posterior Probability') +
  theme_bw() +
  # geom_vline(xintercept = 1, colour = 'grey30') +
  xlab('Diatom-Specific Proteomic Pool Abundance Fold Change (Week 3 / Week 1 Abundance)') +
  # scale_fill_gradient2(low = 'seagreen', mid = 'seagreen1', high = 'seagreen',
  #                      midpoint = 1) +
  # geom_point(data = post_pred_figure_df_2 %>% 
  #              dplyr::rename(coarse_grain = coarse_grains,
  #                            fold_change = relative_change), 
  #                aes(xintercept = fold_change, 
  #                group = coarse_grain),
  #            size = 4, alpha = 0.7,
  #            shape = 'square') +
  geom_vline(data = post_pred_figure_df_2 %>% 
               dplyr::rename(coarse_grain = coarse_grains,
                             fold_change = relative_change),
             aes(xintercept = fold_change, group = coarse_grain), 
             lty = 2, lwd = 1.1, alpha = 0.7) +
  theme(legend.position = 'none', 
        strip.background = element_rect(fill = 'white')) + 
  guides(fill = guide_legend(override.aes = list(shape = c(22, 23), 
                                                 size = 0.1,
                                                 colour = c('blue', 'paleturquoise3'),
                                                 fill = c('blue', 'paleturquoise3'))));post_predi_value

ggsave(post_predi_value, 
       filename = 'figures/metaproteome_diatom_posterior.png',
       width = 9.82, height = 4.19)

## compiled figure

top_right_combined <- ggarrange(growth_rate_plot_2, 
                                fe_transporters_fex, 
                                cost_meta_p1, 
                                avail_meta_p2,
                                labels = c('b', 'c', 'd', 'e'), 
                                align = 'hv', 
                                common.legend = TRUE, legend = 'right');top_right_combined

p1_top <- ggarrange(empty_plot, 
                    top_right_combined,
                    heights = c(1.5, 1), labels = c('a', 
                                                    ''),
                    nrow = 1)
ggsave(p1_top,
       filename = 'figures/intro_to_model_only.png',
       width = 10.8, height = 4.83)

##### plot for fe internal limitation

# plotting underpinnings of growth limitation

# fei_lim <- model_out_mean %>% 
#   mutate(fei_mm = Mni/(10000 + Mni)*Fei/(10000 + Fei)) %>% 
#   ggplot(aes(x = u_trans, 
#              y = fei_mm,
#              fill = Fex)) +
#   # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
#   geom_point(size = 3, 
#              pch = 21,
#              # size = 3, 
#              alpha = 0.7, stroke = 1) +
#   theme_bw() +
#   ylim(0, 0.85) +
#   scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
#   # ylab('Multiplicative Internal Limitation \nof Fe and Mn') +
#   # ylab( expression(paste("Multiplicative Internal Limitation of Fe and Mn ", sigma,",", R^{2},'=',r2.value), italic(Fei))) +
#   ylab( expression(paste("Fe and Mn Internal Limitation", '(', Lambda, ')'))) +
#   annotate(geom = 'label', x = 0.09, y = 0.7, size = 3, alpha = 0.2,
#            label = expression(paste(Lambda, ' = ',
#                                     frac(italic(Fe)[italic(i)], 
#                                          italic(Fe)[italic(i)] + italic(K)[italic(Fe)[italic(i)]]), ~
#                                       frac(italic(Mn)[italic(i)], 
#                                            italic(Mn)[italic(i)] + italic(K)[italic(Mn)[italic(i)]])))) +
#   xlab('Growth Rate (per day)') +
#   # facet_grid(~Mnx) +
#   # geom_line(aes(group = Mnx)) +
#   theme(panel.grid = element_blank(),
#         legend.position = "None");fei_lim
# 
# model_out_mnfe4_model1 %>% 
#   mutate(fei_mm = Mni/(10000 + Mni)*Fei/(10000 + Fei)) %>% 
#   ggplot(aes(x = u_trans, 
#              y = fei_mm,
#              fill = Fex)) +
#   # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
#   geom_point(size = 3, 
#              pch = 21,
#              # size = 3, 
#              alpha = 0.7, stroke = 1) +
#   theme_bw() +
#   ylim(0, 0.85) +
#   scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
#   # ylab('Multiplicative Internal Limitation \nof Fe and Mn') +
#   # ylab( expression(paste("Multiplicative Internal Limitation of Fe and Mn ", sigma,",", R^{2},'=',r2.value), italic(Fei))) +
#   ylab( expression(paste("Fe and Mn Internal\nLimitation", '(', Lambda, ')'))) +
#   annotate(geom = 'label', x = 0.09, y = 0.7, size = 3, alpha = 0.2,
#            label = expression(paste(Lambda, ' = ',
#                                     frac(italic(Fe)[italic(i)], 
#                                          italic(Fe)[italic(i)] + italic(K)[italic(Fe)[italic(i)]]), ~
#                                       frac(italic(Mn)[italic(i)], 
#                                            italic(Mn)[italic(i)] + italic(K)[italic(Mn)[italic(i)]])))) +
#   xlab('Growth Rate (per day)') +
#   # facet_grid(~Mnx) +
#   # geom_line(aes(group = Mnx)) +
#   theme(panel.grid = element_blank(),
#         legend.position = "None")
# 
# 
# translation_lim <- model_out_mean %>% 
#   ggplot(aes(x = u_trans, 
#              y = gamma,
#              fill = Fex)) +
#   # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
#   geom_point(size = 3, 
#              pch = 21,
#              alpha = 0.7, stroke = 1) +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
#   ylab('Total Translational Output \n(amino acids per minute)') +
#   xlab('Growth Rate (per day)') +
#   theme(panel.grid = element_blank(),
#         legend.position = c(0.55, 0.1),
#         legend.background = element_rect(fill=alpha('white', 0)),
#         legend.text = element_text(size = 8),
#         legend.direction = 'horizontal');translation_lim
# 
# model_out_mnfe4_model1 %>% 
#   ggplot(aes(x = u_trans, 
#              y = gamma,
#              fill = Fex)) +
#   # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
#   geom_point(size = 3, 
#              pch = 21,
#              alpha = 0.7, stroke = 1) +
#   theme_bw() +
#   scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
#   ylab('Total Translational Output \n(amino acids per minute)') +
#   xlab('Growth Rate (per day)') +
#   theme(panel.grid = element_blank(),
#         legend.position = c(0.55, 0.1),
#         legend.background = element_rect(fill=alpha('white', 0)),
#         legend.text = element_text(size = 8),
#         legend.direction = 'horizontal')

#  aggregated figure ------------------------------------------------------



# p1_top_right <- ggarrange(growth_rate_plot, 
#                               ggarrange(cost_meta_p1, 
#                                         avail_meta_p2,
#                                         labels = c('d', 'e')),
#                                           nrow = 2, 
#                                           heights = c(1.5, 1),
#                                           labels = c('b', ''))
# 
# 
# p1_top_left <- ggarrange(empty_plot, 
#                          fe_transporters_fex,
#                          heights = c(1.5, 1), labels = c('a', 'c'),
#                          nrow = 2)
# 
# p1_top <- ggarrange(p1_top_left, 
#                          p1_top_right,
#                          heights = c(1.5, 1), labels = c('', ''),
#                          nrow = 1)


# top_panel_growth_transporters <- ggarrange(growth_rate_plot_2, 
#           fe_transporters_fex, 
#           common.legend = TRUE,
#           legend = 'right',
#           labels = c('b', 
#                      'c'))
# 
# bottom_panel_posteriors <- ggarrange(cost_meta_p1, 
#                                      avail_meta_p2,
#                                      align = 'hv',
#                                      labels = c('d', 'e'));bottom_panel_posteriors
# 
# top_right_panel_sub <- ggarrange(top_panel_growth_transporters, 
#                                  bottom_panel_posteriors,
#                                  nrow = 2,
#                                  align = 'hv');top_right_panel_sub

# intro_to_model_plot <- ggarrange(p1_top, post_predi_value, nrow = 2, heights = c(2, 1), labels = c('', 'f'))
# 
# intro_model_plot_with_internal_limitation <- ggarrange(intro_to_model_plot, 
#                                                        ggarrange(fei_lim, translation_lim, nrow = 1,
#                                                                  labels = c('g', 'h')),
#                                                        nrow = 2, heights = c(2, 1))
# 
# ggsave(intro_model_plot_with_internal_limitation, 
#        filename = "figures/intro_to_model.png", 
#        width = 8.89, 
#        height = 9.61)

