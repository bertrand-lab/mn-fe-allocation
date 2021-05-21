# plotting model introduction

library(ggplot2)
library(dplyr)

model_out_mean <- model_out_mnfe4_model1 %>% 
  group_by(Fex, Mnx) %>% 
  summarize_all(mean)# %>% filter(Fex > 49, Mnx > 49)

# empty plot for model schematic ------------------------------------------

empty_plot <- iris %>% ggplot() + 
  theme(plot.background = element_rect(fill = 'white', 
                                       colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'white'));empty_plot


# plotting parameter posteriors -------------------------------------------

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
  ylab('Posterior\nProbability') +
  xlab('MnSOD Efficacy Parameter');eps_meta_p2

## all parameter distributions

par_estimates <- ggarrange(cost_meta_p1, avail_meta_p2, eps_meta_p2, labels = c('A', 'B', 'C'),
                           font.label = list(size = 9))

ggsave(par_estimates, filename = 'figures/parameter_estimate_posteriors.png',
       width = 7.02, height = 6.57)

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
  xlab('Diatom-Specific Proteomic Pool Abundance Fold Change (Week 3 / Week 1 Abundance)') +
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
                                # cost_meta_p1, 
                                # avail_meta_p2,
                                labels = c('B', 'C'),#, 'D', 'E'), 
                                align = 'hv', 
                                common.legend = TRUE, legend = 'top', 
                                nrow = 2, ncol = 1, font.label = list(size = 9));top_right_combined

p1_top <- ggarrange(empty_plot, 
                    top_right_combined,
                    heights = c(2.5, 0.5), widths = c(1.1, 0.5),
                    labels = c('A', ''),
                    nrow = 1, font.label = list(size = 9))
ggsave(p1_top,
       filename = 'figures/intro_to_model_only.png',
       width = 10.8*4/5, height = 4.83)

