## model to metaproteome figure

# looking at the posterior predictive check

post_pred_figure_df <- inner_join(mnfe3_meta_trans_2_relative, 
           mnfe_combined_q_tmn, 
           by = c('cost_par', 'avail_space', 'epsilon_a')) %>% 
  filter(in_set == TRUE, coarse_grain != 'Tmn', day == 3)

levels(post_pred_figure_df$coarse_grain) <- c('Antioxidants', 'Photosystem Units', 
                                              'Iron Transporters', 'Nitrogen Uptake and \nBiosynthesis',
                                              'Ribosomes', 'Manganese Transpoters')

post_pred_figure_df_2 <- coarse_diatoms_specific_norm_relative %>% 
  filter(day == 3) %>% 
  filter(coarse_grains != 'U')

levels(post_pred_figure_df_2$coarse_grains) <- c('Antioxidants', 'Photosystem Units', 'Ribosomes',
                                              'Iron Transporters', 'Nitrogen Uptake and \nBiosynthesis',
                                              'U')

posteriors_fig <- post_pred_figure_df %>% 
  ggplot(aes(x = coarse_grain, y = relative_change))  +
  geom_boxplot(width = 0.2) +
  # facet_grid(~day) +
  # ggtitle('') + 
  theme_bw() +
  xlab('') +
  scale_y_log10() +
  ylab('Proteomic Pool Abundance Fold Change \n(Week 3 / Week 1 Abundance)') +
  geom_point(data = post_pred_figure_df_2, 
             aes(x = coarse_grains, y = relative_change), 
             size = 4,
             colour = 'darkblue',
             alpha = 0.8);posteriors_fig

avail_space_post_fig <- avail_space_post +
  xlab('Proportion of Membrane Space Available') +
  ylab('Density') +
  geom_density(fill = 'darkolivegreen4', alpha = 0.3);avail_space_post_fig


posterior_pars <- mnfe_combined_q_tmn %>%
  filter(in_set == TRUE)

post_model_out <- meta_model_out %>% inner_join(posterior_pars, 
                                                by = c('epsilon_a', 'avail_space', 'cost_par'))

post_model_out$day_form <- ifelse(test = post_model_out$day == 1, yes = 'Week 1', no = 'Week 3')

total_fe_amol_posterior <- post_model_out %>% 
  ggplot(aes(x = factor(day_form), y = total_fe_amol)) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  ylab('Fe Quota\n(aMol per cell)') +
  xlab('')

total_u_trans_posterior <- post_model_out %>% 
  ggplot(aes(x = factor(day_form), y = u_trans)) +
  geom_boxplot(width = 0.2) +
  theme_bw() +
  ylab('Growth Rate\n(per day)') +
  xlab('')


ggarrange_p1 <- ggarrange(avail_space_post_fig, 
          ggarrange(total_fe_amol_posterior, total_u_trans_posterior, 
                    nrow = 2, labels = c('c', 'd')),
          nrow = 1)

ggarrange(posteriors_fig, 
          ggarrange_p1, nrow = 2, labels = c('a', 'b'))
          
