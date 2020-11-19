### making model to metaproteome comparison

# meta_abc_relative.R must be run prior to this

# change in taxonomic composition -----------------------------------------
library(scales)
# base_hues <- hue_pal()(7)
# base_hues[4] <- 'seagreen'
# base_hues[3] <- 'lightgoldenrod1'
# base_hues[7] <- 'cadetblue3'
# 
# tax_data <- read.csv("../ross-sea-meta-omics/data/go-tax-pipeline-output/pipeline_outtfg-all_tax_only.csv")
# 
# levels(tax_data$tax_assign)[20] <- "Haptophyta"                  
# tax_data$week <- paste0('Week ', tax_data$day)                  
# 
# tax_data_figure <- tax_data %>%
#   filter(filter == 3.0, 
#          tax_assign != 'promiscuous-peptide',
#          tax_assign != 'no-assignment-contig',
#          tax_assign %in% c('Cryptophyta',
#                            'Opisthokonta',
#                            'Other Straminopiles',
#                            'Chlorophyta',
#                            'Dinophyta',
#                            'Ciliophora',
#                            'Diatom',
#                            'Haptophyta')) %>% 
#   ggplot(aes(x = tax_assign, y = sum_pep_intensity_tax, fill = tax_assign)) + 
#   theme_bw() +
#   geom_col() + 
#   theme(axis.text.x = element_text(angle = 75, hjust = 1),
#         legend.position = "none") + 
#   facet_grid(~week) +
#   coord_flip() +
#   scale_fill_manual(values = base_hues) +
#   ylab("Sum of Group-Specific Peptide Intensity") +
#   xlab("Taxonomic Group Assignment") +
#   theme(strip.background = element_rect(fill = 'white'));tax_data_figure# +
# theme(axis.title.x = element_text(size = 16),
# axis.title.y = element_text(size = 16),
# axis.text.y = element_text(size = 16),
# strip.text = element_text(size = 16))



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
tfe_post_meta$coarse_grain <- rep('Iron Transporters', nrow(tfe_post_meta))

names(tfe_post_meta)[3] <- 'fold_change'
names(tn_post_meta)[3] <- 'fold_change'
names(p_post_meta)[3] <- 'fold_change'
names(r_post_meta)[3] <- 'fold_change'
names(a_post_meta)[3] <- 'fold_change'

post_meta_prots <- rbind(a_post_meta, p_post_meta, 
                         tn_post_meta, tfe_post_meta, 
                         r_post_meta)

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
                                                 'Iron Transporters', 
                                                 'Nitrogen Uptake and \nBiosynthesis',
                                                 'U')


post_pred_figure_df_2$model <- rep('Observations', nrow(post_pred_figure_df_2))
post_pred_figure_df_2$probability <- rep(0.5, nrow(post_pred_figure_df_2))

post_predi_value <- post_meta_prots %>% 
  # dplyr::filter(fold_change < 1.5) %>% 
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
       width = 9.3, height = 2.85)
  # theme(legend.key.size = unit(0.1, 'in'),
  #       legend.title = element_blank(),
  #       legend.position = c(0.7, 0.16),
  #       legend.text = element_text(size = 7),
  #       legend.box.background = element_rect(colour = 'grey10', fill = alpha('white', 0)));post_predi_value


# looking at the posterior predictive check



# Nunn_growth_rate_subplot ------------------------------------------------

# 
# nunn_growth_data <- data.frame(Fex = c(87, 2876),
#                                fe = c('Low Fe', 'High Fe'),
#                                u_trans = c(0.68, 1.71),
#                                u_trans_sd = c(0.26*2, 0.25*2),
#                                u_trans_se = c(0.26, 0.25), 
#                                probability = c(0.4, 0.4))
# nunn_growth_data$model <- rep('Culture Data', 2)
# 
# nunn_post_high <- read.csv("data/abc_intermediate/u_trans_2876_meta_combined_meta_cohen_nunn_h2_summarized.csv")
# nunn_post_low <- read.csv("data/abc_intermediate/u_trans_87_meta_combined_meta_cohen_nunn_h2_summarized.csv")
# 
# nunn_post_high$fe <- rep('High Fe', nrow(nunn_post_high))
# nunn_post_low$fe <- rep('Low Fe', nrow(nunn_post_low))
# 
# nunn_post_all <- rbind(nunn_post_high %>% rename(u_trans = u_trans_2876), 
#                        nunn_post_low %>% rename(u_trans = u_trans_87))
# nunn_post_all$model <- rep('Model', nrow(nunn_post_all))
# 
# nunn_data_comparison <- nunn_post_all %>%
#   ggplot(aes(x = u_trans, y = probability, fill = model)) +
#   geom_col(colour = 'grey30') +
#   facet_wrap(~fe, nrow = 2) +
#   theme_bw() +
#   ylim(0, 0.5) + 
#   xlab('Growth Rate (per day)')+
#   ylab('Posterior Probability') +
#   theme(strip.background = element_rect(fill = 'white')) +
#   geom_errorbarh(data = nunn_growth_data,
#                  aes(xmin = u_trans - u_trans_sd,
#                      xmax = u_trans + u_trans_sd), 
#                  height = 0.03, colour = 'grey30') +
#   geom_point(data = nunn_growth_data, 
#              aes(x = u_trans, y = probability, fill = model), 
#              shape = 22, size = 3) + #;nunn_data_comparison
#   scale_fill_manual(values = c('paleturquoise3', 'seagreen2'),
#                     labels = c('Cultures', 'Model')) +
#   guides(fill = guide_legend(override.aes = list(shape = 22,
#                                                  size = 0.1,
#                                                  colour = c('paleturquoise3', 'seagreen2'),
#                                                  fill = c('paleturquoise3', 'seagreen2')))) +
#   theme(legend.key.size = unit(0.1, 'in'),
#         legend.title = element_blank(),
#         legend.position = c(0.7, 0.16),
#         legend.text = element_text(size = 7),
#         legend.box.background = element_rect(colour = 'grey10', fill = alpha('white', 0)));nunn_data_comparison
# 
# # inferred growth rate ----------------------------------------------------
# 
# post_growth_out_w1 <- read.csv('data/abc_intermediate/u_trans_1010_meta_combined_meta_cohen_nunn_h2_summarized.csv')
# post_growth_out_w3 <- read.csv('data/abc_intermediate/u_trans_470_meta_combined_meta_cohen_nunn_h2_summarized.csv')
# 
# post_growth_out_w1$week <- rep('Week 1', nrow(post_growth_out_w1))
# post_growth_out_w3$week <- rep('Week 3', nrow(post_growth_out_w3))
# 
# post_model_out <- rbind(post_growth_out_w1 %>% rename(u_trans = u_trans_1010),
#                         post_growth_out_w3 %>% rename(u_trans = u_trans_470))
# 
# total_u_trans_posterior <- post_model_out %>% 
#   ggplot(aes(x = u_trans, y = probability, fill = week)) +
#   geom_col(alpha = 0.9, colour = 'grey30') +
#   facet_wrap(~week, nrow = 2) +
#   theme_bw() +
#   # xlim(0, 0.3) +
#   theme(legend.position = 'none', 
#         strip.background = element_rect(fill = 'White')) +
#   xlab("Inferred Growth Rate (per day)") +
#   ylab("Posterior Probability") +
#   scale_fill_manual(values = c('seagreen', 'seagreen1'));total_u_trans_posterior
# 
# 
# post_fe_out_w1 <- read.csv('data/abc_intermediate/total_fe_amol_1010_meta_summarized.csv')
# post_fe_out_w3 <- read.csv('data/abc_intermediate/total_fe_amol_470_meta_summarized.csv')
# 
# post_fe_out_w1$week <- rep('Week 1', nrow(post_fe_out_w1))
# post_fe_out_w3$week <- rep('Week 3', nrow(post_fe_out_w3))
# 
# post_model_out_fe <- rbind(post_fe_out_w1 %>% rename(total_fe_amol = total_fe_amol_1010),
#                         post_fe_out_w3 %>% rename(total_fe_amol = total_fe_amol_470))
# 
# 
# total_fe_posterior <- post_model_out_fe %>% 
#   ggplot(aes(x = total_fe_amol, y = probability, fill = week)) +
#   geom_col(alpha = 0.9, colour = 'grey30') +
#   facet_wrap(~week, nrow = 2) +
#   theme_bw() +
#   xlim(5, 50) +
#   theme(legend.position = 'none', 
#         strip.background = element_rect(fill = 'White')) +
#   xlab("Inferred Fe Quota (aMol per cell)") +
#   ylab("Posterior Probability") +
#   scale_fill_manual(values = c('seagreen', 'seagreen1'));total_fe_posterior
# 
# 
# 
# # compiling all figures ---------------------------------------------------
# 
# model_to_meta_p <- ggarrange(tax_data_figure, 
#           post_predi_value, 
#           ggarrange(nunn_data_comparison, 
#                     total_u_trans_posterior, 
#                     total_fe_posterior, 
#                     nrow = 1,
#                     labels = c('c', 'd', 'e')),
#           nrow = 3, labels = c('a', 'b'));model_to_meta_p


# ggsave(model_to_meta_p, 
#        filename = 'figures/model_to_meta_pabc.png',
#        width = 8.1, height = 8.45) 





