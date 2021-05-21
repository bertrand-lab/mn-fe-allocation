### inferring in situ rates and biogeochemical metrics

# Nunn_growth_rate_subplot ------------------------------------------------

nunn_growth_data <- data.frame(Fex = c(87, 2876),
                               fe = c('Low Fe', 'High Fe'),
                               u_trans = c(0.68, 1.71),
                               u_trans_sd = c(0.26*2, 0.25*2),
                               u_trans_se = c(0.26, 0.25), 
                               probability = c(0.4, 0.4))
nunn_growth_data$model <- rep('Culture Data', 2)

nunn_post_high <- read.csv("data/abc_intermediate/u_trans_2876_meta_combined_meta_cohen_nunn_h2_summarized.csv")
nunn_post_low <- read.csv("data/abc_intermediate/u_trans_87_meta_combined_meta_cohen_nunn_h2_summarized.csv")

nunn_post_high$fe <- rep('High Fe', nrow(nunn_post_high))
nunn_post_low$fe <- rep('Low Fe', nrow(nunn_post_low))

nunn_post_all <- rbind(nunn_post_high %>% rename(u_trans = u_trans_2876), 
                       nunn_post_low %>% rename(u_trans = u_trans_87))
nunn_post_all$model <- rep('Model', nrow(nunn_post_all))

nunn_data_comparison <- nunn_post_all %>%
  ggplot(aes(x = u_trans, y = probability, fill = model)) +
  geom_col(colour = 'grey30') +
  facet_wrap(~fe, ncol = 2) +
  theme_bw() +
  ylim(0, 0.5) + 
  xlab('Growth Rate (per day)')+
  ylab('Posterior\nProbability') +
  theme(strip.background = element_rect(fill = 'white')) +
  geom_errorbarh(data = nunn_growth_data,
                 aes(xmin = u_trans - u_trans_sd,
                     xmax = u_trans + u_trans_sd), 
                 height = 0.03, colour = 'grey30') +
  geom_point(data = nunn_growth_data, 
             aes(x = u_trans, y = probability, fill = model), 
             shape = 22, size = 3) + #;nunn_data_comparison
  scale_fill_manual(values = c('grey80', '#5f9e60'),
                    labels = c('Cultures', 'Model')) +
  guides(fill = guide_legend(override.aes = list(shape = 22,
                                                 size = 0.1,
                                                 colour = c('grey80', '#5f9e60'),
                                                 fill = c('grey80', '#5f9e60')))) +
  theme(legend.key.size = unit(0.1, 'in'),
        legend.title = element_blank(),
        legend.position = c(0.8, 0.5),
        legend.text = element_text(size = 7),
        legend.box.background = element_rect(colour = 'grey10', fill = alpha('white', 0)));nunn_data_comparison

# inferred fe uptake rate -------------------------------------------------

post_uptake_w1 <- read.csv('data/abc_intermediate/fe_uptake_1010_meta_combined_meta_cohen_nunn_h2_summarized.csv')
post_uptake_w3 <- read.csv('data/abc_intermediate/fe_uptake_470_meta_combined_meta_cohen_nunn_h2_summarized.csv') 

post_uptake_w1$week <- rep('Week 1', nrow(post_uptake_w1))
post_uptake_w3$week <- rep('Week 3', nrow(post_uptake_w3))

post_uptake_out <- rbind(post_uptake_w1 %>% rename(fe_uptake = fe_uptake_1010),
                         post_uptake_w3 %>% rename(fe_uptake = fe_uptake_470))


fe_uptake_posterior <- post_uptake_out %>% 
  ggplot(aes(x = 60*(fe_uptake*1e18)/6.0221409e+23, 
             y = probability, fill = week)) +
  geom_col(alpha = 0.7, colour = 'grey30') +
  facet_wrap(~week, nrow = 2) +
  theme_bw() +
  xlim(0, 0.6) +
  theme(legend.position = 'none', 
        strip.background = element_rect(fill = 'White')) +
  xlab('Inferred Fe Uptake Rate\n(amol Fe / hour) ') +
  ylab("Posterior Probability") +
  scale_fill_manual(values = c('#7db1d6', '#376e9c'));fe_uptake_posterior


# inferred growth rate ----------------------------------------------------

post_growth_out_w1 <- read.csv('data/abc_intermediate/u_trans_1010_meta_combined_meta_cohen_nunn_h2_summarized.csv')
post_growth_out_w3 <- read.csv('data/abc_intermediate/u_trans_470_meta_combined_meta_cohen_nunn_h2_summarized.csv')

post_growth_out_w1$week <- rep('Week 1', nrow(post_growth_out_w1))
post_growth_out_w3$week <- rep('Week 3', nrow(post_growth_out_w3))

post_model_out <- rbind(post_growth_out_w1 %>% rename(u_trans = u_trans_1010),
                        post_growth_out_w3 %>% rename(u_trans = u_trans_470))

mn_fe_concentrations <- data.frame(week = c("Week 1", "Week 3"),
                                   environ_conc = c("0.26 nmol/kg dMn\n 1.01 nmol/kg dFe",
                                                    "0.21 nmol/kg dMn\n 0.47 nmol/kg dFe"),
                                   u_trans = c(0.18, 0.18),
                                   probability = c(0.3, 0.3))

total_u_trans_posterior <- post_model_out %>% 
  ggplot(aes(x = u_trans, y = probability, fill = week)) +
  geom_col(alpha = 0.7, colour = 'grey30') +
  facet_wrap(~week, nrow = 2) +
  theme_bw() +
  xlim(0.03, 0.3) +
  ylim(0, 0.35) +
  ylab("Posterior Probability") +
  theme(legend.position = 'none', 
        strip.background = element_rect(fill = 'White')) +
  xlab("Inferred Growth Rate\n(per day)") +
  geom_label(data = mn_fe_concentrations, 
            aes(label = environ_conc), 
            fill = 'white', size = 3, alpha = 0.6) +
  scale_fill_manual(values = c('#7db1d6', '#376e9c'));total_u_trans_posterior


# fe amol per cell --------------------------------------------------------

post_fe_out_w1 <- read.csv('data/abc_intermediate/total_fe_amol_1010_meta_combined_meta_cohen_nunn_h2_summarized.csv')
post_fe_out_w3 <- read.csv('data/abc_intermediate/total_fe_amol_470_meta_combined_meta_cohen_nunn_h2_summarized.csv')

post_fe_out_w1$week <- rep('Week 1', nrow(post_fe_out_w1))
post_fe_out_w3$week <- rep('Week 3', nrow(post_fe_out_w3))

post_model_out_fe <- rbind(post_fe_out_w1 %>% rename(total_fe_amol = total_fe_amol_1010),
                           post_fe_out_w3 %>% rename(total_fe_amol = total_fe_amol_470))


total_fe_posterior <- post_model_out_fe %>% 
  ggplot(aes(x = total_fe_amol, y = probability, fill = week)) +
  geom_col(alpha = 0.7, colour = 'grey30') +
  facet_wrap(~week, nrow = 2) +
  theme_bw() +
  xlim(5, 50) +
  theme(legend.position = 'none', 
        strip.background = element_rect(fill = 'White')) +
  xlab("Inferred Fe Quota\n(aMol per cell)") +
  ylab("Posterior Probability") +
  scale_fill_manual(values = c('#7db1d6', '#376e9c'));total_fe_posterior

empty_plot <- iris %>% ggplot() + 
  theme(plot.background = element_rect(fill = 'white', 
                                       colour = 'white'),
        panel.background = element_rect(fill = 'white', colour = 'white'));empty_plot


# aggregating figures into one --------------------------------------------

right_posteriors <- ggarrange(nunn_data_comparison,
                              total_u_trans_posterior,
                              total_fe_posterior,
                              fe_uptake_posterior,
                              ncol = 2, nrow = 2,
                              labels = c('B', 'C', 'D', 'E'),
                              align = 'hv', font.label = list(size = 9));right_posteriors

posteriors_metap_only <- ggarrange(total_u_trans_posterior, total_fe_posterior, fe_uptake_posterior,
          ncol = 3, nrow = 1, labels = c('C', 'D', 'E', align = 'hv'), font.label = list(size = 9))

right_posteriors_nunn_top <- ggarrange(empty_plot + 
                                         ggtitle('Cellular model reproduces diatom culture growth rates') + 
                                         # element_text() + 
                                         theme(plot.title = element_text(hjust = 0.5)),
nunn_data_comparison,
                                       empty_plot + ggtitle('Inferences about rates and elemental quotas in the field') + theme(plot.title = element_text(hjust = 0.5)),
                                       posteriors_metap_only, labels = c('B', ''),
                                       ncol = 1, nrow = 4, heights = c(0.2, 1, 0.2, 2.4), 
font.label = list(size = 9))

# top_posterior_panel <- ggarrange(empty_plot, right_posteriors, ncol = 2,
top_posterior_panel <- ggarrange(empty_plot, right_posteriors_nunn_top, ncol = 2,
                                 widths = c(1.3, 2), labels = c('A'), font.label = list(size = 9))

ggsave(top_posterior_panel, 
       filename = 'figures/posterior_panel_w_empty.png',
       width = 10.6, height = 5.33)
