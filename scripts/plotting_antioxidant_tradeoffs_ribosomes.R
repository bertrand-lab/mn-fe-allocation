library(forcats)

# plotting ribosomae antioxidant trade-off

a_r_trade_off <- model_out_mnfe4_model1 %>%
  filter(Fex == 1000) %>% 
  group_by(Mnx, Fex) %>%
  summarize_all(mean) %>%
  ggplot(aes(x = R, y = A/2, 
             fill = Mnx)) +
  theme_bw() +
  geom_point(size = 4, pch = 21) +
  xlab("Ribosomes (per cell)") +
  scale_fill_distiller(name = 'dMn (pM)', 
                       palette = 'GnBu') +
  theme(legend.position = c(0.72, 0.68), axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_blank()) +
  ylab("Antioxidants (MnSOD per cell)");a_r_trade_off


translation_allocation_tradeoff <- model_out_mnfe4_model1 %>% 
  filter(Fex == 1 | Fex == 50 | Fex == 100,
         Mnx == 1000) %>%
  ggplot(aes(x = beta_r, 
             y = beta_a,
             colour = u_trans,
             shape = factor(Fex))) +
  geom_point(size = 4, alpha = 0.8) +
  xlab("Proportion of Ribosomes synthesizing Ribosomes") +
  ylab("Proportion of Ribosomes synthesizing\nAntioxidants") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.74),
        legend.direction = 'horizontal') +
  scale_color_continuous('Growth Rate (per day)') +
  scale_shape_discrete('[dFe] (pM)') +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         shape = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(legend.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(0,0,0,0, unit="cm"));translation_allocation_tradeoff


model_out_mnfe4_model1$Fex_label <- paste('[dFe] =', model_out_mnfe4_model1$Fex, sep = ' ')

# plotting the distribution of Fe quotas
amol_per_cell_tradeoff <- model_out_mnfe4_model1 %>% 
  filter(Fex == 50 | Fex == 100 | Fex == 1,
         Mnx == 1000) %>%
  ggplot() +
  geom_density(aes(total_fe_amol)) +
  facet_wrap(~fct_reorder(Fex_label, Fex), nrow = 3) + 
  xlab('Fe Quota (aMol per cell)') +
  theme_bw() +
  theme(strip.background = element_rect(fill = 'white')) +
  ylab('Density')

## aggregating plots
lower_plots_internal_lim_rearrange <- ggarrange(a_r_trade_off, translation_allocation_tradeoff, 
          amol_per_cell_tradeoff,
          widths = c(1, 2, 1), nrow = 1,
          labels = c('a', 'b', 'c'))

ggsave(lower_plots_internal_lim_rearrange, 
       filename = "figures/internal_lim_consequences_3.png",
       width = 10.7, 
       height = 4.16)

