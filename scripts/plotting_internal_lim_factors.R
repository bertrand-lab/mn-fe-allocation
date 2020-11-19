fei_lim_alone <- model_out_mnfe4_model1 %>% 
  mutate(fei_mm = Fei/(10000 + Fei)) %>% 
  group_by(Mnx, Fex) %>% 
  summarize_all(mean) %>% 
  filter(Mnx == 1000 | Mnx == 1 | Mnx == 50) %>% 
  ggplot(aes(x = u_trans, 
             y = fei_mm,
             fill = Fex,
             shape = as.factor(Mnx))) +
  # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
  geom_point(size = 3, 
             # pch = 21,
             # size = 3, 
             alpha = 0.7, stroke = 1) +
  theme_bw() +
  ylim(0, 1) +
  # ylim(0, 0.85) +
  scale_shape_manual(values = c(21, 22, 23), name = 'dMn (pM)') +
  scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
  # ylab('Multiplicative Internal Limitation \nof Fe and Mn') +
  # ylab( expression(paste("Multiplicative Internal Limitation of Fe and Mn ", sigma,",", R^{2},'=',r2.value), italic(Fei))) +
  ylab( expression(paste("Fe Internal Limitation ", '(', Lambda, ')'))) +
  annotate(geom = 'label', x = 0.20, y = 0.25, size = 4, alpha = 0.2,
           label = expression(paste(Lambda, ' = ',
                                    frac(italic(Fe)[italic(i)], 
                                         italic(Fe)[italic(i)] + italic(K)[italic(Fe)[italic(i)]])))) +
  xlab('Growth Rate (per day)') +
  # facet_grid(~Mnx) +
  # geom_line(aes(group = Mnx)) +
  theme(legend.position = "None");fei_lim_alone

aa_lim_alone <- model_out_mnfe4_model1 %>% 
  mutate(aa_mm = aa/(10000 + aa)) %>% 
  group_by(Mnx, Fex) %>% 
  summarize_all(mean) %>% 
  filter(Mnx == 1000 | Mnx == 1 | Mnx == 50) %>% 
  ggplot(aes(x = u_trans, 
             y = aa_mm,
             fill = Fex,
             shape = as.factor(Mnx))) +
  # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
  geom_point(size = 3, 
             # pch = 21,
             # size = 3, 
             alpha = 0.7, stroke = 1) +
  theme_bw() +
  ylim(0, 1) +
  # ylim(0, 0.85) +
  scale_shape_manual(values = c(21, 22, 23), name = 'dMn (pM)') +
  scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
  # ylab('Multiplicative Internal Limitation \nof Fe and Mn') +
  # ylab( expression(paste("Multiplicative Internal Limitation of Fe and Mn ", sigma,",", R^{2},'=',r2.value), italic(Fei))) +
  ylab( expression(paste("Amino Acid (AA) Internal Limitation ", '(', Lambda, ')'))) +
  annotate(geom = 'label', x = 0.20, y = 0.25, size = 4, alpha = 0.2,
           label = expression(paste(Lambda, ' = ',
                                    frac(italic(AA), 
                                         italic(AA) + italic(K)[italic(AA)])))) +
  xlab('Growth Rate (per day)') +
  # facet_grid(~Mnx) +
  # geom_line(aes(group = Mnx)) +
  theme(legend.position = "None",
        plot.margin = margin(10, 10, 10, 30));aa_lim_alone


e_lim_alone <- model_out_mnfe4_model1 %>% 
  mutate(e_mm = e/(10000 + e)) %>% 
  group_by(Mnx, Fex) %>% 
  summarize_all(mean) %>% 
  filter(Mnx == 1000 | Mnx == 1 | Mnx == 50) %>% 
  ggplot(aes(x = u_trans, 
             y = e_mm,
             fill = Fex,
             shape = as.factor(Mnx))) +
  # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
  geom_point(size = 3, 
             # pch = 21,
             # size = 3, 
             alpha = 0.7, stroke = 1) +
  theme_bw() +
  ylim(0, 1) +
  # ylim(0, 0.85) +
  scale_shape_manual(values = c(21, 22, 23), name = 'dMn (pM)') +
  scale_fill_distiller(palette = 'RdYlBu', name = 'dFe (pM)') +
  # ylab('Multiplicative Internal Limitation \nof Fe and Mn') +
  # ylab( expression(paste("Multiplicative Internal Limitation of Fe and Mn ", sigma,",", R^{2},'=',r2.value), italic(Fei))) +
  ylab( expression(paste("Energetic (e) Internal Limitation", '(', Lambda, ')'))) +
  annotate(geom = 'label', x = 0.20, y = 0.25, size = 4, alpha = 0.2,
           label = expression(paste(Lambda, ' = ',
                                    frac(italic(e)[italic(i)], 
                                         italic(e)[italic(i)] + italic(K)[italic(e)])))) +
  xlab('Growth Rate (per day)') +
  # facet_grid(~Mnx) +
  # geom_line(aes(group = Mnx)) +
  theme(legend.position = "None");e_lim_alone

r_lim_alone<- model_out_mnfe4_model1 %>% 
  group_by(Mnx, Fex) %>% 
  summarize_all(mean) %>% 
  filter(Mnx == 1000 | Mnx == 1 | Mnx == 50) %>% 
  ggplot(aes(x = u_trans, 
             y = R,
             shape = as.factor(Mnx),
             fill = Fex)) +
  # geom_smooth(method = 'lm', colour = 'grey40', fill = 'grey80') +
  geom_point(size = 3, 
             # pch = 21,
             # size = 3, 
             alpha = 0.7, stroke = 1) +
  theme_bw() +
  # ylim(0, 1) +
  # ylim(0, 0.85) +
  scale_shape_manual(values = c(21, 22, 23), 
                     name = 'dMn (pM)') +
  scale_fill_distiller(palette = 'RdYlBu', 
                       name = 'dFe (pM)') +
  guides(shape = guide_legend(title.position="top", 
                              title.hjust = 0.5,
                              order = 1),
         fill = guide_colourbar(title.position="top", title.hjust = 0.5)) +
  # ylab('Multiplicative Internal Limitation \nof Fe and Mn') +
  # ylab( expression(paste("Multiplicative Internal Limitation of Fe and Mn ", sigma,",", R^{2},'=',r2.value), italic(Fei))) +
  ylab( expression(paste("Ribosomes (per cell)"))) +
  xlab('Growth Rate (per day)') +
  # facet_grid(~Mnx) +
  # geom_line(aes(group = Mnx)) +
  theme(legend.position = c(0.7, 0.77),
        legend.direction = "horizontal");r_lim_alone
guide_legend(order = 2)
internal_lim_fig <- ggarrange(empty_plot, fei_lim_alone, 
                              r_lim_alone,
                              aa_lim_alone,
                              nrow = 2, ncol = 2,
                              labels = c('b', 'c', 'd', 'e'),
          align = 'hv',
          common.legend = TRUE,
          legend = 'bottom')

ggsave(internal_lim_fig, filename = "figures/internal_lim_fig.png",
       width = 8.45, 
       height = 6.45)



internal_lim_costs <- ggarrange(growth_experiments_heatmap, 
          internal_lim_fig, 
          ncol = 2, nrow = 1,
          labels = c('a'),
          widths = c(1, 2))

ggsave(internal_lim_costs,
       filename = 'figures/internal_lim_costs.png',
       width = 12.4,
       height = 6.56)


