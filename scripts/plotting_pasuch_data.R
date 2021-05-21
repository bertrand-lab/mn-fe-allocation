## plotting the culture data from Mn and Fe interactions in model organisms and proteomic model

library(ggplot2)
library(magrittr)

pasuch <- read.csv("data/culture_data/Pasuch2019_growth_data.csv")
peers <- read.csv("data/culture_data/peers_and_price2004_data.csv")
peers_oceanica <- read.csv('data/culture_data/peers_price_oceanica.csv')

pasuch_plot <- pasuch %>% 
  ggplot(aes(x = dFe, y = u)) +
  geom_point(size = 4, aes(colour = dMn)) +
  geom_label(aes(label = dMn), nudge_y = 0.005) +
  theme_bw() +
  scale_colour_continuous('[dMn] (nM)') +
  ylab('Growth Rate (per day)') +
  xlab('[dFe] (nM)') +
  labs(title = 'Pasuch et al 2019, PLOS One',
       subtitle = expression(paste(italic(Chaetoceros), ' ', 
                                   italic(debilis))));pasuch_plot

peers_plot <- peers %>% 
  ggplot(aes(x = fe_prime, y = u)) +
  geom_point(size = 4, aes(colour = mn_prime)) +
  geom_label(aes(label = mn_prime), nudge_y = 0.05) +
  scale_colour_continuous('Mn Prime (nM)') +
  theme_bw() +
  ylab('Growth Rate (per day)') +
  xlab('Fe Prime (nM)') +
  labs(title = 'Peers and Price 2004, L&O',
       subtitle = expression(paste(italic(Thalassiosira), ' ', 
                                   italic(pseudonana))));peers_plot

peers_plot_ocean <- peers_oceanica %>% 
  ggplot(aes(x = fe_prime, y = mu)) +
  geom_point(size = 4, aes(colour = mn_prime)) +
  geom_label(aes(label = mn_prime), nudge_y = 0.05) +
  scale_colour_continuous('Mn Prime (nM)') +
  theme_bw() +
  ylab('Growth Rate (per day)') +
  xlab('Fe Prime (nM)') +
  labs(title = 'Peers and Price 2004, L&O',
       subtitle = expression(paste(italic(Thalassiosira), ' ', 
                                   italic(oceanica))));peers_plot_ocean


frag_model <- ros10_out_mean %>% 
  filter(Mnx == 1 | Mnx == 3000, Fex == 1 | Fex == 3000) %>% 
  ggplot(aes(x = Fex, y = u_trans)) +
  geom_point(aes(colour = Mnx),
             size = 4, alpha = 0.7) +
  theme_bw() +
  scale_color_continuous(name = '[dMn] (pM)') +
  xlab('[dFe] (pM)') +
  ylab('Growth Rate (per day)') +
  labs(title = 'Proteomic Allocation Model',
       subtitle = expression(paste(italic(Fragilariopsis), ' ', 
                                   italic(cylindrus))))

growth_rate_comparison <- ggarrange(frag_model,
          pasuch_plot,
          peers_plot,
          peers_plot_ocean,
          labels = c('A', 'B', 'C', 'D'),
          nrow = 2, ncol = 2, font.label = list(size = 9))  

ggsave(growth_rate_comparison,
       filename = 'figures/quant_growth_mn_fe_low_high.png',
       width = 10.2, height = 6.57)  





