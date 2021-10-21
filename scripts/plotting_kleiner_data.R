library(ggplot2)
library(magrittr)
library(dplyr)
library(scales)

kleiner_proc_short <- read.csv('data/kleiner_data/kleiner_data_processed_short.csv')
kleiner_proc_long <- read.csv('data/kleiner_data/kleiner_data_processed_long.csv')

kleiner_proc_long$lc <- rep('Long Chromatographic Run', nrow(kleiner_proc_long))
kleiner_proc_short$lc <- rep('Short Chromatographic Run', nrow(kleiner_proc_short))

kleiner_proc <- rbind(kleiner_proc_short, kleiner_proc_long)

kleiner_data_comparison <- kleiner_proc %>% 
  # filter(protein_amount < 200) %>% 
  ggplot(aes(x = rescale(protein_amount), 
             y = rescale(sum_precursor_intensity),
             size = number_peps)) +
  geom_smooth(method = 'lm') +
  geom_point() +
  facet_grid(~lc) +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0)+
  guides(size=guide_legend(title="Number of Peptides")) +
  scale_x_log10() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  xlab('Rescaled Protein Amount per Organism') +
  ylab('Rescaled Summed Peptide Intensity per Organism') +
  # geom_label(aes(label = Species)) +
  geom_point()

ggsave(kleiner_data_comparison, 
       filename = 'figures/kleiner_data_plotted.png',
       width = 8.86, height = 4.73)
