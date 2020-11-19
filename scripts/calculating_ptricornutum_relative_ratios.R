
library(ggplot2)
library(magrittr)
library(dplyr)

setwd('../Google Drive/Projects/mn-fe-allocation/data/culture_proteomes/')

phaeo <- read.csv("Bertrand2013_S1a_Ptricornutum_proteins.csv")

sum_total <- phaeo %>% dplyr::mutate(felim_avg = (Felim.a + Felim.b + Felim.c)/3,
                                     replete_avg = (Replete.a + Replete.b + Replete.c)/3)
sum(sum_total$replete_avg)
sum(sum_total$felim_avg)

processed_df <- phaeo %>% dplyr::mutate(felim_avg = (Felim.a + Felim.b + Felim.c)/3,
                        replete_avg = (Replete.a + Replete.b + Replete.c)/3) %>% 
  group_by(coarse_grained) %>% 
  dplyr::summarize(sum_pep_int_felim = sum(felim_avg),
            sum_pep_int_replete = sum(replete_avg)) %>% 
  mutate(relative_expression = sum_pep_int_felim/sum_pep_int_replete)


processed_df %>% 
  ggplot(aes(x = coarse_grained, y = relative_expression)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  ggtitle('P. tricornutum proteins') +
  ylab('Relative Expression (Fe Limited / Replete)') +
  theme_bw()







