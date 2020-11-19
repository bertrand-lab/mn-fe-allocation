# quantitative comparison with previously published datasets

# must first run model_output_contour_plots.R to load and format data

library(dplyr)
library(ggplot2)
library(ggforce)
library(readxl)


# uptake rates ------------------------------------------------------------

mcquiad_uptake_rates <- read.csv('data/culture_data/mcquaid2018_isip2a_uptake_rates.csv')

model_uptake_rates_hist <- model_out_mnfe4_model1 %>% 
  mutate(fe_uptake_rate = 60*(fe_uptake*1e18)/6.0221409e+23) %>% 
  ggplot(aes(fe_uptake_rate)) +
  geom_histogram() +
  xlab('Fe Uptake Rate (amol Fe / hour) ') +
  ylab('Count') +
  xlim(0, 2) +
  ggtitle('Model Fe Uptake Rates') +
  theme_bw();model_uptake_rates_hist

mcquaid_uptake_rates_zoomed <- mcquiad_uptake_rates %>% 
  ggplot(aes(x = fe_uptake_rate_amol_cell_hr)) +
  geom_histogram() +
  ggtitle('Compiled Fe Uptake Rates: \nMcQuiad et al 2018, Nature') +
  xlab('Fe Uptake Rate (amol Fe / hour) ') +
  ylab('Count') +
  theme_bw() +
  # scale_x_log10() +
  # scale_x_log10() +
  # xlim(0, 2)
  facet_zoom(xlim = c(0, 2))

model_comparison_uptake_rates <- ggarrange(model_uptake_rates_hist, 
          mcquaid_uptake_rates_zoomed,
          nrow = 1, ncol = 2)

ggsave(model_comparison_uptake_rates, 
       filename = 'figures/uptake_rate_comparison.png',
       width = 9.83, height = 5.64)


# growth rate comparison --------------------------------------------------

loay_data <- read_xlsx('data/culture_data/LJ_EB_Temp_Fe_Interaction_data_LO-Letters.xlsx', 
                       sheet = 1)
culture_loay_data <- loay_data %>%
  ggplot(aes(x = growthrate)) +
  geom_histogram() +
  theme_bw() +
  xlim(0, 0.37) +
  ggtitle('Fragilariopsis cylindrus Culture Growth Rates:\nJabre and Bertrand 2020, L & O') +
  xlab('Growth Rate (per day)') +
  ylab('Count');culture_loay_data

model_growth_rates <- model_out_mnfe4_model1 %>% 
  ggplot(aes(x = u_trans)) +
  geom_histogram() +
  theme_bw() +
  xlim(0, 0.37) +
  ggtitle('Model Growth Rates\n ') +
  xlab('Growth Rate (per day)') +
  ylab('Count');model_growth_rates

growth_rate_comparison <- ggarrange(model_growth_rates, 
                                    culture_loay_data)

ggsave(growth_rate_comparison, 
       filename = 'figures/growth_rate_comparison.png')

## quotas from twining et al 

source('scripts/twining_fe_c_ratios.R')

twining_lower_bound <- umol_fe_mol_c_to_molecules_fe_per_cell(cell_volume = 258, umol_fe_cell = 5.5)
twining_upper_bound <- umol_fe_mol_c_to_molecules_fe_per_cell(cell_volume = 258, umol_fe_cell = 30)

twining_lower_bound*(1/6.022e23)*1e18
twining_upper_bound*(1/6.022e23)*1e18

twining_quota_comparison_plot <- model_out_mnfe4_model1 %>% 
  ggplot(aes(x = total_fe_amol)) +
  geom_histogram() +
  theme_bw() +
  ylab('Count') +
  ggtitle('Model Fe Quota (histogram)\nObserved Fe Quota (Twining et al 2004; dotted lines)') +
  geom_vline(xintercept = twining_lower_bound*(1/6.022e23)*1e18,
             lty = 2) +
  geom_vline(xintercept = twining_upper_bound*(1/6.022e23)*1e18,
             lty = 2) +
  xlab('Fe Quota (amol Fe per cell)')

ggsave(twining_quota_comparison_plot, 
       filename = 'figures/twining_quota_comparison.png',
       width = 7.67, height = 5.64)

####### converting peers and price Mn quota

# converting radius into m cubed, then converted m cubed into L
model_cell_volume <- 1000*(4/3)*pi*0.00000395^3 

convert_model_to_peers_price <- function(amol_per_cell,
                                         model_cell_vol_in = model_cell_volume){
  umol_per_litre <- amol_per_cell/model_cell_vol_in*(1/1e18)*(1e6/1)
  return(umol_per_litre)
}

model_out_mnfe4_model1$umol_fe_per_litre <- convert_model_to_peers_price(model_out_mnfe4_model1$total_fe_amol)
model_out_mnfe4_model1$umol_mn_per_litre <- convert_model_to_peers_price(model_out_mnfe4_model1$total_mn_amol)


peers_price_table1 <- data.frame(umol_mn_per_litre = c(191, 263, 21.5,
                                             78.6, 14.7, 6.9),
                                 umol_fe_per_litre = c(1460, 70, NA,
                                             NA, 1640,
                                             15))

peers_price_fe_quota_plot <- model_out_mnfe4_model1 %>% 
  ggplot(aes(x = umol_fe_per_litre)) +
  geom_histogram() +
  xlab('Fe Quota (uMol Fe per L)') +
  theme_bw() + 
  ylab('Count') +
  ggtitle('Model Fe Quota (histogram)\nObserved Fe Quota (Peers and Price, 2004; \ndotted lines)') +
  geom_vline(data = peers_price_table1,
            aes(xintercept = umol_fe_per_litre),
            lty = 2)

peers_price_mn_quota_plot <- model_out_mnfe4_model1 %>% 
  ggplot(aes(x = umol_mn_per_litre)) +
  geom_histogram() +
  xlab('Mn Quota (uMol Mn per L)') +
  theme_bw() + 
  ylab('Count') +
  ggtitle('Model Mn Quota (histogram)\nObserved Mn Quota (Peers and Price, 2004; \ndotted lines)') +
  geom_vline(data = peers_price_table1,
            aes(xintercept = umol_mn_per_litre),
            lty = 2)

peers_price_quota_comparison <- ggarrange(peers_price_fe_quota_plot, 
          peers_price_mn_quota_plot, nrow = 1, ncol = 2)

ggsave(peers_price_quota_comparison,
       filename = 'figures/peers_price_quota_comparison.png',
       width = 9.83, height = 5.64)


