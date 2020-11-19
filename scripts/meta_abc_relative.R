# metaproteomic ABC analysis

# abc coarse grains from metaproteome
library(reshape2)
library(dplyr)
library(ggplot2)
library(magrittr)

read_in_data_loop_ls <- function(data_name_string, dir_to_look = "../mn-fe-allocation/data/abc_out/"){
  # '''
  # this function reads in the data from the model output in a loop
  # it also parses the file name and adds the replicate number as
  # well as the ROS_sensitivity parameter which is included in the
  # file name
  # '''
  
  model_output_list <- dir(dir_to_look)
  print(model_output_list)
  data_names <- model_output_list[grepl(pattern = data_name_string,
                                        x = model_output_list)]
  
  print(model_output_list)
  # print(data_names)
  master_dataframe <- read.csv(paste0(dir_to_look, data_names[1]))
  # print(master_dataframe)
  # print(get_model_run_replicate(data_names[1]))
  # print(get_ros_sens_par(data_names[1]))
  
  model_rep <- strsplit(data_names[1], split = "_")[[1]][6]
  model_rep_parse <- strsplit(model_rep, 
                              split = ".", 
                              fixed = TRUE)[[1]][1]
  
  master_dataframe$model_rep <- model_rep_parse
  
  # print(get_ros_sens_par(data_names[1]))
  if(length(data_names) > 1){
    for(i in 2:length(data_names)){
      print(data_names[i])
      temp_df <- read.csv(paste0(dir_to_look, data_names[i]))
      
      model_rep <- strsplit(data_names[1], split = "_")[[1]][6]
      model_rep_parse <- strsplit(model_rep, 
                                  split = ".", 
                                  fixed = TRUE)[[1]][1]
      
      temp_df$model_rep <- model_rep_parse
      
      master_dataframe <- rbind(master_dataframe, temp_df)
    }
  }
  
  return(master_dataframe)
}

# read in the model ABC output
# meta_model_out1 <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__april9')
# meta_model_out2 <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__april8')
# meta_model_out3 <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__april13')
# meta_model_out <- rbind(meta_model_out1, meta_model_out2, meta_model_out3)
# 
meta_model_out1 <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__april28')
meta_model_out2 <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__may')

meta_model_out <- rbind(meta_model_out1, meta_model_out2)

# meta_model_out <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__july', dir_to_look = '../mn-fe-allocation/data/abc_out/')

# adding the replicate number and the day comparison for the meta model out
meta_model_out$rep <- rep(1:(nrow(meta_model_out)/2), each=2)
meta_model_out$day <- rep(c(1:2), times = nrow(meta_model_out)/2)

# transforming the model output to be the same as the metaproteomic output
# including Dyn cost here use used because unknown or unassigned proteins will be observed
mnfe3_meta_trans <- meta_model_out %>%
  dplyr::select(A_aa, P_aa, Tfe_aa, Tn_aa, R_aa, Tmn_aa, Fex, 
                cost_par, avail_space,
                epsilon_a) %>% 
  ungroup() %>% 
  melt(id.vars = c("Fex", "cost_par", "avail_space", "epsilon_a"), 
       variable.name = c("coarse_grain"),
       value.name = c("aa_val")) %>% 
  group_by(cost_par, avail_space, epsilon_a, Fex) %>% 
  summarize(sum_aa = sum(aa_val))

# the no dyn cost is used here because that's the best comparison with the observations
mnfe3_meta_trans_2 <- meta_model_out %>%
  dplyr::select(A_aa, P_aa, Tfe_aa, Tn_aa, R_aa, Tmn_aa, Fex, cost_par, avail_space,
                epsilon_a) %>% 
  ungroup() %>% 
  melt(id.vars = c("Fex", "cost_par", "avail_space", "epsilon_a"), 
       variable.name = c("coarse_grain"),
       value.name = c("aa_val")) %>% 
  inner_join(mnfe3_meta_trans, 
             by = c("Fex", "cost_par", "avail_space", "epsilon_a")) %>% 
  mutate(aa_val_frac = aa_val/sum_aa,
         aa_val_frac_80per = aa_val*0.8/sum_aa)

mnfe3_meta_trans_2$day <- rep(c(1, 3), nrow(mnfe3_meta_trans_2)/2)

# renaming factor levels
levels(mnfe3_meta_trans_2$coarse_grain) <- c('A', 'P', 'Tfe', 'Tn', 'R', 'Tmn')

# checking out a uniform prior influences proteome distribution
mnfe3_meta_trans_2 %>%
  group_by(coarse_grain) %>% 
  sample_n(10000) %>% 
  ggplot(aes(x = coarse_grain, y = aa_val_frac)) +
  facet_grid(~Fex) +
  geom_jitter(width = 0.2, alpha = 0.2)

### transform the data to be relative increase or decrease
day_1_vals <- mnfe3_meta_trans_2 %>%
  dplyr::filter(day == 1) %>% 
  dplyr::select(cost_par, avail_space, epsilon_a, coarse_grain, aa_val_frac) %>% 
  dplyr::rename(norm_factor = aa_val_frac)

mnfe3_meta_trans_2_relative <- mnfe3_meta_trans_2 %>% 
  inner_join(day_1_vals, by = c('cost_par', 'avail_space', 
                                'epsilon_a', 'coarse_grain')) %>%
  # dplyr::filter(cost_par > 0.2319033, cost_par < 0.2319035) # checking that the joinging worked as anticipated
  mutate(relative_change = aa_val_frac/norm_factor)


# read in the diatom proteome output from metaproteome
# coarse_diatoms_specific <- read.csv('../ross-sea-meta-omics/data/correction-factor-output/coarse_diatoms_tfg_all_specific_manually_edited.csv')
coarse_diatoms_specific_not_all_four <- read.csv('../ross-sea-meta-omics/data/correction-factor-output/coarse_diatoms_tfg_all_specific.csv')
# coarse_diatoms_specific <- read.csv('../ross-sea-meta-omics/data/correction-factor-output/coarse_diatoms_tfg_all_specific.csv')

# coarse_diatoms_specific <- coarse_diatoms_specific %>% dplyr::filter(grpnorm_assign == 'Fragilariopsis')

peptides_observed_all_days <- coarse_diatoms_specific_not_all_four$peptide %>% 
  table() %>% 
  as.data.frame() %>% 
  dplyr::rename(peptide = '.')

peptides_with_obs_all_four_days <- peptides_observed_all_days[peptides_observed_all_days$Freq == 4, ]$peptide

coarse_diatoms_specific <- coarse_diatoms_specific_not_all_four[coarse_diatoms_specific_not_all_four$peptide %in% peptides_with_obs_all_four_days, ]

# format the metaproteome output to be relevant to model

# normalize peptide abundances to total abundance per day
total_diatom_abundance_per_day <- coarse_diatoms_specific %>% 
  group_by(day, filter) %>% 
  summarize(total_diatom_prot = sum(db_abund_avg))

coarse_diatoms_specific_norm <- coarse_diatoms_specific %>% 
  inner_join(total_diatom_abundance_per_day, by = 'day') %>% 
  mutate(diatom_norm_pep_int = db_abund_avg/total_diatom_prot) %>% 
  group_by(coarse_grains, day) %>% 
  summarize(sum_pep_int = sum(diatom_norm_pep_int))

coarse_diatoms_specific_norm_day_1 <- coarse_diatoms_specific_norm %>% 
  dplyr::filter(day == 1) %>% 
  dplyr::select(coarse_grains, sum_pep_int) %>% 
  dplyr::rename(norm_factor = sum_pep_int)

coarse_diatoms_specific_norm_relative <- coarse_diatoms_specific_norm %>% 
  inner_join(coarse_diatoms_specific_norm_day_1, by = 'coarse_grains') %>% 
  mutate(relative_change = sum_pep_int/norm_factor)

coarse_diatoms_specific_norm_relative %>%
  filter(day != 2,day != 4) %>% 
  ggplot(aes(x = day, y = relative_change)) +
  geom_point() +
  facet_grid(~coarse_grains)

coarse_grain_diatom_meta_plot <- coarse_diatoms_specific_norm %>% 
  # filter(coarse_grains == 'Tn') %>% 
  ggplot(aes(x = day, y = sum_pep_int)) +
  geom_point() +
  facet_grid(~coarse_grains) +
  scale_y_log10() +
  theme_bw() +
  ylab('Sum of Peptide Intensity per Protein Pool') +
  xlab('Day');coarse_grain_diatom_meta_plot

coarse_grain_diatom_meta_plot_rel <- coarse_diatoms_specific_norm_relative %>% 
  # filter(coarse_grains == 'Tn') %>% 
  ggplot(aes(x = day, y = relative_change)) +
  geom_point() +
  facet_grid(~coarse_grains) +
  scale_y_log10() +
  theme_bw() +
  ylab('Relative Change in Sum of Peptide Intensity per Protein Pool') +
  xlab('Day');coarse_grain_diatom_meta_plot_rel


ggsave(coarse_grain_diatom_meta_plot, filename = 'figures/meta_proteomic_plotting_percentage.png', width = 10.4, height = 4.65)
ggsave(coarse_grain_diatom_meta_plot_rel, filename = 'figures/meta_proteomic_plotting_relative.png', width = 10.4, height = 4.65)

# considering the total amount of metaproteome derived diatom protein, but excluding unknown
mnfe_combined_relative_no_tmn_pos <- mnfe3_meta_trans_2_relative %>% 
  dplyr::rename(coarse_grains = coarse_grain,
                model_relative_change = relative_change) %>% 
  inner_join(coarse_diatoms_specific_norm_relative, by = c('day', 'coarse_grains')) %>% 
  mutate(sq_dif = (model_relative_change - relative_change)^2) %>%
  group_by(cost_par, avail_space, epsilon_a) %>%
  summarize(sum_sq_dif = sqrt(sum(sq_dif)))

# par sets that have relative change for Tmn that are positive
tmn_pos_vals <- mnfe3_meta_trans_2_relative %>% 
  filter(coarse_grain == 'Tmn', day == 3) %>% 
  mutate(pos_change = relative_change > 1) %>% 
  dplyr::select(cost_par, avail_space, epsilon_a, pos_change)

# getting the relative change output to get the posteriors (for python pabc_functions)
relative_output_formatted <- dcast(mnfe3_meta_trans_2_relative %>% 
                                     dplyr::filter(day == 3), 
      cost_par + epsilon_a + avail_space ~ coarse_grain, 
      value.var = "relative_change")

# get corresponding model output for a given par set
growth_rate_formatted <- meta_model_out %>% dplyr::select(cost_par, epsilon_a, 
                                                   avail_space, fe_uptake, 
                                                   mn_uptake, u_trans,
                                                   total_fe_amol, 
                                                   total_mn_amol, Fex) %>% 
  melt(id.vars = c("cost_par", "epsilon_a", "avail_space", "Fex"),
       meaure.vars = c("fe_uptake", "mn_uptake", "u_trans", 
                       "total_fe_amol", "total_mn_amol")) %>% 
  dcast(cost_par + epsilon_a + avail_space ~ variable + Fex)


mnfe_combined_relative <- mnfe_combined_relative_no_tmn_pos %>% 
  inner_join(relative_output_formatted, 
             by = c('cost_par', 'avail_space', 'epsilon_a')) %>% 
  inner_join(growth_rate_formatted,
             by = c('cost_par', 'avail_space', 'epsilon_a'))

# writing par set file
write.csv(mnfe_combined_relative, 
          file = '../mn-fe-allocation/data/abc_intermediate/meta_abc_par_sets.csv',
          row.names = FALSE)


######### plotting abc output analysis

percent5_quantile_relative <- quantile(mnfe_combined_relative$sum_sq_dif, 0.01) %>% as.numeric()

mnfe_combined_q <- mnfe_combined_relative %>%
  dplyr::mutate(in_set = sum_sq_dif < percent5_quantile_relative)

# mnfe_combined_q_tmn <- mnfe_combined_relative %>%
  # dplyr::mutate(in_set = sum_sq_dif < percent5_quantile_relative) %>%
  # filter(pos_change == TRUE)

# mnfe_combined_q_tmn_in <- mnfe_combined_q_tmn %>% dplyr::filter(in_set == TRUE)
# 
# mnfe_combined_relative %>% 
#   ggplot(aes(x = sum_sq_dif)) +
#   geom_histogram() +
#   # scale_x_log10() +
#   
#   xlim(0, 2.5) +
#   ggtitle('Meta')
# 
# # plotting quantile cutoff
fe_uptake_post <- mnfe_combined_q %>%
  filter(in_set == TRUE) %>%
  mutate(avail_space_high = avail_space > 0.075) %>%
  ggplot(aes(cost_par)) +
  geom_histogram() +
  ggtitle('Metaproteome-derived Diatom Proteome') +
  # geom_density(fill = 'grey70', kernel = 'gaussian') +
  # facet_grid(avail_space_high~in_set) +
  theme_bw() +
  xlab('Fe Uptake Cost Parameter') +
  ylab('Count');fe_uptake_post
# 
# # plotting quantile cutoff
epsilon_post <- mnfe_combined_q %>%
  filter(in_set == TRUE) %>%
  mutate(avail_space_high = avail_space > 0.075) %>%
  ggplot(aes(epsilon_a)) +
  geom_histogram() +
  ggtitle('') +
  # geom_density(fill = 'grey70') +
  # facet_grid(avail_space_high~in_set) +
  theme_bw() +
  xlab('Epsilon a') +
  ylab('Count');epsilon_post
# 
# # plotting quantile cutoff
avail_space_post <- mnfe_combined_q %>%
  filter(in_set == TRUE) %>%
  ggplot(aes(avail_space)) +
  geom_histogram() +
  # ggtitle('') +
  # geom_density(fill = 'grey70', kernel = 'gaussian') +
  # facet_grid(~in_set) +
  theme_bw() +
  # xlim(0, 0.15) +
  xlab('Avail Space') +
  ylab('Count');avail_space_post
# 
# ggarrange(fe_uptake_post, 
#           epsilon_post, 
#           avail_space_post)
# 
# ######### based on the combined abc, the ss error of 0.56 is a cutoff for parameter sets
# best_par_comb_rel <- mnfe_combined_q %>% 
#   mutate(combined_par_in_set = sum_sq_dif < percent5_quantile_relative) %>% 
#   select(cost_par, 
#          epsilon_a, 
#          avail_space)
# 
# combined_par_sets_meta_model <- meta_model_out %>%
#   inner_join(best_par_comb, by = c('epsilon_a', 'cost_par', 'avail_space'))
# 
# combined_par_sets_meta_model %>% 
#   ggplot(aes(x = Fex %>% as.factor(), y = u)) +
#   geom_boxplot()
# 
# 
# ## distribution of euc distance
# mnfe_combined_q %>% 
#   ggplot(aes(sum_sq_dif)) +
#   geom_histogram() +
#   # xlim(0, 5) +
#   geom_vline(xintercept = percent5_quantile) +
#   ggtitle('Distribution of Euclidean Distances: Metaproteome-derived diatom proteome') +
#   theme_bw() +
#   scale_x_log10() +
#   xlab('Euclidean Distance')
# 
# # looking at the posterior predictive check
# posterior_predictive_check_meta_relative <- inner_join(mnfe3_meta_trans_2_relative, mnfe_combined_q_tmn,
#                                                        by = c('cost_par', 'avail_space', 'epsilon_a')) %>%
#   filter(in_set == TRUE, coarse_grain != 'Tmn', day == 3) %>%
#   ggplot(aes(x = coarse_grain, y = relative_change))  +
#   geom_boxplot(width = 0.3) +
#   facet_grid(~day) +
#   ggtitle('Metaproteome-derived diatom proteome') +
#   theme_bw() +
#   scale_y_log10() +
#   geom_point(data = coarse_diatoms_specific_norm_relative %>%
#                filter(day == 3) %>%
#                filter(coarse_grains != 'U'),
#              aes(x = coarse_grains, y = relative_change),
#              size = 3,
#              colour = 'blue');posterior_predictive_check_meta_relative
# # 
# 
# 
# 
# 
