# nunn_abc_relative

## Nunn data ABC analysis
library(tibble)

read_in_data_loop_ls <- function(data_name_string, dir_to_look = "data/abc_out/"){
  # '''
  # this function reads in the data from the model output in a loop
  # it also parses the file name and adds the replicate number as
  # well as the ROS_sensitivity parameter which is included in the
  # file name
  # '''
  
  model_output_list <- dir(dir_to_look)
  # print(model_output_list)
  data_names <- model_output_list[grepl(pattern = data_name_string,
                                        x = model_output_list)]
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

par_sets_model1 <- read_in_data_loop_ls(data_name_string = 'nunn_par_gen__april28',
                                        dir_to_look = "../mn-fe-allocation/data/abc_out/")
par_sets_model2 <- read_in_data_loop_ls(data_name_string = 'nunn_par_gen__may',
                                        dir_to_look = "../mn-fe-allocation/data/abc_out/")
par_sets_model <- rbind(par_sets_model1, par_sets_model2)

# formatting Nunn model output
# must recalculate the output for the model
nunn_model_par_sets_form <- par_sets_model %>%
  dplyr::select(A_aa, P_aa, Tfe_aa, Tn_aa, R_aa, Tmn_aa, Fex, 
                cost_par, avail_space,
                epsilon_a) %>% 
  ungroup() %>% 
  melt(id.vars = c("Fex", "cost_par", "avail_space", "epsilon_a"), 
       variable.name = c("coarse_grain"),
       value.name = c("aa_val")) %>% 
  group_by(cost_par, avail_space, epsilon_a, Fex) %>% 
  summarize(sum_aa = sum(aa_val))


# Tfe_aa_no_dyn cost is selected here because it was only a ferric reducatse that was identified in this dataset
# FTR1
nunn_model_par_sets_form2 <- par_sets_model %>%
  dplyr::select(A_aa, P_aa, Tfe_aa_no_dyn, Tn_aa, R_aa, Tmn_aa_no_dyn, Fex, cost_par, avail_space,
                epsilon_a) %>% 
  ungroup() %>% 
  melt(id.vars = c("Fex", "cost_par", "avail_space", "epsilon_a"), 
       variable.name = c("coarse_grain"),
       value.name = c("aa_val")) %>% 
  inner_join(nunn_model_par_sets_form, 
             by = c("Fex", "cost_par", "avail_space", "epsilon_a")) %>% 
  mutate(aa_val_frac = aa_val/sum_aa)

# renaming factor levels
levels(nunn_model_par_sets_form2$coarse_grain) <- c('A', 'P', 'Tfe', 'Tn', 'R', 'Tmn')

### transform the data to be relative increase or decrease
high_fe_vals <- nunn_model_par_sets_form2 %>%
  dplyr::filter(Fex > 2000) %>% 
  dplyr::select(cost_par, avail_space, epsilon_a, coarse_grain, aa_val_frac) %>% 
  dplyr::rename(norm_factor = aa_val_frac)

nunn_model_par_sets_form2_relative <- nunn_model_par_sets_form2 %>% 
  inner_join(high_fe_vals, by = c('cost_par', 'avail_space', 'epsilon_a', 'coarse_grain')) %>%
  # dplyr::filter(cost_par > 0.2319033, cost_par < 0.2319035) # checking that the joinging worked as anticipated
  mutate(model_relative_change = aa_val_frac/norm_factor)

# taking the mean across biological replicates from nunn et al 2013
nunn_output_mean <- nunn_output %>% 
  group_by(coarse_grain, Fex) %>% 
  summarize(mean_percentage = mean(percentage))

nunn_output_mean_high_fe <- nunn_output_mean %>% 
  filter(Fex > 2000) %>% 
  dplyr::rename(norm_factor = mean_percentage) %>% 
  dplyr::select(norm_factor, coarse_grain)

nunn_output_mean_relative <- nunn_output_mean %>%
  inner_join(nunn_output_mean_high_fe, by = c('coarse_grain')) %>% 
  mutate(relative_change = mean_percentage/norm_factor)

joined_nunn_par_out_relative <- inner_join(nunn_model_par_sets_form2_relative, 
                                           nunn_output_mean_relative, 
                                           by = c("Fex", "coarse_grain")) %>% 
  # filter(coarse_grain != 'Tfe') %>%
  mutate(squared_percentage_dif = (model_relative_change - relative_change)^2) %>% # mutate a new category to look at the squared difference between percentages
  group_by(cost_par, avail_space, epsilon_a) %>% # group by cost_par and space_par
  summarize(sum_sq_dif = sqrt(sum(squared_percentage_dif)))


nunn_relative_output_formatted <- dcast(nunn_model_par_sets_form2_relative %>% 
                                     dplyr::filter(Fex == 87), 
                                   cost_par + epsilon_a + avail_space ~ coarse_grain, 
                                   value.var = "model_relative_change")

# get corresponding model output for a given par set
nunn_growth_rate_formatted <- par_sets_model %>% dplyr::select(cost_par, epsilon_a, 
                                                          avail_space, fe_uptake, 
                                                          mn_uptake, u_trans,
                                                          total_fe_amol, 
                                                          total_mn_amol, Fex) %>% 
  melt(id.vars = c("cost_par", "epsilon_a", "avail_space", "Fex"),
       meaure.vars = c("fe_uptake", "mn_uptake", "u_trans", 
                       "total_fe_amol", "total_mn_amol")) %>% 
  dcast(cost_par + epsilon_a + avail_space ~ variable + Fex)


joined_nunn_par_out_relative_mod <- joined_nunn_par_out_relative %>% 
  inner_join(nunn_relative_output_formatted, 
             by = c('cost_par', 'avail_space', 'epsilon_a')) %>% 
  inner_join(nunn_growth_rate_formatted,
             by = c('cost_par', 'avail_space', 'epsilon_a'))

# writing abc output summary
write.csv(joined_nunn_par_out_relative_mod,
          file = '../mn-fe-allocation/data/abc_intermediate/nunn_abc_par_sets_relative.csv', 
          row.names = FALSE)
