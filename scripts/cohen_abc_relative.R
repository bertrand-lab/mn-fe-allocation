# cohen abc relative

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

# cohen_par_model1 <- read_in_data_loop_ls(data_name_string = 'hen_par_gen__april9', 
                                        # dir_to_look = "../mn-fe-allocation/data/abc_out/")
# cohen_par_model2 <- read_in_data_loop_ls(data_name_string = 'meta_par_gen__april13', 
                                        # dir_to_look = "../mn-fe-allocation/data/abc_out/")
# cohen_par_model <- rbind(cohen_par_model1, cohen_par_model2)
cohen_par_model1 <- read_in_data_loop_ls(data_name_string = 'cohen_par_gen__april28',
dir_to_look = "../mn-fe-allocation/data/abc_out/")
cohen_par_model2 <- read_in_data_loop_ls(data_name_string = 'cohen_par_gen__may',
dir_to_look = "../mn-fe-allocation/data/abc_out/")
cohen_par_model <- rbind(cohen_par_model1, cohen_par_model2)
# 
# cohen_par_model <- read_in_data_loop_ls(data_name_string = 'cohen_par_gen__july', 
#                                         dir_to_look = "../mn-fe-allocation/data/abc_from_cc/")
# must recalculate the output for the model
# aa free pool is not included because the protein extraction procedure would not pick on this for proteomics MS
cohen_model_par_sets_form <- cohen_par_model %>%
  dplyr::select(A_aa, P_aa, Tfe_aa, Tn_aa, R_aa, Tmn_aa, Fex, 
                cost_par, avail_space,
                epsilon_a) %>% 
  ungroup() %>% 
  melt(id.vars = c("Fex", "cost_par", "avail_space", "epsilon_a"), 
       variable.name = c("coarse_grain"),
       value.name = c("aa_val")) %>% 
  group_by(cost_par, avail_space, epsilon_a, Fex) %>% 
  summarize(sum_aa = sum(aa_val))

cohen_model_par_sets_form2 <- cohen_par_model %>%
  dplyr::select(A_aa, P_aa, Tfe_aa, Tn_aa, R_aa, Tmn_aa, Fex, cost_par, avail_space,
                epsilon_a) %>% 
  ungroup() %>% 
  melt(id.vars = c("Fex", "cost_par", "avail_space", "epsilon_a"), 
       variable.name = c("coarse_grain"),
       value.name = c("aa_val")) %>% 
  inner_join(cohen_model_par_sets_form, 
             by = c("Fex", "cost_par", "avail_space", "epsilon_a")) %>% 
  mutate(aa_val_frac = aa_val/sum_aa)

# renaming factor levels
levels(cohen_model_par_sets_form2$coarse_grain) <- c('A', 'P', 'Tfe', 'Tn', 'R', 'Tmn')

### transform the data to be relative increase or decrease
cohen_high_fe_vals <- cohen_model_par_sets_form2 %>%
  dplyr::filter(Fex > 2000) %>% 
  dplyr::select(cost_par, avail_space, epsilon_a, coarse_grain, aa_val_frac) %>% 
  dplyr::rename(norm_factor = aa_val_frac)

cohen_model_par_sets_form2_relative <- cohen_model_par_sets_form2 %>% 
  inner_join(cohen_high_fe_vals, by = c('cost_par', 'avail_space', 'epsilon_a', 'coarse_grain')) %>%
  # dplyr::filter(cost_par > 0.2319033, cost_par < 0.2319035) # checking that the joinging worked as anticipated
  mutate(model_relative_change = aa_val_frac/norm_factor)


# taking the mean across biological replicates from cohen 2018
cohen_prot_fe <- read.csv('data/culture_proteomes/Cohen2018_protein_combined_tables4_s5_manual.csv')

# cohen_prot_fe[cohen_prot_fe$coarse_grain == 'Tfe' & cohen_prot_fe$kegg_id == 'K02016', ]$coarse_grain <- 'U'


# aggregate the Fe lim treatments 
cohen_prot_fe_sum <- cohen_prot_fe %>% 
  dplyr::mutate(felim = rowMeans(dplyr::select(., dplyr::starts_with("Low.Fe.Repl"))),
         fereplete = rowMeans(dplyr::select(., dplyr::starts_with("High.Fe.Repl")))) %>% 
  group_by(coarse_grain) %>% 
  summarize(sum_pep_nsaf_felim = sum(felim),
            sum_pep_nsaf_fereplte = sum(fereplete)) %>% 
  mutate(relative_change = sum_pep_nsaf_felim/sum_pep_nsaf_fereplte)

## format for the joining below
cohen_prot_fe_sum$Fex <- rep(20, nrow(cohen_prot_fe_sum))
cohen_prot_fe_sum_high_fe <- cohen_prot_fe_sum

cohen_prot_fe_sum_high_fe$Fex <- rep(2700, nrow(cohen_prot_fe_sum))
cohen_prot_fe_sum_high_fe$relative_change <- rep(1, nrow(cohen_prot_fe_sum))

cohen_prot_obs_form <- rbind(cohen_prot_fe_sum,
                             cohen_prot_fe_sum_high_fe)

# cohen_prot_obs_form[cohen_prot_obs_form$coarse_grain == 'A' & cohen_prot_obs_form$Fex == 20, ]$relative_change <- 1.2

# cohen_prot_obs_form <- cohen_prot_obs_form %>% filter(coarse_grain != 'Tfe')

joined_cohen_par_out_relative <- inner_join(cohen_model_par_sets_form2_relative, 
                                            cohen_prot_obs_form, 
                                           by = c("Fex", "coarse_grain")) %>% 
  # filter(coarse_grain != 'Tfe') %>%
  mutate(squared_percentage_dif = (model_relative_change - relative_change)^2) %>% # mutate a new category to look at the squared difference between percentages
  group_by(cost_par, avail_space, epsilon_a) %>% # group by cost_par and space_par
  summarize(sum_sq_dif = sqrt(sum(squared_percentage_dif)))


# writing abc output summary
write.csv(joined_cohen_par_out_relative,
          file = '../mn-fe-allocation/data/abc_intermediate/cohen_abc_par_sets.csv')


# percent5_quantile_relative_cohen <- quantile(joined_cohen_par_out_relative$sum_sq_dif, 0.05) %>% as.numeric()

# gettin gquantiles
# joined_cohen_par_out_q_relative <- joined_cohen_par_out_relative %>%
  # dplyr::mutate(in_set = sum_sq_dif < percent5_quantile_relative_cohen)

# joined_cohen_par_out_relative %>% 
  # ggplot(aes(x = sum_sq_dif)) +
  # geom_histogram(binwidth = 10) +
  # xlim(0, 50) + 
  # scale_x_log10() +
  # ggtitle('Cohen')

# plotting quantile cutoff
# cohen_post_cost <- joined_cohen_par_out_q_relative %>% 
#   filter(in_set == TRUE) %>%
#   ggplot(aes(cost_par)) +
#   geom_histogram() +
#   # geom_density(fill = 'grey70') +
#   # facet_grid(~in_set) +
#   theme_bw() +
#   xlab('Fe Uptake Cost Parameter') +
#   ylab('Count') +
#   ggtitle('Cohen et al (2018)')#;cohen_post_cost

# plotting quantile cutoff
# cohen_post_eps <- joined_cohen_par_out_q_relative %>% 
#   filter(in_set == TRUE) %>%
#   ggplot(aes(epsilon_a)) +
#   ggtitle('') +
#   geom_histogram() +
#   # geom_density(fill = 'grey70') +
#   # facet_grid(~in_set) +
#   theme_bw() +
#   xlab('Epsilon a') +
#   ylab('Count')#;cohen_post_eps
# 
# # plotting quantile cutoff
# cohen_post_avail <- joined_cohen_par_out_q_relative %>% 
#   filter(in_set == TRUE) %>%
#   ggplot(aes(avail_space)) +
#   ggtitle('') +
#   geom_histogram() +
#   # geom_density(fill = 'grey70') +
#   # facet_grid(~in_set) +
#   theme_bw() +
#   xlab('Avail Space') +
#   ylab('Count')#;cohen_post_avail

# ggarrange(cohen_post_cost, cohen_post_eps, cohen_post_avail)
# 
# joined_cohen_par_out_q_relative %>% 
#   ggplot(aes(sum_sq_dif)) +
#   geom_histogram()+
#   # xlim(0, 5) +
#   ggtitle('Distribution of Euclidean Distances: Cohen et al 2018 diatom proteome') +
#   theme_bw()

# predictive posterior check ----------------------------------------------
# 
# posterior_predictive_check_cohen <- inner_join(cohen_model_par_sets_form2_relative, cohen_prot_obs_form, 
#                                               by = c("Fex", "coarse_grain")) %>%
#   inner_join(joined_cohen_par_out_q_relative, 
#              by = c('cost_par', 'avail_space', 'epsilon_a')) %>% 
#   filter(in_set == TRUE,
#          coarse_grain != 'U',
#          Fex < 2000) %>%
#   ggplot(aes(x = coarse_grain, 
#              y = model_relative_change)) +
#   ggtitle('Cohen et al (2018)') + 
#   geom_boxplot(width = 0.3) +
#   facet_grid(~Fex) +
#   # geom_violin() +
#   geom_point(data = cohen_prot_obs_form %>% filter(coarse_grain != 'U', Fex < 2000), 
#              aes(x = coarse_grain, y = relative_change),
#              size = 3, colour = 'blue') +
#   theme_bw() +
#   scale_y_log10();posterior_predictive_check_cohen
