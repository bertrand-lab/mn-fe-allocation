## making heatmap of the model growth rate experiments
library(ggpubr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(ggdendro)
library(readr)
library(ggrepel)
library(forcats)
library(scales)

# read in the par name key

name_key <- read_delim('data/par_changed_name_key_with_full.csv', 
           delim=',', 
           escape_double=FALSE, 
           escape_backslash=TRUE)

# read in experimental conditions

# while reading in get the variable modified to be in the dataframe
read_in_data_loop_growth_experiments <- function(data_name_string, 
                                                 dir_to_look = "data/model_output/cost_experiments_july14/", 
                                                 string_to_not_include = '_exp_'){
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
  
  data_names <- data_names[!grepl(pattern = string_to_not_include, x = data_names)]
  # print(data_names)
  master_dataframe <- read.csv(paste0(dir_to_look, data_names[1]))
  # print(master_dataframe)
  # print(get_model_run_replicate(data_names[1]))
  # print(get_ros_sens_par(data_names[1]))
  
  # master_dataframe$replicate <- get_model_run_replicate(data_names[1])
  # master_dataframe$ros_par <- get_ros_sens_par(data_names[1])
  
  master_dataframe$par_changed <- gsub(pattern = '.csv', 
                                       replacement = '', 
                                       x = strsplit(data_names[1], 
                                                    split = "__")[[1]][2], 
                                       fixed = TRUE)

  
  # print(get_ros_sens_par(data_names[1]))
  if(length(data_names) > 1){
    for(i in 2:length(data_names)){
      print(data_names[i])
      temp_df <- read.csv(paste0(dir_to_look, data_names[i]))
      temp_df$par_changed <- gsub(pattern = '.csv', 
                                  replacement = '', 
                                  x = strsplit(data_names[i], 
                                               split = "__")[[1]][2], 
                                  fixed = TRUE)
      # temp_df$ros_par <- get_ros_sens_par(data_names[i])
      
      master_dataframe <- rbind(master_dataframe, temp_df)
    }
  }
  
  return(master_dataframe)
}

growth_experiments_july14_control <- read_in_data_loop_growth_experiments('baseline', 
                                                                          dir_to_look = 'data/model_output/')

# first caluclate the expected growth rate under the different growth conditions
control_mean <- growth_experiments_july14_control %>% 
  group_by(Fex, Mnx) %>% 
  summarize_all(mean) %>% 
  dplyr::select(Mnx, Fex, u_trans) %>% 
  dplyr::rename(u_trans_control = u_trans) # rename so that you can merge with the perturbation growth rates

# read in the experimental conditions
growth_experiments_july14 <- read_in_data_loop_growth_experiments('july14') %>% group_by(Fex, Mnx, par_changed) %>% summarize_all(mean)

# join these expected growth rates to the experimental conditions, calculate the fold change in growth rate
growth_experiments_july14_control_appended <- growth_experiments_july14 %>% 
  inner_join(control_mean, by = c('Mnx', 'Fex')) %>% 
  mutate(u_trans_ratio = u_trans/u_trans_control)

# Making nice labels for the Mn and Fe concentrations
growth_experiments_july14_control_appended$mn_experimental_condition <- ifelse(growth_experiments_july14_control_appended$Mnx > 1000, 'High Mn', no = 'Low Mn')
growth_experiments_july14_control_appended$fe_experimental_condition <- ifelse(growth_experiments_july14_control_appended$Fex > 1000, 'High Fe', no = 'Low Fe')
  
growth_experiments_july14_control_appended$experimental_condition <- paste(growth_experiments_july14_control_appended$fe_experimental_condition, growth_experiments_july14_control_appended$mn_experimental_condition, sep = '\n')

# calculating the mean growth rate ratio (to the control) for plotting purposes
mean_u_trans_ratio_df <- growth_experiments_july14_control_appended %>% 
  group_by(par_changed) %>% 
  summarize(mean_u_trans_ratio = mean(u_trans_ratio))

quant_utransratio_09 <- quantile(mean_u_trans_ratio_df$mean_u_trans_ratio, 0.75) %>% as.numeric()
quant_utransratio_01 <- quantile(mean_u_trans_ratio_df$mean_u_trans_ratio, 0.25) %>% as.numeric()

par_subset <- mean_u_trans_ratio_df[mean_u_trans_ratio_df$mean_u_trans_ratio > quant_utransratio_09 | mean_u_trans_ratio_df$mean_u_trans_ratio < quant_utransratio_01, ]$par_changed

growth_experiments_july14_control_appended_full_names <- growth_experiments_july14_control_appended %>% 
  inner_join(name_key, by = 'par_changed')

growth_experiments_july14_control_appended_full_names$full_name_2 <- as.character(growth_experiments_july14_control_appended_full_names$full_name)


# make heatmap of all fold changes
growth_experiments_heatmap <- growth_experiments_july14_control_appended_full_names %>% 
  filter(par_changed %in% par_subset) %>% 
  ggplot(aes(x = experimental_condition, 
             y = fct_reorder(full_name_2, u_trans_ratio))) +
  geom_tile(aes(fill = u_trans_ratio)) +
  scale_fill_gradient2(midpoint = 1,
                       name = 'Growth Rate\nFold Change',
                       high = muted("red"),
                       mid = "white",
                       low = muted("blue")) +
  theme_bw() +
  xlab('') +
  theme(legend.position = 'bottom') +
  ylab('Model Parameter');growth_experiments_heatmap

growth_experiments_heatmap_full <- growth_experiments_july14_control_appended_full_names %>% 
  filter(full_name != 'old_k_cat',
         full_name != 'Fe per Antioxidant') %>% # these parameters are filtered out because they are not used in the current model
  ggplot(aes(x = experimental_condition, 
             y = fct_reorder(full_name, u_trans_ratio))) +
  geom_tile(aes(fill = u_trans_ratio)) +
  scale_fill_gradient2(midpoint = 1,
                       name = 'Growth Rate\nFold Change') +
  theme_bw() +
  xlab('') +
  theme(legend.position = 'bottom') +
  ylab('Model Parameter');growth_experiments_heatmap_full


ggsave(growth_experiments_heatmap_full, 
       filename = 'figures/full_heatmap_growth_experiments.png',
       width = 10.2, height = 9.61)


# plotting the interaction index

mean_mid_growth <- growth_experiments_july14_control_appended_full_names %>% 
  filter(experimental_condition == 'High Fe\nLow Mn' |
           experimental_condition == 'Low Fe\nHigh Mn') %>% 
  group_by(par_changed) %>% 
  summarize(mean_growth_rate = mean(u_trans),
            min_growth_rate = min(u_trans))

mean_low_growth <- growth_experiments_july14_control_appended_full_names %>% 
  filter(experimental_condition == 'Low Fe\nLow Mn')

interaction_df <- mean_low_growth %>% 
  inner_join(mean_mid_growth, by = 'par_changed') %>% 
  mutate(int_index_mean = mean_growth_rate - u_trans,
         int_index_min = min_growth_rate - u_trans)

interaction_df_names <- interaction_df[interaction_df$int_index_min > quantile(interaction_df$int_index_min, 0.8) %>% as.numeric(), ]$full_name

# plotting the parameters with the highest interaction index
highest_interaction_parameters <- interaction_df %>%
  dplyr::filter(full_name %in% interaction_df_names) %>%
  ggplot(aes(x = fct_reorder(full_name, int_index_min, .desc = TRUE), 
             y = int_index_min)) +
  geom_col(fill = 'darkblue', alpha = 0.7) +
  theme_bw() +
  coord_flip() +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1)) +
  ylab('Interaction Index') +
  # coord_flip() +
  xlab('Model Parameter');highest_interaction_parameters

ggsave(highest_interaction_parameters, 
       filename = 'figures/interaction_index_proteomic_allocation.png',
       width = 7.95, height = 6.57)


interaction_long_view <- ggarrange(empty_plot, 
          highest_interaction_parameters,
          nrow = 2, heights = c(0.5, 2),
          labels = c('A', ''), font.label = list(size = 9))


all_phenom_out_plots <- ggarrange(phi05_high_k_p, phi05_low_k_p,
  phi05_p, phi09_p, nrow = 2, ncol = 2,
  labels = c('C', 'D', 'E', 'F'),
  common.legend = TRUE, legend = 'top', font.label = list(size = 9))

top_phenom_out_plots <- ggarrange(phi05_high_k_p, phi05_low_k_p, 
                                  nrow = 1, ncol = 2,
  labels = c('B', 'C'),
  common.legend = TRUE, legend = 'none', font.label = list(size = 9))

bottom_phenom_out_plots <- ggarrange(phi05_p, 
                                     phi09_p, 
                                  nrow = 1, ncol = 2,
  labels = c('E', 'F'),
  legend = 'bottom',
  common.legend = TRUE, font.label = list(size = 9))

all_phenom_out_plots_blank_mid <- ggarrange(empty_plot, 
                                            top_phenom_out_plots,
                                            empty_plot,
                                            bottom_phenom_out_plots,
                                            heights = c(1, 2, 1.3, 2.7),
                                            labels = c('A', 'B', 'D'),
                                            ncol = 1,
                                            common.legend = TRUE, font.label = list(size = 9))


ggsave(all_phenom_out_plots_blank_mid, 
       filename = 'figures/phenom_plots_only.png',
       width = 8.42*(4/5)*0.63, height = 8.45*0.75, dpi = 1000)

