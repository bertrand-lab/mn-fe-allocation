## making heatmap of the model growth rate experiments

library(ggplot2)
library(magrittr)
library(dplyr)
library(ggdendro)
library(readr)
library(ggrepel)
library(forcats)

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
  dplyr::rename(u_trans_control = u_trans)

# read in the experimental conditions
growth_experiments_july14 <- read_in_data_loop_growth_experiments('july14') %>% group_by(Fex, Mnx, par_changed) %>% summarize_all(mean)

# join these expected growth rates to the experimental conditions, calculate the fold change in growth rate
growth_experiments_july14_control_appended <- growth_experiments_july14 %>% 
  inner_join(control_mean, by = c('Mnx', 'Fex')) %>% 
  mutate(u_trans_ratio = u_trans/u_trans_control)

growth_experiments_july14_control_appended$mn_experimental_condition <- ifelse(growth_experiments_july14_control_appended$Mnx > 1000, 'High Mn', no = 'Low Mn')
growth_experiments_july14_control_appended$fe_experimental_condition <- ifelse(growth_experiments_july14_control_appended$Fex > 1000, 'High Fe', no = 'Low Fe')
  
growth_experiments_july14_control_appended$experimental_condition <- paste(growth_experiments_july14_control_appended$fe_experimental_condition, growth_experiments_july14_control_appended$mn_experimental_condition, sep = '\n')

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
                       name = 'Growth Rate\nFold Change') +
  theme_bw() +
  xlab('') +
  theme(legend.position = 'bottom') +
  ylab('Model Parameter');growth_experiments_heatmap

growth_experiments_heatmap_full <- growth_experiments_july14_control_appended_full_names %>% 
  filter(full_name != 'old_k_cat',
         full_name != 'Fe per Antioxidant') %>% 
  ggplot(aes(x = experimental_condition, 
             y = fct_reorder(full_name, u_trans_ratio))) +
  geom_tile(aes(fill = u_trans_ratio)) +
  scale_fill_gradient2(midpoint = 1,
                       name = 'Growth Rate\nFold Change') +
  theme_bw() +
  xlab('') +
  theme(legend.position = 'bottom') +
  ylab('Model Parameter');growth_experiments_heatmap_full
# growth_experiments_july14_control_appended %>% 
#   ggplot(aes(x = Fex, u_trans_ratio)) +
#   facet_wrap(~experimental_condition, nrow = 4)
# 
# 
# growth_experiments_july14_control_appended$congruence <- ifelgrowth_experiments_july14_control_appended
# 



external_mn_ribo <- model_out_mean %>%
  # filter(Mnx < 1001, Fex > 1) %>% 
  filter(Mnx < 1001) %>% 
  ggplot(aes(x = Mnx, y = R, fill = Fex)) +
  geom_point(pch = 21, size = 3, alpha = 0.7, stroke = 1) +
  theme_bw() +
  # geom_line(aes(group = Fex)) +
  xlab('dMn (pM)') +
  ylab('Ribosomes (per cell)') +
  scale_fill_distiller(name = 'dFe (pM)', palette = 'RdYlBu') +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8, 0.7));external_mn_ribo

external_mn_a <- model_out_mean %>%
  filter(Mnx < 1001) %>% 
  # filter(Mnx < 1001, Fex > 1) %>% 
  ggplot(aes(x = Mnx, y = A/2, fill = Fex)) +
  geom_point(pch = 21, size = 3, alpha = 0.7, stroke = 1) +
  theme_bw() +
  xlab('dMn (pM)') +
  ylab('Antioxidants (MnSOD per cell)') +
  theme(legend.position = c(0.8, 0.3),
         panel.grid = element_blank()) +
  scale_fill_distiller(name = 'dFe (pM)', palette = 'RdYlBu');external_mn_a


top_mn_anti <- ggarrange(external_mn_a, external_mn_ribo, 
                         nrow = 1, align = 'hv',
                         labels = c('a', 'b'))

# heatmap_with_ribo_anti <- ggarrange(top_mn_anti, 
#           growth_experiments_heatmap,
#           nrow = 2,
#           labels = c('', 'c'))

# growth_experiments_ribo_anti <- ggarrange(growth_experiments_heatmap, 
          # ggarrange(external_mn_a, external_mn_ribo, 
                    # nrow = 2, align = 'hv', 
                    # labels = c('b', 'c')), 
          # nrow = 1, widths = c(2, 1), labels = c('a', ''));growth_experiments_ribo_anti


# ggsave(heatmap_with_ribo_anti, 
#        filename = 'figures/heatmap_with_ribo_anti.png',
#        width = 10.2, height = 9.61)

ggsave(growth_experiments_heatmap_full, 
       filename = 'figures/full_heatmap_growth_experiments.png',
       width = 10.2, height = 9.61)

# growth_experiments_july14_control_appended_full_names %>% 
#   filter(full_name != 'old_k_cat',
#          full_name != 'Fe per Antioxidant') %>% 
#   ggplot(aes(x = experimental_condition, 
#              y = u_trans)) +
#   geom_col() +
#   theme_bw() +
#   geom_point(aes(x = experimental_condition, 
#                  y = u_trans_control)) +
#   xlab('') +
#   facet_wrap(~full_name) +
#   theme(legend.position = 'right') +
#   ylab('Model Parameter Modified')


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

highest_interaction_parameters <- interaction_df %>%
  dplyr::filter(full_name %in% interaction_df_names) %>% 
  ggplot(aes(x = fct_reorder(full_name, int_index_min, .desc = TRUE), 
             y = int_index_min)) +
  geom_col(fill = 'darkblue', alpha = 0.7) +
  theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust = 1)) +
  ylab('Interaction Index') +
  # coord_flip() +
  xlab('Model Parameter');highest_interaction_parameters


interaction_long_view <- ggarrange(empty_plot, 
          highest_interaction_parameters,
          nrow = 2, heights = c(0.5, 2),
          labels = c('a', ''))


all_phenom_out_plots <- ggarrange(phi05_high_k_p, phi05_low_k_p,
  phi05_p, phi09_p, nrow = 2, ncol = 2,
  labels = c('c', 'd', 'e', 'f'),
  common.legend = TRUE, legend = 'top')

top_phenom_out_plots <- ggarrange(phi05_high_k_p, phi05_low_k_p, 
                                  nrow = 1, ncol = 2,
  labels = c('c', 'd'),
  common.legend = TRUE, legend = 'none')

bottom_phenom_out_plots <- ggarrange(phi05_p, 
                                     phi09_p, 
                                  nrow = 1, ncol = 2,
  labels = c('f', 'g'),
  legend = 'bottom',
  common.legend = TRUE)

all_phenom_out_plots_blank_mid <- ggarrange(empty_plot, 
                                            top_phenom_out_plots,
                                            empty_plot,
                                            bottom_phenom_out_plots,
                                            heights = c(1.75, 2, 1, 2.5),
                                            labels = c('b', 'c', 'e'),
                                            ncol = 1,
                                            common.legend = TRUE)

all_phenom_out_plots_blank <- ggarrange(empty_plot,
                                        all_phenom_out_plots, 
                                        heights = c(1, 4),
                                        labels = c('b', ''),
                                        nrow = 2,
                                        ncol = 1)

interaction_index_no_heatmap <- ggarrange(interaction_long_view,
          all_phenom_out_plots_blank, ncol = 2, 
          widths = c(1, 3))

interaction_index_no_heatmap_mid <- ggarrange(interaction_long_view,
                                          all_phenom_out_plots_blank_mid, 
          ncol = 2, 
          widths = c(1.2, 3))

ggsave(interaction_index_no_heatmap_mid, 
       filename = 'figures/interaction_index_no_heatmap_mid.png',
       width = 9.61, height = 7.57)

