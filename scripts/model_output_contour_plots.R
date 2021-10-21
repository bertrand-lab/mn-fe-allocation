## drafting publication ready figures

# setwd("~/../../Google Drive/Projects/mn-fe-allocation/")

# plotting model results from version MnFe14

library(ggplot2)
library(dplyr)
library(magrittr)
library(akima)
library(spatialkernel)
library(fields)
library(gridExtra)
library(readr)
library(ggpubr)
library(stringi)

# functions for plotting --------------------------------------------------

get_model_run_replicate <- function(model_name_i){
  model_run_rep_number <- strsplit(model_name_i, split = "_")[[1]][4]
  
  if(model_run_rep_number == "google"){
    model_run_rep_number <- as.numeric(strsplit(model_name_i, split = "_")[[1]][5])
    model_run_rep_number <- model_run_rep_number*100
  }
  
  return(model_run_rep_number)
}

get_ros_sens_par <- function(model_name_i){
  model_ros_par <- strsplit(model_name_i, split = "_")[[1]][9]
  if(class(model_ros_par) == 'numeric'){
    extracted_number <- parse_number(model_ros_par)
  } else {
    extracted_number <- 'no_sens_par'
  }
  return(extracted_number)
}

# function to read in data
read_in_data_loop <- function(data_name_string, dir_to_look = "data/model_output/", string_to_not_include = '_exp_'){
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
  
  master_dataframe$replicate <- get_model_run_replicate(data_names[1])
  master_dataframe$ros_par <- get_ros_sens_par(data_names[1])
  
  # print(get_ros_sens_par(data_names[1]))
  if(length(data_names) > 1){
    for(i in 2:length(data_names)){
      print(data_names[i])
      temp_df <- read.csv(paste0(dir_to_look, data_names[i]))
      temp_df$replicate <- get_model_run_replicate(data_names[i])
      temp_df$ros_par <- get_ros_sens_par(data_names[i])
      
      master_dataframe <- rbind(master_dataframe, temp_df)
    }
  }
  
  return(master_dataframe)
}

interpolate_model <- function(x_val, y_val, z_val, x_name, y_name, z_name, df){
  
  
  # unique light levels
  unique_light_levels <- df$I %>% unique()
  df_master <- data.frame(x_name = numeric(), y_name = numeric(),
                          z_name = numeric(), 
                          light_level = numeric())
  
  for(light_level in 1:length(unique_light_levels)){
    x_val_filt <- x_val[df$I == unique_light_levels[light_level]]
    y_val_filt <- y_val[df$I == unique_light_levels[light_level]]
    z_val_filt <- z_val[df$I == unique_light_levels[light_level]]
    # print(mean(z_val_filt))
    
    # Create grid
    grd <- expand.grid(x = seq(min(x_val_filt),
                               max(x_val_filt),
                               length=80), 
                       y = seq(min(y_val_filt), 
                               max(y_val_filt), 
                               length=80))
    
    # Find points in convex hull
    mask <- pinpoly(cbind(x_val_filt, 
                          y_val_filt)[chull(x_val_filt, 
                                            y_val_filt),], 
                    as.matrix(grd))
    
    # Crop grid to convex hull
    grd <- grd[mask == 2,]
    
    # Interpolate to points in convex hull
    res <- interpp(x_val_filt, 
                   y_val_filt, 
                   z = z_val_filt, 
                   xo = grd$x, yo = grd$y)
    
    new_df <- data.frame(name1 = res$x, name2 = res$y, 
                         name3 = res$z, 
                         light_level = rep(unique_light_levels[light_level], 
                                           length(res$z)))
    
    names(new_df) <- c(x_name, y_name, z_name, "light_level")
    # add on this to the master
    df_master <- rbind(df_master, new_df)
    
    
  }  
  
  return(df_master)
  
}

# model output interpolation
quota_growth_interpolate <- function(dataset){
  
  # dataset <- dataset %>% dplyr::filter(I == 50 | I == 100 | I == 500)
  
  mni2 <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                            x_val = dataset$Fex, x_name = 'Fex',
                            z_val = dataset$total_mn_amol, z_name = 'total_mn_amol', df = dataset)
  fei2 <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                            x_val = dataset$Fex, x_name = 'Fex',
                            z_val = dataset$total_fe_amol, z_name = 'total_fe_amol', df = dataset)
  u_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                               x_val = dataset$Fex, x_name = 'Fex',
                               z_val = dataset$u_trans, z_name = 'u_trans', df = dataset)
  
  # alpha_inter_mn <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
  #                              x_val = dataset$Fex, x_name = 'Fex',
  #                              z_val = 1e18*dataset$mn_uptake_aff/dataset$total_mn_amol, 
  #                              z_name = 'mn_uptake_aff', df = dataset)
  
  # alpha_inter_fe <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
  #                                  x_val = dataset$Fex, x_name = 'Fex',
  #                                  z_val = 1e18*dataset$fe_uptake_aff/dataset$total_fe_amol, 
  #                                  z_name = 'fe_uptake_aff', df = dataset)
  
  print('finished interpolation')
  
  a_frac_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                                    x_val = dataset$Fex, x_name = 'Fex',
                                    z_val = dataset$A_frac, 
                                    z_name = 'A_frac', df = dataset)
  p_frac_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                                    x_val = dataset$Fex, x_name = 'Fex',
                                    z_val = dataset$P_frac, 
                                    z_name = 'P_frac', df = dataset)
  r_frac_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                                    x_val = dataset$Fex, x_name = 'Fex',
                                    z_val = dataset$R_frac, 
                                    z_name = 'R_frac', df = dataset)
  
  tmn_frac_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                                      x_val = dataset$Fex, x_name = 'Fex',
                                      z_val = dataset$Tmn_no_dyn_frac, 
                                      z_name = 'Tmn_frac', df = dataset)
  tfe_frac_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                                      x_val = dataset$Fex, x_name = 'Fex',
                                      z_val = dataset$Tfe_no_dyn_frac, 
                                      z_name = 'Tfe_frac', df = dataset)
  tn_frac_inter <- interpolate_model(y_val = dataset$Mnx, y_name = 'Mnx',
                                     x_val = dataset$Fex, x_name = 'Fex',
                                     z_val = dataset$Tn_frac, 
                                     z_name = 'Tn_frac', df = dataset)
  
  
  return(list(mni2, fei2, u_inter, #alpha_inter_mn, alpha_inter_fe,
              a_frac_inter, p_frac_inter, r_frac_inter, tmn_frac_inter, 
              tfe_frac_inter, tn_frac_inter))
  
}

quota_growth_plot <- function(interpolate_out, colour_scheme){
  
  # making panel plots for each
  mn_quota_plot <- interpolate_out[[1]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = total_mn_amol)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Cell Mn [aM]") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = total_mn_amol),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  fe_quota_plot <- interpolate_out[[2]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = total_fe_amol)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Cell Fe [aM]") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = total_fe_amol),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10)) #+
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  u_plot <- interpolate_out[[3]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = u_trans)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Growth Rate (per day)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = u_trans),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10)) #+
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  # mn_aff <- interpolate_out[[4]] %>% 
  # ggplot(aes(x = Fex, y = Mnx, fill = mn_uptake_aff)) +
  # theme_bw() +
  # geom_raster(interpolate = TRUE) +
  # # ggtitle(bquote('Mn Uptake Affinity '~(m^3 s^-1))) +
  # ggtitle(expression(paste("Mn Uptake Affinity (", m^3, s^-1, ")"))) +
  # facet_grid(~light_level) +
  # geom_contour(aes(z = mn_uptake_aff),
  #              colour = 'black', lwd = 0.4,
  #              alpha = 0.7) +
  # scale_fill_distiller(palette = colour_scheme, name = '') +
  # labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
  # theme(axis.title = element_text(size = 12),
  #       axis.text = element_text(size = 14),
  #       legend.title = element_text(size = 12)) #+
  # # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  # 
  # fe_aff <- interpolate_out[[5]] %>% 
  # ggplot(aes(x = Fex, y = Mnx, fill = fe_uptake_aff)) +
  # theme_bw() +
  # geom_raster(interpolate = TRUE) +
  # # ggtitle(bquote('Fe Uptake Affinity '~(m^3 s^-1))) +
  # ggtitle(expression(paste("Fe Uptake Affinity (", m^3, s^-1, ")"))) +
  # facet_grid(~light_level) +
  # geom_contour(aes(z = fe_uptake_aff),
  #              colour = 'black', lwd = 0.4,
  #              alpha = 0.7) +
  # scale_fill_distiller(palette = colour_scheme, name = '') +
  # labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
  # theme(axis.title = element_text(size = 12),
  #       axis.text = element_text(size = 14),
  #       legend.title = element_text(size = 12))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  a_frac <- interpolate_out[[4]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = A_frac*100)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Antioxidant (%)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = A_frac),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  p_frac <- interpolate_out[[5]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = P_frac*100)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Photosystem Units (%)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = P_frac),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  r_frac <- interpolate_out[[6]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = R_frac*100)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Ribosomes (%)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = R_frac),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  tmn_frac <- interpolate_out[[7]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = Tmn_frac*100)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Mn Transporters (%)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = Tmn_frac),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  tfe_frac <- interpolate_out[[8]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = Tfe_frac*100)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("Fe Transporters (%)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = Tfe_frac),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  tn_frac <- interpolate_out[[9]] %>% 
    ggplot(aes(x = Fex, y = Mnx, fill = Tn_frac*100)) +
    theme_bw() +
    geom_raster(interpolate = TRUE) +
    ggtitle("AA Biosynthesis (%)") +
    # facet_grid(~light_level) +
    geom_contour(aes(z = Tn_frac),
                 colour = 'black', lwd = 0.4,
                 alpha = 0.7) +
    scale_fill_distiller(palette = colour_scheme, name = '') +
    labs(y = '[Mn] (pM)', x = '[Fe] (pM)') +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10))# +
  # coord_cartesian(ylim = c(0, 2000), expand = FALSE)
  
  return(list(mn_quota_plot, fe_quota_plot, u_plot, 
              # mn_aff, fe_aff, 
              a_frac, p_frac, r_frac, tmn_frac, tfe_frac, tn_frac))
}

# making model plots ------------------------------------------------------

model_out_mnfe4_model1 <- read_in_data_loop(data_name_string = 'mnfe4_model4_july14', 
                                            dir_to_look = 'data/model_output/', 
                                            string_to_not_include = 'baseline')

ros10_out_mean <- model_out_mnfe4_model1 %>% 
  group_by(Fex, Mnx) %>% 
  summarize_all(mean)

model_out_mnfe4_model1_inter10_mean <- quota_growth_interpolate(ros10_out_mean)
plots10 <- quota_growth_plot(interpolate_out = model_out_mnfe4_model1_inter10_mean, 
                             colour_scheme = 'Spectral')

margin_change <- c(0.1, 0.1, 0.1, 0.1)
margin_shift <- 2
margin_box <- c(-1, -1, -1, -1)

plot_a <- plots10[[4]] + ggtitle('Antioxidants (%)') + xlab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_b <- plots10[[5]] + ggtitle('Photosystem Units (%)') + xlab('') + ylab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_c <- plots10[[6]] + ggtitle('Ribosomes (%)') + xlab('') + ylab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_d <- plots10[[7]] + ggtitle('Mn Transporters (%)') + xlab ('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_e <- plots10[[8]] + ggtitle('Fe Transporters (%)') + xlab('') + ylab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_f <- plots10[[9]] + ggtitle('AA Biosynthesis (%)') + xlab('') + ylab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_g <- plots10[[3]] + ggtitle('Growth Rate (per day)') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_h <- plots10[[2]] + ggtitle('Cell Fe (aMol)') + ylab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))
plot_i <- plots10[[1]] + ggtitle('Cell Mn (aMol)') + ylab('') + 
  theme(plot.margin = unit(margin_change, "lines"),
        legend.margin = margin(margin_change + margin_shift),
        legend.box.margin = margin(margin_box))

full_sweep_out <- ggarrange(plot_a, plot_b, plot_c, plot_d, plot_e, 
          plot_f, plot_g, plot_h, plot_i,
          nrow = 3, ncol = 3, align = 'v', labels = 'AUTO', font.label = list(size = 9))

ggsave(full_sweep_out, filename = 'figures/full_model_sweep_base.png',
       width = 10.7, height = 9.61)
