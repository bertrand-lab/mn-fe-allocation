library(dplyr)
library(ggplot2)
library(magrittr)

# plots and equations for phenomenological model in supplementary materials

growth_func <- function(theta_fe, fex, mnx, kfe = 5, kmn = 5){
  out_val <- theta_fe*((fex/(kfe + fex)) - (mnx/(kmn + mnx))) + (mnx/(kmn + mnx))
  
  return(out_val)
}

growth_func_full <- function(theta_fe, fex, mnx, kfe = 5, kmn = 5){
  out_val <- fex/(kfe + fex) - mnx/(kmn + mnx)
  
  return(out_val)
}

growth_mn_deriv <- function(theta_fe, fex, mnx, kfe = 5, kmn = 5){
  out_val <- -1*kfe*(theta_fe - 1)/(mnx + kfe)^2
  return(out_val)
}
growth_fe_deriv <- function(theta_fe, fex, mnx, kfe = 5, kmn = 5){
  out_val <- kfe*(theta_fe)/(fex + kfe)^2
  return(out_val)
}

interpolate_phenom_model <- function(x_val, y_val, z_val, x_name, y_name, z_name, df){
  
  
  # unique light levels
  df_master <- data.frame(x_name = numeric(), 
                          y_name = numeric(),
                          z_name = numeric())
  
  x_val_filt <- df$x_val
  y_val_filt <- df$y_val
  z_val_filt <- df$z_val
  
  # Create grid
  grd <- expand.grid(x = seq(min(x_val_filt),
                             max(x_val_filt),
                             length=50), 
                     y = seq(min(y_val_filt), 
                             max(y_val_filt), 
                             length=50))
  
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
                 xo = grd$x, 
                 yo = grd$y)
  
  new_df <- data.frame(name1 = res$x, 
                       name2 = res$y, 
                       name3 = res$z)
  
  
  return(new_df)
  
}

get_mn_fe_df <- function(fex_seq, theta_fe, k_vec = c(5, 5)){
  
  all_growth_vals <- c()
  all_growth_vals_deriv <- c()
  all_growth_mn_deriv <- c()
  fex_seq_vals <- c()
  mnx_seq_vals <- c()
  
  for(i in 1:length(fex_seq)){
    
    growth_vec <- growth_func(theta_fe = theta_fe, 
                              fex = fex_seq,
                              mnx = fex_seq[i],
                              kfe = k_vec[1],
                              kmn = k_vec[2])
    growth_vec_deriv <- growth_func_full(theta_fe = theta_fe, 
                                         fex = fex_seq,
                                         mnx = fex_seq[i],
                                         kfe = k_vec[1],
                                         kmn = k_vec[2])
    
    all_growth_vals <- c(growth_vec, all_growth_vals)
    all_growth_vals_deriv <- c(growth_vec_deriv, all_growth_vals_deriv)
    fex_seq_vals <- c(fex_seq, fex_seq_vals)
    mnx_seq_vals <- c(rep(fex_seq[i], length(fex_seq)), mnx_seq_vals)
  }
  
  df_out <- data.frame(u_trans = all_growth_vals,
                       fex = fex_seq_vals,
                       mnx = mnx_seq_vals,
                       u_trans_deriv = all_growth_vals_deriv)
  return(df_out)
}

list_of_theta_fe <- function(fex_seq, theta_fe_seq = c(0.1, 0.5, 0.9),
                             k_vec = c(5, 5)){
  
  df_out <- data.frame(u_trans = numeric(),
                       fex = numeric(),
                       mnx = numeric(),
                       cost_inter = numeric())
  
  for(j in 1:length(theta_fe_seq)){
    df_temp <- get_mn_fe_df(fex_seq = fex_seq, 
                            theta_fe = theta_fe_seq[j],
                            k_vec = k_vec)
    df_temp$cost_inter <- rep(theta_fe_seq[j], nrow(df_temp))
    
    df_out <- rbind(df_temp, df_out)
  }
  
  df_out$cost_inter <- factor(df_out$cost_inter)
  
  levels(df_out$cost_inter) <- paste('Fe Protein Cost: ', levels(df_out$cost_inter)) 
  
  return(df_out)
}

## plotting phenom model with background costs and ratios of costs

fex_conc <- seq(0, 50, by = 0.5)

test <- get_mn_fe_df(fex_seq = fex_conc, theta_fe = 0.5)

phi05 <- list_of_theta_fe(fex_seq = fex_conc, theta_fe_seq = c(0.5))
phi05_p <- phi05 %>% 
  ggplot(aes(x = fex, y = mnx, fill = u_trans)) +
  theme_bw() +
  # ggtitle(expression(paste('Fe Cost Parameter (', psi[italic(Fe)], ") = 0.5"))) +
  ggtitle('Similar Cellular Costs\nacross Micronutrients') +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = u_trans),
               colour = 'black', lwd = 0.4,
               alpha = 0.7) +
  scale_fill_distiller(palette = 'Greens', 
                       name = 'Relative\nGrowth Rate') +
  labs(y = 'Relative [Mn]', x = 'Relative [Fe]') +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        title = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill = 'white'));phi05_p


phi09 <- list_of_theta_fe(fex_seq = fex_conc, theta_fe_seq = c(0.9),
                          c(5, 5))
phi09_p <- phi09 %>% 
  ggplot(aes(x = fex, y = mnx, fill = u_trans)) +
  theme_bw() +
  # facet_grid(~cost_inter) +
  # ggtitle(expression(paste('Fe Cost Parameter (', psi[italic(Fe)], ") = 0.9"))) +
  ggtitle('Dissimilar Cellular Costs\nacross Micronutrients') +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = u_trans),
               colour = 'black', lwd = 0.4,
               alpha = 0.7) +
  scale_fill_distiller(palette = 'Greens', 
                       name = 'Relative\nGrowth Rate') +
  labs(y = 'Relative [Mn]', x = 'Relative [Fe]') +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        title = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill = 'white'));phi09_p



phi05_high_k <- list_of_theta_fe(fex_seq = fex_conc, theta_fe_seq = c(0.5),
                          c(20, 20))
phi05_high_k_p <- phi05_high_k %>% 
  ggplot(aes(x = fex, y = mnx, fill = u_trans)) +
  theme_bw() +
  # facet_grid(~cost_inter) +
  # ggtitle(expression(paste('Half Saturation Constant (', italic(K), ") = 20"))) +
  ggtitle('High Background Cost') +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = u_trans),
               colour = 'black', lwd = 0.4,
               alpha = 0.7) +
  scale_fill_distiller(palette = 'Greens', 
                       name = 'Relative\nGrowth Rate') +
  labs(y = 'Relative [Mn]', x = 'Relative [Fe]') +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        title = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill = 'white'));phi05_high_k_p
# 
phi05_low_k <- list_of_theta_fe(fex_seq = fex_conc, theta_fe_seq = c(0.5),
                          c(1, 1))
phi05_low_k_p <- phi05_low_k %>% 
  ggplot(aes(x = fex, y = mnx, fill = u_trans)) +
  theme_bw() +
  # facet_grid(~cost_inter) +
  # ggtitle(expression(paste('Half Saturation Constant (', italic(K), ") = 1"))) +
  ggtitle('Low Background Cost') +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = u_trans),
               colour = 'black', lwd = 0.4,
               alpha = 0.7) +
  scale_fill_distiller(palette = 'Greens', 
                       name = 'Relative\nGrowth Rate') +
  labs(y = 'Relative [Mn]', x = 'Relative [Fe]') +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        title = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.background = element_rect(fill = 'white'));phi05_low_k_p


