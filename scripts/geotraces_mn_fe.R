# this script plots geotraces dissolved mn and fe data
# light levels and MLD were determined from another script (****.py)

# this script also imports data from E. Achterberg group from FISH

library(readxl)
library(dplyr)
library(ggplot2)
library(oce)
library(gridExtra)
library(scales)
library(gtable)
library(grid)

###### ACHTERBERG GROUP FISH DATA

fish <- read.csv("data/oceanographic_data/FISH_trace_metal_data_JR271_JR274.csv")

fish_no_contam <- fish %>% filter(fe_nm < 10000, cruise == 'JR274')

fish_plot <- fish_no_contam %>% 
  filter(cruise == 'JR274', type == 'dissolved') %>%
  ggplot(aes(x = fe_nm, y = mn_nm)) + 
  geom_point(colour = 'white', size = 2) +
  geom_errorbarh(aes(xmin = fe_nm - fe_sd, xmax = fe_nm + fe_sd), alpha = 0.3) +
  geom_errorbar(aes(ymin = mn_nm - mn_sd, ymax = mn_nm + mn_sd),  alpha = 0.3) +
  geom_point(alpha = 0.3, size = 2) +
  scale_y_log10() + 
  scale_x_log10() +
  theme_bw() +
  # facet_grid(~cruise) +
  xlab(expression(paste('[Fe] nmol ', L^-1))) +
  ylab(expression(paste('[Mn] nmol ', L^-1))) +
  theme(text = element_text(size = 14));fish_plot

###### GEOTRACES DATA

# reading in the geotraces discrete data
data <- read.csv('data/oceanographic_data/GEOTRACES_IDP2017_v2/discrete_sample_data/excel/GEOTRACES_IDP2017_v2_Discrete_Sample_Data.csv')

data$cruise_station <- paste(data$Cruise, data$Station, sep = '_')

data2 <- data[grepl(pattern = 'GIPY05', x = data$Cruise) |
                grepl(pattern = 'GI04', x = data$Cruise),]
data3 <- data[, c(1:10)]

# reading in the mld data
station_mld <- read.csv("data/oceanographic_data/geotraces_mld_light_by_station_discrete_data.csv")
station_mld$mld <- station_mld$mld %>% as.character() %>% as.numeric()
station_mld$cruise_station <- paste(station_mld$cruise, station_mld$station, sep = '_')
station_mld$mld2 <- ifelse(station_mld$mld < 0, yes = 0, no = station_mld$mld)

unique_disc <- data$Station %>% unique() %>% as.character()
unique_ctd <- station_mld$station %>% unique() %>% as.character()

## merge the discrete data MLDs
merged_data <- inner_join(x = data, y = station_mld, by = 'cruise_station')

# filtering the data for southern ocean cruises
so_data <- merged_data %>% 
  # filter(DEPTH..m. < 300, 
         filter(Cruise == 'GIPY05' | Cruise == 'GIPY06' | Cruise == 'GIPY04' | Cruise == 'GI04') %>%
  select('Cruise','Station', 'yyyy.mm.ddThh.mm.ss.sss', 'Longitude..degrees_east.', 'Latitude..degrees_north.', 'DEPTH..m.', 'Mn_D_CONC_BOTTLE..nmol.kg.', 'STANDARD_DEV.26', 'Fe_D_CONC_BOTTLE..nmol.kg.', 'STANDARD_DEV.21', 'Mn_D_CONC_FISH..nmol.kg.', 'PRESSURE..dbar.', 'mld2', 'par', 'kd_490', 'cruise_station')
# %>% 
  # mutate(median_light_level = par*exp(-kd_490*mld2/2))

so_light <- so_data %>% 
  filter(Latitude..degrees_north. < -50) %>%
  mutate(in_mld = ifelse(DEPTH..m. > mld2, TRUE, FALSE)) %>% 
  filter(in_mld == TRUE) %>% 
  group_by(cruise_station) %>% 
  # summarize(median_light_level = (1e6/(60*1140))*mean(kd_490)*exp(-mean(par)*mean(mld2)/2),
  summarize(median_light_level = median(par)*1e6/(24*60*60)*exp(-median(kd_490)*median(mld2)/2),
            median_mn_level = median(Mn_D_CONC_BOTTLE..nmol.kg., na.rm = TRUE),
            median_fe_level = median(Fe_D_CONC_BOTTLE..nmol.kg., na.rm = TRUE)) %>%
  ggplot(aes(y = median_mn_level, 
             x = median_fe_level)) + 
  geom_point(size = 4, alpha = 0.9, aes(colour = median_light_level)) +
  scale_colour_gradient(name = expression(paste('Median Mixed Layer Light (uEin', m^-2, sec^-1, ")")),
                        low = 'blue', high = 'yellow') +
  # xlim(0, 1.5) +
  # ylim(0, 1.1) +
  scale_y_log10() + 
  scale_x_log10() +
  xlab(expression(paste('[Fe] nmol ', kg^-1))) +
  ylab(expression(paste('[Mn] nmol ', kg^-1))) +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "top",
        legend.text = element_text(angle = 45, size = 14, hjust = 1));so_light

so_light_df <- so_data %>% 
  filter(Latitude..degrees_north. < -50) %>%
  mutate(in_mld = ifelse(DEPTH..m. > mld2, TRUE, FALSE)) %>% 
  filter(in_mld == TRUE) %>% 
  group_by(cruise_station) %>% 
  summarize(median_light_level = median(par)*1e6/(24*60*60)*exp(-median(kd_490)*median(mld2)/2),
            median_mn_level = median(Mn_D_CONC_BOTTLE..nmol.kg., na.rm = TRUE),
            median_fe_level = median(Fe_D_CONC_BOTTLE..nmol.kg., na.rm = TRUE))

options(scipen=10000)


ggplot(data = so_light_df, aes(y = median_mn_level, 
           x = median_fe_level)) + 
  geom_point(size = 4, alpha = 0.8, aes(colour = median_light_level), 
             shape = 'square') +
  geom_point(data = fish_no_contam, aes(x = fe_nm, y = mn_nm),
             alpha = 0.3) +
  scale_colour_gradient(name = expression(paste('Median Mixed Layer Light (uE', m^-2, sec^-1, ")")),
                        low = 'blue', high = 'yellow') +
  # xlim(0, 1.5) +
  # ylim(0, 1.1) +
  scale_y_log10() + 
  scale_x_log10() +
  xlab(expression(paste('[Fe] nmol ', kg^-1))) +
  ylab(expression(paste('[Mn] nmol ', kg^-1))) +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "top",
        legend.text = element_text(angle = 45, size = 14, hjust = 1)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.3) +
  scale_shape_manual(name="")


### modifying dataframe so they can be together in one
names(so_light_df)
names(fish_no_contam)

so_light_df$method <- rep("GEOTRACES Median Mixed Layer", nrow(so_light_df))
fish_no_contam$method <- rep("Tow-FISH Surface", nrow(fish_no_contam))

fish_no_contam_sub <- fish_no_contam %>% 
  filter(type == 'dissolved') %>% 
  select(fe_nm, mn_nm, method) %>% 
  rename(median_mn_level = mn_nm,
         median_fe_level = fe_nm) 

empty_fish_df <- data.frame(matrix(nrow = nrow(fish_no_contam_sub),
                                   ncol = 2))

names(empty_fish_df) <- names(so_light_df)[1:2]

combined_fish <- cbind(empty_fish_df, fish_no_contam_sub)

combined_fish_geo <- rbind(combined_fish, so_light_df)

# SOME annoying plotting to make ggplot legend locations different on the same figure

#https://stackoverflow.com/questions/13143894/how-do-i-position-two-legends-independently-in-ggplot

# plot with just the light legend
p1 <- ggplot(data = combined_fish_geo, aes(y = median_mn_level, 
                                           x = median_fe_level,
                                           shape = method,
                                           size = method)) + 
  geom_point(alpha = 0.7, aes(colour = median_light_level)) +
  scale_colour_gradient(name = expression(paste('Median Mixed Layer Light (uE', m^-2, sec^-1, ")")),
                        low = 'blue', high = 'yellow') +
  # xlim(0, 1.5) +
  # ylim(0, 1.1) +
  scale_y_log10() + 
  scale_x_log10() +
  xlab(expression(paste('[Fe] nmol ', kg^-1))) +
  ylab(expression(paste('[Mn] nmol ', kg^-1))) +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "top",
        legend.text = element_text(angle = 45, size = 14, hjust = 1)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.3) +
  scale_shape_manual(name = "Sampling Method", values = c("square", "circle")) +
  scale_size_manual(values = c(4, 2)) +
  # scale_shape_manual(name="", values = c("square", "circle"))
  guides(shape = FALSE, size = FALSE);p1

# plot with just the shape legend
p2 <- ggplot(data = combined_fish_geo, aes(y = median_mn_level, 
                               x = median_fe_level,
                               shape = method,
                               size = method)) + 
  geom_point(alpha = 0.7, aes(colour = median_light_level)) +
  scale_colour_gradient(name = expression(paste('Median Mixed Layer Light (uE', m^-2, sec^-1, ")")),
                        low = 'blue', high = 'yellow') +
  # xlim(0, 1.5) +
  # ylim(0, 1.1) +
  scale_y_log10() + 
  scale_x_log10() +
  xlab(expression(paste('[Fe] nmol ', kg^-1))) +
  ylab(expression(paste('[Mn] nmol ', kg^-1))) +
  theme_bw() +
  # theme(text = element_text(size = 14), legend.position = "top",
  #       legend.text = element_text(angle = 45, size = 14, hjust = 1)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.3) +
  # scale_shape_manual(name="", values = c("square", "circle"))
  guides(colour = FALSE, size = FALSE) +
  scale_shape_manual(name = "Sampling Method", values = c("square", "circle")) +
  scale_size_manual(values = c(4, 3)) +
  theme(legend.position = c(0.55, 0.1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"));p2

p3 <- ggplot(data = combined_fish_geo, aes(y = median_mn_level, 
                                           x = median_fe_level,
                                           shape = method,
                                           size = method)) + 
  geom_point(alpha = 0.7, aes(colour = median_light_level)) +
  scale_colour_gradient(name = expression(paste('Median Mixed Layer Light (uE', m^-2, sec^-1, ")")),
                        low = 'blue', high = 'yellow') +
  # xlim(0, 1.5) +
  # ylim(0, 1.1) +
  scale_y_log10() + 
  scale_x_log10() +
  xlab(expression(paste('[Fe] nmol ', kg^-1))) +
  ylab(expression(paste('[Mn] nmol ', kg^-1))) +
  theme_bw() +
  # theme(text = element_text(size = 14), legend.position = "top",
  #       legend.text = element_text(angle = 45, size = 14, hjust = 1)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.3) +
  # scale_shape_manual(name="", values = c("square", "circle"))
  guides(colour = FALSE, size = FALSE, shape = FALSE) +
  scale_shape_manual(name = "Sampling Method", values = c("square", "circle")) +
  scale_size_manual(values = c(4, 2)) +
  theme(legend.position = c(0.55, 0.1),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"));p3

## GETTING legends from the plots above

leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box") 
leg2 <- gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box") 

# adding legends back to plots

plotNew <- p3 + 
  annotation_custom(grob = leg2, xmin = 1, xmax = 4, ymin = -1, ymax = 0.2)

# add the light legend
plotNew <- arrangeGrob(leg1, plotNew,
                       heights = unit.c(leg1$height, unit(1, "npc") -  leg1$height), ncol = 1)



grid.newpage()
grid.draw(plotNew)

mn_fe_observed_concentrations_plot <- ggplot(data = combined_fish_geo, aes(y = median_mn_level, 
                                     x = median_fe_level,
                                     shape = method,
                                     size = method)) + 
  geom_point(alpha = 0.7, aes(colour = median_light_level)) +
  # xlim(0, 1.5) +
  # ylim(0, 1.1) +
  scale_y_log10() + 
  scale_x_log10() +
  xlab(expression(paste('[Fe] nmol ', kg^-1))) +
  ylab(expression(paste('[Mn] nmol ', kg^-1))) +
  theme_bw() +
  # theme(text = element_text(size = 14), legend.position = "top",
  #       legend.text = element_text(angle = 45, size = 14, hjust = 1)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.3) +
  # scale_shape_manual(name="", values = c("square", "circle"))
  guides(size = FALSE) +
  scale_shape_manual(name = "Sampling Method", values = c("square", "circle")) +
  scale_size_manual(values = c(4, 2)) +
  scale_colour_gradient(name = expression(paste('Median Mixed \nLayer Light (uEin', m^-2, sec^-1, ")")),
                        low = 'darkblue', high = 'yellow',
                        guide = guide_colourbar(direction = "horizontal")) +
  theme(legend.position = c(0.3, 0.8),
        legend.spacing = unit(0.02, 'cm'),
        legend.box.background = element_rect(colour = 'black', fill = alpha('white', 0)),
        legend.background = element_rect(fill=alpha('white', 0.7))) +
  guides(shape = guide_legend(override.aes = list(size = 5,
                                                  colour = c("darkblue", "grey30")), 
                              order = 1));mn_fe_observed_concentrations_plot

ggsave(mn_fe_observed_concentrations_plot, filename = 'figures/fig1_mn_fe.png',
       width = 8.01, height = 6.33)
