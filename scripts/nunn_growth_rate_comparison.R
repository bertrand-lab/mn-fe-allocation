posterior_pars_nunn <- joined_nunn_par_out_q_relative #%>%
  # filter(in_set == TRUE)

post_model_out_nunn <- par_sets_model %>% inner_join(posterior_pars_nunn, 
                                                by = c('epsilon_a', 'avail_space', 'cost_par'))

nunn_growth_data <- data.frame(Fex = c(87, 2876),
                               u_trans_mean = c(0.68, 1.71),
                               u_trans_sd = c(0.26*2, 0.25*2),
                               u_trans_se = c(0.26, 0.25))

post_model_out_nunn %>% 
  ggplot(aes(x = factor(Fex), y = u_trans)) +
  geom_boxplot(width = 0.3) +
  facet_grid(~in_set)

summaries_nunn_model_post <- post_model_out_nunn %>% 
  filter(in_set == TRUE) %>% 
  group_by(Fex) %>% 
  summarize(u_trans_mean = mean(u_trans),
            u_trans_sd = sd(u_trans),
            u_trans_se = u_trans_sd/sqrt(n()))

summaries_nunn_model_prior <- post_model_out_nunn %>% 
  group_by(Fex) %>% 
  summarize(u_trans_mean = mean(u_trans),
            u_trans_sd = sd(u_trans),
            u_trans_se = u_trans_sd/sqrt(n()))

summaries_nunn_model_prior$model <- rep('Allocation Model Prior', 2)
summaries_nunn_model_post$model <- rep('Allocation Model', 2)
nunn_growth_data$model <- rep('Culture Data', 2)

all_nunn_data <- rbind(nunn_growth_data, summaries_nunn_model_post)

all_nunn_data %>% 
  ggplot(aes(x = factor(Fex), y = u_trans_mean, shape = model,
             colour = model)) +
  geom_errorbar(aes(ymin = u_trans_mean - u_trans_sd, 
                    ymax = u_trans_mean + u_trans_sd,
                    x = factor(Fex)), 
                position = position_dodge(width = 0.3), width = 0.1) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  theme_bw() + 
  xlab('Fe Prime') +
  ylab('Growth Rate (per day)') +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 14))






