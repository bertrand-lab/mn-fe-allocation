## merging abc results

library(dplyr)

cohen_abc_out <- read.csv("data/abc_intermediate/cohen_abc_par_sets.csv")
meta_abc_out <- read.csv("data/abc_intermediate/meta_abc_par_sets.csv")
nunn_abc_out <- read.csv("data/abc_intermediate/nunn_abc_par_sets_relative.csv")

dim(cohen_abc_out)
dim(meta_abc_out)
dim(nunn_abc_out)

meta_abc_subset <- meta_abc_out[meta_abc_out$cost_par %in% cohen_abc_out$cost_par, ]
nunn_abc_subset <- nunn_abc_out[nunn_abc_out$cost_par %in% cohen_abc_out$cost_par, ]

combined_sum_sq <- meta_abc_subset %>% 
  dplyr::rename(sum_sq_dif_meta = sum_sq_dif) %>% 
  inner_join(nunn_abc_subset %>% 
               dplyr::select(sum_sq_dif, cost_par, avail_space, epsilon_a,
                             u_trans_87, u_trans_2876) %>% 
               dplyr::rename(sum_sq_dif_nunn = sum_sq_dif), 
             by = c("cost_par", "avail_space", "epsilon_a")) %>% 
  inner_join(cohen_abc_out %>% 
               dplyr::select(sum_sq_dif, cost_par, avail_space, epsilon_a) %>% 
               dplyr::rename(sum_sq_dif_cohen = sum_sq_dif), by = c("cost_par", "avail_space", "epsilon_a"))


combined_sum_sq_added_meta_cohen <- combined_sum_sq %>% 
  mutate(sum_sq_dif = sum_sq_dif_meta + sum_sq_dif_cohen)
combined_sum_sq_added_meta_cohen_nunn <- combined_sum_sq %>% 
  mutate(sum_sq_dif = sum_sq_dif_meta + sum_sq_dif_cohen + sum_sq_dif_nunn)
combined_sum_sq_added_meta <- combined_sum_sq %>% 
  mutate(sum_sq_dif = sum_sq_dif_meta)

write.csv(combined_sum_sq_added_meta_cohen, file = "data/abc_intermediate/combined_sum_sq_meta_cohen.csv")
write.csv(combined_sum_sq_added_meta_cohen_nunn, file = "data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv")
write.csv(combined_sum_sq_added_meta, file = "data/abc_intermediate/combined_sum_sq_meta.csv")


cutoff_val <- quantile(combined_sum_sq_added_meta_cohen_nunn$sum_sq_dif, 0.01) %>% as.numeric()

combined_sum_sq_added_meta_cohen_nunn %>% 
  filter(sum_sq_dif < cutoff_val) %>% 
  ggplot(aes(avail_space)) +
  geom_histogram()
combined_sum_sq_added_meta_cohen_nunn %>% 
  filter(sum_sq_dif < cutoff_val) %>% 
  ggplot(aes(cost_par)) +
  geom_histogram()
combined_sum_sq_added_meta_cohen_nunn %>% 
  filter(sum_sq_dif < cutoff_val) %>% 
  ggplot(aes(epsilon_a)) +
  geom_histogram()


# 
# 
# tester <- combined_sum_sq %>% 
#   mutate(sum_sq_dif = sum_sq_dif_meta + sum_sq_dif_cohen + sum_sq_dif_nunn,
#          weights_sq = exp(-(sum_sq_dif^2)/(2*30^2)))
# 
# 
# cutoff_qu <- tester$weights_sq %>% quantile(0.95) %>% as.numeric()
# 
# 
# tester %>% filter(weights_sq > cutoff_qu) %>% 
#   dplyr::select(sum_sq_dif_meta, sum_sq_dif_cohen, sum_sq_dif_nunn) %>% 
#   plot()
# 
# combined_sum_sq_added %>% 
#   ggplot(aes(sum_sq_dif)) +
#   geom_histogram() +
#   scale_x_log10()
# 
# 
h_val <- 2
quantile(exp(-combined_sum_sq_added_meta_cohen_nunn$sum_sq_dif^2/(2*h_val^2)), 0.95)
