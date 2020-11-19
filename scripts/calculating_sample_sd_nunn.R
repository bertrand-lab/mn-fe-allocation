# estimating h (error term in ABC) from Nunn and Cohen data

# https://www.sagepub.com/sites/default/files/upm-binaries/6427_Chapter_4__Lee_(Analyzing)_I_PDF_6.pdf

# Var(r) = Var(X/Y) = r^2 * (Var(y)/y^2 + Var(x)/x^2 - 2*Cov(x,y)/xy)
# where r, x, and y are the means of their distributions

# taking the mean across biological replicates from Nunn et al 2013


nunn_output_high_fe <- nunn_output %>% 
  filter(Fex > 2000) %>% 
  dplyr::rename(norm_factor = percentage) %>% 
  dplyr::select(norm_factor, coarse_grain)

nunn_output_relative <- nunn_output %>%
  inner_join(nunn_output_high_fe, by = c('coarse_grain')) %>% 
  mutate(relative_change = percentage/norm_factor)

nunn_output_relative %>% +
  ggplot(aes(x = coarse_grain, y = relative_change)) +
  geom_point()


## subset out the low fe
nunn_output_low_fe <- nunn_output %>% 
  filter(Fex < 2000)

## subset out the high fe
nunn_output_high_fe <- nunn_output %>% 
  filter(Fex > 2000) %>% 
  rename(percentage_high = percentage)

## assign a fixed joiner id
nunn_output_high_fe$joiner_id <- rep(1:4, 6)

number_monte <- 10000

monte_estimates <- c()
for(i in 1:number_monte){
  ## assign a random joiner id
  random_joiner <- sample(x = c(1, 2, 3, 4), size = 4, replace = FALSE)
  nunn_output_low_fe$joiner_id <- rep(random_joiner, 6)
  
  
  ## calculate the mean difference for all values
  nunn_output_monte <- nunn_output_high_fe %>%
    ungroup() %>% 
    dplyr::select(percentage_high, joiner_id,
                  coarse_grain) %>% 
    inner_join(nunn_output_low_fe, 
               by = c('joiner_id', 'coarse_grain')) %>% 
    ## calculate the difference from high to low fe, square it
    mutate(sq_difs = (percentage_high - percentage)^2)
  
  mean_monte <- mean(nunn_output_monte$sq_difs)
  
  ## sum the squared differences and divide by (n-1)
  n_minus_1 <- length(nunn_output_monte$sq_difs) - 1
  
  monte_estimate_i <- sqrt(sum((nunn_output_monte$sq_difs - mean_monte)^2)/n_minus_1)
  
  ## append the sample standard deviation
  monte_estimates <- c(monte_estimates, monte_estimate_i)
}

monte_estimates %>% mean()







