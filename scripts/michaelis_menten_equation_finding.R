## getting concentrations at half maximal growth

model_out_mnfe4_model1 <- read_in_data_loop(data_name_string = 'mnfe4_model4_july14', 
                                            dir_to_look = 'data/model_output/', 
                                            string_to_not_include = 'exp')

high_mn_values <- model_out_mnfe4_model1 %>% filter(Mnx == 3000)
high_fe_values <- model_out_mnfe4_model1 %>% filter(Fex == 3000)

## making figure and getting kinetic constants for growth for Fe and Mn
dev.off()

png(filename = 'figures/monod_kinetics.png', 
    width = 20, height = 20, units = 'cm', res = 600)

par(mfrow = c(2, 1))

fex_nls <- nls(u_trans ~ vmax*(Fex)/(kfe + Fex), 
    data = high_mn_values,
    start = list(vmax = 0.25, kfe = 750))

# getting kinetic constants after using NLS
fex_nls_sum <- summary(fex_nls)
fex_nls_sum_se <- fex_nls_sum$coefficients[,2][2] %>% as.numeric()
fex_nls_sum_vmax_se <- fex_nls_sum$coefficients[,2][1] %>% as.numeric()

vmaxfe_val <- coef(fex_nls)[1] %>% as.numeric()
kfe_val <- coef(fex_nls)[2] %>% as.numeric()

plot(x = c(1:3000), y = vmaxfe_val*c(1:3000)/(kfe_val + c(1:3000)), 
     type = 'l', ylim = c(0, 0.3),
     ylab = 'Growth Rate', xlab = 'dFe (pM)',
     main = 'A         Monod-type Kinetics for Fe (dMn = 3000pM)')
points(x = high_mn_values$Fex, y = high_mn_values$u_trans)
text(paste('K = ',round(kfe_val, 3), ' pM ', '(SE = ', 
           round(fex_nls_sum_se, 3), ')', sep = ''), x = 500, y = 0.26)
text(paste('V_max = ',round(vmaxfe_val, 3), ' pM ', '(SE = ', 
           round(fex_nls_sum_vmax_se, 3), ')', sep = ''), x = 2000, y = 0.26)

# getting kinetic constants after using NLS
mnx_nls <- nls(u_trans ~ vmax*(Mnx)/(kmn + Mnx), 
    data = high_fe_values,
    start = list(vmax = 0.25, kmn = 75))

vmaxmn_val <- coef(mnx_nls)[1] %>% as.numeric()
kmn_val <- coef(mnx_nls)[2] %>% as.numeric()

mnx_nls_sum <- summary(mnx_nls)
mnx_nls_sum_se <- mnx_nls_sum$coefficients[,2][2] %>% as.numeric()
mnx_nls_sum_vmax_se <- mnx_nls_sum$coefficients[,2][1] %>% as.numeric()

plot(x = c(1:3000), 
     y = vmaxmn_val*c(1:3000)/(kmn_val + c(1:3000)), 
     main = 'B         Monod-type Kinetics for Mn (dFe = 3000pM)',
     type = 'l',
     ylim = c(0, 0.3),
     ylab = 'Growth Rate', xlab = 'dMn (pM)')
points(x = high_fe_values$Mnx, y = high_fe_values$u_trans)
text(paste('K = ',round(kmn_val, 2), ' pM ', '(SE = ', 
           round(mnx_nls_sum_se,digits = 3), ')', sep = ''), x = 500, y = 0.28)
text(paste('V_max = ',round(vmaxmn_val, 3), ' pM ', '(SE = ', 
           round(mnx_nls_sum_vmax_se, 3), ')', sep = ''), x = 2000, y = 0.28)

dev.off()
