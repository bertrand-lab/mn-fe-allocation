### script for showing the equation relationships
R <- 10
kcat_ros <- 10000
eps_p <- 0.2
eps_a <- 0.00001

first_function <- function(R, kcat_ros, eps_p, eps_a, v_e, A = 775051){
  
  v_ros <- kcat_ros*A

  w_u <- (eps_p*v_e - eps_a*v_ros)/(eps_p*v_e + eps_a*v_ros)
  
  return(w_u)
  }

png(filename = 'figures/ros_equation_output.png', 
    width = 20, height = 20, units = 'cm', res = 600)

y_val <- first_function(R = 10, kcat_ros = kcat_ros, eps_p = 0.2, eps_a = eps_a, v_e = seq(from = 1, to = 1e7, by = 10000))
x_val <- seq(from = 1, to = 1e7, by = 10000)
y_val_adjusted <- y_val[which(y_val > 0)]
x_val_adjusted <- x_val[which(y_val > 0)]

p_w <- 2*R^-y_val_adjusted/(R^(-2*y_val_adjusted) + 1)


par(mfrow = c(2, 1), mar = c(4, 6, 3, 3))
plot(x = x_val_adjusted, y = y_val_adjusted, type = 'l',
     xlab = 'Electron flux and antioxidant mismatch (Equations 17, 18, 19)', ylab = expression(paste('Scaled value representing \nelectron flux and antioxidant mismatch (', omega, ')')),
     main = 'A                                                                                        ')

plot(x = y_val_adjusted, 
     y = p_w, 
     type = 'l',
     ylab = expression(paste('Protein Synthesis Penalty Multiplier (p', omega, ')')),
     xlab = expression(paste('Scaled value representing electron flux and antioxidant mismatch (', omega, ')')),
     main = 'B                                                                                        ')
dev.off()
