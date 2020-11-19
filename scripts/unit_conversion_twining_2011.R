tm_quota_umol_per_litre <- 60

# uM * L^-1 into aMol * cell ^-1

###########################
# go to: uM * m^-3
# mutliply by: 1000L / 1m^3
c1 <- 1000/1

# go to: uM * um^-3
# multiply by: 1 m^3 1e6 um^-3
c2 <- (1/1e6)^3

# go to: uM per cell
# multiply by cell volume (*90 um^3 cell ^-1)
c3 <- 90

# go to: Mol per cell
# mulitply by: 1Mol 1e6uMol^-1 
c4 <- (1/1e6)

# go to: Mol to aMol
# multiply by: 1e18 aMol 1Mol^-1
c5 <- (1e18/1)

###########################

tm_quota_umol_per_litre*c1*c2*c3*c4*c5



