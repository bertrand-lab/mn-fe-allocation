### converting Twining Fe/C ratios to upper and lower bound for diatom Fe:C ratios

library(ggplot2)

pgc_cell <- function(vol_um3){
  0.288*vol_um3^0.811
}

molc_cell <- function(pgc_by_cell){
  pgc_by_cell*(1e-12)*(1/12.011)
}

umolfe_cell <- function(molc_by_cell, umol_fe_cell){
  molc_by_cell*umol_fe_cell
}

molecules_fe_cell <- function(umolfe_cell_out){
  umolfe_cell_out*1e-6*6.022e23
}

umol_fe_mol_c_to_molecules_fe_per_cell <- function(cell_volume, umol_fe_cell){
  pgc_cell_out <- pgc_cell(vol_um3 = cell_volume)
  molc_cell_out <- molc_cell(pgc_by_cell = pgc_cell_out)
  umolfe_cell_out <- umolfe_cell(molc_by_cell = molc_cell_out, umol_fe_cell = umol_fe_cell)
  molecules_fe_cell_out <- molecules_fe_cell(umolfe_cell_out)
  
  return(molecules_fe_cell_out)
}

umol_fe_mol_c_to_molecules_fe_per_cell(cell_volume = 258, umol_fe_cell = 5.5)
umol_fe_mol_c_to_molecules_fe_per_cell(cell_volume = 258, umol_fe_cell = 30)

df1 <- data.frame(epsilon_p = c(0.3, 0.05),
           total_fe = c(umol_fe_mol_c_to_molecules_fe_per_cell(cell_volume = 258, umol_fe_cell = 5.5),
                        umol_fe_mol_c_to_molecules_fe_per_cell(cell_volume = 258, umol_fe_cell = 30)))
lm(epsilon_p ~ total_fe, data = df1)