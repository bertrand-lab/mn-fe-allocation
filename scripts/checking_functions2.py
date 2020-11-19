############
# this script outlines functions that go through the output to make sure everything looks good. these are all defensive functions.
############

from __future__ import division

import csv
from numpy import genfromtxt
import numpy as np
import scipy.optimize
import pandas as pd
from scipy.integrate import odeint
# import matplotlib.pyplot as plt
import types
from numba import jit
import math
import time
import scipy.integrate as spi
import argparse
from operator import truediv

from phyto_allocation2 import *

def model_mass_balance(general_model_output, p):
    '''
    get all the sources and sinks from the model and make sure they line up
    '''
    print('general model out')
    print(general_model_output)
    ## convert the list of lists into a new list of lists, but transposed so that it can be input into below
#     general_model_output_trans = map(list, zip(*general_model_output))
#     print(general_model_output_trans)

    all_rel_derivs = []

    for env_condition in general_model_output:
        print('env condition state vars')
        print(env_condition[11:])
        print('betas')
        print(env_condition[5:11])
        derivs = phyto_allocation(S = env_condition[11:], # internal species
                             t = 1, # time
                             beta = env_condition[5:11], # optimized set of parameters for protein synthesis
                             I = env_condition[0], # light
                             mnx = env_condition[1], # external manganese
                             fex = env_condition[2], # external iron
                             no3x = env_condition[3], # external nitrate
                             p = p) # parameter inputs
        #rel_derivs = list(map(truediv, derivs, env_condition[11:]))
        print('dxdt = ')
        print(derivs)
        #print('dxdt/x = ')
        #print(rel_derivs)
        all_rel_derivs.append(derivs)
    return(all_rel_derivs)

def add_model_rates(general_model_output, p):
    '''
    calculate the model rates and add them to the pandas dataframe
    '''
    print('adding model rates')
    print(general_model_output)
    all_model_rates = []
    ## convert the list of lists into a new list of lists, but transposed so that it can be input into below
    for env_condition in general_model_output:
        model_rates = phyto_allocation(S = env_condition[11:], # internal species
                             t = 1, # time
                             beta = env_condition[5:11], # optimized set of parameters for protein synthesis
                             I = env_condition[0], # light
                             mnx = env_condition[1], # external manganese
                             fex = env_condition[2], # external iron
                             no3x = env_condition[3], # external nitrate
                             p = p, rates = True) # parameter inputs
        all_model_rates.append(model_rates)
    print(all_model_rates)
    return(all_model_rates)

# format model output
def add_columns_formatting(model_output_general, p):
    # function for including the formatted model output, not just the state variables
   # computing the uptake rates of Fe and Mn
#     fe_uptake_rates = aksnes_cao_uptake_no_numba(kcat_uptake = p[25],
#                                         n_transporters = model_output_general['Tfe'],
#                                         radius_c = p[0],
#                                         diffusion_coef = 0.9e-9*60,
#                                         bulk_substrate = model_output_general['Fex'],
#                                         trans_size = p[1], affinity_return = False)
    mm_internal_energy = (model_output_general['e'])/(p[33] + model_output_general['e'])

    fe_uptake_rates = aksnes_cao_fe_speciation(percent_species_list = [p[59], p[60]],
                            kcat_uptake_list = [p[57], p[58]],
                            n_transporters = model_output_general['Tfe'],
                            radius_c = p[0],
                            diffusion_coef = 0.9e-9*60,
                            bulk_substrate = model_output_general['Fex'],
                            trans_size = p[1],
                            other_n_transporters = model_output_general['Tmn'],
                            avail_space = p[64])
    
    model_output_general['fe_uptake'] = fe_uptake_rates*mm_internal_energy

    mn_uptake_rates = aksnes_cao_uptake_comp(kcat_uptake = p[24],
                                        n_transporters = model_output_general['Tmn'],
                                        radius_c = p[0],
                                        diffusion_coef = 0.9e-9*60,
                                        bulk_substrate = model_output_general['Mnx'],
                                        trans_size = p[1],
                                        other_n_transporters = model_output_general['Tfe'],
                                        avail_space = p[64])
    
    model_output_general['mn_uptake'] = mn_uptake_rates*mm_internal_energy
    
    model_output_general['total_mn_amol'] = 1e18*(model_output_general['Mni'] + model_output_general['A']*p[14] + model_output_general['P']*p[13])/6.0221409e+23
    model_output_general['total_fe_amol'] = 1e18*(model_output_general['Fei'] + model_output_general['Tn']*p[17] + model_output_general['P']*p[15])/6.0221409e+23
    model_output_general['mn_to_fe'] = model_output_general['total_mn_amol']/model_output_general['total_fe_amol']
#    model_output_general['subs_mm_mni'] = model_output_general['Mni']/(model_output_general['Mni'] + p[29] + model_output_general['Fei'])
#    model_output_general['subs_mm_fei'] = model_output_general['Fei']/(model_output_general['Mni'] + p[28] + model_output_general['Fei'])
#    model_output_general['mni_anti_coef'] = model_output_general['subs_mm_mni']/(model_output_general['subs_mm_mni'] + model_output_general['subs_mm_fei'])
#    model_output_general['fei_anti_coef'] = 1 - model_output_general['mni_anti_coef']
#    model_output_general['A_fei'] = model_output_general['fei_anti_coef']*model_output_general['A']
#    model_output_general['A_mni'] = model_output_general['mni_anti_coef']*model_output_general['A']
#    model_output_general['total_mn_amol_subs'] = 1e18*(model_output_general['Mni'] + model_output_general['A_mni']*p[14] + model_output_general['P']*p[13])/6.0221409e+23
#    model_output_general['total_fe_amol_subs'] = 1e18*(model_output_general['Fei'] + model_output_general['A_fei']*p[14] + model_output_general['P']*p[15])/6.0221409e+23
    model_output_general['u_trans'] = 60*24*model_output_general['u']**2

#     fe_affinity = aksnes_cao_uptake_no_numba(kcat_uptake = p[25],
#                                         n_transporters = model_output_general['Tfe'],
#                                         radius_c = p[0],
#                                         diffusion_coef = 0.9e-9*60,
#                                         bulk_substrate = model_output_general['Fex'],
#                                         trans_size = p[1], affinity_return = True)
#     mn_affinity = aksnes_cao_uptake_no_numba(kcat_uptake = p[24],
#                                         n_transporters = model_output_general['Tmn'],
#                                         radius_c = p[0],
#                                         diffusion_coef = 0.9e-9*60,
#                                         bulk_substrate = model_output_general['Mnx'],
#                                         trans_size = p[1], affinity_return = True)

#     if p[40] == 1:# if there is substitutable SOD

#         model_output_general['fe_affinity_quota'] = fe_affinity/model_output_general['total_fe_amol_subs']
#         model_output_general['mn_affinity_quota'] = mn_affinity/model_output_general['total_mn_amol_subs']

#     if p[40] == 0:# if there is no substitutable SOD

#         model_output_general['fe_affinity_quota'] = fe_affinity/model_output_general['total_fe_amol']
#         model_output_general['mn_affinity_quota'] = mn_affinity/model_output_general['total_mn_amol']

    #if p[48] == 1:
    #    fe_uptake_cost_par = p[49]
    #if p[48] == 0:
    #    fe_uptake_cost_par = p[50]
    fe_uptake_cost_par = p[42]
    print(fe_uptake_cost_par)
    model_output_general['A_aa'] = model_output_general['A']*p[8]
    model_output_general['P_aa'] = model_output_general['P']*p[10]
    model_output_general['R_aa'] = model_output_general['R']*p[5]
    model_output_general['Tmn_aa'] = model_output_general['Tmn']*(p[6] + (model_output_general['mn_uptake']/model_output_general['u']**2)*p[43])
    model_output_general['Tfe_aa'] = model_output_general['Tfe']*(p[7] + (model_output_general['fe_uptake']/model_output_general['u']**2)*fe_uptake_cost_par)
    model_output_general['Tmn_aa_no_dyn'] = model_output_general['Tmn']*(p[6])
    model_output_general['Tfe_aa_no_dyn'] = model_output_general['Tfe']*(p[7])
    model_output_general['Tn_aa'] = model_output_general['Tn']*p[9]
    model_output_general['total_aa'] = model_output_general['A_aa'] + model_output_general['P_aa'] + model_output_general['R_aa'] + model_output_general['Tmn_aa'] + model_output_general['Tfe_aa'] + model_output_general['Tn_aa'] + model_output_general['aa']
    model_output_general['total_aa_no_aa'] = model_output_general['A_aa'] + model_output_general['P_aa'] + model_output_general['R_aa'] + model_output_general['Tmn_aa'] + model_output_general['Tfe_aa'] + model_output_general['Tn_aa']
    model_output_general['A_frac'] = model_output_general['A_aa']/model_output_general['total_aa']
    model_output_general['P_frac'] = model_output_general['P_aa']/model_output_general['total_aa']
    model_output_general['R_frac'] = model_output_general['R_aa']/model_output_general['total_aa']
    model_output_general['Tmn_frac'] = model_output_general['Tmn_aa']/model_output_general['total_aa']
    model_output_general['Tmn_no_dyn_frac'] = model_output_general['Tmn_aa_no_dyn']/model_output_general['total_aa']
    model_output_general['Tfe_frac'] = model_output_general['Tfe_aa']/model_output_general['total_aa']
    model_output_general['Tfe_no_dyn_frac'] = model_output_general['Tfe_aa_no_dyn']/model_output_general['total_aa']
    model_output_general['Tn_frac'] = model_output_general['Tn_aa']/model_output_general['total_aa']
    # model_output_general['ind_frac'] = model_output_general['ind']/model_output_general['total_aa']

    return(model_output_general)
