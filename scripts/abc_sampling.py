## parameter inference for Approximate Bayesian Computation

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
import random

### import the allocation model functions
from phyto_allocation2 import *
from multi_step_opt2 import *
from checking_functions2 import *

def add_columns_formatting_abc(model_output_general, p):
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
    model_output_general['cost_par'] = p[42]
    model_output_general['avail_space'] = p[64]
    model_output_general['epsilon_a'] = p[35]
    
    return(model_output_general)

# generate the parameter sets
def gen_par_sets(eps_low, cost_low, space_low, eps_high, cost_high, space_high):
    '''
    generating parameter sets using bounds (prior distribution)
    '''
    print(eps_low)
    epsilon_a = random.uniform(eps_low, eps_high)
    fe_internal_cost = random.uniform(cost_low, cost_high)
    space_avail = random.uniform(space_low, space_high)
    
    return([epsilon_a, fe_internal_cost, space_avail])

#gener_par = gen_par_sets(cost_low = 0.01, eps_low = 0.0001, space_low = 0.001, 
#                         eps_high = 0.1, cost_high = 1, space_high = 0.15)

def mod_par_sets(para_set, gen_par_set_out):
    '''
    go into the parameter set and alter it with the generated value
    '''
    #epsilon_a = p[35]
    #avail_space = p[64]
    para_set[35] = gen_par_set_out[0]
    para_set[42] = gen_par_set_out[1] # both the Fe and Mn uptake rate costs are changed simulataneously
    para_set[43] = gen_par_set_out[1]
    para_set[64] = gen_par_set_out[2]
    
    return(para_set)

def gen_and_mod_par_set(para_set1, para_set2, para_set3):
    '''
    generate parameter sets and modify par set input wrapper function
    '''
    gener_par_set_i = gen_par_sets(eps_low = 0.00001, cost_low = 0.001, space_low = 0.001, 
                                   eps_high = 0.1, cost_high = 16, space_high = 0.15)
    mod_para_set1 = mod_par_sets(para_set = para_set1,
                               gen_par_set_out = gener_par_set_i)
    mod_para_set2 = mod_par_sets(para_set = para_set2,
                               gen_par_set_out = gener_par_set_i)
    mod_para_set3 = mod_par_sets(para_set = para_set3,
                               gen_par_set_out = gener_par_set_i)

    return([mod_para_set1, mod_para_set2, mod_para_set3])


def generate_and_run_model(number_abc_loops, name_for_out, paras_meta, paras_nunn, paras_cohen,
                           t, S_init, cons, bnds,
                           total_loops = 10, 
                           gpr_eval_number = 10, 
                           number_trials = 10):
    # print counter
    printcounter = 0
    
    # make two new master dataframes
    meta_master_df = []
    nunn_master_df = []
    cohen_master_df = []
    
    # for i loops run model
    # add to
    for abc in range(number_abc_loops):
        # generate parameters
        paras_mod_i = gen_and_mod_par_set(para_set1 = paras_meta, para_set2 = paras_nunn, para_set3 = paras_cohen)

        paras_meta_i = paras_mod_i[0]
        paras_nunn_i = paras_mod_i[1]
        paras_cohen_i = paras_mod_i[2]
        
        # run entire model
        model_output_meta_day1 = run_ode_opt_range(I = [75], 
                                         mnx = [260], 
                                         fex = [1010],
                                         no3x = [10000000], time_step = t, p = paras_meta_i,
                                         x_init = None, S_init = S_init, cons = cons, bnds = bnds, 
                                         total_loops = total_loops,
                                         gpr_eval_number = gpr_eval_number,
                                         number_trials = number_trials)
        model_output_meta_day3 = run_ode_opt_range(I = [75], 
                                         mnx = [210], 
                                         fex = [470],
                                         no3x = [10000000], time_step = t, p = paras_meta_i,
                                         x_init = None, S_init = S_init, cons = cons, bnds = bnds, 
                                         total_loops = total_loops,
                                         gpr_eval_number = gpr_eval_number,
                                         number_trials = number_trials)
        model_output_nunn_felow = run_ode_opt_range(I = [150], 
                                         mnx = [3000], 
                                         fex = [87],
                                         no3x = [10000000], time_step = t, p = paras_nunn_i,
                                         x_init = None, S_init = S_init, cons = cons, bnds = bnds, 
                                         total_loops = total_loops,
                                         gpr_eval_number = gpr_eval_number,
                                         number_trials = number_trials)
        model_output_nunn_fehigh = run_ode_opt_range(I = [150], 
                                         mnx = [3000], 
                                         fex = [2876],
                                         no3x = [10000000], time_step = t, p = paras_nunn_i,
                                         x_init = None, S_init = S_init, cons = cons, bnds = bnds, 
                                         total_loops = total_loops,
                                         gpr_eval_number = gpr_eval_number,
                                         number_trials = number_trials)
        model_output_cohen_felow = run_ode_opt_range(I = [110],
                                         mnx = [3000],
                                         fex = [20],
                                         no3x = [10000000], time_step = t, p = paras_cohen_i,
                                         x_init = None, S_init = S_init, cons = cons, bnds = bnds,
                                         total_loops = total_loops,
                                         gpr_eval_number = gpr_eval_number,
                                         number_trials = number_trials)
        model_output_cohen_fehigh = run_ode_opt_range(I = [110],
                                         mnx = [3000],
                                         fex = [2700],
                                         no3x = [10000000], time_step = t, p = paras_cohen_i,
                                         x_init = None, S_init = S_init, cons = cons, bnds = bnds,
                                         total_loops = total_loops,
                                         gpr_eval_number = gpr_eval_number,
                                         number_trials = number_trials)

        # format model outputs
        model_output_meta_day1_form = add_columns_formatting_abc(model_output_meta_day1, paras_meta_i)
        model_output_meta_day3_form = add_columns_formatting_abc(model_output_meta_day3, paras_meta_i)        

        model_output_nunn_felow_form = add_columns_formatting_abc(model_output_nunn_felow, paras_nunn_i)        
        model_output_nunn_fehigh_form = add_columns_formatting_abc(model_output_nunn_fehigh, paras_nunn_i)

        model_output_cohen_felow_form = add_columns_formatting_abc(model_output_cohen_felow, paras_cohen_i)
        model_output_cohen_fehigh_form = add_columns_formatting_abc(model_output_cohen_fehigh, paras_cohen_i)


        # append to dataframe
        meta_master_df.append(model_output_meta_day1_form)
        meta_master_df.append(model_output_meta_day3_form)

        nunn_master_df.append(model_output_nunn_felow_form)
        nunn_master_df.append(model_output_nunn_fehigh_form)

        cohen_master_df.append(model_output_cohen_felow_form)
        cohen_master_df.append(model_output_cohen_fehigh_form)

        if printcounter == 10:
            printcounter = 0
            # format lists of df to a pd df
            meta_df = pd.concat(meta_master_df)
            nunn_df = pd.concat(nunn_master_df)
            cohen_df = pd.concat(cohen_master_df)
            # write output to a csv
            meta_df.to_csv('../data/abc_out/meta_par_gen_' + name_for_out + '.csv')
            nunn_df.to_csv('../data/abc_out/nunn_par_gen_' + name_for_out + '.csv')
            cohen_df.to_csv('../data/abc_out/cohen_par_gen_' + name_for_out + '.csv')

        printcounter += 1
        
    meta_df = pd.concat(meta_master_df)
    nunn_df = pd.concat(nunn_master_df)
    cohen_df = pd.concat(cohen_master_df)
    
    # write out the final df
    meta_df.to_csv('../data/abc_out/meta_par_gen_' + name_for_out + '.csv')
    nunn_df.to_csv('../data/abc_out/nunn_par_gen_' + name_for_out + '.csv')
    cohen_df.to_csv('../data/abc_out/cohen_par_gen_' + name_for_out + '.csv')

