from __future__ import division

import pickle
import csv
from numpy import genfromtxt
import numpy as np
import scipy.optimize
import pandas as pd
from scipy.integrate import odeint
#import matplotlib.pyplot as plt
import types
from numba import jit
import math
import time
import scipy.integrate as spi
from scipy.optimize import * 
import argparse
import random
from sklearn.cluster import KMeans
from sklearn.datasets import make_friedman2
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
import sys

from phyto_allocation2 import *
from checking_functions2 import *

#@jit
#def calc_ss(beta, I, mnx, fex, no3x, S_init, time_step, p, short_time, number_int_loops, sum_abs_cutoff, report_val = False):
#    '''
#    function for integrating the system of ODES. note that if short_time is False, then the function loops over
#    number_int_loops until a cutoff is reached for the sum of absolute derivatives
#    '''
#    loop_counter = 0
#    for int_loop in range(number_int_loops):
#        S = odeint(phyto_allocation,
#                       y0 = S_init,
#                       t = time_step,
#                       args = (np.abs(beta),
#                               I,
#                               mnx,
#                               fex,
#                               no3x,
#                               p), mxstep = 1000000)
#        if short_time == True:
#            break
#        elif short_time == False:
#            derivs = phyto_allocation(S = S[-1:].tolist()[0],
#                                      t = time_step,
#                                      beta = np.abs(beta),
#                                      I = I,
#                                      mnx = mnx,
#                                      fex = fex,
#                                      no3x = no3x,
#                                      p = p)
            # checks to see if the sum of abs derivs is decreasing to indicate convergence towards a steady state

#            if loop_counter > 0:
#                if sum([abs(val) for val in derivs]) > sum_abs_derivs:
                    # artificially decrease the apparent growth rate for this trial. This is because with
                    # the given system of ODEs (a function of a set of betas), it seems like it's unstable
                    # and integration isn't going to converge.
#                    print(sum([abs(val) for val in derivs]))
#                    print(sum_abs_derivs)

                   # S[-1, 8] = 100
                   # S[-1, 3] = 100
                   # print('derivatives increasing, convergence not reached')
                   # break

#            print('derivative_vals:')
#            print(derivs)
            # calculates the sum of the absolute derivatives to use as a steady state approximation
#            sum_abs_derivs = sum([abs(val) for val in derivs])
#            print('sum derivs good:')
#            print(sum_abs_derivs)
#            sys.stdout.flush() # flush that shiz

#            if sum_abs_derivs < sum_abs_cutoff:
#                break
#            elif sum_abs_derivs > sum_abs_cutoff:
#                # if this sum is above our threshold, use the end of the last integration and re-run the loop.
#                S_init = S[-1:].tolist()[0]
#                loop_counter += 1
#                #print('extending integration time')
#    if report_val == True and short_time == False:
#        print('sum abs derivs:')
#        print(sum_abs_derivs)
#        print(loop_counter)
#    S_sub = np.array([S[-1, 8], S[-1, 3]])
#    return(S_sub)


@jit #(nopython=True)
def optimize_growth_ak(beta, I, mnx, fex, no3x, S_init, time_step, p, short_time): #, number_int_loops, sum_abs_cutoff):
    # run model to steady state with a given set of gene expression controls (beta)
    # and a given set of (fixed) environmental conditions (light, Mn, Fe, NO3)
    # returns the analytical expression for growth rate (Faizi et al 2018)
    # note that the negative of the growth rate is returned so that the minimize function can be
    # used for optimization
    # check that p is a list:

    if short_time:
        time_step_mod = time_step[0:500000]
    if not short_time:
        time_step_mod = time_step

 #   S = calc_ss(beta = np.abs(beta), I = I, mnx = mnx, fex = fex, no3x = no3x, S_init = S_init, time_step = time_step_mod,
 #               p = p, short_time = short_time, number_int_loops = number_int_loops, sum_abs_cutoff = sum_abs_cutoff)

    S = odeint(phyto_allocation,
               y0 = S_init,
               t = time_step_mod,
               args = (np.abs(beta),
                       I,
                       mnx,
                       fex,
                       no3x,
#                        p), mxstep = 10000000)
                       p), mxstep = 1000000)

#                        p))
    #tn_steady_state = S[0]
    #tn_steady_state = S[1]
    tn_steady_state = S[-1, 8] # indexing the last array (ie at steady state)
#     tn_steady_state = S.x[8] # indexing the last array (ie at steady state)
#     tn_steady_state = abs(S[8]) # indexing the last array (ie at steady state)

    e_steady_state = S[-1, 3]
#     e_steady_state = S.x[3]
#     e_steady_state = abs(S[3])
    mm_internal_energy = (e_steady_state)/(p[33] + e_steady_state)
    # uptake function for nutrients
    vn_ak = aksnes_cao_uptake(kcat_uptake = p[26],
                              n_transporters = tn_steady_state,
                              radius_c = p[0], trans_size = p[1],
                              diffusion_coef = 1.17e-8*60,
                              bulk_substrate = no3x)

    vn_ak_energy = vn_ak*mm_internal_energy
    # growth rate solved analytically
    growth_rate = (p[23]*vn_ak_energy)/(p[4]*(1 - p[45]))
    return(-np.sqrt(growth_rate))

@jit #(nopython=True)
def ode_opt(I, mnx, fex, no3x, time_step, p, x_init, S_init, short_time, cons, bnds, 
            #number_int_loops, sum_abs_cutoff, 
            ftol_var=1e-4):

    # run the optimization of protein allocation (gene expression beta's)
    # constraints are that each beta must be between 0 and 1
    # once optimized, then run the model to steady state given the optimized betas,
    # and determine steady state
    opt = scipy.optimize.minimize(optimize_growth_ak,
                                  x_init,
                                  args=(I,
                                        mnx,
                                        fex,
                                        no3x,
                                        S_init,
                                        time_step,
                                        p,
                                        short_time),
#                                        number_int_loops,
#                                        sum_abs_cutoff),
                                  method = 'SLSQP',
                                  constraints = cons,
                                  bounds = bnds,
                                  options={'maxiter':200,'ftol':ftol_var})
    optimized_beta = np.abs(opt.x)
    optimal_growth_rate = -opt.fun
    # once we've found the optimized beta, we need to determine steady state
    opt_ss = odeint(phyto_allocation,  
#     opt_ss = least_squares(fun = phyto_allocation,
#     opt_ss = fsolve(func = phyto_allocation,
                    y0 = S_init,
#                    x0 = S_init,
                    t = time_step, 
#                    bounds = (0, np.inf),
                    args = (np.abs(opt.x), 
                            I, 
                            mnx, 
                            fex, 
                            no3x, 
                            p), mxstep = 1000000)

#    opt_ss  = calc_ss(beta = np.abs(opt.x), I = I, mnx = mnx, fex = fex, no3x = no3x, S_init = S_init, time_step = time_step,
#                p = p, short_time = short_time, number_int_loops = number_int_loops, sum_abs_cutoff = sum_abs_cutoff, report_val = True)


#     opt_ss = opt_ss[0]
#     opt_ss_2 = [abs(i) for i in opt_ss] 
#    print(opt_ss[100,:])

    #print('opt growth rate', optimal_growth_rate)
    return([optimized_beta,
            optimal_growth_rate,
#             opt_ss.x])
            opt_ss[-1,:].tolist()]) # include the steady state solution


def cluster_init_starts(vector_of_betas, num_clust = 50):
    '''
    using sklearn kmeans to reduce the number of initial betas to the most unique starts
    '''
    kmeans_clusters = KMeans(n_clusters = num_clust).fit(vector_of_betas)
    kmeans_centroids = kmeans_clusters.cluster_centers_
    return(kmeans_centroids)

@jit
def make_random_starts(gpr_eval_number):
    '''
    make a np array of random starts
    '''
    random_starts_for_gpr = []
    for random_start_i in range(gpr_eval_number):
        random_start = np.random.random(6)
        random_start /= random_start.sum()
        random_starts_for_gpr.append(random_start)
    random_starts_arr = np.array(random_starts_for_gpr)
    return(random_starts_arr)

@jit
def get_top_growth_betas(random_starts_arr, gpr):
    '''
    get predicted values from gpr
    return the best betas to then run using the objective function
    '''
    # evaluate random starts using gpr
    gpr_growths = gpr.predict(random_starts_arr[:,:])
    # figure out how many runs to return
    how_many_gpr_top = int(round(0.2*len(random_starts_arr)))
    # find which indices have the highest growth
    top_growth_indices = gpr_growths.argsort()[-how_many_gpr_top:][::-1]
    # get the betas that correspond to this
    top_growth_random_starts = random_starts_arr[top_growth_indices]
    return(top_growth_random_starts)

@jit
def get_best_betas_overall(beta_training, growth_training, percentage_of_top = 0.1):
    '''
    function to return the best betas from the bayesian opt, to then feed into the kmeans clustering
    '''
    how_many_gpr_top = int(round(percentage_of_top*len(beta_training)))
    # find which indices have the highest growth
    top_growth_indices = growth_training.argsort()[-how_many_gpr_top:][::-1]
    # get the betas that correspond to this
    top_growth_betas = beta_training[top_growth_indices]
    return(top_growth_betas)


#@jit
def bayes_opt(total_loops, initial_betas, initial_growth, gpr_eval_number, I, mnx, fex, no3x, time_step, p, S_init, cons, bnds): #, number_int_loops, sum_abs_cutoff):
    '''
    bayesian optimization using gaussian process regression to do a guided search
    '''

    # convert the list of lists into np arrays
    initial_betas_arr = np.array(initial_betas)
    initial_growth_arr = np.array(initial_growth)
    # kernel for GPR model
    kernel = DotProduct() + WhiteKernel()

    print(initial_betas_arr.shape)
    print(initial_growth_arr.shape)

    # go through the bayesian opt
    for i in range(total_loops):

        #train the GPR
        print('re-training the GPR')
        gpr = GaussianProcessRegressor(kernel = kernel,
                                       n_restarts_optimizer = 10, alpha = 1e-4,
                                       random_state=0).fit(initial_betas_arr, initial_growth_arr)

        # making an array of random starts to be evaluated by the GPR
        random_starts_arr = make_random_starts(gpr_eval_number)

        # evaluate random starts using gpr
        best_betas_gpr = get_top_growth_betas(random_starts_arr, gpr)
        for gpr_beta in best_betas_gpr:
            try:
                #sys.stdout.flush()
                # run optimization (local)
                sys.stdout.flush()
                trial_i = ode_opt(I = I,
                                  mnx = mnx,
                                  fex = fex,
                                  no3x = no3x,
                                  time_step = time_step,
                                  p = p,
                                  S_init = S_init,
                                  cons = cons,
                                  bnds = bnds,
                                  x_init = gpr_beta,
                                  short_time = True)
                             #number_int_loops = number_int_loops,
                             #sum_abs_cutoff = sum_abs_cutoff)
            except:
                print('odeint failed for some reason')
                continue

            # making the trial_i the correct dimensions for appending
            zero_array = np.zeros((1, 6))
            zero_array[:,] = np.array(trial_i[0])

            # append the trial to the betas array, which will then be used to retrain the GPR
            initial_betas_arr = np.append(initial_betas_arr, zero_array, axis = 0)
            initial_growth_arr = np.append(initial_growth_arr, np.array(trial_i[1]))

            #print(initial_betas_arr.shape)
            #print(initial_growth_arr.shape)

    # writing out the GPR
#    gpr_filename = str(I) + str(fex) + str(mnx) + 'gpr_model.sav'
#    pickle.dump(gpr, open(gpr_filename, 'wb'))

    best_betas = get_best_betas_overall(initial_betas_arr, initial_growth_arr)

    return(best_betas)


def random_start(I, mnx, fex, no3x, time_step, p, S_init, number_trials, cons, bnds, total_loops, gpr_eval_number): # , #number_int_loops, sum_abs_cutoff):
    # function for doing a bunch of random start trials to better approximate global optimum
    # the second part of the function searches with a lower error tolerance to make the optimum found
    # more accurate

    number_of_kmeans_clusters = 30

    number_sub_trials = int(round(0.05*number_trials))
    if number_sub_trials < 1:
        number_sub_trials = 1
    if number_sub_trials < number_of_kmeans_clusters:
        number_of_kmeans_clusters = number_sub_trials
    all_trials = []
    growth_rates = []
    optimized_betas = []

    # initial random guesses
    for trial in range(number_trials):
        # select random values constrained to sum to 1
        print('drunkards walk ', trial)

        random_start = np.random.random(6)
        random_start /= random_start.sum()
        #print(random_start)

        try:
            # run optimization (local)
            trial_i = ode_opt(I = I,
                             mnx = mnx,
                             fex = fex,
                             no3x = no3x,
                             time_step = time_step,
                             p = p,
                             S_init = S_init,
                             cons = cons,
                             bnds = bnds,
                             x_init = np.array(random_start),
                             short_time = True) # , number_int_loops = number_int_loops, sum_abs_cutoff = sum_abs_cutoff)
        #    sys.stdout.flush()
        except:
            print('odeint failed for some reason')
            continue

#        print(trial_i)
        # add to list
        optimized_betas.append(trial_i[0])
        all_trials.append(trial_i)
        growth_rates.append(trial_i[1])

    print('finished initial random starts, now onto Bayesian opt')
    # bayesian optimization
    sys.stdout.flush()
    bayes_opt_estimates = bayes_opt(total_loops = total_loops,
                                    initial_betas = optimized_betas,
                                    initial_growth = growth_rates,
                                    gpr_eval_number = gpr_eval_number,
                                    I = I, mnx = mnx, fex = fex, no3x = no3x,
                                    time_step = time_step,
                                    p = p,
                                    S_init = S_init, cons = cons, bnds = bnds) # ,
                                    #number_int_loops = number_int_loops, sum_abs_cutoff = sum_abs_cutoff)


    # this part is potentially deprecated as the bayes_opt_estimates spits out the top trials anyways
    # find the top ten highest growth rates
    #highest_growth = sorted(range(len(growth_rates)), key=lambda i: growth_rates[i])[-number_sub_trials:]

    # looking at a subset of the best trials
    #best_trials = []
    #for trial in range(number_sub_trials):
    #    x_init_i = all_trials[highest_growth[trial]][0]
    #    best_trials.append(x_init_i)

    # clustering the betas and then using the most unique ones as an initial value
    clustered_betas = cluster_init_starts(vector_of_betas = bayes_opt_estimates, #best_trials,
                                          num_clust = number_of_kmeans_clusters)

    # rerun the optimization but with a lower error tolerance
    # taking the top ten best trials from the previous and running a more precise optimization
    low_tolerance_trials = []
    low_tolerance_growth = []
    low_tolerance_betas = []

    print(clustered_betas)
    sys.stdout.flush()
    for trial in clustered_betas:
        print(trial)
        try:
            sys.stdout.flush()
            precise_trial = ode_opt(I = I,
                                mnx = mnx,
                                fex = fex,
                                no3x = no3x,
                                time_step = time_step,
                                p = p,
                                cons = cons,
                                bnds = bnds,
                                S_init = S_init,
                                x_init = trial,
                                ftol_var = 1e-6,
                                short_time = False)
                            #number_int_loops = number_int_loops,
                            #sum_abs_cutoff = sum_abs_cutoff)
        except:
            print('ode int failed in precise')
            continue

        low_tolerance_trials.append(precise_trial)
        low_tolerance_growth.append(precise_trial[1])
        low_tolerance_betas.append(precise_trial[0])
        print('sober walk ', trial)

    # finding the best trial from the second round
    highest_growth_round2 = low_tolerance_growth.index(max(low_tolerance_growth))
    highest_growth_precise = low_tolerance_trials[highest_growth_round2]
    #print('highest growth:')
    #print(highest_growth_precise)
    return(highest_growth_precise)

def format_ode_opt(I, mnx, fex, no3x, time_step, p, x_init, S_init, number_trials, cons, bnds, total_loops, gpr_eval_number): # , number_int_loops, sum_abs_cutoff):

    if x_init == None:
        print('Initializing random starts for optimized parameters...')
        opt_answer = random_start(I = I, mnx = mnx, fex = fex,
                                 no3x = no3x, time_step = time_step, cons = cons, bnds = bnds,
                                 p = p, S_init = S_init, number_trials = number_trials, total_loops = total_loops,
                                 gpr_eval_number = gpr_eval_number) # , # number_int_loops = number_int_loops, sum_abs_cutoff = sum_abs_cutoff)
    else:
        print('Using designated parameter inputs as start value...')
        # run optimization and format the output nicely
        opt_answer = ode_opt(I = I, mnx = mnx, fex = fex,
                            no3x = no3x, time_step = time_step, p = p, cons = cons, bnds = bnds,
                            x_init = x_init, S_init = S_init)
        #print(opt_answer)
    # aggregate steady state solution
    output = [I, mnx, fex, no3x, opt_answer[1]] + np.abs(opt_answer[0]).tolist() + opt_answer[2]
    return(output)

def run_ode_opt_range(I, mnx, fex, no3x, time_step, p, x_init, S_init, cons, bnds, total_loops, gpr_eval_number,
                      #number_int_loops, sum_abs_cutoff, 
                      number_trials = 10):
    # run optimization across a range of gradients for each environmental variable

    # check that each input is a list
    if not type(I) == type(mnx) == type(fex) == type(no3x):
        ValueError('External input forcing values should be in list form, ie [1, 2, 3]')

    # there probably is a better way of doing this....
    # but this loops through each parameter list and then does the ode opt function above
    all_lists = []
    percentage_done = 0
    for light_val in range(len(I)):
        for mn_val in range(len(mnx)):
            for fe_val in range(len(fex)):
                for no3_val in range(len(no3x)):
                    print('Parameter set: light:', I[light_val],
                          'Mn: ', mnx[mn_val],
                          'Fe: ', fex[fe_val],
                          'NO3: ', no3x[no3_val])
                    sub_list = format_ode_opt(I = I[light_val],
                                             mnx = mnx[mn_val],
                                             fex = fex[fe_val],
                                             no3x = no3x[no3_val],
                                             time_step = time_step,
                                             p = p, x_init = x_init,
                                             cons = cons, bnds = bnds,
                                             S_init = S_init, total_loops = total_loops,
                                             gpr_eval_number = gpr_eval_number,
                                             #number_int_loops = number_int_loops)
                                             #sum_abs_cutoff = sum_abs_cutoff,
                                             number_trials = number_trials)
                    all_lists.append(sub_list)
#                    percentage_done += 1
#                    print('Percent done: ', percentage_done/expected_time*100)

    rel_derivs = model_mass_balance(general_model_output = all_lists, p = p)
    additional_rates = add_model_rates(general_model_output = all_lists, p = p)

    for par_combo in range(len(all_lists)):
        for rate_i in range(len(additional_rates[0])):
            all_lists[par_combo].append(additional_rates[par_combo][rate_i])
        for rel_deriv_i in range(len(rel_derivs[0])):
            all_lists[par_combo].append(rel_derivs[par_combo][rel_deriv_i])

    print(all_lists)

    # aggregating all results
    new_df = pd.DataFrame(columns=['I', 'Mnx', 'Fex', 'NO3x', 'u',
                                'beta_r', 'beta_tmn', 'beta_tfe',
                                'beta_tn',
                                'beta_a',
                                'beta_p',
                                'Mni', 'Fei', 'aa', 'e',
                                 'A',
                                 'P', 'Tmn',
                                'Tfe', 'Tn', 'R', # 'ind',
                                'gr',  'epsilon_p', 'alpha', 've', 'vros', 'omega', 'omega_cost', 'gamma_max_ros', 'gamma',
                                'dmni', 'dfei', 'daa', 'de', 'da', 'dp', 'dtmn', 'dtfe', 'dtn', 'dr'], data = all_lists)


    #return(all_lists)
    return(new_df)
