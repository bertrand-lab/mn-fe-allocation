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
from abc_sampling import add_columns_formatting_abc, gen_par_sets, mod_par_sets, gen_and_mod_par_set, generate_and_run_model

# parsing all the command line Arguments
parser = argparse.ArgumentParser()

parser.add_argument("-n", "--name", required=True, dest = 'name_in_out')

args = parser.parse_args()

name_in_out = args.name_in_out

# run the model
# reading in the parameter file
paras_meta = np.genfromtxt('../data/variable_parameters_mn_fe20.csv', 
                      delimiter=',', skip_header=1)
paras_nunn = np.genfromtxt('../data/variable_parameters_mn_fe20_Nunn_base.csv', 
                      delimiter=',', skip_header=1)
paras_cohen = np.genfromtxt('../data/variable_parameters_mn_fe20_Cohen_base.csv',
                      delimiter=',', skip_header=1)

# setting initial conditions, which are *very* rough guesses from initial model output. Helps with starting
# with somewhat realistic model starts for the optimization.
S_init = np.zeros(10)
S_init[0] = 5000 # Mni
S_init[1] = 5000 # Fei
S_init[2] = 5000000000
S_init[3] = 8000000
S_init[4] = 10000
S_init[5] = 10000
S_init[6] = 500
S_init[7] = 500
S_init[8] = 10000
S_init[9] = 10000

# setting up time step
t = np.linspace(0, int(1000000), int(1000000)/10)

room_left = 1

#if paras[48] == 1:# if ferritin_uptake == 1 then change the fixed proteome amount
#    room_left = 1 - paras[51]# ferritin fixed proteome

#if paras[48] == 0:# if ferritin_uptake == 0 (ie vacuolar uptake) then change the fixed proteome to the regular fixed proteome
#    room_left = 1 - paras[45]

# create some constraints
cons = ({'type':'ineq','fun': lambda x: 1 - x.sum()},
#cons = ({'type':'ineq','fun': lambda x: room_left - x.sum()}, # this is here for testing, the results should be the same
         {'type':'ineq','fun': lambda x: int(all([i >= 0.0 for i in x]))})

bnds = [(0.0, room_left),
        (0.0, room_left),
        (0.0, room_left),
        (0.0, room_left),
        (0.0, room_left),
        (0.0, room_left)]

generate_and_run_model(number_abc_loops = 10000, name_for_out = name_in_out, paras_meta = paras_meta,
                       paras_nunn = paras_nunn, paras_cohen = paras_cohen,
                       t = t, S_init = S_init, cons = cons, bnds = bnds,
                       total_loops = 2,
                       gpr_eval_number = 1000,
                       number_trials = 10)


