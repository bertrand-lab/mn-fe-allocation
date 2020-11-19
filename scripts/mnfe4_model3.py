## Mn/Fe proteomic allocation model
# Original code was written by M. Faizi (from Faizi et al 2018)
# Model and code were changed by Scott McCain (January- 2019)

from __future__ import division

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
import argparse

from datetime import datetime
start = datetime.now()

### import the allocation model functions
from phyto_allocation2 import *
from multi_step_opt2 import *
from checking_functions2 import *

# parsing all the command line Arguments
parser = argparse.ArgumentParser()

parser.add_argument("-I", "--light", nargs='+', required=True, dest = 'light_levels',
                    help="light level in umol Einsteins * m^-2 *s^-1")
parser.add_argument("-m", "--mnx", nargs='+', required=True, dest = 'mn_levels',
                    help="external manganese concentration in pM")
parser.add_argument("-f", "--iron", nargs='+', required=True, dest = 'fe_levels',
                    help="external iron concentration in pM")
parser.add_argument("-n", "--nitrate", nargs='+', required=True, dest = 'no3_levels',
                    help="external nitrate concentration in pM")
parser.add_argument("--iterations", required = True, dest = 'random_trials', help = "an integer designating how many random trials to intialize optimiations")
parser.add_argument("--ss_time", required = True, dest = 'ss_time', help = "maximum time step for calculating steady state growth rates")
parser.add_argument("-p", "--parameters", dest = "parameter_file",
                   help="parameter file for cellular model")
parser.add_argument("--total_loops", required = True, dest = 'total_loops', help = "an integer designating how many times to train the Gaussian Process Regression")
parser.add_argument("--gpr_eval_number", required = True, dest = 'gpr_eval_number', help = "an integer designating how many samples from the GPR should be taken")
parser.add_argument("-o", "--outputname", dest = "output_file_name",
                   help = "model output file name")
#parser.add_argument("--number_int_loops", dest = "number_int_loops",
#                   help = "number of integration loops")
#parser.add_argument("--sum_abs_cutoff", dest = "sum_abs_cutoff",
#                   help = "sum of absolute derivatives cutoff")

args = parser.parse_args()

# collecting the inputs
light_levels = map(int, args.light_levels)
mn_levels = map(int, args.mn_levels)
fe_levels = map(int, args.fe_levels)
no3_levels = map(int, args.no3_levels)
random_trials = args.random_trials
ss_time = args.ss_time
parameter_file = args.parameter_file
output_file_name = args.output_file_name
total_loops = args.total_loops
grp_eval_number = args.gpr_eval_number

#number_int_loops = args.number_int_loops
#sum_abs_cutoff = args.sum_abs_cutoff

# checking parameter inputs are correct
if output_file_name[-4:] != '.csv':
    raise NameError('output file name should be a .csv file type silly goose! :)')
if parameter_file[-4:] != '.csv':
    raise NameError('parameter file name should be a .csv file type silly goose! :)')

# reading in the parameter file
paras = np.genfromtxt(parameter_file, delimiter=',', skip_header=1)

print(paras)

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
#S_init[10] = 36092832672

# setting up time step
t = np.linspace(0, int(ss_time), int(ss_time)/10)

#room_left = 1 - paras[51] # ferritin fixed proteome

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

# run entire model
model_output = run_ode_opt_range(I = light_levels, mnx = mn_levels, fex = fe_levels,
                                 no3x = no3_levels, time_step = t, p = paras,
                                 x_init = None, S_init = S_init, cons = cons, bnds = bnds, total_loops = int(total_loops),
                                 gpr_eval_number = int(grp_eval_number),
                                 number_trials = int(random_trials)) #  , number_int_loops = int(number_int_loops),
                                 #sum_abs_cutoff = float(sum_abs_cutoff))

add_model_output = add_columns_formatting(model_output_general = model_output, p = paras)

print(add_model_output)

add_model_output.to_csv(output_file_name)
#model_output.to_csv(output_file_name)
difference = datetime.now() - start

print('total time', difference)
