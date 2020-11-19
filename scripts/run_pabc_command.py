import numpy as np
import random
from numba import jit
import pandas as pd
import argparse

from pabc_functions import *

parser = argparse.ArgumentParser(description='Run ABC.')

parser.add_argument("-b", nargs='+', required=True, dest = 'b_input',
                    help="number of breaks")
parser.add_argument("-s", nargs='+', required=True, dest = 'h_input',
                    help="h value")
parser.add_argument("-t", nargs='+', required=True, dest = 'i_input',
                    help="iterations")
parser.add_argument("--lower_bound", nargs = "+", required=True, dest = "lower_bound",
                    help = "lower bound of making the histogram bins")
parser.add_argument("--upper_bound", nargs = "+", required = True, dest = "upper_bound",
                    help = "upper bound of making the histogram bins")

parser.add_argument("-p", dest = "parameter_chosen",
                   help = "which parameter to get posterior")
parser.add_argument("-o", "--outputname", dest = "output_file_name",
                   help = "model output file name")
parser.add_argument("-i", "--inputname", dest = "input_file_name",
                   help = "abc input file name")

args = parser.parse_args()


if args.input_file_name[-4:] != '.csv':
    raise NameError('input file name should be a .csv file type silly goose! :)')

# import data
my_data = np.genfromtxt(args.input_file_name, delimiter=',',
                       skip_header = False, names = True)

par_out_cost = run_full_prob_abc(parameter_choice = args.parameter_chosen,
                  data_set = my_data,
                  number_of_breaks = int(args.b_input[0]),
                  iterations_out = int(args.i_input[0]),
                  h_val = float(args.h_input[0]),
                  lower_bound = float(args.lower_bound[0]),
                  upper_bound = float(args.upper_bound[0]))

#out_file_name_long = args.output_file_name + 'all_samples.csv'
out_file_name_short = args.output_file_name + '_summarized.csv'

pd.DataFrame(np.transpose(par_out_cost[0])).to_csv(out_file_name_short, header = ['probability', args.parameter_chosen])
#pd.DataFrame(np.transpose(par_out_cost[0])).to_csv(out_file_name_long)

