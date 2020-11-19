import numpy as np
import random
from numba import jit

def set_histogram_points(parameter_vector, number_of_breaks, 
                         lower_bound = None, 
                         upper_bound = None, 
                         manual_over = True):
    if not manual_over:
        t_seq = np.linspace(start = min(parameter_vector), 
                            stop = max(parameter_vector), 
                            num = number_of_breaks)
    elif manual_over:
        t_seq = np.linspace(start = lower_bound,
                            stop = upper_bound,
                            num = number_of_breaks)
    return(t_seq)

@jit(nopython=True)
def calculate_weights(vector_of_sum_sq, h_val):
    weight_vector = np.exp(-(vector_of_sum_sq**2)/(2*h_val**2)) 
    return(weight_vector)

@jit(nopython=True)
def calculate_expectation_at_intervals(number_of_expect_iter, 
                                       parameter_vector,
                                       t_seq,
                                       calculated_weights):
    '''
    this function calculates the *cumulative* expectation at fixed intervals of a histogram 
    '''

    # initialize an empty dataframe to collect the results in
    par_i_to_B = parameter_vector

    collected_expectations = np.zeros((number_of_expect_iter, len(t_seq)))

    # we also want to collect the weights (1 or 0 vectors) so that we can do a 
    # posterior predictive check afterwards
#     weight_number_rows = number_of_expect_iter*len(t_seq)
#    collected_weights = np.zeros((weight_number_rows, len(par_i_to_B)))

#    weight_counter = 0
    for overall_rep in range(number_of_expect_iter):
#         print(overall_rep)
        # getting the time start for one rep
        t_i_expectation = np.zeros(len(t_seq))
        # we start from the second iteration of t_seq, because below the first one there are no observations
        for i in range(1, len(t_seq)):
            t_i = t_seq[i]
             # probabilistically include parameter less than t_i
            para_accepted_weights = np.zeros(len(par_i_to_B))

            for par_i in range(0, len(par_i_to_B)):
                generated_u = random.uniform(0, 1)
                if generated_u < calculated_weights[par_i]:
                    para_accepted_weights[par_i] = 1
                else:
                    para_accepted_weights[par_i] = 0

            if np.sum(para_accepted_weights) == 0:
                para_accepted_vec_summary = 0
            else:
                w_i_factor = para_accepted_weights/np.sum(para_accepted_weights)
                #### THIS IS WHERE TO LOOP OVER ALL T seqs
                g_i_indicator_f = par_i_to_B <= t_i
                para_accepted_vec_summary = np.sum(g_i_indicator_f*w_i_factor)

            t_i_expectation[i] = para_accepted_vec_summary
#            collected_weights[weight_counter] = para_accepted_weights
            # update the weight counter
#            weight_counter += 1

        collected_expectations[overall_rep] = t_i_expectation
 #   print(collected_weights)
    return(collected_expectations)

def format_expectation_to_histogram(expectation_all_output,
                                    t_seq):
    '''
    formats the expectation to be in a histogram
    '''

    # getting the mean expectation at interal t
    zero_t_i_expectation = np.mean(expectation_all_output, axis = 0)
    
    # OLD LOOPS, NEW CODE BELOW and all from numpy direclty. ugh. read documentation before
    # re-inventing the wheel!
    # transforming the expectation from cumulative to regular prob histogram
#     transformed_expectation = np.zeros(len(zero_t_i_expectation) - 1)
#     transformed_t_seq = np.zeros(len(zero_t_i_expectation) - 1)

#     for k in range(1, len(zero_t_i_expectation) - 1):
#         trans_exp_k = zero_t_i_expectation[k] - zero_t_i_expectation[k - 1]
#         if trans_exp_k < 0:
#             trans_exp_k = 0
#         transformed_expectation[k] = trans_exp_k

#     for i in range(1, len(t_seq) - 1):
#         transformed_t_seq[i] = (t_seq[i] + t_seq[i - 1])/2
    # transforming the histogram break points to get their centroid values
    transformed_t_seq = t_seq[:-1] + np.diff(t_seq)/2
    # transforming expectation from cumulative to regular prob histogram
    transformed_expectation = np.diff(zero_t_i_expectation)
    # removing values that have a negative probability, which are errors in approximation of low prob's
    transformed_expectation[transformed_expectation < 0] = 0

    exp_result = np.array([transformed_expectation, transformed_t_seq])

    return(exp_result)

def run_full_prob_abc(parameter_choice, data_set, 
                      number_of_breaks,
                      iterations_out, h_val, 
                      lower_bound = None, 
                      upper_bound = None, 
                      manual_over = True):

    par_vector = data_set[parameter_choice]
    vector_sum_sq_out = data_set['sum_sq_dif']

    t_seq_out = set_histogram_points(parameter_vector = par_vector,
                                     number_of_breaks = number_of_breaks,
                                     lower_bound = lower_bound,
                                     upper_bound = upper_bound,
                                     manual_over = manual_over)

    calc_weights_out = calculate_weights(vector_of_sum_sq = vector_sum_sq_out,
                                         h_val = h_val)

    expectation_all_out = calculate_expectation_at_intervals(number_of_expect_iter = iterations_out,
                                       parameter_vector = par_vector,
                                       t_seq = t_seq_out,
                                       calculated_weights = calc_weights_out)

    formatted_out = format_expectation_to_histogram(expectation_all_output = expectation_all_out,
                                                    t_seq = t_seq_out)

#    print(expectation_all_out[1])
    return([formatted_out])


def calculate_expectation_at_multipar(number_of_expect_iter,
                                       parameter_vector_list,
                                       t_seq_list,
                                       calculated_weights):
    '''
    this function calculates the *cumulative* expectation at fixed intervals of a histogram
    '''

    # initialize an empty dataframe to collect the results in
    par_i_to_B_1 = parameter_vector[0]
    par_i_to_B_2 = parameter_vector[1]
    par_i_to_B_3 = parameter_vector[2]

    collected_expectations = np.zeros((number_of_expect_iter, len(t_seq)))

    for overall_rep in range(number_of_expect_iter):
        print(overall_rep)
        # getting the time start for one rep
        t_i_expectation_1 = np.zeros(len(t_seq))
        t_i_expectation_2 = np.zeros(len(t_seq))
        t_i_expectation_3 = np.zeros(len(t_seq))
        # we start from the second iteration of t_seq, because below the first one there are no observations

        for i in range(1, len(t_seq)):
            t_i = t_seq[i]
             # probabilistically include parameter less than t_i
            para_accepted_weights = np.zeros(len(par_i_to_B))

            for par_i in range(0, len(par_i_to_B)):
                generated_u = random.uniform(0, 1)
                if generated_u < calculated_weights[par_i]:
                    para_accepted_weights[par_i] = 1
                else:
                    para_accepted_weights[par_i] = 0

            if np.sum(para_accepted_weights) == 0:
                para_accepted_vec_summary = 0
            else:
                w_i_factor = para_accepted_weights/np.sum(para_accepted_weights)
                g_i_indicator_f = np.where(par_i_to_B <= t_i,
                                           1,
                                           0)
                para_accepted_vec_summary = np.sum(g_i_indicator_f*w_i_factor)
            t_i_expectation[i] = para_accepted_vec_summary

        collected_expectations[overall_rep] = t_i_expectation
    return(collected_expectations)

