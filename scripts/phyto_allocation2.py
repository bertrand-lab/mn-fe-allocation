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
# import matplotlib.pyplot as plt
import types
from numba import jit
import math
import time
import scipy.integrate as spi
import argparse

@jit(nopython=True)
def phyto_allocation(S, # internal species (all proteins, internal [Mn] etc.)
                     t, # time
                     beta, # optimized set of parameters for protein synthesis
                     I, # light
                     mnx, # external manganese
                     fex, # external iron
                     no3x, # external nitrate
                     p, #parameter inputs
                     rates = False): # run in checking mode
    ## Define species (initial conditions)
    Mni = S[0]
    Fei = S[1]
    aa = S[2]
    e = S[3]

    # Protein pools
    A = S[4]
    P = S[5]
    Tmn = S[6]
    Tfe = S[7]
    Tn = S[8]
    R = S[9]

    radius_cell = p[0] 
    transporter_size = p[1]
    N = p[2] #avocado constant[2]
    Ki = p[3]
    Mcell = p[4] # total cell mass
    nr = p[5] # ribosome lengthp
    ntmn_unadjusted = p[6] # length of mn transporter
    ntfe_unadjusted = p[7] # length of fe transporter
    nanti = p[8] # length of antioxidant enzyme complexnt
    ntn = p[9] # nitrate uptake enzyme and the precursor synthesis
    nphoto = p[10] # length of photosystem unit

    factor = p[11]  # conversion factor between seconds and min
    micro = p[12] # metric prefix

    # stoichiometric coefficients
    phi_mn_p = p[13] # amount of mn per photosynthetic unit
    phi_mn_a = p[14] # amount of mn per antioxidant enzyme

    phi_fe_p = p[15] # amount of fe per photosynthetic unit
    phi_fe_a = p[16] # amount of fe per antioxidant enzyme
    phi_fe_n = p[17] # amount for converting nitate to ammonium and incorporating into aas

    phi_e = p[18] # amount of energy [e] produced per photosynthetic cycle (equivalent to mphi in faizi model)

    phi_t_mn = p[19] # amount of energy consumed per uptake of Mn
    phi_t_fe = p[20] # amount of energy consumed per uptake of Fe
    phi_t_n = p[21] # amount of energy consumed by nitrate uptake and amino acid biosynthesis (ie. a lot!) (mmu in faizi model)

    mgamma = p[22] # molecules of e for one translational elongation step
    mn = p[23] # average number of aa per nitrogen (the inverse of average number of nitrogens in amino acids)

    k_cat_mn = p[24] # Mn uptake max turnover rate (s^-1)

    k_cat_fe = p[25] # Fe uptake max turnover rate (s^-1)

    k_cat_n = p[26] # NO3 uptake max turnover rate

    k_cat_ros = p[27] # max turnover rate of ROS metabolizing enzyme (min^-1)
    k_fei_p = p[28] # half saturation constant for iron uptake for protein synthesis
    k_mni_p = p[29] # half saturation constant for mn uptake for prot synthesis
    k_fei_a = p[30] # half saturation constant for iron uptake for antioxidant proten synthesis
    k_mni_a = p[31] # half sat constant for mn uptake for antioxidant prot synth
    k_fei_n = p[32] # half sat constant for fe uptake for nitrate uptake and prec. synth

    k_e = p[33] # half saturation constant for energy uptake
    epsilon_a = p[35] # antioxidant ros consumption
    gamma_max = p[36] # maximal translation elongation rate, multiplied by 60 to make it per min
    sigma = p[37] # effective cross-section of photosystem

    k2 = (p[38]) # turnover rate of photosystem
    kd = p[39] # photodamage rate

    subs_mn_anti = p[40]
    temperature = p[41]

    fe_uptake_rate_cost = p[42]
    mn_uptake_rate_cost = p[43]
    ros_sens = p[44]
    fixed_proteome_percentage = p[45]

    k_cat_fe1 = p[57]
    k_cat_fe2 = p[58]
    spec_percent_fe1 = p[59]
    spec_percent_fe2 = p[60]

    avail_space = p[64]

    # gene regulatory parameters
    beta_r = beta[0]
    beta_tmn = beta[1]
    beta_tfe = beta[2]
    beta_tn = beta[3]
    beta_a = beta[4]
    beta_p = beta[5]

    ## Define rate equations

    ### uptake rate equations using Aksnes and Cao 2011 uptake model
    vmn_no_e = aksnes_cao_uptake_comp(kcat_uptake = k_cat_mn,
                            n_transporters = Tmn,
                            radius_c = radius_cell,
                            diffusion_coef = 0.9e-9*60, # convert the molecular diffusion coefficient to m^2 min^-1
                            bulk_substrate = mnx,
                            trans_size = transporter_size,
                            other_n_transporters = [Tfe],
                            avail_space = avail_space)

#     vmn_no_e = aksnes_cao_uptake(kcat_uptake = k_cat_mn, n_transporters = Tmn,
#                            radius_c = radius_cell,
#                            diffusion_coef = 0.9e-9*60,# convert the molecular diffusion coefficient to m^2 min^-1
#                            bulk_substrate = mnx,
#                            trans_size = transporter_size)

    vmn = vmn_no_e*e/(k_e + e)

#     vfe_no_e = aksnes_cao_uptake(kcat_uptake = k_cat_fe,
#                            n_transporters = Tfe,
#                            radius_c = radius_cell,
#                            diffusion_coef = 0.9e-9*60,
#                            bulk_substrate = fex,
#                            trans_size = transporter_size)

    vfe_no_e = aksnes_cao_fe_speciation(percent_species_list = [spec_percent_fe1, spec_percent_fe2],
                            kcat_uptake_list = [k_cat_fe1, k_cat_fe2],
                            n_transporters = Tfe,
                            radius_c = radius_cell,
                            diffusion_coef = 0.9e-9*60,
                            bulk_substrate = fex,
                            trans_size = transporter_size,
                            other_n_transporters = [Tmn],
                            avail_space = avail_space)

    vfe = vfe_no_e*e/(k_e + e)
    # comp uptake model is not used because nitrate is saturating and thus membrane crowding is not important for
    # nitrate uptake
    vno3_no_e = aksnes_cao_uptake(kcat_uptake = k_cat_n,
                            n_transporters = Tn,
                            radius_c = radius_cell,
                            diffusion_coef = 1.17e-8*60,
                            bulk_substrate = no3x,
                            trans_size = transporter_size)
    vn = vno3_no_e*e/(k_e + e)

    total_fe = phi_fe_p*P + Fei + phi_fe_n*Tn

    gr = mn*vn/(Mcell*(1 - fixed_proteome_percentage))
    # dynamic electron leakiness from PSU
    #### linear dynamic leakiness from Twining et al 2004 SOFEX cruise
    if total_fe <= 7173653:
        epsilon_p = 0.3
    if total_fe >= 39129014:
        epsilon_p = 0.05
    if total_fe > 7173653 and total_fe < 39129014:
        epsilon_p = 3.561e-1 - 7.823e-9*total_fe

    # multiply by 60 to get to uE * m^-2 * min^-1
    alpha = (sigma*I*60)/(sigma*I*60 + k2 + gr)
    ve = (1 - epsilon_p)*k2*alpha*P
    vros = k_cat_ros*A

    # adjusting the protein cost of iron uptake machinery - at higher uptake rates, higher cost
    ntfe = ntfe_unadjusted + (vfe/gr)*fe_uptake_rate_cost
    ntmn = ntmn_unadjusted + (vmn/gr)*mn_uptake_rate_cost

    # ROS cost function
    omega = (epsilon_p*ve - epsilon_a*vros)/(epsilon_p*ve + epsilon_a*vros)

    if omega < 0:
        omega = 0

    omega_cost = 2*ros_sens**(-1*omega)/(ros_sens**(-2*omega) + 1)

    ######################### Kept here for testing purposes only
#    omega_cost = 1
    #########################

    # function for adjusting the rate of protein synthesis as a function of temperature
    gamma_max_temp = gamma_max*2**((temperature - 20)/10)
    gamma_max_ros = gamma_max_temp*omega_cost
    gamma = R*gamma_max_ros*(e/(k_e + e))*(aa/(k_e + aa))

    gamma_Tmn = beta_tmn*(gamma/ntmn)

    # substitutable antioxidants ## not used in paper but kept for future work
    if subs_mn_anti == 1:
        subs_mm_mni = Mni/(k_mni_a + Mni + Fei)
        subs_mm_fei = Fei/(k_fei_a + Mni + Fei)
        gamma_A = beta_a*(gamma/nanti)*(subs_mm_mni + subs_mm_fei)
        mni_anti_coef = subs_mm_mni/(subs_mm_mni + subs_mm_fei)
        fei_anti_coef = 1 - subs_mm_mni/(subs_mm_mni + subs_mm_fei)
    else:
        gamma_A = beta_a*(gamma/nanti)*Mni/(k_mni_a + Mni)

    gamma_Tfe = beta_tfe*(gamma/ntfe)
    gamma_Tn = beta_tn*(gamma/ntn)*Fei/(k_fei_n + Fei)
    gamma_R = beta_r*(gamma/nr)
    gamma_P = beta_p*(gamma/nphoto)*Fei/(k_fei_p + Fei)*Mni/(k_mni_p + Mni)

    if subs_mn_anti == 1:
        ## fraction of substitutable resources going to antioxidants
        mni_anti_coef = subs_mm_mni/(subs_mm_mni + subs_mm_fei)
        fei_anti_coef = 1 - subs_mm_mni/(subs_mm_mni + subs_mm_fei)

        ## Define differential equations - additional terms for a sink of Fe/Mn is required
        dmnidt = vmn - phi_mn_p*gamma_P - gr*Mni - mni_anti_coef*phi_mn_a*gamma_A
        dfeidt = vfe - phi_fe_p*gamma_P - gr*Fei - fei_anti_coef*phi_fe_a*gamma_A - phi_fe_n*gamma_Tn
    else:
        dmnidt = vmn - phi_mn_p*gamma_P - gr*Mni - phi_mn_a*gamma_A
        dfeidt = vfe - phi_fe_p*gamma_P - gr*Fei - phi_fe_n*gamma_Tn

    dedt = phi_e*ve - phi_t_mn*vmn - phi_t_fe*vfe - mgamma*(gamma_Tmn*ntmn +
                                                           gamma_Tfe*ntfe +
                                                           gamma_Tn*ntn +
                                                           gamma_R*nr +
                                                           gamma_A*nanti +
                                                           gamma_P*nphoto) - gr*e - phi_t_n*vn
    daadt = mn*vn - (gamma_Tmn*ntmn +
                  gamma_Tfe*ntfe +
                  gamma_Tn*ntn +
                  gamma_R*nr +
                  gamma_A*nanti +
                  gamma_P*nphoto) - gr*aa
    dadt = gamma_A - gr*A
    dpdt = gamma_P - gr*P
    dtmndt = gamma_Tmn - gr*Tmn
    dtfedt = gamma_Tfe - gr*Tfe
    dtndt = gamma_Tn - gr*Tn
    drdt = gamma_R - gr*R

    if not rates:
        return_list = [dmnidt, dfeidt, daadt, dedt,
                       dadt,
                       dpdt, dtmndt, dtfedt, dtndt, drdt]
    if rates:
        return_list = [gr, epsilon_p, alpha, ve, vros, omega, omega_cost, gamma_max_ros, gamma]

    return(return_list)

@jit
def aksnes_cao_uptake(kcat_uptake, n_transporters, radius_c, diffusion_coef, bulk_substrate, trans_size):
    '''
    Function for calculating the uptake rate of a nutrient using the aksnes cao model
    '''
    #### nutrient uptake function from Aksnes and Cao (2011) extended from
    #### Aksnes and Egge (1991) and examined in Fiksen et al (2013)
    ### kcat = maximum turnover rate, converted to handling time, in minutes
    ### n = unitless, number of uptake sites
    ### r = cell radius, m
    ### D = molecular diffusion coefficient (m^2 * s^-1)
    ### bulk substrate concentration, S_infinity, in nanomol * Litre
    ###### Note that: the input has to be picomole per litre, which is then converted to molecules per m^3
    ### s = uptake site radius (m)

    bulk_substrate_converted = bulk_substrate*1e3*1e-12*6.0221409e+23
    h = 1/kcat_uptake # handling time is equivalent to 1/kcat
    # the fraction of surface area covered with uptake sites, p
    p_dens = (n_transporters*math.pi*(trans_size**2))/(4*math.pi*(radius_c**2))
    # uptake affinity, alpha_inft
    alpha_inft = 4*math.pi*diffusion_coef*radius_c*((n_transporters*trans_size)/(n_transporters*trans_size + math.pi*radius_c*(1 - p_dens)))
    ## uptake term including handling time, h, diffusion layer, and variable
    ## number of transporter sites (Aksnes and Cao 2011, Fiksen et al 2013)
    c1 = h/(4*n_transporters*math.pi*radius_c*diffusion_coef*bulk_substrate_converted)
    c2 = (1 - (math.pi*radius_c*p_dens)/(n_transporters*trans_size))

    c = c1*c2
    b = (1/(alpha_inft*bulk_substrate_converted)) + (h/n_transporters)
    v1 = b/(2*c)
    v2 = (1 - math.sqrt(1 - (4*c)/(b**2)))
    v_aksnes = v1*v2

    return(v_aksnes)


################# NOTE:

## These two functions below didn't compile for some reason with Numba, so I've added them to below.

@jit
def p_dens_func(n_transporters, trans_size, radius_c):
    '''
    Helper function to calculate the *p* term in the Aksnes Cao model
    '''
    p_dens = (n_transporters*math.pi*(trans_size**2))/(4*math.pi*(radius_c**2))
    return(p_dens)

@jit
def p_dens_func_interfere(n_transporters, trans_size, radius_c, avail_space):
    '''
    Helper function to calculate the proportion of available diffusive flux based on the 
    proportion of area transporters take up vs. available area on the sphere surface
    '''
    p_dens = (n_transporters*math.pi*(trans_size**2))/(avail_space*4*math.pi*(radius_c**2))
    return(p_dens)


@jit
def total_interference(other_n_transporters, trans_size, radius_c, avail_space):
    '''
    Calculating the total amount of interference at the membrane surface
    '''

    # input must be a list for other_n_transporters
    sum_val = 0
    # goes through all other transporters and calculates the available diffusion flux to the cell
    for other_trans in range(len(other_n_transporters)):
        p_dens_out = p_dens_func_interfere(n_transporters = other_n_transporters[other_trans],
                                           trans_size = trans_size,
                                           radius_c = radius_c,
                                           avail_space = avail_space)
        sum_val += p_dens_out

    # required because the correction factor can be negative, giving strange results
    if sum_val > 1:
        sum_val = 0.99999999
    correction_term = 1 - sum_val

    return(correction_term)


@jit
def aksnes_cao_uptake_comp(kcat_uptake, n_transporters, radius_c,
                           diffusion_coef, bulk_substrate, trans_size,
                           other_n_transporters, avail_space):
    '''
    Uptake rates with membrane competition. A new term is added that modifies the available space 
    for uptake.
    '''
    #### nutrient uptake function from Aksnes and Cao (2011) extended from
    #### Aksnes and Egge (1991) and examined in Fiksen et al (2013)
    ### kcat = maximum turnover rate, converted to handling time, in minutes
    ### n = unitless, number of uptake sites
    ### r = cell radius, m
    ### D = molecular diffusion coefficient (m^2 * s^-1)
    ### bulk substrate concentration, S_infinity, in nanomol * Litre
    ###### Note that: the input has to be picomole per litre, which is then converted to molecules per m^3
    ### s = uptake site radius (m)

    bulk_substrate_converted = bulk_substrate*1e3*1e-12*6.0221409e+23
    h = 1/kcat_uptake # handling time is equivalent to 1/kcat
    # the fraction of surface area covered with uptake sites, p
    p_dens = p_dens_func(n_transporters = n_transporters,
                        trans_size = trans_size,
                        radius_c = radius_c)

    #p_dens = (n_transporters*math.pi*(trans_size**2))/(4*math.pi*(radius_c**2))
    # uptake affinity, alpha_inft
    alpha_mod = total_interference(other_n_transporters = other_n_transporters,
                                  trans_size = trans_size,
                                  radius_c = radius_c,
                                  avail_space = avail_space)
    alpha_inft_no_mod = 4*math.pi*diffusion_coef*radius_c*((n_transporters*trans_size)/(n_transporters*trans_size + alpha_mod*avail_space*math.pi*radius_c*(1 - p_dens)))

#    summed_vals = (other_n_transporters*math.pi*(trans_size**2))/(4*math.pi*(radius_c**2))
#    if summed_vals > 4*math.pi*(radius_c**2):
#        summed_vals = 0.9999999
#    alpha_mod = 1 - summed_vals

    # adjusting alpha by the modification term calculated above
    alpha_inft = avail_space*alpha_inft_no_mod*alpha_mod
    ## uptake term including handling time, h, diffusion layer, and variable
    ## number of transporter sites (Aksnes and Cao 2011, Fiksen et al 2013)
    c1 = h/(4*n_transporters*math.pi*radius_c*diffusion_coef*bulk_substrate_converted)
    c2 = (1 - (math.pi*radius_c*p_dens)/(n_transporters*trans_size))

    c = c1*c2
    b = (1/(alpha_inft*bulk_substrate_converted)) + (h/n_transporters)
    v1 = b/(2*c)
    v2 = (1 - np.sqrt(1 - (4*c)/(b**2)))

    v_aksnes = v1*v2
    return(v_aksnes)


@jit
def aksnes_cao_fe_speciation(percent_species_list, kcat_uptake_list, n_transporters, radius_c,
                           diffusion_coef, bulk_substrate, trans_size,
                           other_n_transporters, avail_space):
    '''
    function for including a range of kcats for different fe species
    '''
    fe_rates_total = 0
    # check that the fe rate constants and the fe speciation percentages are the same
#    assert len(percent_species_list) == len(kcat_uptake_list)
    # go through each fe species and calculate uptake rates
    for fe_species in range(len(percent_species_list)):
        fe_species_rates = aksnes_cao_uptake_comp(kcat_uptake = kcat_uptake_list[fe_species],
                               n_transporters = n_transporters,
                               radius_c = radius_c,
                               diffusion_coef = diffusion_coef,
                               bulk_substrate = bulk_substrate*percent_species_list[fe_species],
                               trans_size = trans_size,
                               other_n_transporters = other_n_transporters,
                               avail_space = avail_space)
        fe_rates_total += fe_species_rates
    return(fe_rates_total)


def aksnes_cao_uptake_no_numba(kcat_uptake, n_transporters, radius_c, diffusion_coef, bulk_substrate, trans_size, affinity_return):
    '''
    This function is the same as above, but is used for the output file to 
    calculate the uptake rates. Because that output uses a pandas df, it can't work with numba, 
    so this function is needed.
    '''

    #### nutrient uptake function from Aksnes and Cao (2011) extended from
    #### Aksnes and Egge (1991) and examined in Fiksen et al (2013)
    ### kcat = maximum turnover rate, converted to handling time, in minutes
    ### n = unitless, number of uptake sites
    ### r = cell radius, m
    ### D = molecular diffusion coefficient (m^2 * s^-1)
    ### bulk substrate concentration, S_infinity, in nanomol * Litre
    ###### Note that: the input has to be picomole per litre, which is then converted to molecules per m^3
    ### s = uptake site radius (m)
    bulk_substrate_converted = bulk_substrate*1e3*1e-12*6.0221409e+23
    h = 1/kcat_uptake # handling time is equivalent to 1/kcat
    # the fraction of surface area covered with uptake sites, p
    p_dens = (n_transporters*math.pi*(trans_size**2))/(4*math.pi*(radius_c**2))
    alpha_inft = 4*math.pi*diffusion_coef*radius_c*((n_transporters*trans_size)/(n_transporters*trans_size + math.pi*radius_c*(1 - p_dens)))
    ## uptake term including handling time, h, diffusion layer, and variable
    ## number of transporter sites (Aksnes and Cao 2011, Fiksen et al 2013)
    c1 = h/(4*n_transporters*math.pi*radius_c*diffusion_coef*bulk_substrate_converted)
    c2 = (1 - (math.pi*radius_c*p_dens)/(n_transporters*trans_size))
    c = c1*c2
    b = (1/(alpha_inft*bulk_substrate_converted)) + (h/n_transporters)
    v1 = b/(2*c)
    v2 = (1 - np.sqrt(1 - (4*c)/(b**2)))

    v_aksnes = v1*v2

    if affinity_return == False:
        returned_value = v_aksnes
    if affinity_return == True:
        returned_value = alpha_inft
    return(returned_value)

def range_of_uptakes(kcat_uptake, n_transporters, radius_c, 
                           diffusion_coef, bulk_substrate, trans_size,
                           other_n_transporters_list):

    '''
    This function is for plotting purposes. It takes a list of transporter amounts,
    'other_n_transporters_list', and then calculates the uptake rate with only the other transporters
    changing.
    '''
    if not isinstance(other_n_transporters_list, list):
        ValueError('input must be a list for this')
        
    list_of_uptakes = []
    for sub_val in range(len(other_n_transporters_list)):
        uptake_val = aksnes_cao_uptake_comp(kcat_uptake = kcat_uptake, 
                               n_transporters = n_transporters, 
                               radius_c = radius_c, 
                           diffusion_coef = diffusion_coef, 
                               bulk_substrate = bulk_substrate, 
                               trans_size = trans_size,
                           other_n_transporters = other_n_transporters_list[sub_val])
        list_of_uptakes.append(uptake_val)
        
    return(list_of_uptakes)
