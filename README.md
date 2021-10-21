This page describes how to run all the code for "Cellular costs underpin micronutrient limitation in phytoplankton", by J.S.P. McCain, A. Tagliabue, E. Susko, A.E. Achterberg, A.E. Allen, and E.M. Bertrand, published in Science Advances. [Here](https://www.science.org/doi/10.1126/sciadv.abg6501) is the paper.

`/scripts` contains all the files for running the model and conducting the Approximate Bayesian Computation for estimating model parameters.

To re-run the code, you need the following file structure:
```
\data
    \model_output
    \abc_out
    \abc_intermediate
    \model_growth_parameterizations
    \kleiner_data
    \oceanographic_data
    \culture_proteomes
    \culture_data
\scripts
\figures
```

Also all of the data are deposited [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.xd2547dfs) so you can run any stage of the script.

### Parameter Estimation

To recreate all analyses, the first stage is parameter estimation, which requires a lot of time and computational power. We used Compute Canada and ran the following script:

You first need to create an environment for running on the Compute Canada cluster:

```
environment setup

mkdir mn-fe-allocation
cd mn-fe-allocation
mkdir scripts
mkdir data
cd data
mkdir model_output
mkdir abc_out
cd ..

module load scipy-stack
module load python/2.7

# setting up virtualenv

virtualenv --no-download ~/pyteo_27
source ~/pyteo_27/bin/activate
pip install --no-index --upgrade pip
pip install -Iv numba==0.37.0
```

Then you need to run the generation of models with parameters sampled from your prior distribution (warning this step takes a while -- ~1-2 months using 500-2000 cores). If you don't want to run the inference step and just play around with the model parameters, go straight to "Running Base Model" section below.
```
mnfe4_model3_abc2_cedar_sbatch.sh

```

This batch script accesses two python scripts: `abc_sampling_mnfe_setup.py`, which sets up the ABC sampling, and `abc_sampling.py`, which is the central part of the generation of model runs.

These scripts also require the main components of the model, which are in the following scripts: `phyto_allocation2.py` (main model equations); `multi_step_opt2.py` (optimization methods and ODE integration); and `checking_functions2.py`.

You also need three parameter files in \data, `variable_parameters_mn_fe20.csv`; `variable_parameters_mn_fe20_Nunn_base.csv`; and `variable_parameters_mn_fe20_Cohen_base.csv`.

### Completing Parameter Estimation with ABC

This step requires some interactive work. You first need to run the following R scripts to generate a file that has the combined Euclidean distances for a given parameter vector:

```
nunn_abc_relative.R
meta_abc_relative.R
cohen_abc_relative.R
```

Which are then all combined using: `merging_abc_results.R`.

Once that is all run, you run: `pabc_run_meta_combinedssq.sh` and `pabc_run_nunn_combinedssq.sh` to get the posterior distributions. These scripts call on two scripts: `run_pabc_command.py` (which conducts the ABC analysis using a command line function) and `pabc_functions.py` (which is the heart of the ABC analysis).

### Running Base Model

To run the base model from the command line, you use the python script `mnfe4_model3.py`, which has command-line inputs for environmental parameters and outputs a csv file. To run the model as shown in the manuscript you can use the bash script `mnfe4_model4_post_inference.sh`, which loops through multiple runs of the model under various environmental conditions.

This requires the parameter file `variable_parameters_mn_fe22_inference.csv`. (which is the same as the above parameter file, except it uses the posterior modes from the ABC analysis).

To run the model experiments (from the manuscript), you need to run: `mnfe4_model4_post_inference_growth_exp.sh`; `mnfe4_model4_post_inference_growth_exp_baseline.sh`; `mnfe4_model4_post_inference_light_experiments.sh`

### Figures and Smaller Analyses

Once you have all the model output, one script (to rule them all!) runs all the plots and smaller analyses, mostly conducted in R.

```
all_figures_script.R
```

### GEOTRACES Data Analysis

Calculating Mixed Layer Depths from the GEOTRACES intermediate data product, and then figuring out the mixed layer concentrations of Mn and Fe is done using: `geotraces_and_par.py`.To get all the data for that, you need to download data from the 2017 IDP, as well as data from the Ocean Colour Database, which can be done using: `get_par_atten_coef.sh` and then some naming formatting using: `name_nc_files.sh` 

### Kleiner et al (2017) Nature Comms Analysis

To re-run the analysis using this awesome dataset, you first need to download the data from their paper. (We have used the same naming to help with consistency). Then you map peptides to taxa uniquely using: `processing_kleiner_data.py` and plot the data using `plotting_kleiner_data.R`

### Miscellaneous Analysis Scripts

Nunn_dFe_calculations.R

- Script written by Loay Jabre to calculate Fe' (Jabre and Bertrand, 2020, L&O). Specifically, we wanted to calculate the amount of Fe' from Nunn et al 2013 (PLOS ONE).

calculating_sample_sd_nunn.R

- The ABC method requires an input for the expected standard error. But for the standard error of the ratio of two time points, it's not clear which replicates you should pair with each other. To get past this, we sampled all possible combinations of pairs. That is what this script does.

formatting_culture_proteomes.R

- old version of culture proteome formatting

nunn2013_proteome_vals.R

- old script examining Nunn et al 2013 protein expression vals

plotting_ros_equation.R

- script that creates one of the supplementary figures to illustrate the reactive oxygen species compartment of the model.
