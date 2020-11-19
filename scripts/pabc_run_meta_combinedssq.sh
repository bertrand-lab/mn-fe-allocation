python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 70 -s 2 -t 400000 --lower_bound 0 --upper_bound 0.1 -p 'epsilon_a' -o 'epsilon_a_meta_combined_meta_cohen_nunn_h2' &
python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 70 -s 2 -t 400000 --lower_bound 0 --upper_bound 16 -p 'cost_par' -o 'cost_par_meta_combined_meta_cohen_nunn_h2' &
python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 70 -s 2 -t 400000 --lower_bound 0 --upper_bound 0.15 -p 'avail_space' -o 'avail_space_meta_combined_meta_cohen_nunn_h2' &


# Declare all the variables from the model to get the posteriors
declare -a PARS=("A" "P" "Tfe" "Tn" "R" "Tmn")

## now loop through the above array to get each posterior histogram
for par_i in "${PARS[@]}"
do
   python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 40 -s 2 -t 400000 --lower_bound 0 --upper_bound 3 -p $par_i -o $par_i'_meta_combined_meta_cohen_nunn_h2' &
done

wait

declare -a UP_RATES=("mn_uptake_470" "mn_uptake_1010" "fe_uptake_470" "fe_uptake_1010")

# loop through each  value of uptake rates
for up_rate_i in "${UP_RATES[@]}"
do
   python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 30 -s 2 -t 400000 --lower_bound 0 --upper_bound 12000 -p $up_rate_i -o $up_rate_i'_meta_combined_meta_cohen_nunn_h2' &
done

wait
declare -a QUOTAS=("total_mn_amol_470" "total_mn_amol_1010" "total_fe_amol_470" "total_fe_amol_1010")

# loop through each  value of uptake rates
for quota_i in "${QUOTAS[@]}"
do
   python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 40 -s 2 -t 400000 --lower_bound 0 --upper_bound 150 -p $quota_i -o $quota_i'_meta_combined_meta_cohen_nunn_h2' &
done

wait
declare -a G_RATES=("u_trans_470" "u_trans_1010")

# loop through each  value of uptake rates
for g_rate_i in "${G_RATES[@]}"
do
   python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 40 -s 2 -t 400000 --lower_bound 0 --upper_bound 0.75 -p $g_rate_i -o $g_rate_i'_meta_combined_meta_cohen_nunn_h2' &
done

