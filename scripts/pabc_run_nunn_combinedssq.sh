declare -a G_RATES=("u_trans_87" "u_trans_2876")

# loop through each  value of uptake rates
for g_rate_i in "${G_RATES[@]}"
do
   python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 40 -s 2 -t 400000 --lower_bound 0 --upper_bound 3.0 -p $g_rate_i -o $g_rate_i'_meta_combined_meta_cohen_nunn_h2' &
   python run_pabc_command.py -i "../data/abc_intermediate/combined_sum_sq_meta_cohen_nunn.csv" -b 40 -s 2 -t 25000 --lower_bound 0 --upper_bound 3.0 -p $g_rate_i -o $g_rate_i'_meta_combined_meta_cohen_nunn_h2_25000' &
done

