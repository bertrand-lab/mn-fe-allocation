iron_levels=(78 3000)
mang_levels=(2 3000)
light_levels=(50)
par_string="../data/model_growth_parameterizations/variable_parameters_mn_fe22_inference*"
model_string='mnfe4_model4_july14_post_inference_base'

DIR='../data/'

# multiple replicates of the same model
for i in {1..10}
do
        echo "$i"

for iron_i in "${iron_levels[@]}"
do
#       echo $FILE
#       echo $iron_i
#       echo ${iron_levels[@]}

#       python mnfe3_model7.py -I ${light_levels[@]} -m $iron_i -f ${iron_levels[@]} -n 10000000 -p $FILE -o '../data/model_output/'$model_string$iron_i$temp_string --iterations 3500 --ss_time 100000000 &

for FILE in $par_string
do
        echo "Processing $FILE..."
        temp_string=${FILE#$'../data/model_growth_parameterizations/'}
        echo $temp_string

for mang_i in "${mang_levels[@]}"
do
        echo $mang_i
#        echo '../data/model_output/'$model_string'_'$counter'_'$iron_i$mang_i$temp_string
        python mnfe4_model3.py -I ${light_levels[@]} -m $mang_i -f $iron_i -n 10000000 -p $FILE -o '../data/model_output/'$model_string'_'$i'_'$iron_i$mang_i$temp_string --iterations 20 --ss_time 3000000 --total_loops 10 --gpr_eval_number 300 &
        echo '../data/model_output/'$model_string'_'$i'_'$iron_i$mang_i$temp_string
done
done
wait
done
done
