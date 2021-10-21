iron_levels=(1 50 100 500 1000 2000 3000)
light_levels=(50)
par_string="../data/variable_parameters_mn_fe22_inference.csv"
model_string='mnfe4_model4_july14_post_inference_base'

DIR='../data/'

# multiple replicates of the same model

for iron_i in "${iron_levels[@]}"
do

for FILE in $par_string
do
        echo "Processing $FILE..."
        temp_string=${FILE#$'../data/'}
        echo $temp_string

for mang_i in "${iron_levels[@]}"
do
        echo $mang_i
#        echo '../data/model_output/'$model_string'_'$counter'_'$iron_i$mang_i$temp_string
        python mnfe4_model3.py -I ${light_levels[@]} -m $mang_i -f $iron_i -n 10000000 -p $FILE -o '../data/model_output/'$model_string'_5_'$iron_i$mang_i$temp_string --iterations 20 --ss_time 3000000 --total_loops 10 --gpr_eval_number 300 &
        python mnfe4_model3.py -I ${light_levels[@]} -m $mang_i -f $iron_i -n 10000000 -p $FILE -o '../data/model_output/'$model_string'_17_'$iron_i$mang_i$temp_string --iterations 20 --ss_time 3000000 --total_loops 10 --gpr_eval_number 300 &

        echo '../data/model_output/'$model_string'_'$i'_'$iron_i$mang_i$temp_string
done
done
wait
done
