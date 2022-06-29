#!/bin/bash

set -e

echo "Performing optimizations"

# Generate many different versions of the simulated data

batch_col=Batch

unadjusted_file=/data/simulated_expression/unadjusted.csv
results_file=/outputs/optimizations/simulated_expression.tsv

tmp_dir=/tmp/confounded
mkdir -p ${tmp_dir}/adjusted ${tmp_dir}/results
rm -f ${tmp_dir}/adjusted/* ${tmp_dir}/results/*

python /scripts/optimize/calculate_metrics.py -i ${unadjusted_file} -o ${tmp_dir}/results/0 -b Batch -p Class -y "minibatch_size,ds_layers,ae_layers,code_size,scaling,loss_weight,minibatch_iterations,learning_rate"

for minibatch_size in 100 50 200
#for minibatch_size in 100
do
    #for ds_layers in 10 5 20
    for ds_layers in 10
    do
        #for ae_layers in 2 5 10
        for ae_layers in 2
        do
            #for code_size in 20 50 10
            for code_size in 50
            do
                #for scaling in linear sigmoid
                for scaling in sigmoid
                do
                    #for loss_weight in 1.0 1.5 2.0
                    for loss_weight in 1.0
                    do
                        #for minibatch_iterations in 10000 1000 20000
                        for minibatch_iterations in 10000
                        do
                            #for learning_rate in 0.0001 0.001 0.01
                            for learning_rate in 0.0001
                            do
                                tmp_adjusted_file=${tmp_dir}/adjusted/${minibatch_size}_${ds_layers}_${ae_layers}_${code_size}_${scaling}_${loss_weight}_${minibatch_iterations}_${learning_rate}.tsv
                                tmp_results_file=${tmp_dir}/results/${minibatch_size}_${ds_layers}_${ae_layers}_${code_size}_${scaling}_${loss_weight}_${minibatch_iterations}_${learning_rate}.tsv

                                #log_file=/outputs/optimizations/confounded_log.csv
                                #  -f ${log_file} \
                                # Leaving this out for now. The default is None.
                                #  -e ${early_stopping} \

                                # https://github.com/jdayton3/Confounded/blob/master/confounded/__main__.py
                                confounded ${unadjusted_file} \
                                  -o ${tmp_adjusted_file} \
                                  -m ${minibatch_size} \
                                  -l ${ds_layers} \
                                  -a ${ae_layers} \
                                  -b ${batch_col} \
                                  -c ${code_size} \
                                  -s ${scaling} \
                                  -w ${loss_weight} \
                                  -i ${minibatch_iterations} \
                                  -g ${learning_rate}

                                python /scripts/optimize/calculate_metrics.py -i ${tmp_adjusted_file} -o ${tmp_results_file} -b Batch -p Class -z "${minibatch_size},${ds_layers},${ae_layers},${code_size},${scaling},${loss_weight},${minibatch_iterations},${learning_rate}"
                            done
                        done
                    done
                done
            done
        done
    done
done

python /scripts/optimize/cat_results.py "${tmp_dir}/results/*" ${results_file}

echo Results are saved in ${results_file}
#/scripts/metrics/single_metric.sh
