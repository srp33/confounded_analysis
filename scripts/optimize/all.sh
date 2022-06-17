#!/bin/bash

set -e

printf "\033[0;32mPerforming optimizations\033[0m\n"

# Generate many different versions of the simulated data

batch_col=Batch

results_file=/outputs/optimizations/simulated_expression.tsv

#for minibatch_size in 100 50 200
for minibatch_size in 100
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
                        #for iterations in 10000 1000 20000
                        for iterations in 10000
                        do
                            #for learning_rate in 0.0001 0.001 0.01
                            for learning_rate in 0.0001
                            do
                                adjusted_file=/data/simulated_expression/optimizations/confounded.csv
                                #log_file=/outputs/optimizations/confounded_log.csv
                                #  -f ${log_file} \

                                # Leaving this out for now. The default is None.
                                #  -e ${early_stopping} \

                                # https://github.com/jdayton3/Confounded/blob/master/confounded/__main__.py
                                #confounded /data/simulated_expression/unadjusted.csv \
                                #  -o ${adjusted_file} \
                                #  -m ${minibatch_size} \
                                #  -l ${layers} \
                                #  -a ${ae_layers} \
                                #  -b ${batch_col} \
                                #  -c ${code_size} \
                                #  -s ${scaling} \
                                #  -w ${loss_weight} \
                                #  -i ${iterations} \
                                #  -g ${learning_rate}

                                python /scripts/optimize/calculate_metrics.py -i ${adjusted_file} -o ${results_file} -b Batch -p Class -y "minibatch_size,ds_layers,ae_layers,code_size,scaling,loss_weight,iterations,learning_rate" -z "${minibatch_size},${ds_layers},${ae_layers},${code_size},${scaling},${loss_weight},${iterations},${learning_rate}"
#/scripts/metrics/single_metric.sh
                            done
                        done
                    done
                done
            done
        done
    done
done

# Calculate metrics for every combination of the simulated data.

#out_file_path="/outputs/optimizations/classification_results.csv"
#
#rm -f ${out_file_path}
#
#script_path="$(dirname $0)/classify.py"
#
#python "${script_path}" -i /data/simulated_expression -o ${batch_out_path} -c Batch 
#python "${script_path}" -i /data/simulated_expression -o ${true_out_path} -c Class
