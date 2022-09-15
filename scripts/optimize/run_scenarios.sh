#!/bin/bash

set -e

dataset=$1
unadjusted_file_name=$2
batch_col=$3
class_col=$4

echo "Performing optimizations for ${dataset}"

unadjusted_file=/data/${dataset}/${unadjusted_file_name}
results_file=/outputs/optimizations/${dataset}.tsv

tasks_file1=/tmp/confounded_tasks1.sh
tasks_file2=/tmp/confounded_tasks2.sh
scale=False

tmp_dir=/tmp/confounded
mkdir -p ${tmp_dir}/adjusted/${dataset} ${tmp_dir}/results/${dataset}
#rm -f ${tmp_dir}/adjusted/* ${tmp_dir}/results/*
rm -f ${tasks_file1} ${tasks_file2}

#for minibatch_size in 100 50 200
for minibatch_size in 100 50
do
    #for ds_layers in 10 5 20
    for ds_layers in 10 5
    do
        #for ae_layers in 2 5 10
        for ae_layers in 2 5
        do
            #for code_size in 20 50 10
            #for code_size in 20 50
            for code_size in 20
            do
                for scaling in linear sigmoid
                do
                    #for loss_weight in 1.0 1.5 2.0
                    for loss_weight in 1.0 2.0
                    do
                        #for minibatch_iterations in 10000 1000 20000
                        #for minibatch_iterations in 10000 1000
                        for minibatch_iterations in 10000
                        do
                            #for learning_rate in 0.0001 0.001 0.01
                            for learning_rate in 0.0001 0.001
                            do
                                tmp_adjusted_file=${tmp_dir}/adjusted/${dataset}/${minibatch_size}_${ds_layers}_${ae_layers}_${code_size}_${scaling}_${loss_weight}_${minibatch_iterations}_${learning_rate}.tsv
                                tmp_results_file=${tmp_dir}/results/${dataset}/${minibatch_size}_${ds_layers}_${ae_layers}_${code_size}_${scaling}_${loss_weight}_${minibatch_iterations}_${learning_rate}.tsv

                                #log_file=/outputs/optimizations/confounded_log.csv
                                #  -f ${log_file} \
                                # Leaving this out for now. The default is None.
                                #  -e ${early_stopping} \
                                # https://github.com/jdayton3/Confounded/blob/master/confounded/__main__.py

                                if [ ! -f ${tmp_adjusted_file} ]
                                then
                                    echo "confounded ${unadjusted_file} -o ${tmp_adjusted_file} -m ${minibatch_size} -l ${ds_layers} -a ${ae_layers} -b ${batch_col} -c ${code_size} -s ${scaling} -w ${loss_weight} -i ${minibatch_iterations} -g ${learning_rate}" >> ${tasks_file1}
                                fi

                                if [ ! -f ${tmp_results_file} ]
                                then
                                    echo "python /scripts/optimize/calculate_metrics.py -i ${tmp_adjusted_file} -o ${tmp_results_file} -b ${batch_col} -p ${class_col} -z ${minibatch_size},${ds_layers},${ae_layers},${code_size},${scaling},${loss_weight},${minibatch_iterations},${learning_rate} -s ${scale}" >> ${tasks_file2}
                                fi
                            done
                        done
                    done
                done
            done
        done
    done
done

if [ -f ${tasks_file1} ]
then
    parallel --jobs 12 --retries 0 --progress --eta -- < ${tasks_file1}
fi

if [ -f ${tasks_file2} ]
then
    parallel --jobs 20 --retries 0 --progress --eta -- < ${tasks_file2}
fi

if [ ! -f ${tmp_dir}/results/${dataset}/0 ]
then
    python /scripts/optimize/calculate_metrics.py -i ${unadjusted_file} -o ${tmp_dir}/results/${dataset}/0 -b ${batch_col} -p ${class_col} -y "minibatch_size,ds_layers,ae_layers,code_size,scaling,loss_weight,minibatch_iterations,learning_rate" -s ${scale}
fi
