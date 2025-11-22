#!/usr/bin/env nextflow

/*
Nextflow workflow for breast cancer batch-effect analysis:
1. Combine gold unadjusted datasets
2. Subset combined data into 2–14 dataset groups
3. Apply batch effect adjusters (gmm, log_combat, min_mean, mnn, log_transformed)
4. Run classifier
5. Aggregate classifier outputs
6. Generate performance plots
*/

nextflow.enable.dsl=2

params.gold_dir = '/data/gold'
params.adjusted_dir = '/data/adjusted_datasets'
params.subset_dir = '/data/all_combined_subsets'
params.output_dir = '/outputs/classify_er_all'
params.figures_dir = '/outputs/classify_all_figures'
params.scripts_dir = '/scripts'

adjusters = ['gmm', 'log_combat', 'min_mean', 'mnn', 'log_transformed']

// Helper to run in apptainer
def run_in_apptainer(script_name, args='') {
    return "bash ~/confounded_analysis/run_in_apptainer.sh ${script_name} ${args}"
}
/************************************************
 * 1. COMBINE
 ************************************************/
process COMBINE {
    input:
        val gold_dir

    output:
        path "all_combined.csv", emit: combined_csv

    publishDir "/data/all_combined_data", mode: 'copy'

    script:
    def cmd = run_in_apptainer("/scripts/evaluations/classify_er_all_datasets/combine_all.py", "--input-dir ${gold_dir} --output-file all_combined.csv")
    """
    ${cmd}
    """

}

/************************************************
 * 2. ORDER
 ************************************************/
process ORDER {
    // input combined CSV file
    input:
        path combined_csv

    // output order files
    output:
        path "*_order.txt", emit: order_files

    script:
    """
    mkdir -p gse_order_files

    python /scripts/evaluations/classify_er_all_datasets/make_order_files.py \
        --input "\$combined_csv" \
        --output-dir gse_order_files

    cp gse_order_files/*_order.txt .
    """
}


/************************************************
 * 3. SUBSET
 ************************************************/
process SUBSET {

    input:
        tuple(path(order_file), val(k), val(test_source), path(combined_csv))

    output:
        path "${k}studies_test_${test_source}.csv"

    script:
    """
    Rscript /scripts/evaluations/classify_er_all_datasets/subset_prep.R \
        --input "\$combined_csv" \
        --order-file "\$order_file" \
        --k \$k \
        --test-source \$test_source \
        --output "${k}studies_test_${test_source}.csv"
    """
}


// /************************************************
//  * 3. ADJUST
//  ************************************************/
// process ADJUST {
//     input:
//         tuple val(adjuster), path(subset_csv)
//     output:
//         path "${adjuster}_subset*test_*.csv", emit: adjusted

//     script:
//     def subset_csv_str = subset_csv.toString()
//     """
//     basename_file=\$(basename $subset_csv_str)
//     subset_idx=\$(echo \$basename_file | sed -E 's/subset_([0-9]+)studies\\.csv/\\1/')

//     Rscript /scripts/evaluations/classify_er_all_datasets/run_scaling_experiment.R $adjuster \$subset_idx
//     """
// }

// /************************************************
//  * 4. CLASSIFIER
//  ************************************************/
// process CLASSIFIER {
//     container = '/apptainer/remove-batch-effects.sif'
//     publishDir params.output_dir, mode: 'copy'
//     input:
//         tuple val(adjuster), path(adjusted_csv)
//     output:
//         path "${adjuster}.csv"
//     script:
//     """
//     # Define the paths
//     CLASSIFIER_SCRIPT="/scripts/evaluations/classify_er_all_datasets/run_classifier.py"
//     OUTPUT_DIR="${params.output_dir}/${adjuster}"
//     mkdir -p \$OUTPUT_DIR

//     # Extract the test source from the filename
//     BASENAME=\$(basename \$adjusted_csv)
//     TEST_SOURCE=\$(echo \$BASENAME | sed -E "s/.*_test_(.*)\\.csv/\\1/")

//     # Run the python classifier
//     python \$CLASSIFIER_SCRIPT \$adjusted_csv \$OUTPUT_DIR

//     # Check if the metrics file was created
//     METRICS_FILE="\$OUTPUT_DIR/\${TEST_SOURCE}_metrics.csv"
//     if [ ! -f "\$METRICS_FILE" ]; then
//         echo "ERROR: No metrics file generated for \$adjusted_csv"
//         exit 1
//     fi

//     echo "Classifier finished for \$adjusted_csv, metrics saved to \$METRICS_FILE"
//     """
// }

// /************************************************
//  * 5. Aggregate per adjuster
//  ************************************************/
// process AGGREGATE_PER_ADJUSTER {
//     container = '/apptainer/remove-batch-effects.sif'
//     publishDir params.output_dir, mode: 'copy'

//     input:
//         tuple val(adjuster), path(adjuster_dir)

//     output:
//         path "${adjuster}.csv"

//     script:
//     """
//     python /scripts/aggregate_adjuster_results.py \
//         --input-dir $adjuster_dir \
//         --output ${adjuster}.csv
//     """
// }


// /************************************************
//  * 6. Aggregate ALL adjusters
//  ************************************************/
// process AGGREGATE_ALL {
//     container = '/apptainer/remove-batch-effects.sif'
//     publishDir params.output_dir, mode: 'copy'

//     input:
//         path csv_files
//         val output_file

//     output:
//         path output_file

//     script:
//     """
//     python /scripts/aggregate_all_adjusters.py \
//         --input-files ${csv_files.join(' ')} \
//         --output $output_file
//     """
// }


// /************************************************
//  * 7. Plot
//  ************************************************/
// process PLOT_PERFORMANCE {
//     container = '/apptainer/remove-batch-effects.sif'
//     publishDir params.figures_dir, mode: 'copy'

//     input:
//         path metrics_csv
//         val output_dir

//     output:
//         path "*"

//     script:
//     """
//     mkdir -p plots_tmp
//     Rscript /scripts/plot_classifier_performance.R \
//         --metrics $metrics_csv \
//         --output plots_tmp
//     cp plots_tmp/* .
//     """
// }


/************************************************
 * WORKFLOW DEFINITION
 ************************************************/
workflow {

    combined_ch = Channel.of(params.gold_dir) | COMBINE
    order_ch = combined_ch | ORDER

    subsets_ch = order_ch.flatMap { order_file ->
        (2..15).collect { k ->
            def test_source = file(order_file).getName().replaceFirst('_order\\.txt$', '')
            tuple(file(order_file), k, test_source)
        }
    }

    subsets_ch_with_csv = subsets_ch
        .cross(combined_ch)
        .map { subset_tuple, combined_csv ->
            def (order_file, k, test_source) = subset_tuple
            tuple(order_file, k, test_source, combined_csv)
        }

    subsets_ch_with_csv | SUBSET
}



// workflow {

//     /* 1–2 */
//     combined_ch = COMBINE(params.gold_dir)
//     subsets_ch  = SUBSET(combined_ch)

//     /* 3. Make 65 jobs (5 adjusters × 13 subset_csv files) */
//     adjust_pairs_ch = Channel
//         .from(adjusters)
//         .cross(subsets_ch)
//         .map { adj, subset -> tuple(adj, subset) }

//     adjusted_ch = ADJUST(adjust_pairs_ch)

//     /* 4. classification */
//     classified_ch = adjusted_ch
//         .map { file -> 
//             def adj = file.name.tokenize('_')[0]   // extract adjuster prefix
//             tuple(adj, file)
//         } \
//         | CLASSIFIER

//     /* 5. aggregate by adjuster */
//     adjuster_dirs_ch = Channel
//         .from(adjusters)
//         .map { adj ->
//             tuple(adj, file("${params.adjusted_dir}/${adj}"))
//         }

//     adjuster_dirs_ch | AGGREGATE_PER_ADJUSTER

//     /* 6 */
//     all_csvs = Channel.of(
//         file("${params.output_dir}/mnn.csv"),
//         file("${params.output_dir}/min_mean.csv"),
//         file("${params.output_dir}/log_combat.csv"),
//         file("${params.output_dir}/gmm.csv"),
//         file("${params.output_dir}/log_transformed.csv")
//     )

//     AGGREGATE_ALL(all_csvs, "${params.output_dir}/all_adjusters_metrics.csv")

//     /* 7 */
//     PLOT_PERFORMANCE(
//         file("${params.output_dir}/all_adjusters_metrics.csv"),
//         "${params.output_dir}/plots"
//     )
// }