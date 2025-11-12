#!/usr/bin/env nextflow

/*
Nextflow workflow for breast cancer batch-effect analysis:
1. Combine gold unadjusted datasets
2. Subset combined data into 2–14 dataset groups
3. Apply batch effect adjusters (gmm, log_combat, min_mean, mnn) + unadjusted
4. Run classifier
5. Aggregate classifier outputs
6. Generate performance plots
*/

nextflow.enable.dls = 1   // safe DLS mode

process {
    executor = 'slurm'
    queue = 'normal'
    cpus = 8
    memory = '16 GB'
    time = '10h'

    container = '~/confounded_analysis/apptainer/remove-batch-effects.sif'
    singularity.enabled = true

    scratch = true
    stageInMode = 'copy'
}

params.gold_dir = '/data/gold'
params.all_combined = '/data/all_combined_data/all_combined.csv'
params.adjusted_dir = '/data/adjusted_datasets'
params.output_dir = '/outputs/classify_er_all'

workflow {

    // 1️⃣ Combine gold datasets
    combined_ch = COMBINE(params.gold_dir)

    // 2️⃣ Subset combined data
    subsets_ch = SUBSET(combined_ch)

    // 3️⃣ Run classifiers for each adjuster
    adjusters = ['gmm', 'log_combat', 'min_mean', 'mnn']
    adjuster_results_ch = Channel.empty()

    adjusters.each { adjuster ->
        csv_files_ch = Channel.fromPath("${params.adjusted_dir}/${adjuster}/*.csv")
        classifier_ch = csv_files_ch.map { file -> tuple(file, adjuster) } | CLASSIFIER
        AGGREGATE_PER_ADJUSTER(classifier_ch, adjuster)
        adjuster_results_ch = adjuster_results_ch.merge(classifier_ch)
    }

    // 4️⃣ Run classifier on unadjusted datasets
    unadjusted_ch = Channel.fromPath("${params.gold_dir}/*/unadjusted.csv")
        .map { file -> tuple(file, 'unadjusted') } | RUN_CLASSIFIER_UNADJUSTED
    AGGREGATE_PER_ADJUSTER(unadjusted_ch, 'unadjusted')

    // 5️⃣ Aggregate all adjusters into a single CSV
    all_csvs = Channel.fromPath([
        "${params.output_dir}/gmm.csv",
        "${params.output_dir}/log_combat.csv",
        "${params.output_dir}/min_mean.csv",
        "${params.output_dir}/mnn.csv",
        "${params.output_dir}/unadjusted.csv"
    ])
    AGGREGATE_ALL(all_csvs, "${params.output_dir}/all_adjusters_metrics.csv")

    // 6️⃣ Plot performance
    PLOT_PERFORMANCE("${params.output_dir}/all_adjusters_metrics.csv", "${params.output_dir}/plots")
}


// ========================================================
// PROCESS DEFINITIONS
// ========================================================

process COMBINE {
    input:
        path gold_dir
    output:
        path "data/all_combined_data/all_combined.csv", emit: combined_csv
    script:
        """
        mkdir -p data/all_combined_data
        python scripts/combine_gold_unadjusted.py \
            --input-dir $gold_dir \
            --output-file data/all_combined_data/all_combined.csv
        """
}

process SUBSET {
    input:
        path combined_csv
    output:
        path "data/adjusted_datasets/*", emit: adjusted_subsets
    script:
        """
        mkdir -p data/adjusted_datasets
        Rscript scripts/evaluations/classify_er_all_datasets/subset_combined_data.R \
            --input $combined_csv \
            --output data/adjusted_datasets
        """
}

process CLASSIFIER {
    input:
        tuple path(csv_file), val(adjuster)
    output:
        path "${params.output_dir}/${adjuster}.csv", emit: classifier_output
    script:
        """
        mkdir -p ${params.output_dir}
        bash scripts/run_in_apptainer.sh $csv_file
        """
}

process RUN_CLASSIFIER_UNADJUSTED {
    input:
        tuple path(csv_file), val(adjuster)
    output:
        path "${params.output_dir}/unadjusted.csv", emit: unadjusted_output
    script:
        """
        mkdir -p ${params.output_dir}
        bash scripts/run_in_apptainer.sh $csv_file
        """
}

process AGGREGATE_PER_ADJUSTER {
    input:
        path results
        val adjuster
    output:
        path "${params.output_dir}/${adjuster}.csv"
    script:
        """
        python scripts/aggregate_adjuster_results.py \
            --input-dir $results \
            --output ${params.output_dir}/${adjuster}.csv
        """
}

process AGGREGATE_ALL {
    input:
        path csv_files
        val output_file
    output:
        path output_file
    script:
        """
        python scripts/aggregate_all_adjusters.py \
            --input-files ${csv_files.join(' ')} \
            --output $output_file
        """
}

process PLOT_PERFORMANCE {
    input:
        path metrics_csv
        val output_dir
    output:
        path "${output_dir}/*"
    script:
        """
        mkdir -p ${output_dir}
        Rscript scripts/plot_classifier_performance.R \
            --metrics $metrics_csv \
            --output ${output_dir}
        """
}
