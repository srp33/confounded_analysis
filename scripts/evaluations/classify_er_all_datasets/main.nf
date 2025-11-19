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
    adjusters = ['gmm', 'log_combat', 'min_mean', 'mnn', 'log_transformed']

    adjuster_dirs_ch = Channel.from(adjusters).map { adj ->
        tuple(adj, file("${params.adjusted_dir}/${adj}"))
    }
    
    adjusters.each { adjuster ->
        csv_files_ch = Channel.fromPath("${params.adjusted_dir}/${adjuster}/*.csv")
        classifier_ch = csv_files_ch.map { file -> tuple(file, adjuster) } | CLASSIFIER
    }

    adjuster_dirs_ch | AGGREGATE_PER_ADJUSTER

    // 5️⃣ Aggregate all adjusters into a single CSV
    all_csvs = Channel.of(
        file("${params.output_dir}/gmm.csv"),
        file("${params.output_dir}/log_combat.csv"),
        file("${params.output_dir}/min_mean.csv"),
        file("${params.output_dir}/mnn.csv"),
        file("${params.output_dir}/log_transformed.csv")
    )
    AGGREGATE_ALL(all_csvs, "${params.output_dir}/all_adjusters_metrics.csv")

    // 6️⃣ Plot performance
    PLOT_PERFORMANCE(file("${params.output_dir}/all_adjusters_metrics.csv"), "${params.output_dir}/plots")
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
        path combined_csv  // ignored by your script but required for chaining

    output:
        path "data/all_combined_subsets/subset_*.csv", emit: adjusted_subsets

    script:
        """
        mkdir -p data/all_combined_subsets

        # Loop over target subset sizes (2–14)
        for k in \$(seq 2 14); do
            out_file="data/all_combined_subsets/subset_\${k}studies.csv"
            if [ ! -f "\$out_file" ]; then
                echo "Creating missing subset file: \$out_file"
                Rscript /scripts/evaluations/classify_er_all_datasets/subset_combined_data.R
            else
                echo "Skipping existing subset file: \$out_file"
            fi
        done
        """
}


process CLASSIFIER {
    publishDir params.output_dir, mode: 'copyNoOverwrite'
    input:
        tuple path(csv_file), val(adjuster)
    when: 
        !file("${params.output_dir}/${adjuster}.csv").exists()
    output:
        path "${adjuster}.csv"
    script:
        """
        bash scripts/run_in_apptainer.sh $csv_file > ${adjuster}.csv
        """
}

process AGGREGATE_PER_ADJUSTER {
    publishDir params.output_dir, mode: 'copyNoOverwrite'

    input:
        tuple val(adjuster), path(adjuster_dir)
    when: 
        !file("${params.output_dir}/${adjuster}.csv").exists()
    output:
        path "${adjuster}.csv"
    script:
        """
        python scripts/aggregate_adjuster_results.py \
            --input-dir $adjuster_dir \
            --output ${adjuster}.csv
        """
}

process AGGREGATE_ALL {
    publishDir params.output_dir, mode: 'copyNoOverwrite'
    input:
        path csv_files
        val output_file
    when:
        !file(output_file).exists()
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
    publishDir params.output_dir + "/plots", mode: 'copyNoOverwrite'
    input:
        path metrics_csv
        val output_dir
    when: 
        !file("${output_dir}/performance_overview.png").exists()
    output:
        path "*"
    script:
        """
        tmpdir="plots_tmp"
        mkdir -p $tmpdir

        Rscript scripts/plot_classifier_performance.R \
            --metrics $metrics_csv \
            --output $tmpdir

        cp \$tmpdir/* .
        """
}
