process COMBINE {
    input:
        path gold_dir
    output:
        path "data/all_combined.tsv"
    script:
        """
        python scripts/combine_gold_unadjusted.py --input-dir $gold_dir --output-file data/all_combined.tsv
        """
}

process SUBSET {
    input:
        path "data/all_combined.tsv"
    output:
        path "data/subsets/"
    script:
        """
        Rscript scripts/subset_combined.R
        """
}

process MODEL {
    input:
        path subset from SUBSET
    output:
        path "results/"
    script:
        """
        Rscript scripts/run_model_pipeline.R $subset
        """
}

workflow {
    COMBINE(gold_dir)
    SUBSET(COMBINE.out)
    MODEL(SUBSET.out)
}
