import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from argparse import ArgumentParser
import os
import sys
import itertools
from pathlib import Path

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(project_root)

from adjust.invert_autoclass import BatchCorrectImpute, take_norm
from metrics.classify import score_adjusters
from metrics.util import DataFrameCache


def split_and_normalize(df, batch_column):
    """
    Splits the DataFrame into gene expression data and batch information,
    normalizes the gene expression data if necessary, and returns the normalized data and batch information.
    """
    # Check for negative values. Assuming that microarray data is mean 0, while RNA-seq data is counts
    # Make sure to filter out the metadata columns if they exist. Metadata columns start with 'meta_'
    meta_cols = [col for col in df.columns if col.startswith('meta_')]
    print(f"Found metadata columns: {meta_cols}")
    gene_cols = [col for col in df.columns if col not in meta_cols]

    genes = df[gene_cols]
    meta_df = df[meta_cols]
    normalized = genes.min().min() < 0

    # Normalize and logp1-transform the data if necessary
    if not normalized:
        print("Data is not normalized. Normalizing and log1p-transforming the data.")
        X_norm = take_norm(genes)

    batches = df[batch_column].values
    return genes, batches, meta_df, gene_cols, normalized


def save(result, gene_cols, meta_cols, output_file):
    """
    Adjusts the gene expression data for batch effects and saves the result to a CSV file.
    """

    # Convert the result to a DataFrame
    result = pd.DataFrame(result, columns=gene_cols)

    # Add the metadata columns back to the result
    if meta_cols:
        result = pd.concat([df[meta_cols], result], axis=1)

    # Save the imputed data to a CSV file
    try:
        result.to_csv(output_file, index=False)
        print(f"Imputed data saved to {output_file}")
    except Exception as e:
        raise ValueError(f"AutoClass Error saving the output file: {output_file}\n{e}")





def score_result(result, prediction_column, input_dir, cache):
    """
    Score the imputed data using the specified prediction column and input directory.
    Returns the ROC AUC score for the imputed data.
    """
    dataset = os.path.basename(input_dir)
    key = ("HistGradientBoostingClassifier", "autoclass", dataset, prediction_column)

    cache.set_dataframe(result, file_path=os.path.join(input_dir, f"autoclass.csv"))

    scores = score_adjusters(
        input_dir=input_dir,
        prediction_column=prediction_column,
        adjusters=["autoclass"],
        learner_names=["HistGradientBoostingClassifier"],
        write_over=True,
        cache=cache
    ).popitem()[1]["roc_auc_score"]
    return np.array(scores)

def get_score_stats(df, column, input_dir, cache):
    score = score_result(df, column, input_dir, cache)
    mean = score.mean(axis=0)
    std = score.std(axis=0)
    return mean, std


def stats_from_run(df, batch_column, true_column, input_dir, cache):
    batch_mean, batch_std = get_score_stats(df, batch_column, input_dir, cache)
    true_mean, true_std = get_score_stats(df, true_column, input_dir, cache)
    return batch_mean, batch_std, true_mean, true_std


def experiment(genes, batches, meta_df, hyperparams, batch_column, true_column, input_dir, cache):
    """
    Run the AutoClass experiment with the given parameters and return the results.
    """
    print(f"Running AutoClass with hyperparameters: {hyperparams}")
    res = run(genes, batches, hyperparams)
    res = pd.DataFrame(res, columns=genes.columns)
    res_with_meta = pd.concat([meta_df, res], axis=1)

    # Score the result
    batch_mean, batch_std, true_mean, true_std = stats_from_run(
        res_with_meta, batch_column, true_column, input_dir, cache
    )
    
    print(f"True ROC AUC score mean : {true_mean}, std: {true_std}")
    print(f"Batch ROC AUC score mean : {batch_mean}, std: {batch_std}")

    return batch_mean, batch_std, true_mean, true_std


def run(genes, batches, hyperparams):
    """
    Perform batch correction and imputation using the BatchCorrectImpute class.
    Returns the imputed gene expression data.
    """
    res = BatchCorrectImpute(
        genes,
        batches,
        cellwise_norm=hyperparams.get('cellwise_norm', False),
        log1p=hyperparams.get('log1p', False),
        verbose=hyperparams.get('verbose', True),
        encoder_layer_size=hyperparams.get('encoder_layer_size', [128]),
        dropout_rate=hyperparams.get('dropout_rate', 0.1),
        epochs=hyperparams.get('epochs', 300),
        adversarial_weight=hyperparams.get('adversarial_weight', 0.6),
        reg=hyperparams.get('reg', 0.000),
        batch_size=hyperparams.get('batch_size', 32),
        es=hyperparams.get('es', 30),
        lr=hyperparams.get('lr', 15)
    )
    new_genes = res['imp']
    return new_genes


def main(input_dir, output_file, batch_column, true_column):
    # Load the dataset --------------------------------
    input_file = Path(input_dir) / "unadjusted.csv"

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        raise ValueError(f"AutoClass Error reading the input file: {input_file}\n{e}")

    genes, batches, meta_df, gene_cols, normalized = split_and_normalize(df, batch_column)

    cache = DataFrameCache(folder=input_dir)
    cache.set_dataframe(df, file_path=input_file)

    dataset = os.path.basename(input_dir)

    adv_weights = [0.002, 0.02, 0.005]
    encoder_layer_sizes = [[128], [256]]
    epochs = [400]
    dropout_rates = [0.1, 0.2, 0.5]
    regularizations = [0.0001, 0.001, 0.000, 0.1]
    lr = [5, 10, 15]

    default_hyperparams = {
        'cellwise_norm': False,
        'log1p': False,
        'adversarial_weight': 0.6,
        'encoder_layer_size': [128],
        'dropout_rate': 0.1,
        'epochs': 300,
        'reg': 0.000,
        'batch_size': 32,
        'verbose': False,
        'es': 30,
        'lr': 15
    }

    # with open(output_file, 'w') as f:
    #     f.write("encoder_layer_size,adversarial_weight,epochs,lr,reg,dropout_rate,"
    #             "batch_mean,std_batch,true_mean,std_true\n")

    hyperparam_combinations = itertools.product(regularizations, lr, epochs, dropout_rates, adv_weights,
                                                encoder_layer_sizes)

    for reg, learning_rate, epochs, dropout_rate, adversarial_weight, encoder_layer_size in hyperparam_combinations:
        hyperparams = default_hyperparams.copy()
        hyperparams['dropout_rate'] = dropout_rate
        hyperparams['reg'] = reg
        hyperparams['lr'] = learning_rate
        hyperparams['adversarial_weight'] = adversarial_weight
        hyperparams['epochs'] = epochs
        hyperparams['encoder_layer_size'] =  encoder_layer_size

        batch_mean, batch_std, true_mean, true_std = experiment(
            genes, batches, meta_df, hyperparams, batch_column, true_column, input_dir, cache
        )
        encoder_layer_size = str(encoder_layer_size).replace(", ", "_").replace("[", "").replace("]", "")
        criteria = 5 * true_mean - batch_mean
        print()
        print(f"Criteria: {criteria:.6f} for encoder_layer_size: {encoder_layer_size}, "
              f"adversarial_weight: {adversarial_weight}, epochs: {epochs}, "
              f"learning_rate: {learning_rate},"
              f"reg: {reg}, dropout_rate: {dropout_rate}")
        print()
        with open(output_file, 'a') as f:
            f.write(f"{encoder_layer_size},{adversarial_weight},{epochs},{learning_rate},"
                    f"{reg},{dropout_rate},"
                    f"{batch_mean:.6f},{batch_std:.6f},{true_mean:.6f},{true_std:.6f},{criteria:.6f}\n")


if __name__ == "__main__":
    # Parse command line args --------------------------
    parser = ArgumentParser(description="AutoClass Imputation Example")
    parser.add_argument("-i", "--input-dir", help="Input directory containing the dataset", required=True)
    parser.add_argument("-o", "--output-file", help="Path to output results file", required=True)
    parser.add_argument("-b", "--batch-column", help="Column name for batch information", required=True)
    parser.add_argument("-c", "--true-column", help="Column name for true labels", required=True)
    args = parser.parse_args()
    main(args.input_dir, args.output_file, args.batch_column, args.true_column)


