
import yaml
import os

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

OUTPUT_FOLDER = config["output_folder"]
ALL_STUDIES = ["GSE37250_SA", "USA", "India", "GSE37250_M", "Africa", "GSE39941_M"]

def get_test_studies_for_n(n):
    return ALL_STUDIES[:n]

classifiers = config["classifiers"]
num_datasets = config["num_datasets"]

combinations = []
for classifier in classifiers:
    for n in num_datasets:
        for test_study in get_test_studies_for_n(n):
            path = f"{OUTPUT_FOLDER}/results/within_study_cv/individual/{classifier}_{n}_{test_study}.csv"
            exists = os.path.exists(path)
            print(f"{path}: {'EXISTS' if exists else 'MISSING'}")
            if not exists:
                print(f"  Check parent: {os.path.dirname(path)}: {'EXISTS' if os.path.exists(os.path.dirname(path)) else 'MISSING'}")
