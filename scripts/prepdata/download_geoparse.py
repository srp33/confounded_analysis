import GEOparse
from collections import Counter

gse = GEOparse.get_GEO(filepath="/home/aw998/groups/grp_batch_effects/data/raw_data/gse96058_hiseq/grp_batch_effects/data/raw_download/gse96058_hiseq/analysis_ready_metadata_GSE96058_HiSeq.tsv", how="brief")

platform_counts = Counter(gsm.metadata["platform_id"][0] for gsm in gse.gsms.values())

for platform, count in platform_counts.items():
    print(f"{platform}: {count} samples")