#!/bin/bash

set -e

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194
# https://pubmed.ncbi.nlm.nih.gov/20064235/
bash /scripts/prepdata/gse20194.sh

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24080
# https://pubmed.ncbi.nlm.nih.gov/20064235/
bash /scripts/prepdata/gse24080.sh

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711
# https://pubmed.ncbi.nlm.nih.gov/25150839/
bash /scripts/prepdata/gse49711.sh

#Other possibilities:
#bash /scripts/prepdata/bladderbatch.sh
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47792 (SEQC superseries)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49711 (same as GSE49711 but uses Agilent microarrays)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54275 (specifically, the samples for GPL15932)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25507 (Affymetrix Human Genome U133 Plus 2.0)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65204 (Agilent)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58979 (Affymetrix PrimeView)
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19750 (Affymetrix Human Genome U133 Plus 2.0)
#bash /scripts/prepdata/tcga.sh
