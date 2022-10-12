library(doParallel)
library(dplyr)
library(readr)
library(readxl)
library(SCAN.UPC)
library(stringr)
library(tibble)

registerDoParallel(cores=16)

CEL_file_pattern = "/tmp/GSE24080/*.CEL.gz"

eSet = SCAN(CEL_file_pattern, probeSummaryPackage="hgu133plus2hsentrezgprobe")

eData = t(data.matrix(exprs(eSet))) %>%
  as.data.frame() %>%
  rownames_to_column("CEL_file") %>%
  as_tibble() %>%
  mutate(CEL_file = str_replace_all(CEL_file, "\\.gz", "")) %>%
  mutate(CEL_file = str_replace_all(CEL_file, "^GSM\\d+_", ""))

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE24080&format=file&file=GSE24080%5FMM%5FUAMS565%5FClinInfo%5F27Jun2008%5FLS%5Fclean%2Exls%2Egz", "/tmp/gse24080_meta.xls.gz")
gunzip("/tmp/gse24080_meta.xls.gz", overwrite = TRUE)

# It has ArrayScanDate, but it's unclear how these correspond to batches.
pData = read_excel("/tmp/gse24080_meta.xls") %>%
    filter(`MAQC_Distribution_Status` %in% c("Training", "Validation")) %>%
    dplyr::rename(batch = `MAQC_Distribution_Status`) %>%
    dplyr::mutate(batch = str_replace_all(batch, "Training", "1")) %>%
    dplyr::mutate(batch = str_replace_all(batch, "Validation", "2")) %>%
    dplyr::rename(Sample = PATID) %>%
    dplyr::rename(CEL_file = `CELfilename`) %>%
    dplyr::rename(cytogenetic_abnormality = `Cyto Abn`) %>%
    dplyr::rename(age = `AGE`) %>%
    dplyr::rename(race = `RACE`) %>%
    dplyr::rename(efs_outcome_label = `EFS_MO JUN2008`) %>%
    dplyr::rename(os_outcome_label = `OS_MO JUN2008`) %>%
    dplyr::rename(sex_label = `CPS1`) %>%
    dplyr::rename(random_label = `CPR1`) %>%
    dplyr::select(batch, Sample, CEL_file, cytogenetic_abnormality, age, race, efs_outcome_label, os_outcome_label, sex_label, random_label)

if (!dir.exists("/data/gse24080"))
    dir.create("/data/gse24080")

inner_join(eData, pData) %>%
    dplyr::select(-CEL_file) %>%
    dplyr::select(batch, Sample, cytogenetic_abnormality, age, race, efs_outcome_label, os_outcome_label, sex_label, random_label, matches("^\\d.+")) %>%
    write_csv("/data/gse24080/unadjusted.csv")
