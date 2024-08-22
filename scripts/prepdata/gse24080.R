library(dplyr)
library(readr)
library(readxl)
library(SCAN.UPC)
library(stringr)
library(tibble)

CEL_file_pattern = "/tmp/GSE24080/*.CEL.gz"

normalized_dir_path = "/tmp/GSE24080"

for (CEL_file_path in list.files(path = normalized_dir_path, pattern = ".CEL.gz", full.names=TRUE)) {
    CEL_file_name = basename(CEL_file_path)
    CEL_file_name = str_replace_all(CEL_file_name, "\\.gz", "")
    CEL_file_name = str_replace_all(CEL_file_name, "^GSM\\d+_", "")

    out_file_path = paste0(normalized_dir_path, "/", CEL_file_name)

    if (!file.exists(out_file_path)) {
        SCAN(CEL_file_path, outFilePath=out_file_path, probeSummaryPackage="hgu133plus2hsentrezgprobe")
    }
}

eData = NULL
#for (f in list.files(path = normalized_dir_path, pattern = ".CEL$", full.names=TRUE)[1:3]) {
for (f in list.files(path = normalized_dir_path, pattern = ".CEL$", full.names=TRUE)) {
    print(paste0("Parsing from ", f))
    fData = read.table(f, sep="\t", col.names=TRUE)
    colnames(fData) = basename(f)

    if (is.null(eData)) {
        eData = tibble(Gene = rownames(fData)) %>%
            bind_cols(fData)
    } else {
        eData = inner_join(eData, tibble(Gene = rownames(fData)) %>%
            bind_cols(fData), by = "Gene")
    }
}

genes = pull(eData, Gene)
genes = sub("_at", "", genes)
eData = select(eData, -Gene) %>%
    as.matrix() %>%
    t()
colnames(eData) = genes
cel_files = rownames(eData)
eData = bind_cols(tibble(CEL_file = cel_files), as_tibble(eData))

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE24080&format=file&file=GSE24080%5FMM%5FUAMS565%5FClinInfo%5F27Jun2008%5FLS%5Fclean%2Exls%2Egz", "/tmp/gse24080_meta.xls.gz")
gunzip("/tmp/gse24080_meta.xls.gz", overwrite = TRUE)

# It has ArrayScanDate, but it's unclear how these correspond to batches.
# The multiple myeloma (MM) data set (endpoints F, G, H and I) was contributed by the Myeloma Institute for Research and Therapy at the University of Arkansas for Medical Sciences. Gene expression profiling of highly purified bone marrow plasma cells was performed in newly diagnosed patients with MM57,58,59. The training set consisted of 340 cases enrolled in total therapy 2 (TT2) and the validation set comprised 214 patients enrolled in total therapy 3 (TT3)59.  https://www.nature.com/articles/nbt.1665
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
