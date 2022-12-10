library(dplyr)
library(GEOquery)
library(readr)
library(stringr)
library(tibble)

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE49711&format=file&file=GSE49711%5FSEQC%5FNB%5FTAV%5FG%5Flog2%2Efinal%2Etxt%2Egz", "/tmp/gse49711.expr.tsv.gz")

expr_data = read_tsv("/tmp/gse49711.expr.tsv.gz") %>%
    filter(Gene_set == "Gene_AceView") %>%
    filter(RefSeq_transcript_ID != ".") %>%
    select(Gene, starts_with("SEQC_")) %>%
    as.data.frame()

rownames(expr_data) = pull(expr_data, Gene)

expr_data = data.matrix(expr_data[,2:ncol(expr_data)])

# Remove genes with lots of zero values
num_zero = apply(expr_data, 1, function(x) { sum(x==0) })
expr_data = expr_data[-which(num_zero > (ncol(expr_data) / 2)),]

# This returned zero.
#print(sum(is.na(expr_data)))

expr_data = t(expr_data) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample_ID")

metadata = as_tibble(as.data.frame(getGEO("GSE49711"))) %>%
    dplyr::rename(Sample_ID = `GSE49711_series_matrix.txt.gz.title`) %>%
    dplyr::rename(Class = `GSE49711_series_matrix.txt.gz.class.label.ch1`) %>%
    #dplyr::rename(Age_at_Diagnosis = `GSE49711_series_matrix.txt.gz.age.at.diagnosis.ch1`) %>%
    #dplyr::rename(Death_from_Disease = `GSE49711_series_matrix.txt.gz.death.from.disease.ch1`) %>%
    #dplyr::rename(High_Risk = `GSE49711_series_matrix.txt.gz.high.risk.ch1`) %>%
    dplyr::rename(INSS_Stage = `GSE49711_series_matrix.txt.gz.inss.stage.ch1`) %>%
    #dplyr::rename(MYCN_Status = `GSE49711_series_matrix.txt.gz.mycn.status.ch1`) %>%
    #dplyr::rename(Progression = `GSE49711_series_matrix.txt.gz.progression.ch1`) %>%
    dplyr::rename(Sex = `GSE49711_series_matrix.txt.gz.Sex.ch1`) %>%
    select(!starts_with("GSE49711_")) %>%
    filter(Class %in% c("0", "1"))

data = inner_join(metadata, expr_data, by="Sample_ID")

#print(table(data$Sex))
#print(table(data$MYCN_Status))

if (!dir.exists("/data/gse49711"))
    dir.create("/data/gse49711")

write_csv(data, "/data/gse49711/unadjusted.csv")
