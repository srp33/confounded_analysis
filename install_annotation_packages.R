library(GEOquery)

tmpDir = tempdir()

#getGEOSuppFiles("GSM505327", makeDirectory=FALSE, baseDir=tmpDir)
#celFilePath = file.path(tmpDir, "GSM555237.CEL.gz")

pkgUrl = "http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133ahsentrezgprobe_25.0.0.tar.gz"
pkgFilePath = paste0(tmpDir, "hgu133ahsentrezgprobe_25.0.0.tar.gz")
download.file(pkgUrl, pkgFilePath)
install.packages(pkgFilePath, repos=NULL, type="source")


pkgUrl = "http://mbni.org/customcdf/25.0.0/entrezg.download/hgu133plus2hsentrezgprobe_25.0.0.tar.gz"
pkgFilePath = paste0(tmpDir, "hgu133plus2hsentrezgprobe_25.0.0.tar.gz")
download.file(pkgUrl, pkgFilePath)
install.packages(pkgFilePath, repos=NULL, type="source")

BiocManager::install(c("pd.hg.u133a", "pd.hg.u133.plus.2"))
