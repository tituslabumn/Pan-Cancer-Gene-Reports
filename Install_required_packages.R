################# Install_required_packages.R ########################
#checks for installed packages, installs if absent

cat("#####################################################################\n")
cat("############# Checking for/installing required R packages ###########\n")
cat("#####################################################################\n\n\n")

# packages
#   cgdsr
#   rtracklayer
#   biomaRt
#   trackViewer
#   BSgenome.Hsapiens.NCBI.GRCh38
#   knitr

#cgdsr: cBioPortal API
cat("cgdsr:\t")
if("cgdsr" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('cgdsr')
  cat("installation complete\n\n")
}

#BiocManager: for installing some of the below packages
cat("BiocManager:\t")
if("rtracklayer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages("BiocManager")
  cat("installation complete\n\n")
}

library(BiocManager)

#rtracklayer: for converting genomic cordinates to different genome versions
cat("rtracklayer:\t")
if("rtracklayer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rtracklayer")
  cat("installation complete\n\n")
}

#biomaRt: ensembl database API
cat("biomaRt:\t")
if("biomaRt" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt")
  cat("installation complete\n\n")
}

#trackViewer: lolliplot visualization package
cat("trackViewer:\t")
if("trackViewer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("trackViewer")
  cat("installation complete\n\n")
}

#BSgenome.Hsapiens.NCBI.GRCh38: provides this creates Hsapiens object for getting seequences for GRCh38 genome coordinsates
cat("BSgenome.Hsapiens.NCBI.GRCh38:\t")
if("BSgenome.Hsapiens.NCBI.GRCh38" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  source("https://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
  cat("installation complete\n\n")
}

#markdown: required for pdf report generation
cat("rmarkdown:\t")
if("rmarkdown" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages("rmarkdown")
  cat("installation complete\n\n")
}

#knitr: required for pdf report generation
cat("knitr:\t")
if("knitr" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('knitr', dependencies = TRUE)
  cat("installation complete\n\n")
}

#knitr: required for reading in ExAC restAPI pages from GOI queries
cat("rjson:\t")
if("rjson" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('rjson')
  cat("installation complete\n\n")
}

#for codon to AA symbol lookup in relative maping
cat("Biostrings:\t")
if("Biostrings" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
  cat("installation complete\n\n")
}

cat("############# All packages installed ##################\n\n\n") #unsure if failed installations would still progress???
