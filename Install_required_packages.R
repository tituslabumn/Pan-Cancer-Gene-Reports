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

#RCurl (for web queries)
cat("RCurl:\t")
if("RCurl" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('RCurl', dependencies = TRUE)
  if ("RCurl" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}
#XML (for web queries)
cat("XML:\t")
if("XML" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('XML', dependencies = TRUE)
  if ("XML" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#cgdsr: cBioPortal API
cat("cgdsr:\t")
if("cgdsr" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('cgdsr', dependencies = TRUE)
  if ("cgdsr" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#BiocManager: for installing some of the below packages
cat("BiocManager:\t")
if("BiocManager" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages("BiocManager", dependencies = TRUE)
  if ("BiocManager" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}
library(BiocManager)

#rtracklayer: for converting genomic cordinates to different genome versions
cat("rtracklayer:\t")
if("rtracklayer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  BiocManager::install("rtracklayer")
  if ("rtracklayer" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#biomaRt: ensembl database API
cat("biomaRt:\t")
if("biomaRt" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  BiocManager::install("biomaRt")
  if ("biomaRt" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#trackViewer: lolliplot visualization package
cat("trackViewer:\t")
if("trackViewer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  BiocManager::install("trackViewer")
  if ("trackViewer" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#BSgenome.Hsapiens.NCBI.GRCh38: provides this creates Hsapiens object for getting seequences for GRCh38 genome coordinsates
cat("BSgenome.Hsapiens.NCBI.GRCh38:\t")
if("BSgenome.Hsapiens.NCBI.GRCh38" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
  if ("BSgenome.Hsapiens.NCBI.GRCh38" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#markdown: required for pdf report generation
cat("rmarkdown:\t")
if("rmarkdown" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages("rmarkdown", dependencies = TRUE)
  if ("rmarkdown" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#knitr: required for pdf report generation
cat("knitr:\t")
if("knitr" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('knitr', dependencies = TRUE)
  if ("knitr" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#knitr: required for reading in ExAC restAPI pages from GOI queries
cat("rjson:\t")
if("rjson" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  install.packages('rjson', dependencies = TRUE)
  if ("rjson" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

#for codon to AA symbol lookup in relative maping
cat("Biostrings:\t")
if("Biostrings" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  cat("installing...\n")
  BiocManager::install("Biostrings")
  if ("Biostrings" %in% row.names(installed.packages())) cat("installation complete\n\n") else stop("instalation FAILED\n")
}

cat("############# All packages installed ##################\n\n\n") #unsure if failed installations would still progress???
