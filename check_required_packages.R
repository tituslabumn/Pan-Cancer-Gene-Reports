################# check_required_packages.R ########################
#checks for installed packages, stops if absent

cat("#####################################################################\n")
cat("############# Checking for required R packages ######################\n")
cat("#####################################################################\n\n\n")

# packages
#   BiocManager (version set to 'devl')
#   RCurl
#   XML
#   xml2
#   cgdsr
#   rtracklayer
#   biomaRt
#   trackViewer
#   BSgenome.Hsapiens.NCBI.GRCh38
#   knitr
#   rjson
#   Biostrings
#   tinytex (run 'tinytex::install_tinytex()' to instal distribution of host)

#RCurl (for web queries)
cat("RCurl:\t")
if("RCurl" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}
#XML (for web queries)
cat("XML:\t")
if("XML" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}
#XML (for web queries)
cat("xml2:\t")
if("xml2" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#cgdsr: cBioPortal API
cat("cgdsr:\t")
if("cgdsr" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#BiocManager: for installing some of the below packages
cat("BiocManager:\t")
if("BiocManager" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}
library(BiocManager)

BiocManager::install(version='devel')

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

#rtracklayer: for converting genomic cordinates to different genome versions
cat("rtracklayer:\t")
if("rtracklayer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#biomaRt: ensembl database API
cat("biomaRt:\t")
if("biomaRt" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#trackViewer: lolliplot visualization package
cat("trackViewer:\t")
if("trackViewer" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#BSgenome.Hsapiens.NCBI.GRCh38: provides this creates Hsapiens object for getting seequences for GRCh38 genome coordinsates
cat("BSgenome.Hsapiens.NCBI.GRCh38:\t")
if("BSgenome.Hsapiens.NCBI.GRCh38" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#markdown: required for pdf report generation
cat("rmarkdown:\t")
if("rmarkdown" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#knitr: required for pdf report generation
cat("knitr:\t")
if("knitr" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#knitr: required for reading in ExAC restAPI pages from GOI queries
cat("rjson:\t")
if("rjson" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#for codon to AA symbol lookup in relative maping
cat("Biostrings:\t")
if("Biostrings" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

#for codon to AA symbol lookup in relative maping
cat("tinytex:\t")
if("tinytex" %in% row.names(installed.packages())){
  cat("\talready installed\n")
}else{
  stop("\tPackage missing")
}

cat("############# All packages installed ##################\n\n\n") #unsure if failed installations would still progress???
