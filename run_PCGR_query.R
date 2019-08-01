#!/usr/bin/env Rscript

cat("\n#######################################################","\n")
cat("################# Initialize PCGR #####################","\n")
cat("#######################################################","\n","\n")


#check that required packages are installed and install if needed
source("check_required_packages.R")

initial_directory <- getwd()

#load previously generated initialization workspace, if in output directory load from there preferentially
cat("Loading existing cBP workspace","\n")
setwd("/OUTPUT/")
if ("initialized_workspace_cBP.RData" %in% list.files()){
  cat("...","\t")
  load("initialized_workspace_cBP.RData")
  cat("load successful","\n")
  cat("Objects in workspace: ","\n")
  print(ls())
}else if ("initialized_workspace_cBP.RData" %in% list.files(path = initial_directory)){
  setwd(initial_directory)
  cat("...","\t")
  load("initialized_workspace_cBP.RData")
  cat("load successful","\n")
  cat("Objects in workspace: ","\n")
  print(ls())
  setwd("/OUTPUT/")
}else{
  stop("initialized_workspace_cBP.RData not found")
}
cat("done","\n\n\n")
# use getCancerStudies(mycgds) to check if there are new studies not included in old workspace

bash_args = commandArgs(trailingOnly=TRUE)
if(is.na(bash_args[1])){
  stop("\t","no gene of interest supplied as argument")
} else {
  GOI <- bash_args[1]
  cat("GOI supplied: ",GOI,"\n")
}
#create output dir for final output, dont show warning if it already exists, will eventually overwrite by default
output_directory <- paste0("/OUTPUT/",GOI)
dir.create(output_directory,showWarnings = FALSE)

####################################################################################
#query cBP mutations for GOI using custom function declared in initialized_workspace
cat("mutation query initiated\n")

cat("loading CGDSR API","\n")
library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/") #connection object used for all cBP queries
cat("\n")

cat("loading rtracklayer package","\n")
library(rtracklayer)
cat("\n")

source("/PCGR/Pan-Cancer-Gene-Reports/cBP_mutations_query_func.R")
cBP.mutation.query.2(GOI,all.mut.studies.filtered)
cat("mutation query complete\n\n\n")
####################################################################################
#ensemble gene/transcript/protein database annotation API
library("biomaRt")
  #set homo sapian mart for GOI annotations
  hs_ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

cat("UniProtKB query initiated\n")
source("/PCGR/Pan-Cancer-Gene-Reports/UniProtKB_query.R")

####################################################################################
cat("Ensembl query initiated\n")

#keys for PFAM IDs
# library("PFAM.db")
source("/PCGR/Pan-Cancer-Gene-Reports/Ensemble_query.R")
cat("Ensembl query complete\n\n\n")

save.image("troubleshooting_workspace.RData") #####################

#retrieve ExAC variant data for GOI
source("/PCGR/Pan-Cancer-Gene-Reports/ExAC_query.R")

save.image("troubleshooting_workspace.RData") #####################

#relative position mapping for visualization
setwd(initial_directory)
source("/PCGR/Pan-Cancer-Gene-Reports/Relative_position_mapping.R")

save.image("troubleshooting_workspace.RData") #####################

#analasis and figure data initialization
source("/PCGR/Pan-Cancer-Gene-Reports/Generate_figures.R")

#Knit data to output PDF
cat("############### Knitting report to PDF #################\n\n\n")
Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") #must specify if running in bash rather than Rstudio
rmarkdown::render("/PCGR/Pan-Cancer-Gene-Reports/output_markdown_format.Rmd",output_file = paste0(GOI," query output report (",Sys.Date(),").pdf"),output_dir = output_directory)

#save final workspace
cat("saving final workspace to GOI output folder\n\n\n")
setwd(output_directory)
save.image("final_workspace.RData") #####################

cat("################################################ Session info #################################################\n\n\n")
sessionInfo()

cat("\n\n\n############################################# Pipeline finished :] ##########################################################\n\n\n")
