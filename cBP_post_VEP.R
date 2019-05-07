#!/usr/bin/env Rscript

#######################################################
########## Run cBioPortal Query: post VEP #############
#######################################################

#Read in ensembl VEP output and complete analysis


bash_args_2 = commandArgs(trailingOnly=TRUE)
setwd(paste0("./output/",bash_args_2[1]))
cat("recovering temp workspace image:...")
load("temp_workspace_cBP.RData")
cat("load successful","\n\n\n")
setwd(initial_directory)

save.image("troubleshooting_workspace.RData") #####################

read.ensemblVEP.output.GOI <- function(output_file){
  setwd(output_directory)
  cat("\n","reading in file","\n")
  df <- read.table(output_file , header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  setwd(initial_directory)
  cat("\n","filtering for cononical transcript","\n")
  #filter for only canonical transcript
  df <- df[df$Feature == GOI_TRANSCRIPT,]
  
  #protein_position is a range in some cases, change to numeric as new col that changes range to first AA position
  cat("\n","adding first amino acid position col as numeric","\n")
  #expect coercion warning
  df$first_amino_acid_position <- sapply(as.character(df$Protein_position), function(x) as.numeric(unlist(strsplit(x,"-")))[1])
  
  cat("\n","changing col variable types","\n")
  #change variable type for many cols
  cols_to_make_num <- c(
    "AF",
    "MAX_AF"
  )
  #expect coercion warnings
  for (x in cols_to_make_num) {
    cat("\t",x,"\n")
    df[,x]<-as.numeric(df[,x])
  }
  #assign
  cat("\n","assigning","\n")
  assign(output_file,df , envir = .GlobalEnv)
}

save.image("troubleshooting_workspace.RData") #####################

#parse back in
read.ensemblVEP.output.GOI(paste0(GOI,"_cBP_mutations_VEP_OUT"))
read.ensemblVEP.output.GOI(paste0(GOI,"_biomaRt_SNPs_filtered_VEP_OUT"))
cat("read in succesful","\n")

save.image("troubleshooting_workspace.RData") #####################

#relative position mapping for visualization
setwd(initial_directory)
source("Relative_position_mapping.R")

save.image("troubleshooting_workspace.RData") #####################

#analasis and figure data initialization
source("Generate_figures.R")

#Knit data to output PDF
cat("############### Knitting report to PDF #################\n\n\n")
Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc") #must specify if running in bash rather than Rstudio
rmarkdown::render("output_markdown_format.Rmd",output_file = paste0(GOI," query output report (",Sys.Date(),").pdf"),output_dir = output_directory)

#save final workspace
cat("saving final workspace to GOI output folder\n\n\n")
setwd(output_directory)
save.image("final_workspace.RData") #####################

cat("################################################ Session info #################################################\n\n\n")
sessionInfo()

cat("\n\n\n############################################# Pipeline finished :] ##########################################################\n\n\n")

