#!/usr/bin/env Rscript

#######################################################
################# Run cBioPortal Query ################
#######################################################

#run in bash:   $ Rscript cBP_single_gene_query.R [initialize("y"/"n")] ["GOI"]
  #run from full pipeline directory

#first argument passed from bash ("y" or "n") is whether to initialize/re-initialize query workspace
  #data not specific to gene of interest "GOI"
  #will not do so by default

#second argument is Gene of interest gene symbol
  #eg "TP53"

cat("\n#######################################################","\n")
cat("################# Run cBioPortal Query ################","\n")
cat("#######################################################","\n","\n")

bash_args = commandArgs(trailingOnly=TRUE)
#default to initializing if no arg supplied
if(is.na(bash_args[1])){
  initialize_cBP <- "y"
  no_args <- TRUE
  cat("\t","no first argument supplied.","\n\t","Initializing workspace by default","\n\n")
} else {
  initialize_cBP <- bash_args[1]
}

#check that required packages are installed and install if needed
source("Install_required_packages.R")

cat("loading CGDSR API","\n")
library(cgdsr)
mycgds <- CGDS("http://www.cbioportal.org/") #connection object used for all cBP queries
cat("\n")

cat("loading rtracklayer package","\n")
library(rtracklayer)
cat("\n")

#initialize or load cBP pre-query workspace
    if(initialize_cBP == "y" | initialize_cBP == "Y"){
      #run initialization R script to generate a workspace that can be loaded next time rather than re-querying database for each GOI
      source("initialize_cBP.R")
      cat("initialization complete ##########################################","\n")
      cat("Objects in workspace: ","\n")
      print(ls())
      cat("##################################################################\n")
    }else if (initialize_cBP == "n" | initialize_cBP == "N"){
      #load previously generated initialization workspace
      cat("Loading existing cBP workspace","\n")
      if ("initialized_workspace_cBP.RData" %in% list.files()){
        cat("...","\t")
        bash_args_temp <- bash_args
        load("initialized_workspace_cBP.RData")
        bash_args <- bash_args_temp
        cat("load successful","\n")
        cat("Objects in workspace: ","\n")
        print(ls())
      }else{
        stop("initialized_workspace_cBP.RData not found")
      }
    }else{
      stop("wrong first argument (initialize cBP?)")
    }
cat("done","\n\n\n")
# use getCancerStudies(mycgds) to check if there are new studies not included in old workspace

stop("############################################################### stop here for now")

if(is.na(bash_args[2])){
  stop("\t","no gene of interest supplied as second argument")
} else {
  GOI <- bash_args[2]
  cat("GOI supplied: ",GOI,"\n")
}
#create output dir for final output, dont show warning if it already exists, will eventually overwrite by default
initial_directory <- getwd()
output_directory <- paste0("./output/",GOI)
dir.create(output_directory,showWarnings = FALSE)

####################################################################################
#query cBP mutations for GOI using custom function declared in initialized_workspace
cat("mutation query initiated\n")
source("cBP_mutations_query_func.R")
cBP.mutation.query.2(GOI,all.mut.studies.filtered)
cat("mutation query complete\n\n\n")
####################################################################################
#ensemble gene/transcript/protein database annotation API
library("biomaRt")
  #set homo sapian mart for GOI annotations
  hs_ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

cat("UniProtKB query initiated\n")
source("UniProtKB_query.R")

####################################################################################
cat("Ensembl query initiated\n")

#keys for PFAM IDs
# library("PFAM.db")
source("Ensemble_query.R")
cat("Ensembl query complete\n\n\n")

#saving temp workspace
cat("saving temp workspace image\n")
setwd(output_directory)
save.image("temp_workspace_cBP.RData")

setwd(initial_directory)
cat("\n\n\n############################ stage one complete #########################################\n")

