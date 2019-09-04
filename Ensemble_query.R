cat("\n#######################################################","\n")
cat("################# Ensembl query #######################","\n")
cat("#######################################################","\n","\n")

#use biomaRt package to query GOI annotations 

################################################################################
############### retrieve GOI ensembl annotations ###############################
################################################################################
cat("#################### Fetching ensemble annotations for GOI #########################\n")

cat("GOI: ",GOI,"\n")

cat("assigning canonical transcript per key used by cBP:","\n")
GOI_TRANSCRIPT <- cBP_canonical_transcripts[GOI,"ensembl_transcript_id"]
cat(GOI_TRANSCRIPT,"\n\n")
if(is.na(GOI_TRANSCRIPT)) cat("\t!!!!! no canonicical transcript found, default to uniprot sequence; resolve with relative mapping below\n\n")

save.image("troubleshooting_workspace.RData") #####################

#function for getting annotations for a single gene
#assigns df of transcript level annotation
#asssigns df of metadata (e.g. number of transcripts, domains etc.)
BM_GOI_annotation <- function(filter_type = "hgnc_symbol", value = GOI) {
  cat("Calling BM_GOI_annotoation()")

  cat("\nretrieving ENSG id","\n")
  #assign GOI ENSG id
  GOI_ENSG <<- getBM(
    attributes = c(
      "ensembl_gene_id"
    ),
    filters =  "ensembl_transcript_id", 
    values = GOI_TRANSCRIPT,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  print(GOI_ENSG)
  # some genes (e.g. MYH11) return two ENSG ids. 
  if(length(unlist(GOI_ENSG)) > 1){
    cat("\n\n\t#### WARNING! More than one ENSG id returned ################\n")
    GOI_ENSG <<- unlist(GOI_ENSG) # make compatible for gnomAD query iteration
    #there may no longer be an issue now that ENSG is being colected with ENST filter
  }
  
  #get accurate CHR (cBP returns 23 for X)
  cat("\nretrieving chromosome name","\n")
  GOI_CHR <- getBM(
    attributes = c(
      "chromosome_name"
    ),
    filters =  "ensembl_gene_id", 
    values = GOI_ENSG,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  # convert data frame to character
  GOI_CHR <<- as.character(unlist(GOI_CHR))[1]
  print(GOI_CHR)
  
  #exon info
  cat("\nretrieving exon annotations","\n")
  annotation_df_exon <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "ensembl_exon_id",
      "rank",
      "exon_chrom_start",
      "exon_chrom_end",
      "cds_start",
      "cds_end",
      "cds_length",
      "genomic_coding_start",
      "genomic_coding_end",
      "transcript_start",
      "transcript_end",
      "transcript_length",
      "strand"
    ),
    filters =  "ensembl_gene_id", 
    values = GOI_ENSG,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  #filter for selected transcript; save as seperate df
  annotation_df_exon_main_transcript <- annotation_df_exon[annotation_df_exon$ensembl_transcript_id == GOI_TRANSCRIPT,]
  print(annotation_df_exon[!duplicated(annotation_df_exon$ensembl_transcript_id),c("ensembl_transcript_id","cds_length")])
  
  cat("\nretrieving transcript support annotations and sorting/filtering","\n")
  # retrive support information for returned transcripts
  #   https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html
  #   TSL
  #   GENCODE BASIC
  #   APPRIS
  GOI_transcript_support <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "ensembl_transcript_id",
      "transcript_tsl",
      "transcript_gencode_basic",
      "transcript_appris"
    ),
    filters =  "ensembl_gene_id", 
    values = GOI_ENSG,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  # remove cases that are tslNA
  GOI_transcript_support <- GOI_transcript_support[!grepl("NA",GOI_transcript_support$transcript_tsl),]
  # sort based on level of support, GOI_TRANSCRIPT always first even if not highest support level (often the case)
  GOI_transcript_support$sort <- 0
  # 'main cBP transcript' set to 1000
  GOI_transcript_support[GOI_transcript_support$ensembl_transcript_id == GOI_TRANSCRIPT,"sort"] <- 1000
  # (6 - tsl rating) * 10 (0 if missing)
  # if no tsl returned substr returns '' ; convert to 6 with gsub
  GOI_transcript_support$sort <- GOI_transcript_support$sort + 10*(6-as.numeric(gsub("^$","6",substr(gsub("tsl","",GOI_transcript_support$transcript_tsl),1,1))))
  # add 60 for GENCODE BASIC
  GOI_transcript_support$sort <- GOI_transcript_support$sort + 60 * !(GOI_transcript_support$transcript_gencode_basic == "")
  # principal APRIS add 500
  GOI_transcript_support$sort <- GOI_transcript_support$sort + 500 * grepl("principal",GOI_transcript_support$transcript_appris)
  # alternative APPris add 100
  GOI_transcript_support$sort <- GOI_transcript_support$sort + 100 * grepl("alternative",GOI_transcript_support$transcript_appris)
  # sort
  GOI_transcript_support <- GOI_transcript_support[rev(order(GOI_transcript_support$sort)),]
  # remove values under 40 (non appris/basic with tsl3 and greater)
  cat("\tFiltering out ",sum(GOI_transcript_support$sort <50),"transcripts with low support\n\n")
  GOI_transcript_support <<- GOI_transcript_support[GOI_transcript_support$sort >=50 ,]
  cat("GOI_TRANSCRIPT ('canonical transcript' used by cBP): ",GOI_TRANSCRIPT,"\n\n")
  print(GOI_transcript_support[,2:5])
  
  #get peptide sequences for all transcripts
  cat("\nretrieving peptide sequences","\n")
  peptide_sequences <- getBM(
    attributes = c(
      "ensembl_transcript_id",
      "peptide"
    ),
    filters =  "ensembl_gene_id", 
    values = GOI_ENSG,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  #add peptide length col
  peptide_sequences$peptide_length <- sapply(peptide_sequences$peptide,
                                             function(x){
                                               if(x == "Sequence unavailable") return(0) else return(nchar(x))
                                             })
  print(peptide_sequences[peptide_sequences$ensembl_transcript_id == GOI_TRANSCRIPT,"peptide"])
  
  
  #assign GOI ENSG from GRCh37 (must be used for ExAC querries)
  GOI_ENSG_GRCh37 <<- getBM(
    attributes = c(
      "ensembl_gene_id"
    ),
    filters =  "ensembl_transcript_id", 
    values = GOI_TRANSCRIPT ,
    mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", dataset="hsapiens_gene_ensembl")
  )
  
  #assign global data frames for both
  cat("\n\nassigning and writing","\n\n")
  assign(paste(value,"exon_annotation",sep = "_"),annotation_df_exon,envir = .GlobalEnv)
  assign(paste(value,"exon_annotation_main_transcript",sep = "_"),annotation_df_exon_main_transcript,envir = .GlobalEnv)
  assign(paste(value,"peptide_sequences",sep = "_"),peptide_sequences,envir = .GlobalEnv)
}

#get/assign annotations
BM_GOI_annotation()
cat("GOI annotation complete","\n\n\n")

save.image("troubleshooting_workspace.RData") #####################
