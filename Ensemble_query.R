cat("\n#######################################################","\n")
cat("################# Ensembl query #######################","\n")
cat("#######################################################","\n","\n")

#use biomaRt package to query GOI annotations 

################################################################################
############### retrieve GOI ensembl annotations ###############################
################################################################################
cat("#################### Fetching ensemble annotations for GOI #########################\n")


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

  #exon info
  cat("\nretrieving exon annotations","\n")
  annotation_df_exon <- getBM(
    attributes = c(
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
    filters =  filter_type, 
    values = value ,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  #filter for selected transcript; save as seperate df
  annotation_df_exon_main_transcript <- annotation_df_exon[annotation_df_exon$ensembl_transcript_id == GOI_TRANSCRIPT,]
  
  #get peptide sequences for all transcripts
  cat("\nretrieving peptide sequences","\n")
  peptide_sequences <- getBM(
    attributes = c(
      "ensembl_transcript_id",
      "peptide"
    ),
    filters =  filter_type, 
    values = value ,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  #add peptide length col
  peptide_sequences$peptide_length <- sapply(peptide_sequences$peptide,
                                             function(x){
                                               if(x == "Sequence unavailable") return(0) else return(nchar(x))
                                             })
  #assign GOI ENSG id
  GOI_ENSG <<- getBM(
    attributes = c(
      "ensembl_gene_id"
    ),
    filters =  filter_type, 
    values = value ,
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  )
  # some genes (e.g. MYH11) return two ENSG ids. 
  if(length(unlist(GOI_ENSG)) > 1){
    cat("\n\n\t#### WARNING! More than one ENSG id returned ################\n")
    print(GOI_ENSG)
    GOI_ENSG <<- unlist(GOI_ENSG) # make compatible for gnomAD query iteration
  }
  
  #assign GOI ENSG from GRCh37 (must be used for ExAC querries)
  GOI_ENSG_GRCh37 <<- getBM(
    attributes = c(
      "ensembl_gene_id"
    ),
    filters =  filter_type, 
    values = value ,
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
