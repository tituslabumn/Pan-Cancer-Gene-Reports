cat("\n#######################################################","\n")
cat("################# Ensembl query #######################","\n")
cat("#######################################################","\n","\n")

#use biomaRt package to query GOI annotations and to fetch snps

################################################################################
############### retrieve GOI ensembl annotations ###############################
################################################################################
cat("#################### Fetching ensemble annotations for GOI #########################\n")

#get key lists for PFAM ids
#important note - there seems to be an issue with the select() function related to biomaRt
#need to map keys the long way :[

cat("assigning canonical transcript per key used by cBP:","\n")
GOI_TRANSCRIPT <- cBP_canonical_transcripts[GOI,"ensembl_transcript_id"]
cat(GOI_TRANSCRIPT,"\n\n")

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

################################################################################
############### retrieve ensembl SNPs ##########################################
################################################################################

cat("#################### Fetching ensemble SNPs for GOI #########################\n")

cat("Loading BSgenome.Hsapiens.NCBI.GRCh38 for sequence retrival (filling missing ref alleles)","\n")
library(BSgenome.Hsapiens.NCBI.GRCh38) #for retrieving missing ref alleles

#set mart for SNP query
snp_mart <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

get.ensembl.SNPs <-function(gene = GOI){
  cat("calling get.ensembl.SNPs")
  #there is not hgnc_symbol filter in snp dataset so get ensemble stable ID to query with
  cat("\n","getting ensemble gene ID","\n")
  ensembl_gene <- getBM(attributes = "ensembl_gene_id", filters = "hgnc_symbol", values = gene, mart = hs_ensembl)
  
  cat("\n","getting all SNPs in biomaRt for that ensemble gene ID.","\n")
  snps_BM <- getBM(
    attributes = c(
      "ensembl_gene_stable_id",
      "ensembl_transcript_stable_id",
      "refsnp_id",
      "chrom_start",
      "chrom_end",
      "chrom_strand",
      "allele",
      "translation_start",
      "translation_end",
      "cdna_start",
      "cdna_end",
      "consequence_type_tv",
      "refsnp_source",
      "refsnp_source_description"
    ),
    filters = "ensembl_gene",
    values = ensembl_gene,
    mart = snp_mart
  )
  troubleshooting_snp.df <<- snps_BM
  #returns duplicates for different transcripts (total >300k)
  cat("\n","Unfiltered SNPs returned:","\n")
  cat("\t","-> ",length(snps_BM[,1]),"\n")
  #filter for desired transcript
  cat("\n","filtering for desired transcript","\n")
  snps_BM <-snps_BM[snps_BM$ensembl_transcript_stable_id == GOI_TRANSCRIPT,]
  cat("\t","-> ",length(snps_BM[,1]),"\n")
  #filter for desired mut types
  cat("\n","filtering out intronic, synonymous and UTR varients","\n")
  snps_filtered_BM <- snps_BM[!(snps_BM$consequence_type_tv %in% c("intron_variant","3_prime_UTR_variant","5_prime_UTR_variant","synonymous_variant")),]
  cat("\t","-> ",length(snps_filtered_BM[,1]),"\n")
  #fix instances of missing ref alleles
  cat("\n","Fixing missing ref alleles if any are returned")
  na_ref_index <- (startsWith(snps_filtered_BM$allele,"/"))
  cat("\t","instances:",sum(na_ref_index),"\n")
  if(sum(na_ref_index)>0){
    snps_filtered_BM[na_ref_index,"allele"] <- sapply(which(na_ref_index), function(x) paste0(as.character(getSeq(Hsapiens, as.character(GOI_CHR), start = snps_filtered_BM$chrom_start[x], end = snps_filtered_BM$chrom_end[x])),snps_filtered_BM$allele[x])) #use BSgenome to get sequence
  }
  #remove HGMD variants; cannot obtain alleles with any high thruput API (only available on web page to registered members)
  cat("\n","Removing HGMD variants")
  HGMD_index <- (snps_filtered_BM$allele != "HGMD_MUTATION")
  cat("\t","instances:",sum(!HGMD_index),"\n")
  snps_filtered_BM <- snps_filtered_BM[HGMD_index,]
  #change '-/PhenCode_variation*' alleles on -1 strand to include N's for cariant allele numerated by refsnpid and change strand to +
  cat("\n","Fixing PhenCode variant insertions to 'N's")
  phencodevar_ins_index <- (grepl("PhenCode_variation",snps_filtered_BM$allele) & grepl("ins",snps_filtered_BM$refsnp_id)) #should all be insertions but filter just in case
  cat("\t","instances:",sum(phencodevar_ins_index),"\n")
  if(sum(phencodevar_ins_index,na.rm = TRUE) >0){
    snps_filtered_BM[phencodevar_ins_index,"allele"] <- sapply(which(phencodevar_ins_index), function(x) sub("PhenCode_variation",paste(rep("N", as.numeric(gsub(".*ins","\\1",snps_filtered_BM$refsnp_id[x]))),collapse = ""),snps_filtered_BM$allele[x]))
    snps_filtered_BM[phencodevar_ins_index,"chrom_strand"] <- 1
  }
  
  #assign global data frames for both
  cat("\n","assigning","\n")
  
  assign("biomaRt_SNPs_df_name",paste0(GOI,"_biomaRt_SNPs"), envir = .GlobalEnv)
  assign("biomaRt_SNPs_df_filtered_name",paste0(GOI,"_biomaRt_SNPs_filtered"), envir = .GlobalEnv)
  
  assign(paste0(GOI,"_biomaRt_SNPs"),snps_BM, envir = .GlobalEnv)
  assign(paste0(GOI,"_biomaRt_SNPs_filtered"), snps_filtered_BM, envir = .GlobalEnv)
}

#et/write/assign biomaRt SNPs
get.ensembl.SNPs()

save.image("troubleshooting_workspace.RData") #####################

#assign GOI chromosome
cat("assigning GOI chromosome number","\n")
GOI_CHR<- as.character(getBM(
  attributes = c('chromosome_name'),
  filters =  "hgnc_symbol", 
  values = GOI ,
  mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
))

cat("SNP fetching complete","\n\n\n")

################################################################################
############### prepare data for ensemble VEP annotation #######################
################################################################################
cat("#################### Writing out data for Ensemble VEP inport via PERL API #########################\n")

#function to write tab delim file that will serve as input for ensemblVEP API query
#datatype; 1 = cBP (IDs are patient samples IDs pasted to AA change) ; 0 = dbSNP or other data base output (IDs are rs### ids)
#chr is chromosome number as numeric
#strand.direction; "+" or "-" NOTE - based on ensembl dbSNP and cBP outputs allele and chrom_start/end should always be '+' (regardless of gene strand)
#assign.df is TRUE or FALSE; can return df.out for diagnostic purposes
write.ensemblVEP.input.GOI <- function(df.in, data.type = 1, chr = GOI_CHR, strand.direction = "+", write_name, assign.df = TRUE) {
  VEPin_name <- paste0(write_name,"_VEP_IN")
  cat("\n","VEPin_name: ",VEPin_name,"\n")
  if(data.type == 1){ #cBP data
    cat("\n","prepping df for cBP variant data","\n")
    df.out <- data.frame(
      #don't assign rownames so line nubers can be associated with warnings returened by VEP API query
      "chromosome" = chr,
      "start" = df.in$start_position,
      "end" = df.in$end_position,
      "allele" = paste(df.in$reference_allele, df.in$variant_allele, sep = "/"),
      "strand" = strand.direction,
      "identifyer" = df.in$unique_id,
      stringsAsFactors = FALSE
    )
    
  }else if (data.type == 0) { #BM SNP data
    cat("\n","prepping df for biomaRt SNP data","\n")
    df.out <- data.frame(
      #dont assign row.names, there are duplicate rsIDs
      "chromosome" = chr,
      "start" = df.in$chrom_start,
      "end" = df.in$chrom_end,
      "allele" = df.in$allele,
      "strand" = strand.direction,
      "identifyer" = df.in$refsnp_id,
      stringsAsFactors = FALSE
    )
  }else{
    cat("\n","! wrong input for data.type arg !","\n")
  }
  cat("\n","writing out and assigning","\n")
  setwd(output_directory)
  write.table(df.out,paste0(VEPin_name,".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  if (assign.df == TRUE){
    assign(VEPin_name, df.out, envir = .GlobalEnv)
  }
  setwd(initial_directory)
}



cat("#################### Writing VEP input files #########################\n")
#write cBP data
write.ensemblVEP.input.GOI(get(GOI_mut_df_name), data.type = 1, write_name = GOI_mut_df_name)
#write BM SNP data all test
write.ensemblVEP.input.GOI(get(biomaRt_SNPs_df_filtered_name), data.type = 0, write_name = biomaRt_SNPs_df_filtered_name)

cat("#################### Writing VEP run script #########################\n")
#write bash script to execute Ensembl VEP PERL API annotations
setwd(output_directory)
fileConn <- file("VEP_run.sh")
lines <- c("#!/bin/bash",
        paste0("./../.././ensembl_VEP_query.sh ./",GOI,"_cBP_mutations_VEP_IN.txt '",GOI,"_cBP_mutations_VEP_OUT' /mnt/DATA/.vep"),
        paste0("./../.././ensembl_VEP_query.sh ./",GOI,"_biomaRt_SNPs_filtered_VEP_IN.txt '",GOI,"_biomaRt_SNPs_filtered_VEP_OUT' /mnt/DATA/.vep")
        )
writeLines(lines,fileConn)
close(fileConn)
rm(lines,fileConn)
setwd(initial_directory)