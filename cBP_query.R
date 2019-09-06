################# cBP_query.R ########################

  cat("##############################################################\n")
  cat("############### Running cBP mutation query ###################\n")
  cat("##############################################################\n\n\n")
  #declare empty df with colnames for return of getMutationData()
  cat("Creating blank data frame with colnames of cBP output","\n")
  mut.df <- data.frame(
    "entrez_gene_id" = integer(),
    "gene_symbol" = character(),
    "case_id" = character(),
    "sequencing_center" = character(),
    "mutation_status" = logical(),
    "mutation_type" = character(),
    "validation_status" = character(),
    "amino_acid_change" = character(),
    "functional_impact_score" = character(),
    "xvar_link" = character(),
    "xvar_link_pdb" = character(),
    "xvar_link_msa" = character(),
    "chr" = integer(),
    "start_position" = integer(),
    "end_position" = integer(),
    "reference_allele" = character(),
    "variant_allele" = character(),
    "reference_read_count_tumor" = integer(),
    "variant_read_count_tumor" = integer(),
    "reference_read_count_normal" = logical(),
    "variant_read_count_normal" = logical(),
    "genetic_profile_id" = character(),
    stringsAsFactors = FALSE
  )
  
  save.image("troubleshooting_workspace.RData") #####################
  
  cat("query mutations and clinical data from each mutation-data-containing study:","\n\n")
  #get mut data for all relavent studies and append all hits to the df
  
  count <- 1
  for (x in all.mut.studies.filtered) {
    cat("[",count,"/",length(all.mut.studies.filtered),"]\t\t")
    mut.study <- getMutationData(mycgds,paste(x,"all",sep = "_"),paste(x,"mutations",sep = "_"),GOI)
    cat(length(mut.study[,1]),"\t\t",x,"\n")
    mut.df <- rbind(mut.df,mut.study)
    count <- count + 1
  }
  rm(x,count)
  
  assign("pre_filter_diagotstic_df",mut.df)
  save.image("troubleshooting_workspace.RData") #####################
  
  #make fusion.df and remove fusions from mut.df
  cat("splitting fusions and other muts","\n")
  fusion.df <- mut.df[mut.df$start_position == -1,] #fusions are the only types to return -1
  cat("\n","total fusions returned:",length(fusion.df[,1]),"\n")
  mut.df <- mut.df[mut.df$start_position != -1,]
  cat("\n","total alterations so far:",length(mut.df[,1]),"\n")
  
  
  cat("\n\n\nGOI is on chromosome #:",GOI_CHR,"\n")
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #convert cordinates
  cat("\n","converting chr locations from hg19/GRCh37 to GRCh38.p12 via  rtracklayer/liftOver","\n")
  #convert chr locations from hg19/GRCh37 to GRCh38.p12 via  rtracklayer/liftOver
  #chain imported above
  #in case variant positions return empty converted value (found errors in end positions in DACH1); set those values to 0 and remove after
  
  mut.df <- tryCatch({
    mut.df$start_position <- unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = mut.df$start_position, width = 1)),Chain_19to38)))
    mut.df #return if no errors
    },
    error=function(err){
      cat("\terror in  start position:\n")
      print(err)
      not_empty_index <- sapply(sapply(mut.df$start_position, function(x) unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = x, width = 1)),Chain_19to38)))),length)
      cat("\tempty conversion positions:", sum(!(not_empty_index)),"\n")
      #remove problematic variants
      mut.df <- mut.df[as.logical(not_empty_index),]
      #re-try
      mut.df$start_position <- unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = mut.df$start_position, width = 1)),Chain_19to38)))
      return(mut.df)
    }
  )
  
  mut.df <- tryCatch({
    mut.df$end_position <- unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = mut.df$end_position, width = 1)),Chain_19to38)))
    mut.df #return if no errors
  },
  error=function(err){
    cat("\terror in  end position:\n")
    print(err)
    not_empty_index <- sapply(sapply(mut.df$end_position, function(x) unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = x, width = 1)),Chain_19to38)))),length)
    cat("\tempty conversion positions:", sum(!(not_empty_index)),"\n")
    #remove problematic variants
    mut.df <- mut.df[as.logical(not_empty_index),]
    #re-try
    mut.df$end_position <- unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = mut.df$end_position, width = 1)),Chain_19to38)))
    return(mut.df)
  }
  )
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #filtering
  cat("\n","####### Filtering ########","\n")
  #add unique_id column that pastes case_ID and amino acid change sep with "@@"
  #used to check for duplicates etc
  cat("adding $unique_id","\n")
  mut.df$unique_id <- paste(mut.df$case_id,mut.df$amino_acid_change, sep = "@@")
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #df will have duplicate studies
  cat("total alterations so far:","\n")
  print(length(mut.df[,1]))
  cat("\n","checking for duplicate unique ids:","\n")
  dup.cases <- mut.df[mut.df$unique_id %in% unique(mut.df$unique_id[duplicated(mut.df$unique_id)]),]
  cat("\n","number of duplicated unique ids:","\t")
  cat(length(unique(dup.cases$unique_id)),"\t")
  cat("within total numbr of duplicates","\t")
  cat(length(dup.cases$unique_id), "\n")
  if(length(dup.cases$unique_id) > 0){
    cat("\n","removing duplications","\n")
    mut.df <- mut.df[-as.numeric(row.names(dup.cases[duplicated(dup.cases$unique_id),])),]
    cat("\n","total alterations so far:",length(mut.df[,1]),"\n")
  }
  save.image("troubleshooting_workspace.RData") #####################
  
  #'sclc_cancercell_gardner_2017' appears to have normal and relapse samples for each patient leading to false 140% of cases annoatation
  #not all querys will return mutations in this study!
  #fix this
  #upon futher inspection most of these samples, not just the normal and resistant pairs, harbor the exact same mutation(s)
  #reading the source paper [PMID: 28196596] these are supposed to be PDX mouse model pairs derived from ~10 independant cases
  #it seeems essentially impossible that all of these patients would harbor the exact same mutation(s)
  #assume this is some sort of upload/annotation error (could check other mutations later)
  #these are all identified via the sequencing_center = 'mskcc.org'
  #remove all duplicates so only one of each is represented in this df
  cat("\n","filtering mskcc.org studies","\n")
  #only filter if there are mutations from this study returned
  if ('mskcc.org' %in% mut.df$sequencing_center) {
    mskcc.org_index <- which(mut.df$sequencing_center == 'mskcc.org')###########################################################!!!!!!!!!!!!!!!!!!
    cat("total instances: ","\t",length(mskcc.org_index),"\n")
    #only filter if there are duplicates
    if(sum(duplicated(mut.df$amino_acid_change[mskcc.org_index]))>0){
      mut.df <- mut.df[-(mskcc.org_index[ duplicated(mut.df$amino_acid_change[mskcc.org_index]) ] ),]
    }else{
      cat("No duplicates; no filtering","\n")
    }
  }else{
    cat("There were no mutations from mskcc.org study in this gene","\n")
  }
  cat("\n","total alterations so far:",length(mut.df[,1]),"\n")
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #some bases in ref and alt column are "TRUE" rather than "T"
  #I confirmed this is correct based on the ref base from genome browser and AAchange listed
  #substitute all "TRUE" values to "T" in these columns
  cat("fixing TRUE > T alleles","\n")
  mut.df$reference_allele <- gsub("TRUE","T",mut.df$reference_allele)
  mut.df$variant_allele <- gsub("TRUE","T",mut.df$variant_allele)
  cat("total alterations so far:","\n")
  print(length(mut.df[,1]))
  
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #some ref alleles can be marked as NA - these appear to be for insertions where NA should be "-"
  if(sum(is.na(mut.df$reference_allele)) > 0){
    cat("fixing NA ref alleles","\n")
    refNAindex <- which(is.na(mut.df$reference_allele))
    cat("\t","mut-types with ref-allele NAs:","\n")
    print(unique(mut.df$mutation_type[refNAindex]))
    mut.df[refNAindex,"reference_allele"] <- "-"
  }
  #same as above for some deletions where NA should be "-"
  if(sum(is.na(mut.df$variant_allele)) > 0){
    cat("fixing NA var alleles","\n")
    varNAindex <- which(is.na(mut.df$variant_allele))
    cat("\t","mut-types with var-allele NAs:","\n")
    print(unique(mut.df$mutation_type[varNAindex]))
    mut.df[varNAindex,"variant_allele"] <- "-"
  }
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #convert insertion chr cordinates to ensembl-style (start = end + 1)   
  cat("converting insertion chr cordinates to ensembl-style (start = end + 1)","\n")
  mut.df[mut.df$reference_allele == "-" & !is.na(mut.df$reference_allele),"start_position"] <- mut.df[mut.df$reference_allele == "-"& !is.na(mut.df$reference_allele),"start_position"] + 1
  mut.df[mut.df$reference_allele == "-" & !is.na(mut.df$reference_allele),"end_position"] <- mut.df[mut.df$reference_allele == "-" & !is.na(mut.df$reference_allele),"end_position"] - 1
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #remove never-used cols
  mut.df <- mut.df[,!(colnames(mut.df) %in% c("sequencing_center","mutation_status","validation_status","functional_impact_score","xvar_link","xvar_link_pdb","xvar_link_msa","chr"))]
  
  #add columns for later use
  #add column to mut.df that has the study name w/o profile suffex
  cat("adding $study","\n")
  mut.df$study <- sub("_mutations","", mut.df$genetic_profile_id)
  #add $AA_change_freq
  cat("adding $AA_change_freq","\n")
  mut.df$AA_change_freq <- sapply(mut.df$amino_acid_change, function(x) sum(mut.df$amino_acid_change == x))
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #add $AA
  mut.df$AA <- sapply(mut.df$amino_acid_change, function(x) as.numeric(sub(".*?([0-9]+).*", "\\1", x))) #keeps "MUTATED" aa changes for now
  #add $AA_freq
  mut.df$AA_freq <- sapply(mut.df$AA, function(x) sum(mut.df$AA == x,na.rm = TRUE))
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #add $altered_case_id
  cat("adding $altered_case_id","\n")
  mut.df$altered_case_id <- gsub("_",".",mut.df$case_id)
  mut.df$altered_case_id <- gsub("-",".",mut.df$altered_case_id)
  cat("unique case_ids:",length(unique(mut.df$case_id)),"\n")
  cat("unique altered_case_ids:",length(unique(mut.df$altered_case_id)),"\n")
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #add $altered_unique_case_id
  cat("adding $altered_unique_case_id","\n")
  mut.df$altered_unique_case_id <- paste(mut.df$amino_acid_change,mut.df$altered_case_id,sep = "@@")
  #if any duplicates after converting to altered case ids remove them
  if(any(duplicated(mut.df$altered_unique_case_id))) {
    cat("\n####### duplicates detected with altered case ids!","\n")
    mut.df <- mut.df[!(duplicated(mut.df$altered_unique_case_id)),]
    cat("new alterations total:",length(mut.df[,1]),"\n\n")
  }
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #add $case_ID_freq
  cat("adding $case_ID_freq","\n")
  mut.df$case_ID_freq <- sapply(mut.df$altered_case_id, function(x) sum(mut.df$altered_case_id == x))
  #add cancer_type
  cat("adding $cancer_type","\n")
  save.image("troubleshooting_workspace.RData") #####################
  mut.df$cancer_type <- sapply(1:length(mut.df[,1]), function(x) master_case_df$cancer_type[master_case_df$study == mut.df$study[x] & master_case_df$altered_case_id == mut.df$altered_case_id[x]])
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #create mut.study.seqtotals
  cat("creating mut.study.seqtotals","\n")
  mut.study.seqtotals <- data.frame(
    study = all.mut.studies.filtered,
    n_sequenced = as.numeric(sapply(all.mut.studies.filtered, function(x) master_case_list_df$n_cases[master_case_list_df$case_list_id == paste0(x,"_sequenced")])),
    n_altered = sapply(all.mut.studies.filtered, function(x) length(unique(mut.df$altered_case_id[mut.df$study == x]))),
    stringsAsFactors = FALSE
  )
  
  save.image("troubleshooting_workspace.RData") #####################
  
  mut.study.seqtotals$percent_altered <- (mut.study.seqtotals$n_altered / mut.study.seqtotals$n_sequenced )*100
  mut.study.seqtotals <- mut.study.seqtotals[order(-mut.study.seqtotals$percent_altered),]
  #create cancer.type.totals
  cat("creating cancer.type.totals","\n")
  cancer.type.totals <- data.frame(
    cancer_type = unique(master_case_df$cancer_type),
    stringsAsFactors = FALSE 
  )
  
  save.image("troubleshooting_workspace.RData") #####################
  
  cancer.type.totals$n_studies <- sapply(cancer.type.totals$cancer_type, function(x) length(unique(master_case_df$study[master_case_df$study %in% all.mut.studies.filtered & master_case_df$cancer_type == x])))
  cancer.type.totals$n_cases <- sapply(cancer.type.totals$cancer_type, function(x) sum((master_case_df$study %in% all.mut.studies.filtered) & (master_case_df$cancer_type == x),na.rm = TRUE))
  cancer.type.totals$n_altered <- sapply(cancer.type.totals$cancer_type, function(x) length(unique(mut.df$altered_case_id[mut.df$cancer_type == x])))
  cancer.type.totals$percent_altered <- (cancer.type.totals$n_altered / cancer.type.totals$n_cases )*100
  cancer.type.totals <- cancer.type.totals[order(-cancer.type.totals$percent_altered),]
  
  save.image("troubleshooting_workspace.RData") #####################
  
  #print out metadata info
  cat("\n","################### METADATA #####################","\n")
  cat("\n","total alterations so far:",length(mut.df[,1]),"\n")
  #grand total
  cat("\n","total number of cases:","\n")
  print(sum(mut.study.seqtotals$n_sequenced, na.rm = TRUE)) 
  cat("total number of studies with altered GOI:","\n")
  print(length(unique(mut.df$case_id)))
  cat("percent of cases with altered GOI","\n")
  print( 100 * (length(unique(mut.df$case_id))) / (sum(mut.study.seqtotals$n_sequenced, na.rm = TRUE)))    
  
  cat("\n","table of mutation types","\n")
  print(table(mut.df$mutation_type))
  cat("\n")
  
  cat("Number of mutation-data-containing studies (filtered) containing GOI mutations:","\n")
  print(sum(all.mut.studies.filtered %in% sub("_mutations", "",mut.df$genetic_profile_id)))
  cat("out of:","\n")
  print(length(mut.study.seqtotals[,1]))
  cat("\n")
  
  cat("Samples with multiple alterations:","\n")
  print(table(as.data.frame(table(mut.df$case_id))$Freq))
  cat("\n")
  
  cat("Most frequently occuring specific alterations:","\n")
  AAC_freq <- as.data.frame(table(mut.df$amino_acid_change))
  print(head(AAC_freq[order(-AAC_freq$Freq) ,],20))
  cat("\n\n\n")
  
  #assign out final df and save the df name for referencing below
  GOI_cBP_mutations <<- mut.df
  assign("GOI_mut_df_name",paste0(GOI,"_cBP_mutations"))  #for accessing via get() later
  assign(paste0(GOI,"_cBP_mutations"),mut.df)
  assign("GOI_fusion_df_name",paste0(GOI,"_cBP_fusions"))  #for accessing via get() later
  assign(paste0(GOI,"_cBP_fusions"),fusion.df)
  assign(paste0(GOI,"_cBP_study_mut_totals"),mut.study.seqtotals) # won't reference later, dont need name holder
  assign(paste0(GOI,"_cBP_tissue_totals"),cancer.type.totals) # won't reference later, dont need name holder

  save.image("troubleshooting_workspace.RData") #####################
  
############ expression query #################
  cat("################### gene expression query #######################\n\n")
  
  # master_genetic_profile_df contains genetic profiles needed to query expression data, there are several types
  MRNA_EXPRESSION_list <- unique(master_genetic_profile_df[master_genetic_profile_df$genetic_alteration_type == "MRNA_EXPRESSION","id_suffix"])
  # remove zscore profiles
  MRNA_EXPRESSION_list <- MRNA_EXPRESSION_list[!grepl("z",MRNA_EXPRESSION_list,ignore.case = TRUE)]
  #remove mirna
  MRNA_EXPRESSION_list <- MRNA_EXPRESSION_list[!grepl("mirna",MRNA_EXPRESSION_list)]
  #remove mrna_seq_fpkm_polya
  MRNA_EXPRESSION_list <- MRNA_EXPRESSION_list[!grepl("mrna_seq_fpkm_polya",MRNA_EXPRESSION_list)]
  #remove mrna_outliers
  MRNA_EXPRESSION_list <- MRNA_EXPRESSION_list[!grepl("mrna_outliers",MRNA_EXPRESSION_list)]
  
  # generate a list of profiles to query, scrape all to start
  query_list <- master_genetic_profile_df[master_genetic_profile_df$id_suffix %in% MRNA_EXPRESSION_list,"genetic_profile_id"]
  first <- TRUE
  count <- 1
  for (x in query_list) {
    cat("\t",count,"\t\t",x,"\n")
    genProf <- x
    study <- master_genetic_profile_df[master_genetic_profile_df$genetic_profile_id == x,"cancer_study_name"]
    case_list_id <- paste0(study,"_all")
    expr_data <- getProfileData(mycgds,
                                genes = GOI,
                                geneticProfiles = genProf,
                                caseList = case_list_id)
    data <- data.frame(case_id = row.names(expr_data),
                       expression = expr_data[,1],
                       stringsAsFactors = FALSE)
    data$study <- study
    data$genetic_profile_id <- genProf
    data$id_suffix <- master_genetic_profile_df[master_genetic_profile_df$genetic_profile_id == x,"id_suffix"]
    data$manual_tissue <- master_case_df[master_case_df$study == study,"manual_tissue_annotation"][1]
    if(first){
      Data.df <- data
      first <- FALSE
    } else {
      Data.df <- rbind(Data.df,data)
    }
    count <- count + 1
  }
  Data.df$altered_case_id <- gsub("_",".",Data.df$case_id)
  Data.df$altered_case_id <- gsub("-",".",Data.df$altered_case_id)
  Data.df <- Data.df
  Data.df$duplicated_within_suffix <- FALSE
  for(x in unique(Data.df$id_suffix)){
    Data.df[Data.df$id_suffix == x, "duplicated_within_suffix"] <- duplicated(Data.df[Data.df$id_suffix == x, "altered_case_id"])
  }
  rm(first,count,genProf,case_list_id,expr_data,data,study)
  
  save.image("troubleshooting_workspace.RData") #####################
  
  # cat out nuber duplicated and unique for each suffix
  cat("suffix\t\tn_dublicated\t\tn_unique\n")
  for(x in unique(Data.df$id_suffix)) cat(x,
                                          "\t",
                                          sum(duplicated(Data.df[Data.df$id_suffix == x,
                                          "altered_case_id"])),
                                          "\t",
                                          length(unique(Data.df[Data.df$id_suffix == x,
                                          "altered_case_id"])),
                                          "\n"
                                          )
  
  # add expression values to mutation data frame