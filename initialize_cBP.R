################# initialize_cBP.R ########################

#Initialize cBP workspace

cat("\n#######################################################","\n")
cat("############### Initialize cBP workspace ##############","\n")
cat("#######################################################","\n\n")

#check that required packages are installed and install if needed
source("Install_required_packages.R")

cat("Establish metadata data frame for final output","\n\n")
initial_metadata <- data.frame(
  total_studies = 0,
  total_cases = 0,
  studies_without_genetic_profiles = "",
  studies_without_case_lists = "",
  stringsAsFactors = FALSE
)


cat("Load 'rtacklayer' and check for chain to convert cBP data to GRCh38.p12","\n\n")
#cBioPortal uses genome version hg19/GRCh37 while biomaRt(Ensembl) is currently using GRCh38.p12
#for converting:
#library(rtracklayer)
    #download chain file from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/
    # hg19ToHg38.over.chain.gz
    #place in this directory
    # chain_path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
    Chain_19to38 <- import.chain("hg19ToHg38.over.chain") # used in cBP mutation query function below
cat("\n\n")

cat("Fetching canonical transcript key used by cBioPortal","\n\n\n")
cBP_canonical_transcripts <- read.delim("https://raw.githubusercontent.com/mskcc/vcf2maf/master/data/isoform_overrides_uniprot",stringsAsFactors = FALSE, fileEncoding = "latin1")
row.names(cBP_canonical_transcripts)<- cBP_canonical_transcripts$gene_name
colnames(cBP_canonical_transcripts)[1] <- "ensembl_transcript_id"

save.image("troubleshooting_workspace.RData") #####################

################################################################################
############### retrieve study/sample annotation data from cBioPortal ##########
################################################################################

#Construct a CGDS connection object, required for most CGDSR querries
mycgds <- CGDS("http://www.cbioportal.org/")

#test URL
cat("testing cBP connection object","\n")
test(mycgds)
cat("\n\n")

#get list of all cancer studies in cBioPortal database
#n=240 (7,406 total 'cases')
cat("Retrieving list of all cancer studies","\n")
all.studies <- getCancerStudies(mycgds)
cat("Totoal number of studies returned by cBP:","\t",length(all.studies[,1]),"\n\n")

initial_metadata$total_studies <- length(all.studies[,1])

#old version of this code contained manual annotations for tissue type
#this is automated now downstream via getClinicalData()

cat("Generating master annotation data frame for all studies and their respective genetic-profiles","\n")
#make a master dataframe with all genetic profiles for each study as returned by getGeneticProfiles()
    master_genetic_profile_df <- data.frame(
      cancer_study_name = character(),
      genetic_profile_id = character(),
      genetic_profile_name = character(),
      genetic_profile_description = character(),
      cancer_study_id = character(),
      genetic_alteration_type = character(),
      show_profile_in_analysis_tab = character(),
      stringsAsFactors = FALSE
    )
    
    #for each study append all its genetic profiles to the master df 
    y<-1
    null_studies <- character()
    for (x in all.studies$cancer_study_id) {
      cat("[",y,"/",length(all.studies$cancer_study_id),"]\t",x)
      genetic_profile <- getGeneticProfiles(mycgds,x)
      #if cBP does not return any genetic profiles skip apending this study
      if(length(genetic_profile[,1]) == 0){
        cat("\t","### WARNING: no genetic profiles returned ###")
        null_studies <- c(null_studies,x)
      }else{
        master_genetic_profile_df <- rbind(master_genetic_profile_df,cbind(data.frame(cancer_study_name = x),genetic_profile))
      }
      y<-y+1
      cat("\n")
    }
    #add null studies to metadata as space delim char string
    initial_metadata$studies_without_case_lists <- paste(null_studies,collapse = " ")
    rm(x,y,genetic_profile)
    
    save.image("troubleshooting_workspace.RData") #####################
    
    #re-order colums and remove "show_profile_in_analysis_tab"
    master_genetic_profile_df <- master_genetic_profile_df[,c("cancer_study_id","cancer_study_name","genetic_profile_id","genetic_alteration_type","genetic_profile_name","genetic_profile_description")]
    #set row names to the case list/genetic profile id
    row.names(master_genetic_profile_df) <- master_genetic_profile_df$genetic_profile_id
    #make col of genetic_profile_id suffixes
    master_genetic_profile_df$id_suffix <- sapply(row.names(master_genetic_profile_df), function(x) gsub(paste0(master_genetic_profile_df[x,"cancer_study_name"],"_"),"",master_genetic_profile_df[x,"genetic_profile_id"]))
    #not sure how this happened but fix:
    master_genetic_profile_df$cancer_study_name <- as.character(master_genetic_profile_df$cancer_study_name)
    
    save.image("troubleshooting_workspace.RData") #####################
    
cat("\n\n","Generating master annotation data frame for all studies and their respective case-lists","\n")
#make a master dataframe with all case lists for each study as returned by getCaseLists()    
    master_case_list_df <- data.frame(
      cancer_study_name = character(),
      case_list_id = character(),
      case_list_name = character(),
      cancer_study_id = character(),
      case_ids = character(),
      stringsAsFactors = FALSE
    )
    
    #for each study append all its genetic profiles to the master df 
    y<-1
    null_studies <- character()
    for (x in all.studies$cancer_study_id) {
      cat("[",y,"/",length(all.studies$cancer_study_id),"]\t",x)
      case_list <- getCaseLists(mycgds,x)
      #if cBP does not return any case lists skip apending this study
      if(length(case_list[,1]) == 0){
        cat("\t","### WARNING: no case lists returned ####")
        null_studies <- c(null_studies,x)
      }else{
        master_case_list_df <- rbind(master_case_list_df,cbind(data.frame(cancer_study_name = x),case_list))
      }
      y<-y+1
      cat("\n")
    }
    #add null studies to metadata as space delim char string
    initial_metadata$studies_without_case_lists <- paste(null_studies,collapse = " ")
    rm(x,y,case_list,null_studies)
    
    save.image("troubleshooting_workspace.RData") #####################
    
    #re-order colums and remove "show_profile_in_analysis_tab"
    master_case_list_df <- master_case_list_df[,c("cancer_study_id","cancer_study_name","case_list_id","case_list_name","case_list_description","case_ids")]
    #set row names to the case list/genetic profile id
    row.names(master_case_list_df) <- master_case_list_df$case_list_id
    #add col with suffix of case list id ('genetic profile' type)
    master_case_list_df$id_suffix <- sapply(row.names(master_case_list_df), function(x) gsub(paste0(master_case_list_df[x,"cancer_study_name"],"_"),"",master_case_list_df[x,"case_list_id"]))
    #add col with number of cases
    master_case_list_df$n_cases <- sapply(master_case_list_df$case_ids, function(x) length(unlist(strsplit(x,split = " "))))
    
    initial_metadata$total_cases <- sum(master_case_list_df$n_cases)
    
    cat("\n\n\n")

    save.image("troubleshooting_workspace.RData")  #####################   
    
#make master case id data frame (derived from each 'xxxxxxx_all'  case id list)
cat("########### Making master case id data frame ##############","\n\n")
all_clindata_colnames <- character()
skipped_no_clindata_list <- character()
y<-1
z<- paste0(all.studies$cancer_study_id,"_all")
for (x in z) {
  cat("[",y,"/",length(z),"]","\t",x,"\t")
  
  #in some case_lists getClinicalData seems to be broken. Query from web API in all cases
      # getClinicalData(mycgds,x) no longer used
      # will produce errors in some at least one case if fileEncoding arg is not specified
      clindata <- read.delim(paste0('http://www.cbioportal.org/webservice.do?cmd=getClinicalData&case_set_id=',x),stringsAsFactors = FALSE, fileEncoding = "latin1")
      cat("retrived data frame with the following dimensions",dim(clindata),"\n")
      #remove any duplicated rown names; has yet to be called for full cBP data
      if(any(duplicated(clindata$CASE_ID))){
        cat("removing duplicate case ids")
        clindata <- clindata[-duplicated(clindata$CASE_ID)]
      }
      row.names(clindata) <- clindata$CASE_ID

  #skip if no clinical data is returned (currently true for crc_msk_2017_all)
  if(length(clindata[1,])<1){
      cat("################ No clinical data; skipping... #######################\n")
      assign("skipped_no_clindata_list",c(skipped_no_clindata_list,x), envir = .GlobalEnv)
  }else{
    #append colnames before populating empty cols below
    all_clindata_colnames <- c(all_clindata_colnames,colnames(clindata))
    
    #populate empty lists for instances of missing colnames (data that is not available for all studies)
    for (xx in c("CANCER_TYPE","CANCER_TYPE_DETAILED","MUTATION_COUNT","AGE","SEX","GENDER","RACE","ETHNICITY","FRACTION_GENOME_ALTERED","OS_STATUS","OS_MONTHS","SAMPLE_TYPE")) {
      if(!(xx %in% colnames(clindata))) clindata[,xx] <- NA
    }  
    
    #apend data to master df
    current_case_df <- data.frame(
      case_id = row.names(clindata), 
      study = rep(sub("_all","",x),length(row.names(clindata))),
      sample_count = clindata$SAMPLE_COUNT,
      #below cols: data not available for all
      cancer_type = clindata$CANCER_TYPE, 
      cancer_type_detailed = clindata$CANCER_TYPE_DETAILED, 
      age = clindata$AGE,
      sex = clindata$SEX,
      gender = clindata$GENDER,
      race = clindata$RACE,
      ethnicity = clindata$ETHNICITY,
      fraction_genome_altered = clindata$FRACTION_GENOME_ALTERED,
      mutation_count = clindata$MUTATION_COUNT,
      os_status = clindata$OS_STATUS,
      os_months = clindata$OS_MONTHS,
      sample_type = clindata$SAMPLE_TYPE,
      stringsAsFactors = FALSE
    )
    #for cases missing cancer_type or cancer_type_detailed populate with study for unique manual tissue type annotation
    if(anyNA(current_case_df$cancer_type) | anyNA(current_case_df$cancer_type_detailed)){
      current_case_df$cancer_type <- current_case_df$study[1]
      current_case_df$cancer_type_detailed <- current_case_df$study[1]
    }
    #some instances have empty char for cancer_type and cancer_type_detailed but not for entire study, replace these with study name for manual tissue type annotation
    if(any(current_case_df$cancer_type == "")) current_case_df[current_case_df$cancer_type == "","cancer_type"] <- current_case_df[current_case_df$cancer_type == "","study"]
    if(any(current_case_df$cancer_type_detailed == "")) current_case_df[current_case_df$cancer_type_detailed == "","cancer_type_detailed"] <- current_case_df[current_case_df$cancer_type_detailed == "","study"]
    
    if(y == 1){
      master_case_df <- current_case_df
    } else {
      master_case_df <- rbind(master_case_df,current_case_df)
    }
  }
  y<-y+1
}
rm(x,y,xx,z,clindata,current_case_df)

save.image("troubleshooting_workspace.RData") #####################

cat("\n\n\n")

cat("altering case id dashes and underscores\n\n")
#some studies share case ids but are encoded differently where '.' has been swapped for '_'
  #this may only be for cases where clinical data was obtained from web API
  #add col with case ideas with all '_'s and '-'s changed to '.'
  master_case_df$altered_case_id <- gsub("_",".",master_case_df$case_id)
  master_case_df$altered_case_id <- gsub("-",".",master_case_df$altered_case_id)

  cat("################ Annotating cancer type #######################\n\n")
#get cancer type key (not avaliable through CGDSR, only via web API).
#this returns what is cancer_type_detailed in getClinicalData query
cancer_types_key_web_API <- read.delim("http://www.cbioportal.org/webservice.do?cmd=getTypesOfCancer",stringsAsFactors = FALSE, fileEncoding = "latin1")
#use this key to fill in empty ("") cancer typed in master case df
  master_case_df$webAPI_cancer_type <-  cancer_types_key_web_API$name[match(sapply(master_case_df$study, function(x) unlist(strsplit(x,split = "_"))[1]),cancer_types_key_web_API$type_of_cancer_id)]  

save.image("troubleshooting_workspace.RData")   #####################

cat("creating/writing manual annotation template tab delim file \n\n")
#creatE manual annotation file to be manually annotated
    #keep "study","cancer_type","cancer_type_detailed" cols
    master_case_manual_annotation_out <- master_case_df[,c("study","cancer_type","cancer_type_detailed")]
    #convert $study to prefix only for key matching with cancer_types_key_web_API
    master_case_manual_annotation_out$prefix <- sapply(master_case_manual_annotation_out$study, function(x) unlist(strsplit(x,split = "_"))[1])
    #some prefixes with internal '_' are truncated. These don't overlap anyway, ignore for now
    master_case_manual_annotation_out$prefix_key_match <- cancer_types_key_web_API$name[match(master_case_manual_annotation_out$prefix,cancer_types_key_web_API$type_of_cancer_id)]
    #filter out duplicate cancer_type/cancer_type_detailed pairs
    master_case_manual_annotation_out <- master_case_manual_annotation_out[!duplicated(paste0(master_case_manual_annotation_out$cancer_type,master_case_manual_annotation_out$cancer_type_detailed)),]
    #write out
    write.table(master_case_manual_annotation_out, file = "cancer_type_manual_annotation_out.tab",sep = "\t", quote = FALSE)    
#read in manual tissue type annotations
  # this is a tab delim file created via Manual_cancer_types_annotation.R
  # a manual tissue type (32 total types) is assigned for each unique cancer_type/cancer_type_detailed pair from master_case_df

cat("reading in manual annotation file \n\n")
if(!("cancer_type_manual_annotation_in.tab" %in% list.files())){
  stop("'cancer_type_manual_annotation_in.tab' file not found:\n\tmanually annotate 'cancer_type_manual_annotation_out.tab' file \n\tadd final_tissue column \n\tsace as 'cancer_type_manual_annotation_in.tab'\n\trun again\n\n\n")
}
master_case_manual_annotation_in <- read.delim("cancer_type_manual_annotation_in.tab")
all.studies.old <- read.delim("past_manual_annotation_all_studies.tab")

cat("\tchecking that this file is up to date... \n\n")
if(length(master_case_manual_annotation_in[,1]) != length(master_case_manual_annotation_out[,1])) {
  cat("new studies in database with no manual annotation:\n")
  print(all.studies$cancer_study_id[which(!(all.studies$cancer_study_id %in% all.studies.old$cancer_study_id))])
  stop("manual annotation input and output files do not have equal number of rows:\n\tmanually annotate 'cancer_type_manual_annotation_out.tab' file \n\tadd final_tissue column \n\tsace as 'cancer_type_manual_annotation_in.tab'\n\trun again\n\n\n")
}

#overwrite past_manual_annotation_all_studies.tab if manual annotation passes check
write.table(all.studies, file = "past_manual_annotation_all_studies.tab",sep = "\t", quote = FALSE) 

#map these types to master_case_df
#for cases that are new (not yet annotated) populate with NA
#use _out data frames for pasted keys as manual annotation altered some escape characters and leads to small number of mismatch (i.e. "Head and Neck CancerHead and Neck Squamous Cell Carcinoma\xc3?" vs "Head and Neck CancerHead and Neck Squamous Cell Carcinoma�?")
master_case_manual_annotation_in$match_key <- paste0(master_case_manual_annotation_in$cancer_type,master_case_manual_annotation_in$cancer_type_detailed)
master_case_df$manual_tissue_annotation <- NA
master_case_df$manual_tissue_annotation <- master_case_manual_annotation_in$final_tissue[match(paste0(master_case_df$cancer_type,master_case_df$cancer_type_detailed),master_case_manual_annotation_in$match_key)]

save.image("troubleshooting_workspace.RData")   #####################

#make df table of all clinical data colnames
all_clinical_data_colnames_table <- as.data.frame(table(all_clindata_colnames))
all_clinical_data_colnames_table <- all_clinical_data_colnames_table[order(-all_clinical_data_colnames_table$Freq),]
#make df table of all cancer types returned by getClinicalData() and filled in by web API key
all_cancer_types_table <- as.data.frame(table(master_case_df$cancer_type))
all_cancer_types_table <- all_cancer_types_table[order(-all_cancer_types_table$Freq),]
#make df table of all tissue types from manual annotation
all_tissue_types_table <- as.data.frame(table(master_case_df$manual_tissue_annotation))
all_tissue_types_table <- all_tissue_types_table[order(-all_tissue_types_table$Freq),]

save.image("troubleshooting_workspace.RData") #####################
     
cat("########### Identifying overlaping studies ##############","\n\n")       
    #add col to all.studies signifying if study has overlapping case Ids with any other
    all.studies$unique_ids <- sapply(all.studies$cancer_study_id, function(x) !(sum(master_case_df$altered_case_id[master_case_df$study == x] %in% master_case_df$altered_case_id[master_case_df$study != x])))
    cat("Number of studies that have cases that overlap with any other study: ", sum(!(all.studies$unique_ids))," out of ",length(all.studies$cancer_study_id)," total\n")
    #add col with ':::' delim list of which studies are overlapping
    all.studies$overlap <- NA
    all.studies$overlap <- sapply(all.studies$cancer_study_id, function(x) paste(unique(master_case_df$study[(master_case_df$altered_case_id %in% master_case_df$altered_case_id[master_case_df$study == x]) & (master_case_df$study != x)]),collapse = ":::"))
    #add col with number of cases
    all.studies$n_cases <- sapply(all.studies$cancer_study_id, function(x) sum(master_case_df$study == x))
    #add col for number of cases represented in overlap (are there instances where overlap is incomplete?)
    all.studies$n_overlaped_cases <- sapply(all.studies$cancer_study_id, function(x) sum(master_case_df$altered_case_id[master_case_df$study == x] %in% master_case_df$altered_case_id[master_case_df$study != x]))
    cat("Number of unique case ids: ", length(unique(master_case_df$altered_case_id))," out of ",length(master_case_df$altered_case_id)," total case ids\n")
    initial_metadata$n_unique_case_ids <- length(unique(master_case_df$altered_case_id))
#get df of tabulated genetic profile names    
    genetic.profile.names.table <- as.data.frame(table(master_genetic_profile_df$genetic_profile_name))
    genetic.profile.names.table <- genetic.profile.names.table[order(-genetic.profile.names.table$Freq),c(2,1)]
    colnames(genetic.profile.names.table)[2] <- "Genetic profile name"
    cat(" ########### There are ",length(genetic.profile.names.table$Freq)," unique genetic profile names:","\n")
    print(genetic.profile.names.table,right = FALSE, row.names = FALSE)
    
    cat("\n\n\n")

    save.image("troubleshooting_workspace.RData")   #####################  
    
#get df of tabulated case list IDs    
    case.list.id.suffix.table <- as.data.frame(table(master_case_list_df$id_suffix))
    case.list.id.suffix.table <- case.list.id.suffix.table[order(-case.list.id.suffix.table$Freq),c(2,1)]
    colnames(case.list.id.suffix.table)[2] <- "case list ID suffix"
    cat(" ########### There are ",length(case.list.id.suffix.table$Freq)," unique case list ID suffixes:","\n")
    print(case.list.id.suffix.table,right = FALSE, row.names = FALSE)
    
    cat("\n\n\n")

    save.image("troubleshooting_workspace.RData") #####################
    
#get list of studies with some known overlaps filtered out
    cat("applying some initial filtering (only keeping provisional tcga data)","\n\n")
    #tcga data is known to be reflected in two published and one provisional dataset
    #refer to cBP website FAQ for the distinction
    #keep only provisional (most extensive) datasets for the purposes of this analysis
    filtered_studies <- all.studies$cancer_study_id[!(grepl("tcga_pan",all.studies$cancer_study_id) | grepl("tcga_pub",all.studies$cancer_study_id))]
    
    save.image("troubleshooting_workspace.RData") #####################
        
#get list of all studies that have mutation data (all studies with case id list with 'mutations' suffix)       
all.mut.studies <- master_genetic_profile_df$cancer_study_name[master_genetic_profile_df$id_suffix == 'mutations']
all.mut.studies.filtered <- all.mut.studies[all.mut.studies %in% filtered_studies]
cat("number of studies (filtered) with mutation data available:",length(all.mut.studies.filtered),"\n\n\n")

save.image("troubleshooting_workspace.RData") #####################



save.image("initialized_workspace_cBP.RData")