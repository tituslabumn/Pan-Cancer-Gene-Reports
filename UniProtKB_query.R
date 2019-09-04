################# UniProtKB_query.R ########################
#Obtain structural annotation info for GOI from UniProtKB web API


cat("##############################################################\n")
cat("####### Running UniProtKB structural annotation query ########\n")
cat("##############################################################\n\n\n")

#declare query function
parse_uniprotKB_annotation <- function(gene = GOI){
  cat("############# retriving Uniprot/Swissprot ID from biomaRt ############","\n\n")
  
  GOI_uniprot_id <<- as.character(getBM(attributes = "uniprotswissprot", filters = "ensembl_gene_id", values = GOI_ENSG, mart = hs_ensembl))
  # in rare cases biomaRt returns emty logical for uniprotID. This happens for MYH16 regardless of genome/mart version.
  if(GOI_uniprot_id == "logical(0)"){
    cat("\t########### no UniProtID returned by biomaRt ############\n\t\tAttempting to retrive from bioDBnet instead\n\tIDs mapped:\n")
    # rjson:: returns error likely due to bug? use jsonlite:: to query bioDBnet
    library(jsonlite)
    bioDBnet_UniProtIds <- jsonlite::fromJSON(paste0("https://biodbnet-abcc.ncifcrf.gov/webServices/rest.php/biodbnetRestApi.json?method=db2db&format=row&input=genesymbol&inputValues=",GOI,"&outputs=UniProtAccession&taxonId=9606"))
    print(bioDBnet_UniProtIds$`UniProt Accession`)
    GOI_uniprot_id <<- bioDBnet_UniProtIds$`UniProt Accession`[1]
  }
  
  cat("\tUniProtKB ID:",GOI_uniprot_id,"\n\n")
  cat("Querying UniProtKB GOI page via RCurl","\n\n")
  # retrieve features with some filetering
  
  # old method produced SSL connect errors; use RCurl::getURL
  # still an issue, seems to just require persistant attempts?
  #   read.df <- read.delim(paste0("https://www.uniprot.org/uniprot/",GOI_uniprot_id,".txt"),stringsAsFactors = FALSE) 
  #   read.delim(textConnection(getURL(UniProt_URL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)),stringsAsFactors = FALSE)
  UniProt_URL <- paste0("https://www.uniprot.org/uniprot/",GOI_uniprot_id,".txt") #returns data frame with one col of lines
  cat("\tURL: ",UniProt_URL,"\n\n")
  library(RCurl)
  # attempt to access url, if fail try up to 19 more times every 15s
  attempt<-1
  success <- FALSE
  while(attempt < 20 & success == FALSE){
    read.df <- tryCatch({
      success <- TRUE
      #read.delim(textConnection(getURL(UniProt_URL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)),stringsAsFactors = FALSE)
      read.delim(url(UniProt_URL),stringsAsFactors = FALSE)
      },error=function(cond){
      print(cond)
      cat("\tAttempt",attempt,"failed. Attemting ",20-attempt," more times.\n")
      success <<- FALSE
      attempt<<-attempt+1
      if(attempt == 20) stop("Unable to connect with UniProt; try again later :[")
      Sys.sleep(15)
    })
  }
  rm(attempt,success)
  Uniprot_txt_parsed <<- read.df
  cat("Succesffully retrieved\n\n")

  #colnames of this parsed txt file contains AA length <- extract for future use
  GOI_UNIPROT_AA_LENGTH <<- as.numeric(rev(unlist(strsplit(colnames(read.df),split = "\\.+")))[2])
  cat("\tAA length parsed:",GOI_UNIPROT_AA_LENGTH,"\n\n")
  #subset for feature lines that are not secondary stucture elements
  cat("\t","Filtering out secondary structure elements","\n")
  read.df.ft <- as.data.frame(read.df[grepl("FT {3}\\S",read.df[,1]) & !grepl("PDB:",read.df[,1]),],stringsAsFactors = FALSE) 
  #filter out conflict and variant rows (switch to list from df)
  read.df.ft <- as.character(read.df.ft[!grepl("FT   VARIANT",read.df.ft[,1]) & !grepl("FT   CONFLICT",read.df.ft[,1]),]) 
  ft.col.names <- unique(sapply(read.df.ft, function(x) unlist(strsplit(x,"\\s+"))[2]))
  cat("Feature names:","\n")
  print(ft.col.names)
  #subseting with "FT {3}\\S" omits some continuations of lines that are prefixed by "FT {>3}". 
  #Append continued lines in those cases ( were the line does not end in "." ("\\.$") or contain ". " ("\\.\\s") just in cast there is internal period)
  truncated_line_index <- !grepl("\\.$|\\.\\s",read.df.ft)
  if(sum(truncated_line_index)>0){
    cat("\t","Fixing truncated lines","\n")
    read.df.ft[truncated_line_index] <- sapply(read.df.ft[truncated_line_index],
                                       function(x){
                                         if(length(unlist(strsplit(x,"\\s+"))) > 4){ #skip if no label (also lacks terminal period)
                                           holder <- paste0(x," ",sub("FT\\s+","\\1",read.df[(which(read.df[,1] == x)+1),1]))
                                           if(grepl("\\.",holder)){ #may need to add third line
                                             return(holder)
                                           }else{
                                              holder <- paste0(holder," ",sub("FT\\s+","\\1",read.df[(which(read.df[,1] == x)+2),1]))
                                              if(!grepl("\\.",holder)) cat("\t\tWARNING - no terminal period; still truncated after merging 3 lines","\n") 
                                              return(holder)
                                           }
                                        }else{
                                          cat("\t\tWARNING - Line contains no label:",x,"\n") 
                                          return(x)
                                        }
                                       })
  }
 
  cat("\nCreating feature data frame","\n")
  #create feature_df
  feature_df <- data.frame(
    TYPE = character(length(read.df.ft)),
    AA_start = integer(length(read.df.ft)),
    AA_end = integer(length(read.df.ft)),
    LABEL = character(length(read.df.ft)),
    stringsAsFactors = FALSE
  )
  for (x in 1:length(read.df.ft)) {
    cat("\t",read.df.ft[x],"\n")
    line.split <- unlist(strsplit(read.df.ft[x],split = "\\s+")) #first element is 'FT'
    feature_df[x,"TYPE"] <- line.split[2]
    feature_df[x,"AA_start"] <- as.integer(line.split[3])
    feature_df[x,"AA_end"] <- as.integer(line.split[4])
    if(!is.na(line.split[5])){
      feature_df[x,"LABEL"] <- sub(paste0(".* ",line.split[4]," *(.*?) *\\..*"),"\\1",read.df.ft[x])
    }else{
      feature_df[x,"LABEL"] <-line.split[2] #if no label provided (i.e. DNA_BIND) use the TYPE as the LABEL
    }
  }
  cat("\nCreating feature metadata data frame","\n")
  domain_annotation_metadata <- data.frame(stringsAsFactors = FALSE)
  for (x in ft.col.names) {
    domain_annotation_metadata[1,x]<- sum(feature_df$TYPE == x)
  }
  
  #write out data
  cat("Writing out data","\n")
  assign("GOI_uniprot_id",GOI_uniprot_id, envir = .GlobalEnv)
  assign(paste0(gene,"_protein_feature_annotation"),feature_df, envir = .GlobalEnv)
  assign(paste0(gene,"_protein_feature_annotation_metadata"),domain_annotation_metadata, envir = .GlobalEnv)
  
  cat("######## UniProtKB query complete ##########","\n\n\n")
}

#query GOI
parse_uniprotKB_annotation()

save.image("troubleshooting_workspace.RData") #####################
