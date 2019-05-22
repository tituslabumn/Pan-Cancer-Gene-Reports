################# ExAC_query.R ########################
#Obtain all Exome variants for GOI found in the exome aggregation consortium

cat("#####################################################################\n")
cat("############# Retrieving all ExAC variants for GOI ##################\n")
cat("#####################################################################\n\n\n")

cat("retrieving ExAC variants from RESTful API; returns JSON page\n\n")

library("rjson")
#don't worry about warning about final line
ExAC_json_data <- fromJSON(paste(readLines(paste0("http://exac.hms.harvard.edu/rest/gene/variants_in_gene/",GOI_ENSG_GRCh37)), collapse=""))

save.image("troubleshooting_workspace.RData") #####################
cat("\t\tVariants:",length(ExAC_json_data),"\n\n")

cat("\tConverting complex list to data frame\n")
#parse complex list into data.frame
#note: pos refers to start of variant in case of del
ExAC_df <- data.frame(t(unlist(ExAC_json_data[[1]])),stringsAsFactors = FALSE)
rownames(ExAC_df)<- ExAC_df$variant_id
for(x in 2:length(ExAC_json_data)){
  current_row <- unlist(ExAC_json_data[[x]])
  for (y in names(current_row)) {
    ExAC_df[current_row["variant_id"],y]<- current_row[y]
  }
}
rm(x,y,current_row)

cat("\tFiltering out variants that lack 'PASS' annotation\n")
ExAC_df <- ExAC_df[ExAC_df$filter == "PASS",]
cat("\t\tFiltered variants:",length(ExAC_df[,1]),"\n\n")

save.image("troubleshooting_workspace.RData") #####################

cat("\tConverting coordinates from hg19 to hg38\n")
#this returns coordinates aligned to GRCh37/hg19; change to hg38
ExAC_df$pos <- unlist(start(liftOver(GRanges(paste0('chr',GOI_CHR), IRanges(start = as.numeric(ExAC_df$pos), width = 1)),Chain_19to38)))

save.image("troubleshooting_workspace.RData") #####################

cat("######## ExAC query complete #########\n\n\n")