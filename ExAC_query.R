################# ExAC_query.R ########################
#Obtain all Exome variants for GOI found in the exome aggregation consortium

cat("#####################################################################\n")
cat("############# Retrieving all ExAC variants for GOI ##################\n")
cat("#####################################################################\n\n\n")

library("rjson")
json_data <- fromJSON(paste(readLines(paste0("http://exac.hms.harvard.edu/rest/gene/variants_in_gene/",GOI_ENSG)), collapse=""))

save.image("troubleshooting_workspace.RData") #####################

cat("######## ExAC query complete #########\n\n\n")