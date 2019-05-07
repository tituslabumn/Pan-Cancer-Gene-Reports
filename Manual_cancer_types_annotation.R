################### Manual_cancer_types_annotation.R #####################

# manual annotation of cancer types
# not sourced in pipeline run

#keep "study","cancer_type","cancer_type_detailed" cols
master_case_manual_annotation_out <- master_case_df[,2:4]
#convert $study to prefix only for key matching with cancer_types_key_web_API
master_case_manual_annotation_out$prefix <- sapply(master_case_manual_annotation_out$study, function(x) unlist(strsplit(x,split = "_"))[1])
#some prefixes with internal '_' are truncated. These don't overlap anyway, ignore for now
master_case_manual_annotation_out$prefix_key_match <- cancer_types_key_web_API$name[match(master_case_manual_annotation_out$prefix,cancer_types_key_web_API$type_of_cancer_id)]
#filter out duplicate cancer_type/cancer_type_detailed pairs
master_case_manual_annotation_out <- master_case_manual_annotation_out[!duplicated(paste0(master_case_manual_annotation_out$cancer_type,master_case_manual_annotation_out$cancer_type_detailed)),]
#write out
setwd("myosin_pancancer_project/code/Full_single_gene_cBioPortal_pipeline/")
write.table(master_case_manual_annotation_out, file = "cancer_type_manual_annotation_out.tab",sep = "\t", quote = FALSE)

#read in after manual annotation to check data
master_case_manual_annotation_in <- read.delim("cancer_type_manual_annotation_in.tab")
#incorperate into initialization script seperately


