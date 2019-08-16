########### map_cBP_gnomAD_overlap.R ################

cat("##################### label overlap between cBioPortal data and gnomAD variants #########################")

GOI_cBP_mutations$unified_label <- sapply(1:length(GOI_cBP_mutations[,1]), function(x) paste(GOI_cBP_mutations[x,c("unified_pos","reference_allele","variant_allele")],collapse = "@"))
GOI_gnomAD_df_filtered$unified_label <- sapply(1:length(GOI_gnomAD_df_filtered[,1]), function(x) paste(GOI_gnomAD_df_filtered[x,c("unified_pos","Reference","Alternate")],collapse = "@"))
GOI_cBP_mutations$cBPgnomAD_overlap <- GOI_cBP_mutations$unified_label %in% GOI_gnomAD_df_filtered$unified_label
cat("\tcBP total cases represented in gnomAD variants returned: ", sum(GOI_cBP_mutations$cBPgnomAD_overlap),"\n")
cat("\tcBP unique variants represented in gnomAD variants returned: ", sum(GOI_cBP_mutations[!duplicated(GOI_cBP_mutations$amino_acid_change),"cBPgnomAD_overlap"]),"\n\n\n")
GOI_gnomAD_df_filtered$cBPgnomAD_overlap <- GOI_gnomAD_df_filtered$unified_label %in% GOI_cBP_mutations$unified_label
GOI_cBP_mutations[!duplicated(GOI_cBP_mutations$amino_acid_change),"cBPgnomAD_overlap"]