########### map_cBP_gnomAD_overlap.R ################

cat("##################### label overlap between cBioPortal data and gnomAD variants #########################\n\n")

GOI_cBP_mutations$unified_label <- sapply(1:length(GOI_cBP_mutations[,1]), function(x) paste(GOI_cBP_mutations[x,c("unified_pos","reference_allele","variant_allele")],collapse = "@"))
GOI_gnomAD_df_filtered$unified_label <- sapply(1:length(GOI_gnomAD_df_filtered[,1]), function(x) paste(GOI_gnomAD_df_filtered[x,c("unified_pos","Reference","Alternate")],collapse = "@"))
GOI_cBP_mutations$cBPgnomAD_overlap <- GOI_cBP_mutations$unified_label %in% GOI_gnomAD_df_filtered$unified_label
cat("\tcBP total cases represented in gnomAD variants returned: ", sum(GOI_cBP_mutations$cBPgnomAD_overlap),"\n")
cat("\tcBP unique variants represented in gnomAD variants returned: ", sum(GOI_cBP_mutations[!duplicated(GOI_cBP_mutations$amino_acid_change),"cBPgnomAD_overlap"]),"\n\n\n")
GOI_gnomAD_df_filtered$cBPgnomAD_overlap <- GOI_gnomAD_df_filtered$unified_label %in% GOI_cBP_mutations$unified_label
GOI_cBP_mutations[!duplicated(GOI_cBP_mutations$amino_acid_change),"cBPgnomAD_overlap"]

# add filter threshold to both data.frames
# for now set threshold for filtering (variants more freqent than this value are filtered out of cBP data) at 1 in 10,000 0.0001
GOI_gnomAD_df_filtered$overlap_filter_threshold <- FALSE
GOI_gnomAD_df_filtered[GOI_gnomAD_df_filtered$cBPgnomAD_overlap & GOI_gnomAD_df_filtered$Allele.Frequency >= 0.0001 ,"overlap_filter_threshold"] <- TRUE
GOI_cBP_mutations$overlap_filter_threshold <- FALSE
GOI_cBP_mutations[GOI_cBP_mutations$cBPgnomAD_overlap & GOI_cBP_mutations$unified_label %in% GOI_gnomAD_df_filtered$unified_label[GOI_gnomAD_df_filtered$overlap_filter_threshold],"overlap_filter_threshold"] <- TRUE

