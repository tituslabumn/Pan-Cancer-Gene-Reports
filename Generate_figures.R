################# Generate_figures.R ########################
#generate figures for output report


cat("#########################################################################\n")
cat("############### Generating figure data for output report ################\n")
cat("#########################################################################\n\n\n")

library("trackViewer")
library(markdown)
library(knitr)
library(RColorBrewer)
#make combined qualitative color palette for max of 21 distinct features repeaded once (total = 42)
combined_qualitative_palette <- c(brewer.pal(9,"Set1"),brewer.pal(12,"Set3"),brewer.pal(9,"Set1"),brewer.pal(12,"Set3"))

cat("parse in some data to non-gene-specific names\n")

#parse some data frames to non-gene-specific names
GOI_cBP_fusions <- get(paste0(GOI,"_cBP_fusions"))

cat("declare function to resolve overlapping features\n")
#make function to segregate overlapped features between rows
#acepts filtered data frame with AA_start and AA_end cols ; returns list of layer numbers
sep_overlap_features <- function(df){
  if(length(df[,1]) <= 1){
    return(1)
  }else{
    layer_id <- integer(length = length(df[,1]))
    layer_id[1] <- 1 #first feature is always layer 1
    for(x in 2:length(df[,1])){ 
      layer_check <- 1
      while(any( (df$AA_end[layer_id == layer_check] >= df$AA_start[x]) & (df$AA_start[layer_id == layer_check] <= df$AA_end[x]) )){
        #overlaps with at leat on in that layer
        layer_check <- layer_check + 1
      }
      layer_id[x] <- layer_check
    }
    #swap layer IDs so top layer is on top
    layer_id <- match(layer_id,max(layer_id):1)
    return(layer_id)
  }
}

save.image("troubleshooting_workspace.RData") #####################

cat("declare function to truncate labels\n")
#declare function that accepts full or subseted feature data frame and returns list of labels
#number of characters to limit to as char_number_cap argument
truncate.feature.labels <- function(df,char_number_cap){
  #truncate
  df$LABEL <- sapply(df$LABEL, function(x){
    if(nchar(x) > char_number_cap){
      truncating_pattern <- paste0("(.{",char_number_cap,"}).*$")
      return(paste0(sub(truncating_pattern,"\\1",x),"..."))
    }else{
      return(x)
    }
  })
  #sandwich with TYPE and position range
  return(paste0(df$TYPE,"-",df$LABEL,"[",df$AA_start,"-",df$AA_end,"]"))
}

cat("Figure1\n")
#Figure 1 - features; domains,regions, DNA_bind, MOTIF
F1_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION","DNA_BIND","MOTIF"),]
#F1_site_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("METAL","SITE","MOD_RES","CROSSLNK"),]
if(length(F1_feature_df[,1])>0){
  F1_feature_labels <- truncate.feature.labels(F1_feature_df,20)
  F1_features <- GRanges(seqnames = "chr", IRanges(start = F1_feature_df$AA_start, end = F1_feature_df$AA_end, names = F1_feature_labels))
  F1_features$height <- 0.06
  F1_features$fill <- combined_qualitative_palette[1:length(F1_features)]
  F1_features$featureLayerID <- sep_overlap_features(F1_feature_df)
  
  #F1_variants <- GRanges(seqnames = "chr", IRanges(start = F1_site_df$AA_start, width = 1), label = as.character(1:length(F1_site_df[,1])))
  #F1_variants$label.col <- "black"
  F1_variants <- GRanges()
  
  F1_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
  F1_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
  
  #lolliplot(SNP.gr = F1_variants ,features = F1_features, ranges = F1_ranges, ylab = FALSE, xaxis = F1_x_axis)
}
#seperately plot a legend
# par(mar = c(0,0,0,0))
# plot.new()
# legend("center", legend = F1_feature_labels, fill = F1_features$fill)

cat("Figure2\n")

#Figure 2 - all transcripts plotted to relative transcript position
F2_feature_df <- GOI_exon_annotation
  F2_features <- GRanges(seqnames = "chr", IRanges(start = F2_feature_df$relative_union_start, end = F2_feature_df$relative_union_end, names = NULL))
  F2_features$height <- 0.02
  F2_features$fill <- "grey"
  F2_features$fill[F2_feature_df$ensembl_transcript_id == GOI_TRANSCRIPT] <- "red"
  F2_features$featureLayerID <- F2_feature_df$ensembl_transcript_id
  F2_variants <- GRanges()
  F2_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
  F2_x_axis <- round(seq(from = 1, to = max(c(F2_feature_df$relative_union_start,F2_feature_df$relative_union_end)), length.out = 10),-1) #even split by 5 rounded to nearest 10
  #lolliplot(SNP.gr = F2_variants ,features = F2_features, ylab = FALSE, xaxis = F2_x_axis)


save.image("troubleshooting_workspace.RData") #####################  
  
#get master case df with only studies that survived cBP mut query filtering
all_mut_cases <- master_case_df[master_case_df$study %in% unique(GOI_cBP_mutations$study),]
GOI_cBP_mutations$manual_tissue <- all_mut_cases$manual_tissue_annotation[match(paste0(GOI_cBP_mutations$study,GOI_cBP_mutations$altered_case_id),paste0(all_mut_cases$study,all_mut_cases$altered_case_id))]

cat("Tissue table\n")

#tissue type enrichment
colnames(all_tissue_types_table) <- c("Primary Tissue","Total cases in database")
all_tissue_types_table$`Total cases sequenced` <- sapply(all_tissue_types_table$`Primary Tissue`, function(x) sum(all_mut_cases$manual_tissue_annotation == x,na.rm = TRUE))
GOI_cBP_mutations$manual_tissue <- as.character(GOI_cBP_mutations$manual_tissue)
all_tissue_types_table$`Total altered cases` <- sapply(all_tissue_types_table$`Primary Tissue`, function(x) sum(GOI_cBP_mutations$manual_tissue == x,na.rm = TRUE))
all_tissue_types_table$`Percent altered` <- 100*all_tissue_types_table$`Total altered cases`/all_tissue_types_table$`Total cases sequenced`
all_tissue_types_table$`Percent altered` <- round(all_tissue_types_table$`Percent altered`,2)
all_tissue_types_table <- all_tissue_types_table[order(-all_tissue_types_table$`Percent altered`),]

cat("Multiple var table\n")
#table representing frequency of multiple GOI muts per sample
  multi_mut_table <- as.data.frame(table(GOI_cBP_mutations$case_ID_freq))[,2:1]
  colnames(multi_mut_table) = c("Number of samples","Number of mutations per sample")
cat("map variants to nearest exon junction\n")
#map variants to nearest exon junction
  GOI_cBP_mutations$imaging_AA <- GOI_mapping_key[as.character(GOI_cBP_mutations$unified_pos),"nearest_junction_codon"]
  
  #remove rare out of range cases
  GOI_cBP_mutations <- GOI_cBP_mutations[!is.na(GOI_cBP_mutations$imaging_AA),]
  
  Unique_mutations_plot <- GOI_cBP_mutations[!duplicated(GOI_cBP_mutations$amino_acid_change),]
  mut_type_color_key <- data.frame(
    row.names = c(
      "missense_variant",
      "inframe_insertion",
      "inframe_deletion",
      "stop_gained",
      "frameshift_insertion",
      "frameshift_deletion",
      "splice_site_variant",
      "splice_region_variant",
      "start_lost",
      "stop_lost"
    ),
    color = c(
      "#999999", # grey - missense_variant
      "#8DD3C7", # teal - inframe_insertion
      "#F781BF", # pink - inframe_deletion
      "#E41A1C", #red - stop_gained
      "#FF7F00", #orange - frameshift_insertion
      "#FFFF33", # yellow - frameshift_deletion
      "#377EB8", #light blue - splice_site_variant
      "#984EA3", #purple - splice_region_variant
      "#4DAF4A", #green - start_lost
      "#A65628" # brown - stop_lost
    ),
    stringsAsFactors = FALSE
  )
  
  legend <- list(
    labels = row.names(mut_type_color_key),
    fill = mut_type_color_key$color
  )

cat("Figure3\n")    
#Figure3 - all unique variants
  F3_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION"),]
  if(length(F3_feature_df[,1])>0){
    F3_feature_labels <- truncate.feature.labels(F3_feature_df,20)
    F3_features <- GRanges(seqnames = "chr", IRanges(start = F3_feature_df$AA_start, end = F3_feature_df$AA_end, names = NULL))
    F3_features$height <- 0.01
    F3_features$fill <- combined_qualitative_palette[1:length(F3_features)]
    F3_features$featureLayerID <- sep_overlap_features(F3_feature_df)
    F3_variants <- GRanges(seqnames = "chr", IRanges(start = Unique_mutations_plot$imaging_AA, width = 1,names = NULL))
    F3_variants$score <- Unique_mutations_plot$AA_change_freq
    F3_variants$color <- mut_type_color_key[Unique_mutations_plot$unified_annotation,"color"]
    F3_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
    F3_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
    
  }  
  #seperately plot a legend
  # par(mar = c(0,0,0,0))
  # plot.new()
  # legend("center", legend = F3_feature_labels, fill = F3_features$fill)  
cat("Figure4\n")


# GOI_gnomAD_df_filtered$imaging_AA <- GOI_mapping_key[as.character(GOI_gnomAD_df_filtered$unified_pos),"nearest_junction_codon"]
# GOI_gnomAD_df_filtered$imaging_score <- log2(GOI_gnomAD_df_filtered$Allele.Count)
# 
# Ftest_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION"),]
# 
# if(length(Ftest_feature_df[,1])>0){
#   Ftest_feature_labels <- truncate.feature.labels(Ftest_feature_df,20)
#   Ftest_features <- GRanges(seqnames = "chr", IRanges(start = Ftest_feature_df$AA_start, end = Ftest_feature_df$AA_end, names = NULL))
#   Ftest_features$height <- 0.04
#   Ftest_features$fill <- combined_qualitative_palette[1:length(Ftest_features)]
#   Ftest_features$featureLayerID <- sep_overlap_features(Ftest_feature_df)
#   Ftest_features_2 <- Ftest_features
#   Ftest_features_2$height <- 0.04
#   
#   Ftest_gnomAD <- GRanges(seqnames = "chr", IRanges(start = GOI_gnomAD_df_filtered$imaging_AA, width = 1,names = NULL))
#   Ftest_gnomAD$score <- GOI_gnomAD_df_filtered$imaging_score
#   Ftest_gnomAD$color <- mut_type_color_key[GOI_gnomAD_df_filtered$unified_annotation,"color"]
#   #Ftest_gnomAD$SNPsideID <- "bottom"
#   
#   Ftest_variants <- GRanges(seqnames = "chr", IRanges(start = Unique_mutations_plot$imaging_AA, width = 1,names = NULL))
#   Ftest_variants$score <- Unique_mutations_plot$AA_change_freq
#   Ftest_variants$color <- mut_type_color_key[Unique_mutations_plot$unified_annotation,"color"]
#   
#   Ftest_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
#   Ftest_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
#   
# }  

# lolliplot(SNP.gr = list(C=Ftest_gnomAD,B=GRanges(),A=Ftest_variants),
#           features = list(z=Ftest_features,y=Ftest_features_2,x=Ftest_features),
#           ranges = Ftest_ranges,
#           ylab = c("gnomAD variants [log2(allele count)]","UniProt \nFeatures","cBP cancer variants [# of instances]"),
#           xaxis = Ftest_x_axis,
#           cex = 0.75,
#           jitter = "label",
#           legend = legend)


#top 30 variants
  # F4_plot_ceiling <- 50
  # if(length(Unique_mutations_plot[,1]) < 50) F4_plot_ceiling <- length(Unique_mutations_plot[,1]) #if less than 50 mutations just use max; will be redundant plot :/
  # Top_unique_muts_plot <- Unique_mutations_plot[order(-Unique_mutations_plot$AA_change_freq),]
  # Top_unique_muts_plot <- Top_unique_muts_plot[1:F4_plot_ceiling,]
  # 
  # F4_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION"),]
  # if(length(F4_feature_df[,1])>0){
  #   F4_feature_labels <- paste0(F4_feature_df$TYPE,"-",F4_feature_df$LABEL,"[",F4_feature_df$AA_start,"-",F4_feature_df$AA_end,"]")
  #   F4_features <- GRanges(seqnames = "chr", IRanges(start = F4_feature_df$AA_start, end = F4_feature_df$AA_end, names = NULL))
  #   F4_features$height <- 0.01
  #   F4_features$fill <- combined_qualitative_palette[1:length(F4_features)]
  #   F4_features$featureLayerID <- sep_overlap_features(F4_feature_df)
  #   F4_variants <- GRanges(seqnames = "chr", IRanges(start = Top_unique_muts_plot$imaging_AA, width = 1,names = NULL))
  #   F4_variants$score <- Top_unique_muts_plot$AA_change_freq
  #   F4_variants$color <- mut_type_color_key[Top_unique_muts_plot$mutation_type,"color"]
  #   F4_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
  #   F4_yaxis <- round(seq(from = 1, to = Top_unique_muts_plot$AA_change_freq[1], length.out = 10),-1)
  #   F4_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
  #   
  #   #lolliplot(SNP.gr = F4_variants ,features = F4_features, ranges = F4_ranges, ylab = FALSE, xaxis = F4_x_axis,yaxis = F4_yaxis,cex = 1, jitter = "label")
  # }  
  # 
  # 
  # F5_variant_df <- ExAC_df[!(ExAC_df$major_consequence %in% c("intron_variant","non_coding_transcript_exon_variant")),]
  # F5_variant_df$image_AA <- GOI_mapping_key[as.character(F5_variant_df$pos),"nearest_junction_codon"]
  # exome_mut_type_color_key <- data.frame(row.names = unique(F5_variant_df$major_consequence), color = combined_qualitative_palette[1:length(unique(F5_variant_df$major_consequence))])
  # F5_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION"),]
  # F5_feature_labels <- paste0(F5_feature_df$TYPE,"-",F5_feature_df$LABEL,"[",F5_feature_df$AA_start,"-",F5_feature_df$AA_end,"]")
  # F5_features <- GRanges(seqnames = "chr", IRanges(start = F5_feature_df$AA_start, end = F5_feature_df$AA_end, names = NULL))
  # F5_features$height <- 0.01
  # F5_features$fill <- combined_qualitative_palette[1:length(F5_features)]
  # F5_features$featureLayerID <- sep_overlap_features(F5_feature_df)
  # F5_variants <- GRanges(seqnames = "chr", IRanges(start = as.numeric(F5_variant_df$image_AA), width = 1,names = NULL))
  # F5_variants$score <- as.numeric(F5_variant_df$allele_count)
  # F5_variants$color <- exome_mut_type_color_key[F5_variant_df$major_consequence,"color"]
  # F5_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
  # #F5_yaxis <- round(seq(from = 1, to = Top_unique_muts_plot$AA_change_freq[1], length.out = 10),-1)
  # F5_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1)
  # lolliplot(SNP.gr = F5_variants ,features = F5_features, ranges = F5_ranges, ylab = FALSE,cex = .5, jitter = "label")
  # 
  # 
  # F6_variant_df <- ExAC_df[!(ExAC_df$major_consequence %in% c("intron_variant","non_coding_transcript_exon_variant")) & as.numeric(ExAC_df$allele_count) >= 1,]
  # F6_variant_df$image_AA <- GOI_mapping_key[as.character(F6_variant_df$pos),"nearest_junction_codon"]
  # F6_variant_df <- F6_variant_df[!(is.na(F6_variant_df$image_AA)),]
  # exome_mut_type_color_key <- data.frame(row.names = unique(F6_variant_df$major_consequence), color = combined_qualitative_palette[1:length(unique(F6_variant_df$major_consequence))])
  # F6_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION"),]
  # F6_feature_labels <- paste0(F6_feature_df$TYPE,"-",F6_feature_df$LABEL,"[",F6_feature_df$AA_start,"-",F6_feature_df$AA_end,"]")
  # F6_features <- GRanges(seqnames = "chr", IRanges(start = F6_feature_df$AA_start, end = F6_feature_df$AA_end, names = NULL))
  # F6_features$height <- 0.03
  # F6_features$fill <- combined_qualitative_palette[1:length(F6_features)]
  # F6_features$featureLayerID <- sep_overlap_features(F6_feature_df)
  # F6_variants <- GRanges(seqnames = "chr", IRanges(start = as.numeric(F6_variant_df$image_AA), width = 1,names = NULL))
  # F6_variants$score <- log2(as.numeric(F6_variant_df$allele_count))
  # F6_variants$color <- exome_mut_type_color_key[F6_variant_df$major_consequence,"color"]
  # F6_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
  # #F6_yaxis <- round(seq(from = 1, to = Top_unique_muts_plot$AA_change_freq[1], length.out = 10),-1)
  # F6_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1)
  # # lolliplot(SNP.gr = F6_variants ,features = F6_features, ranges = F6_ranges,cex = .5, jitter = "label",ylab = "Log2(allele count)")
  # # 
  # # par(mar = c(0,0,0,0))
  # # plot.new()
  # # legend("center", legend = row.names(exome_mut_type_color_key), fill = exome_mut_type_color_key$color)
  # 
  
  cat("############### Figure data ready for export #################\n\n\n")

