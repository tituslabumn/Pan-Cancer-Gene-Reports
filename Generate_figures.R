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

#parse some data frames to non-gene-specific names
GOI_cBP_mutations <- get(paste0(GOI,"_cBP_mutations"))
GOI_cBP_fusions <- get(paste0(GOI,"_cBP_fusions"))

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


#Figure 1 - features; domains,regions, DNA_bind, MOTIF
F1_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION","DNA_BIND","MOTIF"),]
#F1_site_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("METAL","SITE","MOD_RES","CROSSLNK"),]
if(length(F1_feature_df[,1])>0){
  F1_feature_labels <- truncate.feature.labels(F1_feature_df,20)
  F1_features <- GRanges(seqnames = "chr", IRanges(start = F1_feature_df$AA_start, end = F1_feature_df$AA_end, names = NULL))
  F1_features$height <- 0.04
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


#get master case df with only studies that survived cBP mut query filtering
all_mut_cases <- master_case_df[master_case_df$study %in% unique(GOI_cBP_mutations$study),]
GOI_cBP_mutations$manual_tissue <- all_mut_cases$manual_tissue_annotation[match(paste0(GOI_cBP_mutations$study,GOI_cBP_mutations$altered_case_id),paste0(all_mut_cases$study,all_mut_cases$altered_case_id))]

#tissue type enrichment
colnames(all_tissue_types_table) <- c("Primary Tissue","Total cases in database")
all_tissue_types_table$`Total cases sequenced` <- sapply(all_tissue_types_table$`Primary Tissue`, function(x) sum(all_mut_cases$manual_tissue_annotation == x,na.rm = TRUE))
GOI_cBP_mutations$manual_tissue <- as.character(GOI_cBP_mutations$manual_tissue)
all_tissue_types_table$`Total altered cases` <- sapply(all_tissue_types_table$`Primary Tissue`, function(x) sum(GOI_cBP_mutations$manual_tissue == x,na.rm = TRUE))
all_tissue_types_table$`Percent altered` <- 100*all_tissue_types_table$`Total altered cases`/all_tissue_types_table$`Total cases sequenced`
all_tissue_types_table$`Percent altered` <- round(all_tissue_types_table$`Percent altered`,2)
all_tissue_types_table <- all_tissue_types_table[order(-all_tissue_types_table$`Percent altered`),]

#table representing frequency of multiple GOI muts per sample
  multi_mut_table <- as.data.frame(table(GOI_cBP_mutations$case_ID_freq))[,2:1]
  colnames(multi_mut_table) = c("Number of samples","Number of mutations per sample")

#map variants to nearest exon junction
  if(GOI_STRAND == -1){ #this only matters for large deletions
    GOI_cBP_mutations$imaging_AA <- GOI_mapping_key[as.character(GOI_cBP_mutations$end_position),"nearest_junction_codon"]
  }else{
    GOI_cBP_mutations$imaging_AA <- GOI_mapping_key[as.character(GOI_cBP_mutations$start_position),"nearest_junction_codon"]
  }
  #remove rare out of range cases
  GOI_cBP_mutations <- GOI_cBP_mutations[!is.na(GOI_cBP_mutations$imaging_AA),]
  
  Unique_mutations_plot <- GOI_cBP_mutations[!duplicated(GOI_cBP_mutations$amino_acid_change),]
  mut_type_color_key <- data.frame(row.names = unique(Unique_mutations_plot$mutation_type), color = combined_qualitative_palette[1:length(unique(Unique_mutations_plot$mutation_type))])
  
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
    F3_variants$color <- mut_type_color_key[Unique_mutations_plot$mutation_type,"color"]
    F3_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
    F3_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
    
    #lolliplot(SNP.gr = F3_variants ,features = F3_features, ranges = F3_ranges, ylab = FALSE, xaxis = F3_x_axis,cex = 1, jitter = "label")
  }  
  #seperately plot a legend
  # par(mar = c(0,0,0,0))
  # plot.new()
  # legend("center", legend = F3_feature_labels, fill = F3_features$fill)  

#top 30 variants
  Top_unique_muts_plot <- Unique_mutations_plot[order(-Unique_mutations_plot$AA_change_freq),]
  Top_unique_muts_plot <- Top_unique_muts_plot[1:50,]
  
  F4_feature_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("DOMAIN","REGION"),]
  if(length(F4_feature_df[,1])>0){
    F4_feature_labels <- paste0(F4_feature_df$TYPE,"-",F4_feature_df$LABEL,"[",F4_feature_df$AA_start,"-",F4_feature_df$AA_end,"]")
    F4_features <- GRanges(seqnames = "chr", IRanges(start = F4_feature_df$AA_start, end = F4_feature_df$AA_end, names = NULL))
    F4_features$height <- 0.01
    F4_features$fill <- combined_qualitative_palette[1:length(F4_features)]
    F4_features$featureLayerID <- sep_overlap_features(F4_feature_df)
    F4_variants <- GRanges(seqnames = "chr", IRanges(start = Top_unique_muts_plot$imaging_AA, width = 1,names = NULL))
    F4_variants$score <- Top_unique_muts_plot$AA_change_freq
    F4_variants$color <- mut_type_color_key[Top_unique_muts_plot$mutation_type,"color"]
    F4_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH))
    F4_yaxis <- round(seq(from = 1, to = Top_unique_muts_plot$AA_change_freq[1], length.out = 10),-1)
    F4_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
    
    #lolliplot(SNP.gr = F4_variants ,features = F4_features, ranges = F4_ranges, ylab = FALSE, xaxis = F4_x_axis,yaxis = F4_yaxis,cex = 1, jitter = "label")
  }  
  
  
cat("############### Figure data ready for export #################\n\n\n")

