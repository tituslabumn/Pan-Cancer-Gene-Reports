################# Generate_figures.R ########################
#generate figures for output report


cat("#########################################################################\n")
cat("############### Generating figure data for output report ################\n")
cat("#########################################################################\n\n\n")

library("trackViewer")
library(markdown)
library(knitr)
library(RColorBrewer)
library(magick)
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
F1_feature_df$fill <- combined_qualitative_palette[1:length(F1_feature_df[,1])]
#F1_site_df <- GOI_protein_feature_annotation[GOI_protein_feature_annotation$TYPE %in% c("METAL","SITE","MOD_RES","CROSSLNK"),]
if(length(F1_feature_df[,1])>0){
  F1_feature_labels <- truncate.feature.labels(F1_feature_df,20)
  F1_features <- GRanges(seqnames = "chr", IRanges(start = F1_feature_df$AA_start, end = F1_feature_df$AA_end, names = F1_feature_labels))
  F1_features$height <- 0.06
  F1_features$fill <- F1_feature_df$fill
  F1_features$featureLayerID <- sep_overlap_features(F1_feature_df)
  
  #F1_variants <- GRanges(seqnames = "chr", IRanges(start = F1_site_df$AA_start, width = 1), label = as.character(1:length(F1_site_df[,1])))
  #F1_variants$label.col <- "black"
  F1_variants <- GRanges()
  
  F1_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH+1)) # +1 for stop
  F1_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10
  
  #lolliplot(SNP.gr = F1_variants ,features = F1_features, ranges = F1_ranges, ylab = FALSE, xaxis = F1_x_axis)
}
#seperately plot a legend
# par(mar = c(0,0,0,0))
# plot.new()
# legend("center", legend = F1_feature_labels, fill = F1_features$fill)

cat("Figure2\n")

#Figure 2 - all transcripts plotted to relative transcript position
F2_feature_df <- GOI_transcript_exons_filtered
  F2_features <- GRanges(seqnames = "chr", IRanges(start = F2_feature_df$relative_union_start, end = F2_feature_df$relative_union_end, names = NULL))
  F2_features$height <- 0.02
  F2_features$fill <- "grey"
  F2_features$fill[F2_feature_df$ensembl_transcript_id == GOI_TRANSCRIPT] <- "red"
  F2_features$featureLayerID <- F2_feature_df$sort
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
      "#4DAF8A", #green - start_lost
      "#A65628" # brown - stop_lost
    ),
    stringsAsFactors = FALSE
  )
  
  legend <- list(
    labels = row.names(mut_type_color_key),
    fill = mut_type_color_key$color,
    cex = 1.5
  )


  
cat("Figure3\n")    
# Figure3 - all unique variants

    F3_feature_df <- F1_feature_df[F1_feature_df$TYPE %in% c("DOMAIN","REGION"),]
    F3_feature_labels <- truncate.feature.labels(F3_feature_df,20)
    F3_features <- GRanges(seqnames = "chr", IRanges(start = F3_feature_df$AA_start, end = F3_feature_df$AA_end, names = F3_feature_labels))
    F3_features$height <- 0.01
    F3_features$fill <- F3_feature_df$fill
    F3_features$featureLayerID <- sep_overlap_features(F3_feature_df)
    
    # make feature object with unnamed features to plot without label legend
    No_name_features <- GRanges(seqnames = "chr", IRanges(start = F3_feature_df$AA_start, end = F3_feature_df$AA_end, names = NULL))
    No_name_features$height <- 0.03
    No_name_features$fill <- F3_feature_df$fill
    No_name_features$featureLayerID <- sep_overlap_features(F3_feature_df)
    # make feature object that will only display legend
    Feature_legend <- F3_features
    Feature_legend$height <- 0
    Feature_legend$featureLayerID <- 1
    
    #seperate into data frames by annotation type
    # F3_missense
    # F3_frameshift
    # F3_InsDel
    
    
    F3_variants <- GRanges(seqnames = "chr", IRanges(start = Unique_mutations_plot$imaging_AA, width = 1,names = NULL))
    F3_variants$score <- Unique_mutations_plot$AA_change_freq
    F3_variants$color <- mut_type_color_key[Unique_mutations_plot$unified_annotation,"color"]
    F3_variants$type <- Unique_mutations_plot$unified_annotation # for filtering later
    F3_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH+1))
    F3_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 10 rounded to nearest 10
    

save.image("troubleshooting_workspace.RData") #####################  
    
cat("Figure4\n")
F4_types <- c("missense_variant")
# break down further into different classes of missense mutation based on nature of AA change (e.g. charge change)
F4_df <- Unique_mutations_plot[Unique_mutations_plot$unified_annotation == "missense_variant",]
# all cBP $amino_acid_change annotations appear to follow [AA]###[AA] format
# however delins follows [AA]###_[AA]###delins[AAAA] format
# add ref_AA and var_AA cols and a flag for delins
F4_df$ref_AA <- ""
F4_df$var_AA <- ""
F4_df$delins <- FALSE
for (x in 1:length(F4_df[,1])) {
  AA_change <- F4_df[x,"amino_acid_change"]
  if(grepl("delins",AA_change)){
    F4_df[x,"delins"] <- TRUE
    ref <- unlist(strsplit(AA_change, split = "delins"))[1]
    ref <- paste(
      substr(
        unlist(strsplit(ref, split = "_")
        ),1,1),
      collapse = ""
    )
    var <- unlist(strsplit(AA_change, split = "delins"))[2]
  } else {
    ref <- substr(AA_change,1,1)
    var <- substr(AA_change,nchar(AA_change),nchar(AA_change))
  }
  F4_df[x,"ref_AA"] <- ref
  F4_df[x,"var_AA"] <- var
}
rm(x,var,ref)
# assign AA_change_class (from https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/amino-acid-physical-properties.html)
# amino_acid_properties data.frame loaded during initialization from supplied .csv
# add color and shape cols to F4_df
F4_df$color <- ""
F4_df$shape <- "circle"
for(x in 1:length(F4_df[,1])){
  #if has multiple changed AAs
  if(F4_df[x,"delins"]){
    F4_df[x,"color"] <- "#FFFF33" #yellow
    ref <- F4_df[x,"ref_AA"]
    var <- F4_df[x,"var_AA"]
    ref_list <- unlist(strsplit(ref,split = ""))
    var_list <- unlist(strsplit(var,split = ""))
    charge_loss <- FALSE
    charge_gain <- FALSE
    charge_rev <- FALSE
    for(y in 1:nchar(ref)){
      ref <- ref_list[y]
      var <- var_list[y]
      ref_category <- amino_acid_properties[amino_acid_properties$Code == ref,"Category"]
      var_category <- amino_acid_properties[amino_acid_properties$Code == var,"Category"]
      # loss of charge
      if( (ref_category == "Negative/Acidic" | ref_category == "Positive/Basic") & (var_category != "Negative/Acidic" & var_category != "Positive/Basic") ){
        charge_loss <- TRUE
      }
      # gain of charge
      if( (ref_category != "Negative/Acidic" & ref_category != "Positive/Basic") & (var_category == "Negative/Acidic" | var_category == "Positive/Basic") ){
        charge_gain <- TRUE
      }
      # charge reversal
      if( (ref_category == "Negative/Acidic" & var_category == "Positive/Basic") | (var_category == "Negative/Acidic" & ref_category == "Positive/Basic") ){
        charge_rev <- TRUE
      }
    }
    if(sum(c(charge_loss,charge_gain,charge_rev))>1){
      F4_df[x,"shape"] <- "square"
    } else if(sum(c(charge_loss,charge_gain,charge_rev))==1){
      F4_df[x,"shape"] <- c("triangle_point_down", "triangle_point_up","diamond")[which(c(charge_loss,charge_gain,charge_rev))]
    }
  }else{
    # if single AA substitution
    ref <- F4_df[x,"ref_AA"]
    var <- F4_df[x,"var_AA"]
    ref_category <- amino_acid_properties[amino_acid_properties$Code == ref,"Category"]
    var_category <- amino_acid_properties[amino_acid_properties$Code == var,"Category"]
    
    F4_df[x,"color"] <- amino_acid_properties[amino_acid_properties$Code == ref,"Color"]
    # loss of charge
    if( (ref_category == "Negative/Acidic" | ref_category == "Positive/Basic") & (var_category != "Negative/Acidic" & var_category != "Positive/Basic") ){
      F4_df[x,"shape"] <- "triangle_point_down"
    }
    # gain of charge
    if( (ref_category != "Negative/Acidic" & ref_category != "Positive/Basic") & (var_category == "Negative/Acidic" | var_category == "Positive/Basic") ){
      F4_df[x,"shape"] <- "triangle_point_up"
    }
    # charge reversal
    if( (ref_category == "Negative/Acidic" & var_category == "Positive/Basic") | (var_category == "Negative/Acidic" & ref_category == "Positive/Basic") ){
      F4_df[x,"shape"] <- "triangle_point_down"
    }
  }
}
rm(x,y,ref,var,ref_category,var_category,ref_list,var_list)
# assign legend for missense mutation plots
amino_acid_legend <- list(
  labels = c("Negative/Acidic",
             "Positive/Basic",
             "Amidic",
             "Hydroxylic",
             "Sulfur-containing",
             "Aliphatic",
             "Aromatic",
             "Multi-amino-acid/-category"
  ),
  fill = c("#E41A1C",
           "#FF7F00",
           "#984EA3",
           "#4DAF8A",
           "#8DD3C7",
           "#999999",
           "#377EB8",
           "#FFFF33"
  ),
  cex = 1.5
)
F4_variants <- GRanges(seqnames = "chr", IRanges(start = F4_df$imaging_AA, width = 1,names = NULL))
F4_variants$score <- F4_df$AA_change_freq
F4_variants$color <- F4_df$color
F4_variants$shape <- F4_df$shape

cat("Figure5\n")
F5_types <- c("frameshift_insertion","frameshift_deletion","stop_gained","stop_lost","start_lost")
F5_variants <- F3_variants[F3_variants$type %in% F5_types]


cat("Figure6\n")
F6_types <- c("inframe_insertion","inframe_deletion")
F6_variants <- F3_variants[F3_variants$type %in% F6_types]


cat("Figure7\n")
F7_types <- c("splice_site_variant","splice_region_variant")
F7_variants <- F3_variants[F3_variants$type %in% F7_types]
# F7_features will be concatenated with F3 features to add exons to top
F7_exons <- GOI_transcript_exons_filtered[GOI_transcript_exons_filtered$ensembl_transcript_id == GOI_TRANSCRIPT,c("rank","exon_chrom_start","exon_chrom_end")]

# make strand agnostic labels
F7_exons$start <- sapply(1:length(F7_exons[,1]),
                         function(x) {
                           min(GOI_mapping_key[as.character(F7_exons[x,"exon_chrom_start"]),"nearest_junction_codon"],
                           GOI_mapping_key[as.character(F7_exons[x,"exon_chrom_end"]),"nearest_junction_codon"]
                           )
                           }
                         )
F7_exons$end <- sapply(1:length(F7_exons[,1]),
                         function(x) {
                           max(GOI_mapping_key[as.character(F7_exons[x,"exon_chrom_start"]),"nearest_junction_codon"],
                               GOI_mapping_key[as.character(F7_exons[x,"exon_chrom_end"]),"nearest_junction_codon"]
                           )
                         }
)
F7_features <- GRanges(seqnames = "chr", IRanges(start = F7_exons$start,
                                                 end = F7_exons$end,
                                                 names = NULL,
                                                 )
                       )
F7_features$height <- 0.01
F7_features$fill <- rep(c("green","red"),100)[1:length(F7_exons[,1])]#colorRampPalette(c("white", "black"))( length(F7_exons[,1]) )
F7_features$featureLayerID <- 7
colnames(F7_exons)<- c("Exon","Genomic Start","Genomic End","Nearest Peptide Start","Nearest Peptide End")
cat("Figure8\n")
# F8 will be identical to F2 except with variants displayed
F8_exons <- F7_exons # note altered colnames
F8_df <- Unique_mutations_plot[Unique_mutations_plot$unified_annotation %in% F7_types,]
# find nearest relative position from union of exonic sequence
exonic_relative_positions <- as.numeric(rownames(union_transcripts_relative_pos_key))
# for each unified_pos what is the minimal abs(list_of_exonic_pos's - unified_pos)
F8_df$nearest_exonic_pos <- sapply(1:length(F8_df[,1]),
                                   function(x){
                                     union_transcripts_relative_pos_key[which.min(abs(exonic_relative_positions - F8_df$unified_pos[x])),"relative_union_transcript_position"]
                                   }
                                   )
F8_variants <- GRanges(seqnames = "chr", IRanges(start = F8_df$nearest_exonic_pos, width = 1,names = NULL))
F8_variants$score <- F8_df$AA_change_freq
F8_variants$color <- mut_type_color_key[F8_df$unified_annotation,"color"]
#use same features except for main transcript adopt same alternating green/red color scheme
F8_features <- F2_features
F8_features$fill[F8_features$fill == "red"] <- rep(c("green","red"),100)[1:length(F7_exons[,1])]

cat("Figure9\n")
GOI_gnomAD_df_filtered$imaging_AA <- GOI_mapping_key[as.character(GOI_gnomAD_df_filtered$unified_pos),"nearest_junction_codon"]
GOI_gnomAD_df_filtered$imaging_score <- log2(GOI_gnomAD_df_filtered$Allele.Count)
# filter
gnomAD_imaging_df <- GOI_gnomAD_df_filtered[!is.na(GOI_gnomAD_df_filtered$exon_intron),]
gnomAD_leftover_df <-gnomAD_imaging_df[!(gnomAD_imaging_df$unified_annotation %in% legend$labels),]
gnomAD_imaging_df <- gnomAD_imaging_df[(gnomAD_imaging_df$unified_annotation %in% legend$labels),]

max_gnomAD_allele_count <- max(gnomAD_imaging_df[,"Allele.Count"], na.rm = TRUE)
F9_features <- F3_features
F9_gnomAD <- GRanges(seqnames = "chr", IRanges(start = gnomAD_imaging_df[,"imaging_AA"], width = 1,names = NULL))
F9_gnomAD$score <- gnomAD_imaging_df$imaging_score
F9_gnomAD$color <- mut_type_color_key[gnomAD_imaging_df$unified_annotation,"color"]
F9_gnomAD$type <- gnomAD_imaging_df$unified_annotation
F9_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH+1))
F9_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10


cat("Figure10\n") # missense gnomAD
F10_df <- gnomAD_imaging_df[gnomAD_imaging_df$unified_annotation == "missense_variant",]
# at the moment assume no variants are recorded in gnomAD under missense_variant that correspond to delins's of equal size (e.g. ref:AG var:CC)
#  I have been unable to find any examples so far
F10_df$ref_AA <- sapply(F10_df$Protein.Consequence, function(x){
  substr(gsub("p\\.","",x),1,3)
})
F10_df$var_AA <- sapply(F10_df$Protein.Consequence, function(x){
  substr(gsub("p\\.","",x),nchar(gsub("p\\.","",x))-2,nchar(gsub("p\\.","",x)))
})
# assign AA_change_class (from https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-biology-learning-center/protein-biology-resource-library/pierce-protein-methods/amino-acid-physical-properties.html)
# amino_acid_properties data.frame loaded during initialization from supplied .csv
# add color and shape cols to F10_df
F10_df$color <- ""
F10_df$shape <- "circle"
for(x in 1:length(F10_df[,1])){
    # if single AA substitution
    ref <- F10_df[x,"ref_AA"]
    var <- F10_df[x,"var_AA"]
    ref_category <- amino_acid_properties[amino_acid_properties$Three_letter == ref,"Category"]
    var_category <- amino_acid_properties[amino_acid_properties$Three_letter == var,"Category"]
    
    F10_df[x,"color"] <- amino_acid_properties[amino_acid_properties$Three_letter == ref,"Color"]
    # loss of charge
    if( (ref_category == "Negative/Acidic" | ref_category == "Positive/Basic") & (var_category != "Negative/Acidic" & var_category != "Positive/Basic") ){
      F10_df[x,"shape"] <- "triangle_point_down"
    }
    # gain of charge
    if( (ref_category != "Negative/Acidic" & ref_category != "Positive/Basic") & (var_category == "Negative/Acidic" | var_category == "Positive/Basic") ){
      F10_df[x,"shape"] <- "triangle_point_up"
    }
    # charge reversal
    if( (ref_category == "Negative/Acidic" & var_category == "Positive/Basic") | (var_category == "Negative/Acidic" & ref_category == "Positive/Basic") ){
      F10_df[x,"shape"] <- "triangle_point_down"
    }
}
rm(x,ref,var,ref_category,var_category)
F10_gnomAD <- GRanges(seqnames = "chr", IRanges(start = F10_df[,"imaging_AA"], width = 1,names = NULL))
F10_gnomAD$score <- F10_df$imaging_score
F10_gnomAD$color <- F10_df$color
F10_gnomAD$shape <- F10_df$shape
F10_max_gnomAD_allele_count <- 2^max(F10_gnomAD$score, na.rm = TRUE)


cat("Figure11\n")
F11_types <- c("frameshift_insertion","frameshift_deletion","stop_gained","stop_lost","start_lost")
F11_gnomAD <- F9_gnomAD[F9_gnomAD$type %in% F11_types]
F11_max_gnomAD_allele_count <- 2^max(F11_gnomAD$score, na.rm = TRUE)

cat("Figure12\n")
F12_types <- c("inframe_insertion","inframe_deletion")
F12_gnomAD <- F9_gnomAD[F9_gnomAD$type %in% F12_types]
F12_max_gnomAD_allele_count <- 2^max(F12_gnomAD$score, na.rm = TRUE)

cat("Figure13\n")
F13_types <- c("splice_site_variant","splice_region_variant")
F13_gnomAD <- F9_gnomAD[F9_gnomAD$type %in% F13_types]
F13_max_gnomAD_allele_count <- 2^max(F13_gnomAD$score, na.rm = TRUE)
# use F7_exons and F7_features

cat("Figure14\n")
F14_df <- gnomAD_imaging_df[gnomAD_imaging_df$unified_annotation %in% F13_types,]
# exonic_relative_positions generated above
# for each unified_pos what is the minimal abs(list_of_exonic_pos's - unified_pos)
F14_df$nearest_exonic_pos <- sapply(1:length(F14_df[,1]),
                                   function(x){
                                     union_transcripts_relative_pos_key[which.min(abs(exonic_relative_positions - F14_df$unified_pos[x])),"relative_union_transcript_position"]
                                   }
)
F14_variants <- GRanges(seqnames = "chr", IRanges(start = F14_df$nearest_exonic_pos, width = 1,names = NULL))
F14_variants$score <- F14_df$imaging_score
F14_variants$color <- mut_type_color_key[F14_df$unified_annotation,"color"]
#use same features except for main transcript adopt same alternating green/red color scheme
F14_features <- F8_features


cat("Figure15\n")
flat_features <- No_name_features
flat_features$height <- 0
flat_features$featureLayerID <- 1

F15_features <- No_name_features
F15_features$height <- 0.01

F15_shape_key <- data.frame(row.names = c(TRUE,FALSE),
                            shape = c("square","circle"),
                            stringsAsFactors = FALSE)

unique_cBP_overlap_df <- Unique_mutations_plot[Unique_mutations_plot$cBPgnomAD_overlap,]
F15_cBP_overlap <- GRanges(seqnames = "chr", IRanges(start = unique_cBP_overlap_df$imaging_AA, width = 1,names = NULL))
F15_cBP_overlap$score <- unique_cBP_overlap_df$AA_change_freq
F15_cBP_overlap$color <- mut_type_color_key[unique_cBP_overlap_df$unified_annotation,"color"]
F15_cBP_overlap$shape <- F15_shape_key[as.character(unique_cBP_overlap_df$overlap_filter_threshold),]

gnomAD_imageing_overlap_df <- gnomAD_imaging_df[gnomAD_imaging_df$cBPgnomAD_overlap,]
F15_gnomAD_overlap <- GRanges(seqnames = "chr", IRanges(start = gnomAD_imageing_overlap_df[,"imaging_AA"], width = 1,names = NULL))
F15_gnomAD_overlap$score <- gnomAD_imageing_overlap_df$imaging_score
F15_gnomAD_overlap$color <- mut_type_color_key[gnomAD_imageing_overlap_df$unified_annotation,"color"]
F15_gnomAD_overlap$SNPsideID <- "bottom"
F15_gnomAD_overlap$shape <- F15_shape_key[as.character(gnomAD_imageing_overlap_df$overlap_filter_threshold),]
F15_max_gnomAD_allele_count <- 2^max(F15_gnomAD_overlap$score, na.rm = TRUE)

F15_ranges <- GRanges(seqnames = "chr", IRanges(start = 1,end = GOI_UNIPROT_AA_LENGTH+1))
F15_x_axis <- round(seq(from = 1, to = GOI_UNIPROT_AA_LENGTH, length.out = 10),-1) #even split by 5 rounded to nearest 10

save.image("troubleshooting_workspace.RData") #####################  

# overlapping variants table
# add select gnomAD data cols
overlap_summary_table <- gnomAD_imageing_overlap_df[,c("unified_label",
                                                       "Position",
                                                       "rsID",
                                                       "Reference",
                                                       "Alternate",
                                                       "Consequence",
                                                       "unified_annotation",
                                                       "Allele.Count",
                                                       "Homozygote.Count",
                                                       "Allele.Frequency",
                                                       "overlap_filter_threshold"
                                                       )
                                                    ]
# sort by allele freq
overlap_summary_table <- overlap_summary_table[rev(order(overlap_summary_table$Allele.Frequency)),]
# add cBP occurances
overlap_summary_table$cBioPortal.Occurances <- sapply(overlap_summary_table$unified_label, function(x) sum(x == GOI_cBP_mutations$unified_label))
overlap_summary_table <- overlap_summary_table[,-1]



  cat("############### Figure data ready for export #################\n\n\n")

