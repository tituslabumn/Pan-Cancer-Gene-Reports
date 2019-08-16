################# Relative_position_mapping.R ########################
# Obtain relative mapping for variatns
    # Relative to main transcript exons/introns
    # Relative to the union of all exons/transcripts
    # Relative mapping to features
    # Chromosome position to AA position (Main transcript)
    # Use map to fix "MUTATED" instances of $amino_acid_change rarely returned by some cBP queries


cat("#####################################################################\n")
cat("#### Relative mapping to main transcript AA, exons and features #####\n")
cat("#####################################################################\n\n\n")

#parse in data frames to non-specific names
GOI_exon_annotation_main_transcript <- get(paste0(GOI,"_exon_annotation_main_transcript"))
GOI_exon_annotation <- get(paste0(GOI,"_exon_annotation"))
GOI_peptide_sequence <- get(paste0(GOI,"_peptide_sequences"))
GOI_protein_feature_annotation <- get(paste0(GOI,"_protein_feature_annotation"))

GOI_STRAND <- GOI_exon_annotation$strand[1]

save.image("troubleshooting_workspace.RData") #####################

#check that main transcript selected is the same cds returned by UnipProtKB
cat("Check that UniProt coding sequence is the same length as main transcript selected:","\n")
peptide_transcript_check <-GOI_UNIPROT_AA_LENGTH == (GOI_exon_annotation_main_transcript$cds_length[1]/3-1)
cat("\t",peptide_transcript_check,"\n\n\n")
#if fail stop
if(is.na(peptide_transcript_check)){
  cat("No canonical cBP transcript; default to UniProt/n/n")
  #select transcript (first that matches) with same cds size as uniprot product
  GOI_TRANSCRIPT <- GOI_exon_annotation$ensembl_transcript_id[which((GOI_exon_annotation$cds_length/3-1) == GOI_UNIPROT_AA_LENGTH)][1]
  GOI_exon_annotation_main_transcript <- GOI_exon_annotation[GOI_exon_annotation$ensembl_transcript_id == GOI_TRANSCRIPT,]
}else if(peptide_transcript_check != TRUE) {
  stop("peptide cds check against chosen GOI transcript FAILED/n/n")
}


#make new data frame that maps all introns based on transcript stat/end and exon boundries
GOI_intron_map <- data.frame(
  rank = 1:(max(GOI_exon_annotation_main_transcript$rank)-1), #intron 1 is b/t exons 1 and 2, intron 2 is b/t exon 2 &3 etc
  stringsAsFactors = FALSE
)
GOI_intron_map$intron_chrom_start <- integer(length(GOI_intron_map$rank))
GOI_intron_map$intron_chrom_end <- integer(length(GOI_intron_map$rank))
for (x in GOI_intron_map$rank) {
  adjacent_exons <- c(GOI_exon_annotation_main_transcript$exon_chrom_start[GOI_exon_annotation_main_transcript$rank == x],
                      GOI_exon_annotation_main_transcript$exon_chrom_end[GOI_exon_annotation_main_transcript$rank == x],
                      GOI_exon_annotation_main_transcript$exon_chrom_start[GOI_exon_annotation_main_transcript$rank == x+1],
                      GOI_exon_annotation_main_transcript$exon_chrom_end[GOI_exon_annotation_main_transcript$rank == x+1]
                      )
  adjacent_exons<-adjacent_exons[order(adjacent_exons)] #lowest @ [1]
  GOI_intron_map[x,"intron_chrom_start"] <- adjacent_exons[2]+1
  GOI_intron_map[x,"intron_chrom_end"] <- adjacent_exons[3]-1
}
rm(x,adjacent_exons)

#retrive and align transcrit sequence with BSgenome
#THIS IS SEQUENCE FOR POSITIVE STRAND (regardless of GOI coding strand)
library(BSgenome.Hsapiens.NCBI.GRCh38)
GOI_mapping_key <- data.frame(
  row.names = (GOI_exon_annotation_main_transcript$transcript_start[1]:GOI_exon_annotation_main_transcript$transcript_end[1]),
  relative_transcript = 1:(1+GOI_exon_annotation_main_transcript$transcript_end[1]-GOI_exon_annotation_main_transcript$transcript_start[1]),
  transcript = unlist(strsplit(as.character(getSeq(Hsapiens, 
                                                   as.character(GOI_CHR), 
                                                   start = GOI_exon_annotation_main_transcript$transcript_start[1], 
                                                   end = GOI_exon_annotation_main_transcript$transcript_end[1]
                                                   )
                                            ),"")
                      ),
  stringsAsFactors = FALSE
)
if(GOI_STRAND == -1) GOI_mapping_key$relative_transcript <- rev(GOI_mapping_key$relative_transcript)
GOI_mapping_key$exon_intron <- character(length(GOI_mapping_key[,1]))

#add col of position (redundant to row.names() but makes easily searchable in rstudio GUI)
GOI_mapping_key$pos <- row.names(GOI_mapping_key)

#map introns
  for (x in GOI_intron_map$rank) {
    GOI_mapping_key$exon_intron[ (rownames(GOI_mapping_key) >= GOI_intron_map$intron_chrom_start[x]) & (rownames(GOI_mapping_key) <= GOI_intron_map$intron_chrom_end[x]) ] <- paste0("Intron_",x)
  }
  rm(x)
  
#map exons
  #re-order exons
  GOI_exon_annotation_main_transcript <- GOI_exon_annotation_main_transcript[order(GOI_exon_annotation_main_transcript$rank),]
  for (x in GOI_exon_annotation_main_transcript$rank) {
    GOI_mapping_key$exon_intron[ (rownames(GOI_mapping_key) >= GOI_exon_annotation_main_transcript$exon_chrom_start[x]) & (rownames(GOI_mapping_key) <= GOI_exon_annotation_main_transcript$exon_chrom_end[x]) ] <- paste0("Exon_",x)
  }
  rm(x)
  
#map "Coding"/"Intronic"/"X'_UTR"
  GOI_mapping_key$coding_segments_key <- character(length(GOI_mapping_key[,1]))
  GOI_mapping_key$coding_segments_key[grepl("Exon",GOI_mapping_key$exon_intron)] <- "Coding"
  #get boundries of coding sequence (strand agnostic)
  cds_genomic_min <- min(c(GOI_exon_annotation_main_transcript$genomic_coding_start,GOI_exon_annotation_main_transcript$genomic_coding_end),na.rm = TRUE)
  cds_genomic_max <- max(c(GOI_exon_annotation_main_transcript$genomic_coding_start,GOI_exon_annotation_main_transcript$genomic_coding_end),na.rm = TRUE)
  transcript_genomic_start <- GOI_exon_annotation_main_transcript$transcript_start[1]
  transcript_genomic_end <- GOI_exon_annotation_main_transcript$transcript_end[1]
  if(GOI_STRAND == 1){ #change labeling of 3' vs 5' UTRs based on GOI_STRAND
    GOI_mapping_key$coding_segments_key[ (row.names(GOI_mapping_key) >= transcript_genomic_start) & (row.names(GOI_mapping_key) < cds_genomic_min) ] <- "5'_UTR"
    GOI_mapping_key$coding_segments_key[ (row.names(GOI_mapping_key) <= transcript_genomic_end) & (row.names(GOI_mapping_key) > cds_genomic_max) ] <- "3'_UTR"
  }else{
    GOI_mapping_key$coding_segments_key[ (row.names(GOI_mapping_key) >= transcript_genomic_start) & (row.names(GOI_mapping_key) < cds_genomic_min) ] <- "3'_UTR"
    GOI_mapping_key$coding_segments_key[ (row.names(GOI_mapping_key) <= transcript_genomic_end) & (row.names(GOI_mapping_key) > cds_genomic_max) ] <- "5'_UTR"
  }
  GOI_mapping_key$coding_segments_key[grepl("Intron",GOI_mapping_key$exon_intron)] <- "Intronic"
  
#map relative_coding_sequence and relative_AA_sequence (codon number)
  #(non-coding sequences are left as 0 for both columns)
  cds_index <- grepl("Coding",GOI_mapping_key$coding_segments_key) #include stop codon for now
  GOI_mapping_key$relative_coding_sequence <- numeric(length(GOI_mapping_key[,1]))
  GOI_mapping_key$relative_AA_position <- numeric(length(GOI_mapping_key[,1]))
  if(GOI_STRAND == 1){
    GOI_mapping_key$relative_coding_sequence[cds_index] <- 1:sum(cds_index)
  }else{
    GOI_mapping_key$relative_coding_sequence[cds_index] <- sum(cds_index):1
  }
  GOI_mapping_key$relative_AA_position <- (GOI_mapping_key$relative_coding_sequence+2)%/%3 #conver cds to codon number
  
#map pptide sequence to relative AA postion
  GOI_mapping_key$peptide <- NA
  GOI_main_peptide_sequence <- GOI_peptide_sequence[GOI_peptide_sequence$ensembl_transcript_id == GOI_TRANSCRIPT,"peptide"]
  if(GOI_STRAND == 1) GOI_mapping_key$peptide[cds_index] <- rep(unlist(strsplit(GOI_main_peptide_sequence,"")), each = 3)
  if(GOI_STRAND == -1) GOI_mapping_key$peptide[cds_index] <- rev(rep(unlist(strsplit(GOI_main_peptide_sequence,"")), each = 3))
  
#map nearest_intron_exon_junction
  #this is for aligning splice region variants to the nearest codon for visualization
  #if equidistant default to upstream on positive strand
  cat("mapping nearest junction for non-coding bases","\n\n")
  #copy over relative_AA_position
  GOI_mapping_key$nearest_junction_codon <- GOI_mapping_key$relative_AA_position
  proximity_df <- data.frame(
    index = 1:length(GOI_mapping_key[,1]),
    stringsAsFactors = FALSE
  )
  proximity_df$junction <- FALSE
  proximity_df$junction_number <- 0
  proximity_df$AA_position <- 0
  junction_num <- 1 #unique number for block of non coding
  for (x in proximity_df$index[2:(length(cds_index)-1)]) {
    if(!cds_index[x] & cds_index[x+1]) { #if next index is transition to coding
      proximity_df$junction[x] <- TRUE
      proximity_df$junction_number[x] <- junction_num
      junction_num <- junction_num + 1
      proximity_df$AA_position[x] <- GOI_mapping_key$relative_AA_position[x+1]
    }else if(!cds_index[x] & cds_index[x-1]){ #if next index is tansition to coding 
      proximity_df$junction[x] <- TRUE
      proximity_df$junction_number[x] <- junction_num
      proximity_df$AA_position[x] <- GOI_mapping_key$relative_AA_position[x-1]
    }
  }
  rm(x,junction_num)
  proximity_df <- proximity_df[proximity_df$junction,]
  #assign nearest AA numbers
    #first UTR
    GOI_mapping_key$nearest_junction_codon[(1:proximity_df$index[1])] <- proximity_df$AA_position[1]
    #last UTR
    GOI_mapping_key$nearest_junction_codon[(proximity_df$index[length(proximity_df[,1])]:length(GOI_mapping_key[,1]))] <- proximity_df$AA_position[length(proximity_df[,1])]
    #introns
    for(x in 2:(max(proximity_df$junction_number)-1)){
      # cat(x,"\n")
      intron_subset <- proximity_df[proximity_df$junction_number == x,]
      cutoff <- intron_subset$index[1] + floor((intron_subset$index[2]-intron_subset$index[1])/2)
      # cat(intron_subset$index[1],"\n")
      # cat(cutoff,"\n")
      # cat(intron_subset$index[2],"\n")
      GOI_mapping_key$nearest_junction_codon[intron_subset$index[1]:cutoff] <- intron_subset$AA_position[1]
      GOI_mapping_key$nearest_junction_codon[(cutoff+1):intron_subset$index[2]] <- intron_subset$AA_position[2]
    }
    rm(x,intron_subset)
  
  
    cat("mapping relative transcript position (union of all transcripts without intronic space)","\n\n")
#make key for relative transcript position (union of all transcripts without intronic space)
  Unique_GOI_ensembl_exons <- GOI_exon_annotation[!(duplicated(GOI_exon_annotation$ensembl_exon_id)),]
  #get union of positions for each range apppended
  union_ranges <- numeric()
  for(x in 1:length(Unique_GOI_ensembl_exons[,1])){
    exon_range <- Unique_GOI_ensembl_exons[x,"exon_chrom_start"]:Unique_GOI_ensembl_exons[x,"exon_chrom_end"]
    union_ranges <- union(union_ranges,exon_range)
  }
  rm(x,exon_range)
  union_ranges <- union_ranges[order(union_ranges)]
  union_transcripts_relative_pos_key <- data.frame(
    row.names = union_ranges,
    relative_union_transcript_position = 1:length(union_ranges),
    stringsAsFactors = FALSE
  )
  if(GOI_STRAND == -1){
    union_transcripts_relative_pos_key$relative_union_transcript_position <- rev(union_transcripts_relative_pos_key$relative_union_transcript_position)
    GOI_exon_annotation$relative_union_start <- union_transcripts_relative_pos_key[as.character(GOI_exon_annotation$exon_chrom_end),"relative_union_transcript_position"]
    GOI_exon_annotation$relative_union_end <- union_transcripts_relative_pos_key[as.character(GOI_exon_annotation$exon_chrom_start),"relative_union_transcript_position"]
  }else{
    GOI_exon_annotation$relative_union_start <- union_transcripts_relative_pos_key[as.character(GOI_exon_annotation$exon_chrom_start),"relative_union_transcript_position"]
    GOI_exon_annotation$relative_union_end <- union_transcripts_relative_pos_key[as.character(GOI_exon_annotation$exon_chrom_end),"relative_union_transcript_position"]
  }

save.image("troubleshooting_workspace.RData") #####################  
  
cat("fixing cBP aa_change labeled as 'MUTATED'\n")
cat("\tinstances:",sum(GOI_cBP_mutations$amino_acid_change == "MUTATED"),"\n\n")
# in some cases cBP returns AA_change as "MUTATED" manually fix annotations in these cases. 
  if(sum(GOI_cBP_mutations$amino_acid_change == "MUTATED") > 0){
    library(Biostrings) #for codon lookup
    # if respective mutation_type is not "1" or NA append position and alt allele to the mutation_type, this will be enough to distinguish
    GOI_cBP_mutations[GOI_cBP_mutations$amino_acid_change == "MUTATED","AA"] <- GOI_mapping_key[as.character(GOI_cBP_mutations[GOI_cBP_mutations$amino_acid_change == "MUTATED","start_position"]),"nearest_junction_codon"]
    known_mutation_types <- c("Frame_Shift_Del","Missense_Mutation","In_Frame_Ins","Frame_Shift_Ins","Nonsense_Mutation","Splice_Site","Splice_Region","In_Frame_Del","Translation_Start_Site","Nonstop_Mutation","Frame_Shift","Targeted_Region")
    cat("\tfixing instances of known mut_types\n")
    known_type_index <- which((GOI_cBP_mutations$amino_acid_change == "MUTATED") & (GOI_cBP_mutations$mutation_type %in% known_mutation_types))
    GOI_cBP_mutations[known_type_index,"amino_acid_change"] <- paste(GOI_cBP_mutations$mutation_type[known_type_index],GOI_cBP_mutations$start_position[known_type_index],GOI_cBP_mutations$variant_allele[known_type_index],sep = "_")
    cat("\t\tinstances:",sum(GOI_cBP_mutations$amino_acid_change == "MUTATED"),"\n\n")
    
    cat("\tfixing other instances\n")
    # if rescpective mutation_type is not a recognized type manually assign one, remove synonymous muts after
    for(x in which(GOI_cBP_mutations$amino_acid_change == "MUTATED")){
      #remove all non-coding variants for now
      if(GOI_mapping_key[as.character(GOI_cBP_mutations$start_position[x]),"coding_segments_key"] != "Coding"){
        GOI_cBP_mutations[x,"amino_acid_change"] <- "remove"
      }else if(nchar(GOI_cBP_mutations$reference_allele[x]) == 1 & nchar(GOI_cBP_mutations$variant_allele[x]) == 1 & GOI_cBP_mutations$reference_allele[x] != "-" & GOI_cBP_mutations$variant_allele[x] != "-"){ #assume all others are single base modifications
        ref_AA <- GOI_mapping_key[as.character(GOI_cBP_mutations$start_position[x]),"peptide"]
        #get list of transcript bases in order of increasing genomic pos (correct order for + strand genes)
        ref_codon <- GOI_mapping_key$transcript[GOI_mapping_key$relative_AA_position == GOI_mapping_key[as.character(GOI_cBP_mutations$start_position[x]),"relative_AA_position"]]
        ref_codon_pos <- which(GOI_mapping_key$pos[GOI_mapping_key$relative_AA_position == GOI_mapping_key[as.character(GOI_cBP_mutations$start_position[x]),"relative_AA_position"]] == as.character(GOI_cBP_mutations$start_position[x]))
        if(GOI_STRAND == -1){ #rev comp
          var_codon <- ref_codon
          var_codon[ref_codon_pos] <- as.character(complement(DNAString(GOI_cBP_mutations$variant_allele[x])))
          var_codon <- as.character(complement(DNAString(paste(rev(var_codon), collapse = ""))))
          var_AA <- GENETIC_CODE[[var_codon]]
          ref_codon <- as.character(complement(DNAString(paste(rev(ref_codon), collapse = ""))))
        }else{
          var_codon <- ref_codon
          var_codon[ref_codon_pos] <- GOI_cBP_mutations$variant_allele[x]
          var_codon <- paste(var_codon, collapse = "")
          var_AA <- GENETIC_CODE[[var_codon]]
          ref_codon <- paste(ref_codon, collapse = "")
        }
        # remove if synonymous
        if(var_AA == ref_AA){
          GOI_cBP_mutations[x,"amino_acid_change"] <- "remove"
        }else{
          GOI_cBP_mutations[x,"amino_acid_change"] <- paste0(ref_AA,GOI_mapping_key[as.character(GOI_cBP_mutations$start_position[x]),"relative_AA_position"],var_AA)
        }
        #assign mut_type
        if(var_AA == "*"){
          GOI_cBP_mutations[x,"mutation_type"] <- "Nonsense_Mutation"
        }else{
          GOI_cBP_mutations[x,"mutation_type"] <- "Missense_Mutation" #label as such even if not correct; those cases are already flagged for removal
        }
      }
    }
    #remove rows flagged for removal
    GOI_cBP_fusions <- GOI_cBP_mutations[GOI_cBP_mutations$amino_acid_change != "remove",]
    cat("\t\tinstances:",sum(GOI_cBP_mutations$amino_acid_change == "MUTATED"),"\n\n")
    #clean out temp variables
    rm(x,known_type_index,ref_AA,var_AA,ref_codon_pos,ref_codon,var_codon)
  }
  
  #update AA_change_freq
  GOI_cBP_mutations$AA_change_freq <- sapply(GOI_cBP_mutations$amino_acid_change, function(x) sum(GOI_cBP_mutations$amino_acid_change == x))
  
  #update AA_freq, this time only count for coding changes (start or end position in cds)
  coding_variant_index <- which((GOI_mapping_key[as.character(GOI_cBP_mutations$start_position),"coding_segments_key"] == "Coding") | (GOI_mapping_key[as.character(GOI_cBP_mutations$end_position),"coding_segments_key"] == "Coding"))
  GOI_cBP_mutations[coding_variant_index,"AA_freq"] <- sapply(GOI_cBP_mutations[coding_variant_index,"AA"], function(x) sum(GOI_cBP_mutations[coding_variant_index,"AA"] == x, na.rm = TRUE))
  
  
  save.image("troubleshooting_workspace.RData") #####################  
  
  # positions for single base changes are not different between the two (always correspond to the actual base pos)
  # insertions and deletions are not different between + and - stranded genes 
  #   ins: gnomAD_pos = cBP_end; cBP_end = base before inserted sequence; cBP_start = pos of base after insertion (always end + 1)
  #   del: gnomAD_pos = cBP_start - 1; cBP start and end match start and end of deletion sequence; gnomAD pos matches base before deleted section
  # make new col in both datasets: $unified_pos
  #   cBP ins = end_pos
  #   gnomAD ins = dont change
  #   cBP del = start_pos
  #   gnomAD del = pos + 1
  #   assume it is possible for ambiguous edge cases where ref and var are both greater than 1 in length and possibly also equal in length 
  cat("Adding $unified_pos and $unified_annotation to both cBP and gnomAD data frames for easy comparison of matching variants\n")
  
  cat("\tgnomAD:\n")
  GOI_gnomAD_df_filtered$unified_pos <- 0
  GOI_gnomAD_df_filtered$unified_annotation <- GOI_gnomAD_df_filtered$Annotation
  #loop thru each var and add unified_pos and unified_annotation
  for(x in 1:length(GOI_gnomAD_df_filtered$unified_pos)){
    ref <- GOI_gnomAD_df_filtered[x,"Reference"]
    var <- GOI_gnomAD_df_filtered[x,"Alternate"]
    pos <- GOI_gnomAD_df_filtered[x,"Position"]
    flag <- ""
    # flag used to assign fs type below
    if(ref == "-" | nchar(ref) < nchar(var)){ #insertion
      GOI_gnomAD_df_filtered[x,"unified_pos"] <- pos #no change
      flag <- "ins"
    } else if(var == "-" | nchar(ref) > nchar(var)){ #deletion
      GOI_gnomAD_df_filtered[x,"unified_pos"] <- pos + 1
      flag <- "del"
    } else if(nchar(ref) == 1 & nchar(var) == 1 & ref != var){ #missense
      GOI_gnomAD_df_filtered[x,"unified_pos"] <- pos #no change
    } else if(nchar(ref) == nchar(var)){ #ins/del???
      cat(paste0("\t\tIns/Del??? - pos: ",pos," ref: ",ref," var: ",var,"\n"))
      GOI_gnomAD_df_filtered[x,"unified_pos"] <- pos #no change
    } else {
      cat(paste0("\t\tEdge case? - pos: ",pos," ref: ",ref," var: ",var,"\n"))
      GOI_gnomAD_df_filtered[x,"unified_pos"] <- pos #no change
      flag <- "other"
    }
    
    if(GOI_gnomAD_df_filtered[x,"Annotation"] == "frameshift_variant" & flag == "ins"){
      GOI_gnomAD_df_filtered[x,"unified_annotation"] <- "frameshift_insertion"
    } else if(GOI_gnomAD_df_filtered[x,"Annotation"] == "frameshift_variant" & flag == "del"){
      GOI_gnomAD_df_filtered[x,"unified_annotation"] <- "frameshift_deletion"
    } else if(GOI_gnomAD_df_filtered[x,"Annotation"] == "splice_acceptor_variant" | GOI_gnomAD_df_filtered[x,"Annotation"] == "splice_donor_variant"){
      GOI_gnomAD_df_filtered[x,"unified_annotation"] <- "splice_site_variant"
    } 
    
  }
  rm(x,var,ref,flag,pos)
  
  save.image("troubleshooting_workspace.RData") #####################  
  
  cat("\tcBP:\n")
  GOI_cBP_mutations$unified_pos <- 0
  GOI_cBP_mutations$unified_annotation <- ""
  #loop thru each var and add unified_pos and unified_annotation
  for(x in 1:length(GOI_cBP_mutations$unified_pos)){
    ref <- GOI_cBP_mutations[x,"reference_allele"]
    var <- GOI_cBP_mutations[x,"variant_allele"]
    start_pos <- GOI_cBP_mutations[x,"start_position"]
    end_pos <- GOI_cBP_mutations[x,"end_position"]
    flag <- ""
    cat(x," ")
    # flag used to assign fs type below
    if(ref == "-" | nchar(ref) < nchar(var)){ #insertion
      GOI_cBP_mutations[x,"unified_pos"] <- end_pos #no change
      flag <- "ins"
    } else if(var == "-" | nchar(ref) > nchar(var)){ #deletion
      GOI_cBP_mutations[x,"unified_pos"] <- start_pos
      flag <- "del"
    } else if(nchar(ref) == 1 & nchar(var) == 1 & ref != var){ #missense
      GOI_cBP_mutations[x,"unified_pos"] <- end_pos
    } else if(nchar(ref) == nchar(var)){ #ins/del???
      cat(paste0("\t\tIns/Del??? - pos: ",pos," ref: ",ref," var: ",var,"\n"))
      GOI_cBP_mutations[x,"unified_pos"] <- start_pos #choose start pos arbitrarily
    } else {
      cat(paste0("\t\tEdge case? - pos: ",start_pos," ref: ",ref," var: ",var,"\n"))
      GOI_cBP_mutations[x,"unified_pos"] <- start_pos #choose start pos arbitrarily
      flag <- "other"
    }
    
    if(GOI_cBP_mutations[x,"mutation_type"] == "Frame_Shift_Ins"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "frameshift_insertion"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Frame_Shift_Del"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "frameshift_deletion"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Missense_Mutation"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "missense_variant"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "In_Frame_Ins"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "inframe_insertion"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Nonsense_Mutation"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "stop_gained"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Splice_Site"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "splice_site_variant"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Splice_Region"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "splice_region_variant"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "In_Frame_Del"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "inframe_deletion"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Translation_Start_Site"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "start_lost"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Nonstop_Mutation"){
      GOI_cBP_mutations[x,"unified_annotation"] <- "stop_lost"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Targeted_Region"){ # found in TP53 I don't know what this is, has ref and var == "-"
      cat("\t\tFound 'Targeted_Region' type; flagging for removal\n")
      GOI_cBP_mutations[x,"unified_annotation"] <- "remove"
    } else if(GOI_cBP_mutations[x,"mutation_type"] == "Frame_Shift"){ # found in TP53, called from incorrect annotation of long ref and long var, is actually single base ins
      cat("\t\tFound 'Frame_Shift' type (not ins or del?) ; flagging for removal\n")
      GOI_cBP_mutations[x,"unified_annotation"] <- "remove"
    } else {
      ann <- GOI_cBP_mutations[x,"mutation_type"]
      cat(paste0("\t\tEdge case? - pos: ",start_pos," ref: ",ref," var: ",var,"\n", "annotation: ", ann))
      GOI_cBP_mutations[x,"unified_annotation"] <- "remove"
    }
    
  }
  rm(x,var,ref,flag,start_pos,end_pos,ann)
  #remove rows flagged for removal
  cat("\tRemoving ",sum(GOI_cBP_mutations$unified_annotation == "remove")," rows flagged for removal\n\n")
  GOI_cBP_mutations <- GOI_cBP_mutations[GOI_cBP_mutations$unified_annotation != "remove",]
  
  save.image("troubleshooting_workspace.RData") #####################  
  
  #map positions of gnomAD variants ; note: some gnomAD variants are outside map key (e.g. where they map to other transcript with larger UTR)
  cat("relative mapping gnomAD variants based on unified position\n\n")
  for (x in c("relative_transcript","exon_intron","coding_segments_key","relative_coding_sequence","relative_AA_position","peptide","nearest_junction_codon")) {
    GOI_gnomAD_df_filtered[,x] <- GOI_mapping_key[as.character(GOI_gnomAD_df_filtered$unified_pos),x]
  }
  #map positions of gnomAD variants ; note: some gnomAD variants are outside map key (e.g. where they map to other transcript with larger UTR)
  cat("relative mapping cBP variants based on unified position\n\n")
  for (x in c("relative_transcript","exon_intron","coding_segments_key","relative_coding_sequence","relative_AA_position","peptide","nearest_junction_codon")) {
    GOI_cBP_mutations[,x] <- GOI_mapping_key[as.character(GOI_cBP_mutations$unified_pos),x]
  }
  
  save.image("troubleshooting_workspace.RData") #####################  
  
cat("\n\n\n############################## Relative mapping complete ####################################\n\n\n")
