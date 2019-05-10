################# Relative_position_mapping.R ########################
#Obtain relative mapping for variatns
    #Relative to main transcript exons/introns
    #Relative to the union of all exons/transcripts
    #Relative mapping to features
    #Chromosome position to AA position (Main transcript)


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
  cat("mapping nearest junction for non-coding bases","\n")
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
  
  
cat("\n\n\n############################## Relative mapping complete ####################################\n\n\n")