# Martina Mijuskovic
# SV pipeline development
# Tiering SVs in the cancer cohort
# Assessment of false positive rate in Domain 1 SVs
# Dec 2017

library(dplyr)
library(ggplot2)
library(data.table)
library(jsonlite)
library(VariantAnnotation)
library(ensembldb)

options(max.print=10000)

####### Prepare a list of SV.VCF paths ####### 
### use FF PCR free with known tumour type in the first instance)

# Read list of suitable FF samples
FF_list <- read.table("./Data/FF_PCRfree_wClinical_forSVtiering.tsv", sep = "\t", header = T)
dim(FF_list)
table(FF_list$TumourType, exclude = NULL)

# Read list on HPC
FF_list <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/FF_PCRfree_wClinical_forSVtiering.tsv", sep = "\t", header = T)

# Add SV.vcf path (execute on HPC)
FF_list$SV_vcf_path <- as.character(sapply(FF_list$Path, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))

# Write out the file with paths
write.table(FF_list, file = "FF_PCRfree_wClinical_forSVtiering_wPaths.tsv", quote = F, row.names = F, col.names = T, sep = "\t")



#######  Prepare the GO list of Domain 1 transcripts ####### 

# Read Alona's file
#actionable_SVs <- read.table("/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Data/GO/GENOMONCOLOGY_SOLID_TUMOUR.SV.v1.6.tsv", sep = "\t", header = T)  # old corrupt file
actionable_SVs <- read.table("/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Data/GO/GENOMONCOLOGY_SOLID_TUMOUR.SV.fixed.v1.6.tsv", sep = "\t", header = T)

head(actionable_SVs)  # looks ok
length(unique(actionable_SVs$gene_name))  # 106 genes (!?)
length(unique(actionable_SVs$transcript_ID))  # 105 unique transcripts
length(unique(actionable_SVs$gene_ID))  # 105 unique gene IDs
dim(actionable_SVs) # 1383 
sum(duplicated(actionable_SVs))  # 0 duplicate rows
#actionable_SVs[duplicated(actionable_SVs),]
# Deduplicate
#actionable_SVs <- actionable_SVs[!duplicated(actionable_SVs),]
#dim(actionable_SVs)

# Remove non-genes
actionable_SVs <- actionable_SVs %>% filter(transcript_ID != "")
dim(actionable_SVs)  # 1380
length(unique(actionable_SVs$gene_name))  # 104
length(unique(actionable_SVs$transcript_ID))  # 104
length(unique(actionable_SVs$gene_ID))  # 104

# Write out the file without non-genes
write.table(actionable_SVs, file = "./Data/GO/GENOMONCOLOGY_SOLID_TUMOUR.SV.fixed.noNonGenes.v1.6.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# Get list of unique transcripts for annotation to Domain 1
GO_transcr <- unique(actionable_SVs$transcript_ID)
length(GO_transcr)  # 104






####### Transcript annotation functions #######

### Annotation function returning VCF info field from Domain 1 variants
### Also returns SAMPLE_WELL_ID, overlapping transript from the list for start and end breaksite, and all related mates/events even if not overlapping a transcript from the list
### Keeps Canvas and Manta calls

# Use on HPC (fixed)
getSVsInTranscripts <- function(vcf_path, transcripts){
  
  require(ensembldb)
  require(VariantAnnotation)
  require(dplyr)
  
  ### Read VCF
  
  vcf <- readVcf(vcf_path)
  # Get VCF info fields
  vcf_info <- as.data.frame(info(vcf))
  vcf_info$FILTER <- rowRanges(vcf)$FILTER
  vcf_info$Application <- ""
  vcf_info[grepl("Canvas", rownames(vcf_info)),]$Application <- "Canvas"
  vcf_info[grepl("Manta", rownames(vcf_info)),]$Application <- "Manta"
  vcf_info$ID <- rownames(vcf_info)
  # Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
  vcf_info$START <- as.data.frame(ranges(vcf))$start
  vcf_info$CHR <- as.character(seqnames(vcf))
  # Add number of supporting paired and reads -----fixed
  vcf_info$PR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][1]
  })
  vcf_info$PR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][2]
  })
  # Add number of supporting split reads  --------fixed
  vcf_info$SR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][1]
  })
  vcf_info$SR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][2]
  })
  # Change chr names in the VCF info table
  vcf_info$CHR <- sub("chr", "", vcf_info$CHR)
  vcf_info$MATEID <- as.character(vcf_info$MATEID)
  
  
  ### Read from Ensembl database
  
  # Read the transcript database
  txdb_pth <- '/home/mmijuskovic/Fusions/Homo_sapiens.GRCh38.84.sqlite'
  if (!file.exists(txdb_pth)) {
    #txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
    print("No Homo_sapiens.GRCh38.84.sqlite found. See https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/ for installation")
  }
  txdb <- EnsDb(txdb_pth)
  
  # Get transcript ranges from Ensembl transcript database
  transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts))
  
  
  ### Add transcript overlap annotations to SV breakpoints
  
  # Flag start
  vcf_info$transcript_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr))
  
  # Flag end
  vcf_info$transcript_related_end <- sapply(1:dim(vcf_info)[1], function(x){
    if(is.na(vcf_info$END[x])){
      return(NA)
    } 
    else {
      if(overlapsAny(GRanges(seqnames=vcf_info$CHR[x], ranges=IRanges(start = vcf_info$END[x], width = 1)), transcripts_gr)) {
        return(1)
      }
      else {
        return(0)
      }
    }
  })
  
  
  # Add overlapping gene to breakpoints
  
  # Initiate fields
  vcf_info$start_ann <- ""
  vcf_info$end_ann <- ""
  
  # Find exact overlaps of start and end position with transcripts
  #overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), transcripts_gr)  # testing
  overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr)
  #overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR[x], ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END[x], width = 1)), transcripts_gr)
  overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR, ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END, width = 1)), transcripts_gr)
  
  # Add start position overlaps to table
  if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
    vcf_info[queryHits(overlaps_start),]$start_ann <- sapply(seq(length(overlaps_start)), function(x){
      i <- queryHits(overlaps_start)[x]
      j <- subjectHits(overlaps_start)[x]
      vcf_info[i,]$start_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  # Add end position overlaps to table --- FIXED
  if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
    vcf_info[!is.na(vcf_info$transcript_related_end),][queryHits(overlaps_end),]$end_ann <- sapply(seq(length(overlaps_end)), function(x){
      i <- queryHits(overlaps_end)[x]
      j <- subjectHits(overlaps_end)[x]
      vcf_info[!is.na(vcf_info$transcript_related_end),][i,]$end_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  
  # Reduce VCF to annotated SVs, keeping all SVs within the same event
  mates <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(MATEID)) %>% pull(MATEID)
  events <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(EVENT)) %>% pull(EVENT)
  extracted_SVs <- vcf_info %>% filter((start_ann != "" | end_ann != "") | ID %in% mates | EVENT %in% events)
  # Deduplicate
  extracted_SVs <- extracted_SVs[!duplicated(extracted_SVs),]
  
  
  # Add sample ID and return table if any results
  if (dim(extracted_SVs)[1] != 0){
    ### Add sample ID
    # Use on HPC
    extracted_SVs$SAMPLE_WELL_ID <- strsplit(strsplit(vcf_path, split = "Cancer")[[1]][2], split = "_Normal")[[1]][1]  # works for HPC path, not local
    # Use locally
    #extracted_SVs$SAMPLE_WELL_ID <- paste0("LP", strsplit(strsplit(samples$local_VCF_path[1], ".somatic")[[1]][1], "_LP")[[1]][2])
    
    return(extracted_SVs)
  }
}

### Annotation function returning VCF info field from Domain 1 variants
### Also returns SAMPLE_WELL_ID, overlapping transript from the list for start and end breaksite, and all related mates/events even if not overlapping a transcript from the list
### Removes non-SOMATIC (Canvas REF)
### Keeps only Manta calls (SOMATIC = TRUE flag is missing in Canvas calls)
### Adds a list of VCF filters to include (default is only "PASS") - NOTE that some might have missing mates/events (if filtered out)

getSVsInTranscripts_v2 <- function(vcf_path, transcripts, filters = "PASS"){
  
  require(ensembldb)
  require(VariantAnnotation)
  require(dplyr)
  
  ### Read VCF
  
  vcf <- readVcf(vcf_path)
  # Get VCF info fields
  vcf_info <- as.data.frame(info(vcf))
  vcf_info$FILTER <- rowRanges(vcf)$FILTER
  vcf_info$Application <- ""
  vcf_info[grepl("Canvas", rownames(vcf_info)),]$Application <- "Canvas"
  vcf_info[grepl("Manta", rownames(vcf_info)),]$Application <- "Manta"
  vcf_info$ID <- rownames(vcf_info)
  # Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
  vcf_info$START <- as.data.frame(ranges(vcf))$start
  vcf_info$CHR <- as.character(seqnames(vcf))
  # Add number of supporting paired and reads -----fixed
  vcf_info$PR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][1]
  })
  vcf_info$PR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][2]
  })
  # Add number of supporting split reads  --------fixed
  vcf_info$SR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][1]
  })
  vcf_info$SR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][2]
  })
  # Change chr names in the VCF info table
  vcf_info$CHR <- sub("chr", "", vcf_info$CHR)
  vcf_info$MATEID <- as.character(vcf_info$MATEID)
  
  
  
  
  ### Subset VCF for SVs with desired FILTER and remove non-SOMATIC
  
  filters_keep <- filters
  
  vcf_info <- vcf_info %>% filter(FILTER %in% filters_keep, SOMATIC == "TRUE")  # NOTE that this removes all Canvas calls bc they don't have "SOMATIC" in the INFO
  
  
  
  ### Read from Ensembl database
  
  # Read the transcript database
  txdb_pth <- '/home/mmijuskovic/Fusions/Homo_sapiens.GRCh38.84.sqlite'
  if (!file.exists(txdb_pth)) {
    #txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
    print("No Homo_sapiens.GRCh38.84.sqlite found. See https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/ for installation")
  }
  txdb <- EnsDb(txdb_pth)
  
  # Get transcript ranges from Ensembl transcript database
  transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts))
  

  
  
  
  ### Add transcript overlap annotations to SV breakpoints
  
  # Flag start
  vcf_info$transcript_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr))
  
  # Flag end
  vcf_info$transcript_related_end <- sapply(1:dim(vcf_info)[1], function(x){
    if(is.na(vcf_info$END[x])){
      return(NA)
    } 
    else {
      if(overlapsAny(GRanges(seqnames=vcf_info$CHR[x], ranges=IRanges(start = vcf_info$END[x], width = 1)), transcripts_gr)) {
        return(1)
      }
      else {
        return(0)
      }
    }
  })
  
  
  # Add overlapping gene to breakpoints
  
  # Initiate fields
  vcf_info$start_ann <- ""
  vcf_info$end_ann <- ""
  
  # Find exact overlaps of start and end position with transcripts
  #overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), transcripts_gr)  # testing
  overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr)
  #overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR[x], ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END[x], width = 1)), transcripts_gr)
  overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR, ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END, width = 1)), transcripts_gr)
  
  # Add start position overlaps to table
  if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
    vcf_info[queryHits(overlaps_start),]$start_ann <- sapply(seq(length(overlaps_start)), function(x){
      i <- queryHits(overlaps_start)[x]
      j <- subjectHits(overlaps_start)[x]
      vcf_info[i,]$start_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  # Add end position overlaps to table --- FIXED
  if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
    vcf_info[!is.na(vcf_info$transcript_related_end),][queryHits(overlaps_end),]$end_ann <- sapply(seq(length(overlaps_end)), function(x){
      i <- queryHits(overlaps_end)[x]
      j <- subjectHits(overlaps_end)[x]
      vcf_info[!is.na(vcf_info$transcript_related_end),][i,]$end_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  
  # Reduce VCF to annotated SVs, keeping all SVs within the same event
  mates <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(MATEID)) %>% pull(MATEID)
  events <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(EVENT)) %>% pull(EVENT)
  extracted_SVs <- vcf_info %>% filter((start_ann != "" | end_ann != "") | ID %in% mates | EVENT %in% events)
  # Deduplicate
  extracted_SVs <- extracted_SVs[!duplicated(extracted_SVs),]
  
  
  # Add sample ID and return table if any results
  if (dim(extracted_SVs)[1] != 0){
    ### Add sample ID
    # Use on HPC
    extracted_SVs$SAMPLE_WELL_ID <- strsplit(strsplit(vcf_path, split = "Cancer")[[1]][2], split = "_Normal")[[1]][1]  # works for HPC path, not local
    # Use locally
    #extracted_SVs$SAMPLE_WELL_ID <- paste0("LP", strsplit(strsplit(samples$local_VCF_path[1], ".somatic")[[1]][1], "_LP")[[1]][2])
    
    return(extracted_SVs)
  }
}


### Improved function that also extracts REF and ALT (essetial for BND directionality)

getSVsInTranscripts_v3 <- function(vcf_path, transcripts, filters = "PASS"){
  
  require(ensembldb)
  require(VariantAnnotation)
  require(dplyr)
  
  ### Read VCF
  
  vcf <- readVcf(vcf_path)
  # Get VCF info fields
  vcf_info <- as.data.frame(info(vcf))
  
  # Get REF and ALT --- added
  vcf_info$REF <- as.data.frame(rowRanges(vcf)$REF)$x
  vcf_info$ALT <- as.data.frame(rowRanges(vcf)$ALT)$value
  
  # Get other fields
  vcf_info$FILTER <- rowRanges(vcf)$FILTER
  vcf_info$Application <- ""
  vcf_info[grepl("Canvas", rownames(vcf_info)),]$Application <- "Canvas"
  vcf_info[grepl("Manta", rownames(vcf_info)),]$Application <- "Manta"
  vcf_info$ID <- rownames(vcf_info)
  # Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
  vcf_info$START <- as.data.frame(ranges(vcf))$start
  vcf_info$CHR <- as.character(seqnames(vcf))
  # Add number of supporting paired and reads -----fixed
  vcf_info$PR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][1]
  })
  vcf_info$PR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][2]
  })
  # Add number of supporting split reads  --------fixed
  vcf_info$SR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][1]
  })
  vcf_info$SR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][2]
  })
  # Change chr names in the VCF info table
  vcf_info$CHR <- sub("chr", "", vcf_info$CHR)
  vcf_info$MATEID <- as.character(vcf_info$MATEID)
  
  
  
  
  ### Subset VCF for SVs with desired FILTER and remove non-SOMATIC
  
  filters_keep <- filters
  
  vcf_info <- vcf_info %>% filter(FILTER %in% filters_keep, SOMATIC == "TRUE")  # NOTE that this removes all Canvas calls bc they don't have "SOMATIC" in the INFO
  
  
  
  ### Read from Ensembl database
  
  # Read the transcript database
  txdb_pth <- '/home/mmijuskovic/Fusions/Homo_sapiens.GRCh38.84.sqlite'
  if (!file.exists(txdb_pth)) {
    #txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
    print("No Homo_sapiens.GRCh38.84.sqlite found. See https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/ for installation")
  }
  txdb <- EnsDb(txdb_pth)
  
  # Get transcript ranges from Ensembl transcript database
  transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts))
  
  
  
  
  
  ### Add transcript overlap annotations to SV breakpoints
  
  # Flag start
  vcf_info$transcript_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr))
  
  # Flag end
  vcf_info$transcript_related_end <- sapply(1:dim(vcf_info)[1], function(x){
    if(is.na(vcf_info$END[x])){
      return(NA)
    } 
    else {
      if(overlapsAny(GRanges(seqnames=vcf_info$CHR[x], ranges=IRanges(start = vcf_info$END[x], width = 1)), transcripts_gr)) {
        return(1)
      }
      else {
        return(0)
      }
    }
  })
  
  
  # Add overlapping gene to breakpoints
  
  # Initiate fields
  vcf_info$start_ann <- ""
  vcf_info$end_ann <- ""
  
  # Find exact overlaps of start and end position with transcripts
  #overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), transcripts_gr)  # testing
  overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr)
  #overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR[x], ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END[x], width = 1)), transcripts_gr)
  overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR, ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END, width = 1)), transcripts_gr)
  
  # Add start position overlaps to table
  if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
    vcf_info[queryHits(overlaps_start),]$start_ann <- sapply(seq(length(overlaps_start)), function(x){
      i <- queryHits(overlaps_start)[x]
      j <- subjectHits(overlaps_start)[x]
      vcf_info[i,]$start_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  # Add end position overlaps to table --- FIXED
  if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
    vcf_info[!is.na(vcf_info$transcript_related_end),][queryHits(overlaps_end),]$end_ann <- sapply(seq(length(overlaps_end)), function(x){
      i <- queryHits(overlaps_end)[x]
      j <- subjectHits(overlaps_end)[x]
      vcf_info[!is.na(vcf_info$transcript_related_end),][i,]$end_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  
  # Reduce VCF to annotated SVs, keeping all SVs within the same event
  mates <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(MATEID)) %>% pull(MATEID)
  events <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(EVENT)) %>% pull(EVENT)
  extracted_SVs <- vcf_info %>% filter((start_ann != "" | end_ann != "") | ID %in% mates | EVENT %in% events)
  # Deduplicate
  extracted_SVs <- extracted_SVs[!duplicated(extracted_SVs),]
  
  
  # Add sample ID and return table if any results
  if (dim(extracted_SVs)[1] != 0){
    ### Add sample ID
    # Use on HPC
    extracted_SVs$SAMPLE_WELL_ID <- strsplit(strsplit(vcf_path, split = "Cancer")[[1]][2], split = "_Normal")[[1]][1]  # works for HPC path, not local
    # Use locally
    #extracted_SVs$SAMPLE_WELL_ID <- paste0("LP", strsplit(strsplit(samples$local_VCF_path[1], ".somatic")[[1]][1], "_LP")[[1]][2])
    
    return(extracted_SVs)
  }
}

getSVsInTranscripts_v3_local <- function(vcf_path, transcripts, filters = "PASS"){
  
  require(ensembldb)
  require(VariantAnnotation)
  require(dplyr)
  
  ### Read VCF
  
  vcf <- readVcf(vcf_path)
  # Get VCF info fields
  vcf_info <- as.data.frame(info(vcf))
  
  # Get REF and ALT --- added
  vcf_info$REF <- as.data.frame(rowRanges(vcf)$REF)$x
  vcf_info$ALT <- as.data.frame(rowRanges(vcf)$ALT)$value
  
  # Get other fields
  vcf_info$FILTER <- rowRanges(vcf)$FILTER
  vcf_info$Application <- ""
  vcf_info[grepl("Canvas", rownames(vcf_info)),]$Application <- "Canvas"
  vcf_info[grepl("Manta", rownames(vcf_info)),]$Application <- "Manta"
  vcf_info$ID <- rownames(vcf_info)
  # Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
  vcf_info$START <- as.data.frame(ranges(vcf))$start
  vcf_info$CHR <- as.character(seqnames(vcf))
  # Add number of supporting paired and reads -----fixed
  vcf_info$PR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][1]
  })
  vcf_info$PR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[1]][,2][x][[1]][2]
  })
  # Add number of supporting split reads  --------fixed
  vcf_info$SR_REF <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][1]
  })
  vcf_info$SR_ALT <- sapply(1:dim(vcf_info)[1], function(x){
    geno(vcf)[[2]][,2][x][[1]][2]
  })
  # Change chr names in the VCF info table
  vcf_info$CHR <- sub("chr", "", vcf_info$CHR)
  vcf_info$MATEID <- as.character(vcf_info$MATEID)
  
  
  
  
  ### Subset VCF for SVs with desired FILTER and remove non-SOMATIC
  
  filters_keep <- filters
  
  vcf_info <- vcf_info %>% filter(FILTER %in% filters_keep, SOMATIC == "TRUE")  # NOTE that this removes all Canvas calls bc they don't have "SOMATIC" in the INFO
  
  
  
  ### Read from Ensembl database
  
  # Read the transcript database
  txdb_pth <- '/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Homo_sapiens.GRCh38.84.sqlite'
  if (!file.exists(txdb_pth)) {
    #txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
    print("No Homo_sapiens.GRCh38.84.sqlite found. See https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/ for installation")
  }
  txdb <- EnsDb(txdb_pth)
  
  # Get transcript ranges from Ensembl transcript database
  transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts))
  
  
  
  
  
  ### Add transcript overlap annotations to SV breakpoints
  
  # Flag start
  vcf_info$transcript_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr))
  
  # Flag end
  vcf_info$transcript_related_end <- sapply(1:dim(vcf_info)[1], function(x){
    if(is.na(vcf_info$END[x])){
      return(NA)
    } 
    else {
      if(overlapsAny(GRanges(seqnames=vcf_info$CHR[x], ranges=IRanges(start = vcf_info$END[x], width = 1)), transcripts_gr)) {
        return(1)
      }
      else {
        return(0)
      }
    }
  })
  
  
  # Add overlapping gene to breakpoints
  
  # Initiate fields
  vcf_info$start_ann <- ""
  vcf_info$end_ann <- ""
  
  # Find exact overlaps of start and end position with transcripts
  #overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), transcripts_gr)  # testing
  overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr)
  #overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR[x], ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END[x], width = 1)), transcripts_gr)
  overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$transcript_related_end),]$CHR, ranges=IRanges(start = vcf_info[!is.na(vcf_info$transcript_related_end),]$END, width = 1)), transcripts_gr)
  
  # Add start position overlaps to table
  if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
    vcf_info[queryHits(overlaps_start),]$start_ann <- sapply(seq(length(overlaps_start)), function(x){
      i <- queryHits(overlaps_start)[x]
      j <- subjectHits(overlaps_start)[x]
      vcf_info[i,]$start_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  # Add end position overlaps to table --- FIXED
  if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
    vcf_info[!is.na(vcf_info$transcript_related_end),][queryHits(overlaps_end),]$end_ann <- sapply(seq(length(overlaps_end)), function(x){
      i <- queryHits(overlaps_end)[x]
      j <- subjectHits(overlaps_end)[x]
      vcf_info[!is.na(vcf_info$transcript_related_end),][i,]$end_ann <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  
  # Reduce VCF to annotated SVs, keeping all SVs within the same event
  mates <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(MATEID)) %>% pull(MATEID)
  events <- vcf_info %>% filter(start_ann != "" | end_ann != "") %>%  filter(Application == "Manta" & !is.na(EVENT)) %>% pull(EVENT)
  extracted_SVs <- vcf_info %>% filter((start_ann != "" | end_ann != "") | ID %in% mates | EVENT %in% events)
  # Deduplicate
  extracted_SVs <- extracted_SVs[!duplicated(extracted_SVs),]
  
  
  # Add sample ID and return table if any results
  if (dim(extracted_SVs)[1] != 0){
    ### Add sample ID
    # Use on HPC
    extracted_SVs$SAMPLE_WELL_ID <- strsplit(strsplit(vcf_path, split = "Cancer")[[1]][2], split = "_Normal")[[1]][1]  # works for HPC path, not local
    # Use locally
    #extracted_SVs$SAMPLE_WELL_ID <- paste0("LP", strsplit(strsplit(samples$local_VCF_path[1], ".somatic")[[1]][1], "_LP")[[1]][2])
    
    return(extracted_SVs)
  }
}

# Read the list of transcripts (HPC)
actionable_SVs <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/GENOMONCOLOGY_SOLID_TUMOUR.SV.fixed.noNonGenes.v1.6.tsv", sep = "\t", header = T)

# Get unique list of transcripts
GO_transcr <- unique(actionable_SVs$transcript_ID)

# Read the list of samples with paths to SV VCF
FF_list <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/FF_PCRfree_wClinical_forSVtiering_wPaths.tsv", sep = "\t", header = T)
FF_list$SV_vcf_path <- as.character(FF_list$SV_vcf_path)

# # Test (works!)
# #getSVsInTranscripts_v2("/genomes/by_date/2016-11-29/CANCP40570/CancerLP3000070-DNA_G08_NormalLP3000069-DNA_E08/SomaticVariations/LP3000069-DNA_E08_LP3000070-DNA_G08.somatic.SV.vcf.gz", GO_transcr)
# ptm <- proc.time()
# #domain1_SVs <- lapply(FF_list$SV_vcf_path[1], getSVsInTranscripts_v2, GO_transcr, filters = c("PASS", "q10"))
# domain1_SVs <- lapply(FF_list$SV_vcf_path[1], getSVsInTranscripts_v2, GO_transcr)
# proc.time() - ptm


##### Test annotation ##### 

# Run (HPC) --- test
setwd("/home/mmijuskovic/SV_dev/SV_tiering/tiered_results/test")

library(dplyr)
library(VariantAnnotation)
library(ensembldb)

# Read the list of samples with paths to SV VCF
FF_list <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/FF_PCRfree_wClinical_forSVtiering_wPaths.tsv", sep = "\t", header = T)
FF_list$SV_vcf_path <- as.character(FF_list$SV_vcf_path)

# Read the list of transcripts (HPC)
actionable_SVs <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/GENOMONCOLOGY_SOLID_TUMOUR.SV.fixed.noNonGenes.v1.6.tsv", sep = "\t", header = T)
# Get unique list of transcripts
GO_transcr <- unique(actionable_SVs$transcript_ID)
GO_transcr <- GO_transcr[GO_transcr != ""]

domain1_SVs <- lapply(FF_list$SV_vcf_path[1:50], getSVsInTranscripts_v2, GO_transcr)  # job 1
#domain1_SVs <- lapply(FF_list$SV_vcf_path[1:5], getSVsInTranscripts_v2, GO_transcr)
#domain1_SVs <- lapply(FF_list$SV_vcf_path[1:5], getSVsInTranscripts_v3, GO_transcr)

domain1_SVs <- bind_rows(domain1_SVs)
domain1_SVs <- domain1_SVs %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
save.image("domain1_SVs_job1.RData")

domain1_SVs <- lapply(FF_list$SV_vcf_path[51:100], getSVsInTranscripts_v2, GO_transcr)  # job 2
domain1_SVs <- bind_rows(domain1_SVs)
domain1_SVs <- domain1_SVs %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
save.image("domain1_SVs_job2.RData")


### Load test files locally

load("./Data/test/domain1_SVs_job1.RData")
domain1_SVs_1 <- domain1_SVs
load("./Data/test/domain1_SVs_job2.RData")
domain1_SVs_2 <- domain1_SVs

# Merge
domain1_SVs <- rbind(domain1_SVs_1, domain1_SVs_2)

# Summarize (NOTE that any DEL >10kb is automatically filtered out and not included here)
dim(domain1_SVs)  # 96 total PASS filter in 100 patients (~1 average per patient, BND counted twice, some may be missing mates)
length(unique(domain1_SVs$SAMPLE_WELL_ID))  # 36/100 patients have PASS SVs in Domain 1
table(domain1_SVs$SVTYPE, domain1_SVs$Application)  # No Canvas - check
table(domain1_SVs$FILTER, domain1_SVs$SOMATIC)
domain1_SVs %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, SVTYPE, SVLEN, CHR, START, END, ColocalizedCanvas, PR_ALT, SR_ALT, start_ann, end_ann)


# Add KEY, calculate recurrent variants, add variant frequency info to the table
domain1_SVs$KEY <- sapply(1:dim(domain1_SVs)[1], function(x){
  paste(domain1_SVs$CHR[x], domain1_SVs$START[x], domain1_SVs$END[x], domain1_SVs$SVTYPE[x], sep = "-")
})

unique_keys <- unique(domain1_SVs$KEY)
length(unique_keys)  # 93
recurrent_SVs <- data.frame(
  KEY = unique_keys,
  NUM_OBS = sapply(unique_keys, function(x){ sum(domain1_SVs$KEY == x)}),
  NUM_PATIENTS = sapply(unique_keys, function(x){ length(unique(domain1_SVs %>% filter(KEY == x) %>% .$SAMPLE_WELL_ID)) })
  )
rownames(recurrent_SVs) <- NULL

# Add variant frequency info to the main table
recurr_keys <- recurrent_SVs %>% filter(NUM_OBS >1) %>% pull(KEY)
domain1_SVs$recurrent <- 0  # 3 total
domain1_SVs[!domain1_SVs$KEY %in% recurr_keys,] <- 1


  
#### Annotate Domain 1 SVs #### 
  
# 4 interactive jobs (190 samples each) on HPC

#setwd("/home/mmijuskovic/SV_dev/SV_tiering/tiered_results")
setwd("/home/mmijuskovic/SV_dev/SV_tiering/tiered_results_fixed")

library(dplyr)
library(VariantAnnotation)
library(ensembldb)

# Read the list of samples with paths to SV VCF
FF_list <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/FF_PCRfree_wClinical_forSVtiering_wPaths.tsv", sep = "\t", header = T)
FF_list$SV_vcf_path <- as.character(FF_list$SV_vcf_path)


# Read the list of transcripts (HPC)
actionable_SVs <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/GENOMONCOLOGY_SOLID_TUMOUR.SV.fixed.noNonGenes.v1.6.tsv", sep = "\t", header = T)
# Get unique list of transcripts
GO_transcr <- unique(actionable_SVs$transcript_ID)


domain1_SVs <- lapply(FF_list$SV_vcf_path[1:190], getSVsInTranscripts_v2, GO_transcr)  # job 1
domain1_SVs <- bind_rows(domain1_SVs)
domain1_SVs <- domain1_SVs %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
save.image("domain1_SVs_job1.RData")

domain1_SVs <- lapply(FF_list$SV_vcf_path[191:380], getSVsInTranscripts_v2, GO_transcr)  # job 2
domain1_SVs <- bind_rows(domain1_SVs)
domain1_SVs <- domain1_SVs %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
save.image("domain1_SVs_job2.RData")

domain1_SVs <- lapply(FF_list$SV_vcf_path[381:570], getSVsInTranscripts_v2, GO_transcr)  # job 3
domain1_SVs <- bind_rows(domain1_SVs)
domain1_SVs <- domain1_SVs %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
save.image("domain1_SVs_job3.RData")

domain1_SVs <- lapply(FF_list$SV_vcf_path[571:759], getSVsInTranscripts_v2, GO_transcr)  # job 4
domain1_SVs <- bind_rows(domain1_SVs)
domain1_SVs <- domain1_SVs %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
save.image("domain1_SVs_job4.RData")   





#### Summary analysis of Domain1 SVs #### 

load("./Data/domain1_SVs_job1.RData")
domain1_SVs_1 <- domain1_SVs
load("./Data/domain1_SVs_job2.RData")
domain1_SVs_2 <- domain1_SVs
load("./Data/domain1_SVs_job3.RData")
domain1_SVs_3 <- domain1_SVs
load("./Data/domain1_SVs_job4.RData")
domain1_SVs_4 <- domain1_SVs

# Merge
domain1_SVs <- rbind(domain1_SVs_1, domain1_SVs_2, domain1_SVs_3, domain1_SVs_4)
rm(domain1_SVs_1, domain1_SVs_2, domain1_SVs_3, domain1_SVs_4)

# Summary ((NOTE that any DEL >10kb is automatically filtered out and not included here)
dim(domain1_SVs)  # 1022 total (TRs listed as 2 BNDs)
table(domain1_SVs$SVTYPE)  # 719 Domain 1 SVs (~1 per patient)
# BND DEL DUP INS INV 
# 606 128  27  18 243

# Proportion of different SV types
table(domain1_SVs$SVTYPE)/719
# BND        DEL        DUP        INS        INV 
# 0.84283727 0.17802503 0.03755216 0.02503477 0.33796940  # transl are 0.4214186

# Median SV per patient
median(as.numeric(table(domain1_SVs$SAMPLE_WELL_ID)))  # 2


# Number of patients with Domain 1 SVs
length(unique(domain1_SVs$SAMPLE_WELL_ID))  # 311/759 (41% patients have Domain 1 SVs)
table(domain1_SVs$SVTYPE, domain1_SVs$Application)  # No Canvas
table(domain1_SVs$FILTER, domain1_SVs$SOMATIC)
table(domain1_SVs$SVTYPE, domain1_SVs$IMPRECISE)
# FALSE TRUE
# BND   464  142
# DEL   127    1
# DUP    27    0
# INS    18    0
# INV   201   42
domain1_SVs %>% dplyr::select(SAMPLE_WELL_ID, ID, IMPRECISE, SVTYPE, SVLEN, CHR, START, END, ColocalizedCanvas, PR_ALT, SR_ALT, start_ann, end_ann)


# Add KEY, calculate recurrent variants, add VAF and variant frequency info to the table
domain1_SVs$KEY <- sapply(1:dim(domain1_SVs)[1], function(x){
  paste(domain1_SVs$CHR[x], domain1_SVs$START[x], domain1_SVs$END[x], domain1_SVs$SVTYPE[x], sep = "-")
})

unique_keys <- unique(domain1_SVs$KEY)
length(unique_keys)  # 1004
recurrent_SVs <- data.frame(
  KEY = unique_keys,
  NUM_OBS = sapply(unique_keys, function(x){ sum(domain1_SVs$KEY == x)})
)
rownames(recurrent_SVs) <- NULL
# Look at recurrent SVs
recurrent_SVs %>% filter(NUM_OBS > 1) %>% arrange(desc(NUM_OBS))

# Add variant frequency info to the main table
recurr_keys <- as.character(recurrent_SVs %>% filter(NUM_OBS >1) %>% pull(KEY))
domain1_SVs$recurrent <- 0  
domain1_SVs[domain1_SVs$KEY %in% recurr_keys,]$recurrent <- 1  # 9 total (3 deletions, 3 transl)
table(domain1_SVs$recurrent, domain1_SVs$SVTYPE)
domain1_SVs %>% filter(recurrent == 1, SVTYPE == "DEL")
domain1_SVs %>% filter(recurrent == 1, KEY == "4-54688371-54688815-DEL")
domain1_SVs %>% filter(recurrent == 1, KEY == "12-68810707-68811111-DEL")
domain1_SVs %>% filter(recurrent == 1, KEY == "5-142596308-142596582-DEL")
domain1_SVs %>% filter(recurrent == 1, SVTYPE == "BND")
domain1_SVs %>% filter(recurrent == 1, KEY == "10-43115587-NA-BND")
domain1_SVs %>% filter(recurrent == 1, ID == "MantaBND:114898:0:1:0:0:0:0") # 2 SVs in same EVENT, sharing one breakopoint, but mate not recurrent
domain1_SVs %>% filter(recurrent == 1, ID == "MantaBND:114898:0:1:1:3:0:1") # 2 SVs in same EVENT, sharing one breakopoint, but mate not recurrent

### Remove BND with missing mates and recalculate summaries
domain1_SVs$missing_mate <- 0
missing_mates <- sapply(unique(domain1_SVs$SAMPLE_WELL_ID), function(x){
                      mates <- domain1_SVs %>% filter(SAMPLE_WELL_ID == x, SVTYPE == "BND") %>% pull(MATEID)
                      # List missing mates
                      mates[!mates %in% (domain1_SVs %>% filter(SAMPLE_WELL_ID == x, SVTYPE == "BND") %>% pull(ID))]
                      })
table(as.character(missing_mates))  # none missing


### Calculate VAF
domain1_SVs$VAF <- sapply(1:dim(domain1_SVs)[1], function(x){
  sum(domain1_SVs$PR_ALT[x], domain1_SVs$SR_ALT[x], na.rm = T)/sum(domain1_SVs$PR_REF[x], domain1_SVs$SR_REF[x], domain1_SVs$PR_ALT[x], domain1_SVs$SR_ALT[x], na.rm = T)
})

# Mean VAF by SV type
domain1_SVs %>% group_by(SVTYPE) %>% summarise(mean(VAF), median(VAF))




#### Annotate Domain 1 SVs with REF/ALT #### 

### Read SV.VCFs from the local mount, local annotation (not HPC)

# Modify VCF paths to read from the local mount
FF_list$SV_vcf_path_local <- paste0("/Users/MartinaMijuskovic", FF_list$SV_vcf_path)

# Test annotation (ok!)
#domain1_SVs_test <- lapply(FF_list$SV_vcf_path_local[1:2], getSVsInTranscripts_v3_local, GO_transcr)

# Full annotation job
domain1_SVs_fix1 <- lapply(FF_list$SV_vcf_path_local[1:10], getSVsInTranscripts_v3_local, GO_transcr)
domain1_SVs_fix1 <- bind_rows(domain1_SVs_fix1)
domain1_SVs_fix1 <- domain1_SVs_fix1 %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))

domain1_SVs_fix2 <- lapply(FF_list$SV_vcf_path_local[11:100], getSVsInTranscripts_v3_local, GO_transcr)
domain1_SVs_fix2 <- bind_rows(domain1_SVs_fix2)
domain1_SVs_fix2 <- domain1_SVs_fix2 %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))

domain1_SVs_fix3 <- lapply(FF_list$SV_vcf_path_local[101:400], getSVsInTranscripts_v3_local, GO_transcr)
domain1_SVs_fix3 <- bind_rows(domain1_SVs_fix3)
domain1_SVs_fix3 <- domain1_SVs_fix3 %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))


domain1_SVs_fix4 <- lapply(FF_list$SV_vcf_path_local[401:759], getSVsInTranscripts_v3_local, GO_transcr)
domain1_SVs_fix4 <- bind_rows(domain1_SVs_fix4)
domain1_SVs_fix4 <- domain1_SVs_fix4 %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))


dim(domain1_SVs_fix1)  # 10
dim(domain1_SVs_fix2)  # 88
dim(domain1_SVs_fix3)  # 422
dim(domain1_SVs_fix4)  # 502

# Merge all SVs
domain1_SVs_fix <- rbind(domain1_SVs_fix1, domain1_SVs_fix2, domain1_SVs_fix3, domain1_SVs_fix4)
rm(domain1_SVs_fix1, domain1_SVs_fix2, domain1_SVs_fix3, domain1_SVs_fix4)

# Summary ((NOTE that any DEL >10kb is automatically filtered out and not included here)
dim(domain1_SVs_fix)  # 1022 total (TRs listed as 2 BNDs)
table(domain1_SVs_fix$SVTYPE)  # 719 Domain 1 SVs (~1 per patient)
# BND DEL DUP INS INV 
# 606 128  27  18 243

# Proportion of different SV types
table(domain1_SVs_fix$SVTYPE)/719
# BND        DEL        DUP        INS        INV 
# 0.84283727 0.17802503 0.03755216 0.02503477 0.33796940  # transl are 0.4214186

# Median SV per patient
median(as.numeric(table(domain1_SVs_fix$SAMPLE_WELL_ID)))  # 2


# Number of patients with Domain 1 SVs
length(unique(domain1_SVs_fix$SAMPLE_WELL_ID))  # 311/759 (41% patients have Domain 1 SVs)
table(domain1_SVs_fix$SVTYPE, domain1_SVs_fix$Application)  # No Canvas
table(domain1_SVs_fix$FILTER, domain1_SVs_fix$SOMATIC)
table(domain1_SVs_fix$SVTYPE, domain1_SVs_fix$IMPRECISE)
# FALSE TRUE
# BND   464  142
# DEL   127    1
# DUP    27    0
# INS    18    0
# INV   201   42
domain1_SVs_fix %>% dplyr::select(SAMPLE_WELL_ID, ID, IMPRECISE, SVTYPE, SVLEN, CHR, START, END, ColocalizedCanvas, PR_ALT, SR_ALT, start_ann, end_ann)


# Add KEY, calculate recurrent variants, add VAF and variant frequency info to the table
domain1_SVs_fix$KEY <- sapply(1:dim(domain1_SVs_fix)[1], function(x){
  paste(domain1_SVs_fix$CHR[x], domain1_SVs_fix$START[x], domain1_SVs_fix$END[x], domain1_SVs_fix$SVTYPE[x], sep = "-")
})

unique_keys <- unique(domain1_SVs_fix$KEY)
length(unique_keys)  # 1004
recurrent_SVs <- data.frame(
  KEY = unique_keys,
  NUM_OBS = sapply(unique_keys, function(x){ sum(domain1_SVs$KEY == x)})
)
rownames(recurrent_SVs) <- NULL
# Look at recurrent SVs
recurrent_SVs %>% filter(NUM_OBS > 1) %>% arrange(desc(NUM_OBS))

# Add variant frequency info to the main table
recurr_keys <- as.character(recurrent_SVs %>% filter(NUM_OBS >1) %>% pull(KEY))
domain1_SVs_fix$recurrent <- 0  
domain1_SVs_fix[domain1_SVs_fix$KEY %in% recurr_keys,]$recurrent <- 1  # 9 total (3 deletions, 3 transl)
table(domain1_SVs_fix$recurrent, domain1_SVs_fix$SVTYPE)
domain1_SVs_fix %>% filter(recurrent == 1, SVTYPE == "DEL")
domain1_SVs_fix %>% filter(recurrent == 1, KEY == "4-54688371-54688815-DEL")
domain1_SVs_fix %>% filter(recurrent == 1, KEY == "12-68810707-68811111-DEL")
domain1_SVs_fix %>% filter(recurrent == 1, KEY == "5-142596308-142596582-DEL")
domain1_SVs_fix %>% filter(recurrent == 1, SVTYPE == "BND")
domain1_SVs_fix %>% filter(recurrent == 1, KEY == "10-43115587-NA-BND")
domain1_SVs_fix %>% filter(recurrent == 1, ID == "MantaBND:114898:0:1:0:0:0:0") # 2 SVs in same EVENT, sharing one breakopoint, but mate not recurrent
domain1_SVs_fix %>% filter(recurrent == 1, ID == "MantaBND:114898:0:1:1:3:0:1") # 2 SVs in same EVENT, sharing one breakopoint, but mate not recurrent

### Remove BND with missing mates and recalculate summaries
domain1_SVs_fix$missing_mate <- 0
missing_mates <- sapply(unique(domain1_SVs_fix$SAMPLE_WELL_ID), function(x){
  mates <- domain1_SVs_fix %>% filter(SAMPLE_WELL_ID == x, SVTYPE == "BND") %>% pull(MATEID)
  # List missing mates
  mates[!mates %in% (domain1_SVs_fix %>% filter(SAMPLE_WELL_ID == x, SVTYPE == "BND") %>% pull(ID))]
})
table(as.character(missing_mates))  # none missing


### Calculate VAF
domain1_SVs_fix$VAF <- sapply(1:dim(domain1_SVs_fix)[1], function(x){
  sum(domain1_SVs_fix$PR_ALT[x], domain1_SVs_fix$SR_ALT[x], na.rm = T)/sum(domain1_SVs_fix$PR_REF[x], domain1_SVs_fix$SR_REF[x], domain1_SVs_fix$PR_ALT[x], domain1_SVs_fix$SR_ALT[x], na.rm = T)
})

# Mean VAF by SV type
domain1_SVs_fix %>% group_by(SVTYPE) %>% summarise(mean(VAF), median(VAF))


# Unlist variables
domain1_SVs_fix$SVINSLEN <- as.numeric(domain1_SVs_fix$SVINSLEN)
domain1_SVs_fix$SVLEN <- as.numeric(domain1_SVs_fix$SVLEN)
domain1_SVs_fix$CIGAR <- as.character(domain1_SVs_fix$CIGAR)
# domain1_SVs_fix$CIPOS <- as.character(domain1_SVs_fix$CIPOS)   # info loss
# domain1_SVs_fix$CIEND <- as.character(domain1_SVs_fix$CIEND)   # info loss
domain1_SVs_fix$HOMLEN <- as.numeric(domain1_SVs_fix$HOMLEN)
domain1_SVs_fix$HOMSEQ <- as.character(domain1_SVs_fix$HOMSEQ)
domain1_SVs_fix$SVINSSEQ <- as.character(domain1_SVs_fix$SVINSSEQ)
domain1_SVs_fix$LEFT_SVINSSEQ <- as.character(domain1_SVs_fix$LEFT_SVINSSEQ)
domain1_SVs_fix$RIGHT_SVINSSEQ <- as.character(domain1_SVs_fix$RIGHT_SVINSSEQ)



#### Add directionality flags ####

# For INV, directionality flags already exist in the table (2 possibilities: INV3, INV5)
# For BND, directionality flags need to be derived from the ALT structure (is relevant region upstream or downstream of POS?)

# Add BND direction flag (equivalent to INV, BND3=chr upstream of POS is included, BND5=chr downstream of POS is included in translocation)
domain1_SVs_fix$BND5 <- grepl("^]|^\\[", domain1_SVs_fix$ALT)   # ALT starts with brackets ([ or ]), indicating that POS is the 5' end of transclocation junction
domain1_SVs_fix$BND3 <- (grepl("]|\\[", domain1_SVs_fix$ALT) & !grepl("^]|^\\[", domain1_SVs_fix$ALT))  # ALT has brackets, but doesn't start with them
# Sanity check
head((domain1_SVs_fix %>% filter(SVTYPE == "BND") %>% dplyr::select(KEY, ALT, BND3, BND5)),30)
head((domain1_SVs_fix %>% filter(SVTYPE != "BND") %>% dplyr::select(KEY, ALT, BND3, BND5)),30)



#### Mappablity and repeat overlap ####

# Average fragment size for FF PCR-free
mean((QC %>% filter(LIBRARY_TYPE == "TruSeq PCR-Free") %>% pull(AV_FRAGMENT_SIZE_BP)), na.rm = T) # 483.1857
median((QC %>% filter(LIBRARY_TYPE == "TruSeq PCR-Free") %>% pull(AV_FRAGMENT_SIZE_BP)), na.rm = T)  # 482.3

# Assess IMPRECISE flag against SR support
table(domain1_SVs_fix$IMPRECISE, domain1_SVs_fix$SR_ALT == 0, exclude = NULL)  # imprecise have SR_ALT = NA; some precise have SR_ALT = 0
#           FALSE TRUE <NA>
# FALSE     833    4    0
# TRUE      0      0   185
# <NA>      0      0    0
domain1_SVs_fix %>% filter(IMPRECISE == FALSE, SR_ALT == 0) %>% dplyr::select(IMPRECISE, SVTYPE, CIPOS, ID, MATEID, HOMLEN, ALT, KEY, PR_REF, PR_ALT, SR_REF, SR_ALT, SAMPLE_WELL_ID)

# CIPOS, CIEND of precise and imprecise variants
table(domain1_SVs_fix$SVTYPE, domain1_SVs_fix$IMPRECISE)
# imprecise (mostly BND and INV, just 1 deletion)
mean(abs(unlist(domain1_SVs_fix[domain1_SVs_fix$IMPRECISE == TRUE,]$CIPOS)))  # 343.6405
# precise
mean(abs(unlist(domain1_SVs_fix[domain1_SVs_fix$IMPRECISE == FALSE,]$CIPOS)), na.rm = T)  #  2.667535
mean(abs(unlist(domain1_SVs_fix[domain1_SVs_fix$IMPRECISE == FALSE,]$CIEND)), na.rm = T)  #  2.18617


# Create a bed file with SVs, start and end separately, remove NAs
# Note that for BED format START has to be adjusted to 0-based (-1 bp) and END position is not included (stays same)
# Add 240 bp to the breakpoint (region defined based on SV type and directionality)

# Change START, END class to numeric
domain1_SVs_fix$START <- as.numeric(domain1_SVs_fix$START)
domain1_SVs_fix$END <- as.numeric(domain1_SVs_fix$END)


# Bed file for START position
start_bed <- domain1_SVs_fix %>% mutate(
  START_bed = case_when(
    SVTYPE == "DEL" ~ START-240,
    SVTYPE == "INS" ~ START-240,
    SVTYPE == "DUP" ~ START,
    (SVTYPE == "INV" & INV3 == TRUE) ~ START-240,
    (SVTYPE == "INV" & INV5 == TRUE) ~ START,
    (SVTYPE == "BND" & BND3 == TRUE) ~ START-240,
    (SVTYPE == "BND" & BND5 == TRUE) ~ START
    ),
  END_bed = case_when(
    SVTYPE == "DEL" ~ START,
    SVTYPE == "INS" ~ START,
    SVTYPE == "DUP" ~ START+240,
    SVTYPE == "INV" & INV3 == TRUE ~ START,
    SVTYPE == "INV" & INV5 == TRUE ~ START+240,
    SVTYPE == "BND" & BND3 == TRUE ~ START,
    SVTYPE == "BND" & BND5 == TRUE ~ START+240
  )
) %>% dplyr::select(CHR, START_bed, END_bed, KEY)

# Bed file for END position
end_bed <- domain1_SVs_fix %>% mutate(
  START_bed = case_when(
    SVTYPE == "DEL" ~ END,
    SVTYPE == "INS" ~ END,
    SVTYPE == "DUP" ~ END-240,
    SVTYPE == "INV" & INV3 == TRUE ~ END-240,
    SVTYPE == "INV" & INV5 == TRUE ~ END,
    SVTYPE == "BND" ~ 1000  # placeholder
    ),
  END_bed = case_when(
    SVTYPE == "DEL" ~ END+240,
    SVTYPE == "INS" ~ END+240,
    SVTYPE == "DUP" ~ END,
    SVTYPE == "INV" & INV3 == TRUE ~ END,
    SVTYPE == "INV" & INV5 == TRUE ~ END+240,
    SVTYPE == "BND" ~ 1300  # placeholder
  )
) %>% dplyr::select(CHR, START_bed, END_bed, KEY)

# Add "chr" to chrom name
start_bed$CHR <- paste0("chr", start_bed$CHR)
end_bed$CHR <- paste0("chr", end_bed$CHR)

# Write bed files
write.table(start_bed, file = "./Data/mappability_repeats/start.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(end_bed, file = "./Data/mappability_repeats/end.bed", quote = F, row.names = F, col.names = F, sep = "\t")



### Call bedtools to find overlaps with Umap36 

# start
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/start.bed", "-b ./Data/mappability_repeats/mappable_k36.umap.bedgraph > ./Data/mappability_repeats/Umap36_start_overlap.bed"), intern = T)
Umap36_start <- read.table("./Data/mappability_repeats/Umap36_start_overlap.bed", sep = "\t")
names(Umap36_start) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/end.bed", "-b ./Data/mappability_repeats/mappable_k36.umap.bedgraph > ./Data/mappability_repeats/Umap36_end_overlap.bed"), intern = T)
Umap36_end <- read.table("./Data/mappability_repeats/Umap36_end_overlap.bed", sep = "\t")
names(Umap36_end) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")




### Call bedtools to find overlaps with WindowMasker 

# start
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/start.bed", "-b ./Data/mappability_repeats/windowmaskerSdust.hg38.bed > ./Data/mappability_repeats/wMasker_start_overlap.bed"), intern = T)
wMasker_start <- read.table("./Data/mappability_repeats/wMasker_start_overlap.bed", sep = "\t")
names(wMasker_start) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/end.bed", "-b ./Data/mappability_repeats/windowmaskerSdust.hg38.bed > ./Data/mappability_repeats/wMasker_end_overlap.bed"), intern = T)
wMasker_end <- read.table("./Data/mappability_repeats/wMasker_end_overlap.bed", sep = "\t")
names(wMasker_end) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")




### Call bedtools to find overlaps with simple repeats

# start
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/start.bed", "-b ./Data/mappability_repeats/simpleRepeat.hg38.bed > ./Data/mappability_repeats/simpleRepeat_start_overlap.bed"), intern = T)
simpleRepeat_start <- read.table("./Data/mappability_repeats/simpleRepeat_start_overlap.bed", sep = "\t")
names(simpleRepeat_start) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/end.bed", "-b ./Data/mappability_repeats/simpleRepeat.hg38.bed > ./Data/mappability_repeats/simpleRepeat_end_overlap.bed"), intern = T)
simpleRepeat_end <- read.table("./Data/mappability_repeats/simpleRepeat_end_overlap.bed", sep = "\t")
names(simpleRepeat_end) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")




### Call bedtools to find overlaps with segmental duplications

# start
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/start.bed", "-b ./Data/mappability_repeats/genomicSuperDups.hg38.bed > ./Data/mappability_repeats/segDups_start_overlap.bed"), intern = T)
segDups_start <- read.table("./Data/mappability_repeats/segDups_start_overlap.bed", sep = "\t")
names(segDups_start) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system(paste("/usr/local/bin/bedtools coverage -a", "./Data/mappability_repeats/end.bed", "-b ./Data/mappability_repeats/genomicSuperDups.hg38.bed > ./Data/mappability_repeats/segDups_end_overlap.bed"), intern = T)
segDups_end <- read.table("./Data/mappability_repeats/segDups_end_overlap.bed", sep = "\t")
names(segDups_end) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")




# Add overlap percentage

domain1_SVs_fix$Umap36_PCT_overlap_START <- Umap36_start[match(domain1_SVs_fix$KEY, Umap36_start$KEY),]$PCT
domain1_SVs_fix$Umap36_PCT_overlap_END <- Umap36_end[match(domain1_SVs_fix$KEY, Umap36_end$KEY),]$PCT

domain1_SVs_fix$wMasker_PCT_overlap_START <- wMasker_start[match(domain1_SVs_fix$KEY, wMasker_start$KEY),]$PCT
domain1_SVs_fix$wMasker_PCT_overlap_END <- wMasker_end[match(domain1_SVs_fix$KEY, wMasker_end$KEY),]$PCT

domain1_SVs_fix$simpleRepeat_PCT_overlap_START <- simpleRepeat_start[match(domain1_SVs_fix$KEY, simpleRepeat_start$KEY),]$PCT
domain1_SVs_fix$simpleRepeat_PCT_overlap_END <- simpleRepeat_end[match(domain1_SVs_fix$KEY, simpleRepeat_end$KEY),]$PCT

domain1_SVs_fix$segDups_PCT_overlap_START <- segDups_start[match(domain1_SVs_fix$KEY, segDups_start$KEY),]$PCT
domain1_SVs_fix$segDups_PCT_overlap_END <- segDups_end[match(domain1_SVs_fix$KEY, segDups_end$KEY),]$PCT

# Set BND END to mappable and repeat overlap as none (placeholder)
domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$Umap36_PCT_overlap_END <- 1
domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$wMasker_PCT_overlap_END <- 0
domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$simpleRepeat_PCT_overlap_END <- 0
domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$segDups_PCT_overlap_END <- 0
  
  
# Add flag for high confidence variants (100% mappability, or 100% mappability and no repeat overlap)
domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Mappable = case_when(Umap36_PCT_overlap_START == 1 & Umap36_PCT_overlap_END == 1 ~ 1, TRUE ~ 0),
                                              Mappable_noRepeats = case_when(Umap36_PCT_overlap_START == 1 & Umap36_PCT_overlap_END == 1 & 
                                                                               wMasker_PCT_overlap_START == 0 & wMasker_PCT_overlap_END == 0 & 
                                                                               simpleRepeat_PCT_overlap_START == 0 & simpleRepeat_PCT_overlap_END == 0 &
                                                                               segDups_PCT_overlap_START == 0 & segDups_PCT_overlap_END ==0 ~ 1, TRUE ~ 0)
)


# Add another flag for mappable+no repeats that excludes windowMasker (windowMasker is over-filtering)
domain1_SVs_fix <- domain1_SVs_fix %>% mutate(
                                              Mappable_noRepeats_2 = case_when(Umap36_PCT_overlap_START == 1 & Umap36_PCT_overlap_END == 1 & 
                                                                               simpleRepeat_PCT_overlap_START == 0 & simpleRepeat_PCT_overlap_END == 0 &
                                                                               segDups_PCT_overlap_START == 0 & segDups_PCT_overlap_END ==0 ~ 1, TRUE ~ 0)
)




### Add mappability for BND mate

# Add MATE_KEY to the table for BND

# Get unmappable BND
unmap_BND <- domain1_SVs_fix %>% filter(SVTYPE == "BND", Mappable == 0) %>% pull(KEY)  # 329
# Check if recurrent also have same directionality (yes)
sum(duplicated(unmap_BND))  # 6
domain1_SVs_fix %>% filter(KEY %in% unmap_BND[duplicated(unmap_BND)]) %>% dplyr::select(KEY, SAMPLE_WELL_ID, VAF,  BND5,  BND3, PR_ALT, SR_ALT) %>% arrange(KEY)

# Get mappable + no repeats BND
map_norep_BND <- domain1_SVs_fix %>% filter(SVTYPE == "BND", Mappable_noRepeats == 1) %>% pull(KEY)  # 57
map_norep_BND_2 <- domain1_SVs_fix %>% filter(SVTYPE == "BND", Mappable_noRepeats_2 == 1) %>% pull(KEY)  # 267


# Add MATE_KEY
domain1_SVs_fix$MATE_KEY <- NA
domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$MATE_KEY <- domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",][match(domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$MATEID, domain1_SVs_fix[domain1_SVs_fix$SVTYPE == "BND",]$ID),]$KEY
# Sanity check
head(domain1_SVs_fix %>% dplyr::select(SAMPLE_WELL_ID, KEY, ID, MATEID, MATE_KEY) %>% arrange(SAMPLE_WELL_ID),20)

# Add MATE_Mappable and MATE_Mappable_noRepeats, MATE_Mappable_noRepeats_2
domain1_SVs_fix <- domain1_SVs_fix %>% mutate(
  MATE_Mappable = case_when(
    SVTYPE == "BND" & (MATE_KEY %in% unmap_BND) ~ 0,
    SVTYPE == "BND" & (!MATE_KEY %in% unmap_BND) ~ 1,
    SVTYPE != "BND" ~ 0
  ),
  MATE_Mappable_noRepeats = case_when(
    SVTYPE == "BND" & (MATE_KEY %in% map_norep_BND) ~ 1,
    SVTYPE == "BND" & (!MATE_KEY %in% map_norep_BND) ~ 0,
    SVTYPE != "BND" ~ 0
  ),
  MATE_Mappable_noRepeats_2 = case_when(
    SVTYPE == "BND" & (MATE_KEY %in% map_norep_BND_2) ~ 1,
    SVTYPE == "BND" & (!MATE_KEY %in% map_norep_BND_2) ~ 0,
    SVTYPE != "BND" ~ 0
  )
)

# Set non-BND to NA
domain1_SVs_fix[domain1_SVs_fix$SVTYPE != "BND",]$MATE_Mappable <- NA
domain1_SVs_fix[domain1_SVs_fix$SVTYPE != "BND",]$MATE_Mappable_noRepeats <- NA
domain1_SVs_fix[domain1_SVs_fix$SVTYPE != "BND",]$MATE_Mappable_noRepeats_2 <- NA



#### High confidence SV analysis #### 

# SVs with 100% mappable regions (both breakpoints)
dim(domain1_SVs_fix %>% filter(Mappable == 1, (is.na(MATE_Mappable) | (MATE_Mappable == 1) )))  # 297, 151 non-BND plus 146 BND (73 TR) with mate mappable
# Mappable by type
table(domain1_SVs_fix$Mappable, domain1_SVs_fix$SVTYPE, domain1_SVs_fix$MATE_Mappable, exclude = NULL)


# SVs with 100% mappable regions and no TRF or segDup overlap
dim(domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, (is.na(MATE_Mappable_noRepeats_2) | (MATE_Mappable_noRepeats_2 == 1))  ))  # 279 
# Mappable and no repeats by type
table(domain1_SVs_fix$Mappable_noRepeats_2, domain1_SVs_fix$SVTYPE, domain1_SVs_fix$MATE_Mappable_noRepeats_2, exclude = NULL)  # 138 BND (69 TR) plus 141


# SVs with 100% mappable regions and with some TRF or segDup overlap
dim(domain1_SVs_fix %>% filter(Mappable == 1, Mappable_noRepeats_2 == 0 , SVTYPE != "BND"))  #  10
domain1_SVs_fix %>% filter(Mappable == 1, MATE_Mappable == 1, (Mappable_noRepeats_2 == 0 | MATE_Mappable_noRepeats_2 == 0), SVTYPE == "BND")  # 4 with both BND mappable and 1 or both overlapping repeats

# SVs NOT 100% mappable but no TRF or segDup overlap

dim(domain1_SVs_fix %>% filter(SVTYPE != "BND", Mappable == 0, simpleRepeat_PCT_overlap_START == 0, simpleRepeat_PCT_overlap_END == 0, segDups_PCT_overlap_START == 0, segDups_PCT_overlap_END == 0))  # 176
table((domain1_SVs_fix %>% filter(SVTYPE != "BND", Mappable == 0, simpleRepeat_PCT_overlap_START == 0, simpleRepeat_PCT_overlap_END == 0, segDups_PCT_overlap_START == 0, segDups_PCT_overlap_END == 0)) %>%  pull(SVTYPE))

dim(domain1_SVs_fix %>% filter(SVTYPE == "BND", (Mappable == 0 | MATE_Mappable == 0), simpleRepeat_PCT_overlap_START == 0, segDups_PCT_overlap_START == 0))  # 164 TR



######### Manual review lists

#### non-BND

# List of mappable with no repeat overlap (no BND): super-high confidence candidates  (6)
domain1_SVs_fix %>% filter(Mappable_noRepeats == 1, SVTYPE != "BND") %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS,CIEND, SOMATICSCORE)
# Mean VAF
mean((domain1_SVs_fix %>% filter(Mappable_noRepeats == 1, SVTYPE != "BND") %>% pull(VAF)))  # 0.1855203



# List of mappable with no TRF, segDup repeat overlap (no BND): high confidence candidates (141)
as.data.frame(domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, SVTYPE != "BND") %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS,CIEND, SOMATICSCORE))[1:76,]
as.data.frame(domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, SVTYPE != "BND") %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS,CIEND, SOMATICSCORE))[77:141,]
# Mean VAF
mean((domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, SVTYPE != "BND") %>% pull(VAF)))  # 0.1634514



#### BND

# List of mappable with no repeat overlap (BND): super-high confidence candidates (2 TR)
domain1_SVs_fix %>% filter(Mappable_noRepeats == 1, MATE_Mappable_noRepeats == 1, SVTYPE == "BND") %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, SOMATICSCORE)
# Mean VAF
mean((domain1_SVs_fix %>% filter(Mappable_noRepeats == 1, MATE_Mappable_noRepeats == 1, SVTYPE == "BND") %>% pull(VAF)))  # 0.02308771


# List of mappable with no TRF, segDup repeat overlap (BND): high confidence candidates (69 TR)
as.data.frame(domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, MATE_Mappable_noRepeats_2 == 1, SVTYPE == "BND") %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, SOMATICSCORE))[1:71,]
as.data.frame(domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, MATE_Mappable_noRepeats_2 == 1, SVTYPE == "BND") %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, SOMATICSCORE))[72:138,]
# Mean VAF
mean((domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, MATE_Mappable_noRepeats_2 == 1, SVTYPE == "BND") %>% pull(VAF)))  # 0.1075102


##### Mean VAF and imprecise summary 

# Mean VAF for all super high conf (adding 2 transl manually)
mean(c((domain1_SVs_fix %>% filter(Mappable_noRepeats == 1, SVTYPE != "BND") %>% pull(VAF)), 0.02614379, 0.02003163))  # 0.1449121

# Mean VAF for all high conf (0.1450707)
mean(c((domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, SVTYPE != "BND") %>% pull(VAF)), ((domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, MATE_Mappable_noRepeats_2 == 1, SVTYPE == "BND") %>% pull(VAF))[!duplicated(domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, MATE_Mappable_noRepeats_2 == 1, SVTYPE == "BND") %>% pull(VAF))])))  

# % imprecise in high conf (wrong, BND mates not taken into account, calc manually) 10/69 = 0.1449275
#table(domain1_SVs_fix[domain1_SVs_fix$Mappable_noRepeats_2 == 1,]$SVTYPE, domain1_SVs_fix[domain1_SVs_fix$Mappable_noRepeats_2 == 1,]$IMPRECISE)

# View imprecise high conf SVs
domain1_SVs_fix %>% filter(Mappable_noRepeats_2 == 1, MATE_Mappable_noRepeats_2 == 1 | is.na(MATE_Mappable_noRepeats_2), IMPRECISE == TRUE) %>% dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, SVINSLEN, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, SOMATICSCORE)



#### Mappable with repeat overlap (TRF, segDup)

# Non-BND
domain1_SVs_fix %>% filter(Mappable == 1, Mappable_noRepeats_2 == 0 , SVTYPE != "BND") %>% 
  dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, CIEND, SOMATICSCORE, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)

# BND
domain1_SVs_fix %>% filter(Mappable == 1, MATE_Mappable == 1, (Mappable_noRepeats_2 == 0 | MATE_Mappable_noRepeats_2 == 0), SVTYPE == "BND")  %>% 
  dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, CIEND, SOMATICSCORE, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)



######### Tiering SVs by confidence #########

# Add tiering flags:

# TIER 1: Super-high confidence
# • both breakpoints 100% mappable
# • no simple repeat, segmental duplication nor windowMasker overlap

# TIER 2: High confidence
# • both breakpoints 100% mappable
# • no simple repeat nor segmental duplication (windowMasker overlap allowed)

# TIER 3: Medium confidence
# • both breakpoints 100% mappable
# • some overlap with simple repeats or segmental duplication in one or both breakpoints

# TIER 4: Low confidence
# • only one breakpoint 100% mappable

# TIER 5: Super-low confidence
# • neither breakpoint 100% mappable


domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Tier1 = case_when(
  (Mappable_noRepeats == 1 & SVTYPE != "BND")  ~ 1,
  (Mappable_noRepeats == 1 & MATE_Mappable_noRepeats == 1 & SVTYPE == "BND") ~ 1,
  TRUE ~ 0
))

domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Tier2 = case_when(
  (Mappable_noRepeats_2 == 1 & SVTYPE != "BND" & Tier1 == 0)  ~ 1,
  (Mappable_noRepeats_2 == 1 & MATE_Mappable_noRepeats_2 == 1 & SVTYPE == "BND" & Tier1 == 0) ~ 1,
  TRUE ~ 0
))

domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Tier3 = case_when(
  (Mappable == 1 & SVTYPE != "BND" & Tier1 == 0 & Tier2 == 0)  ~ 1,
  (Mappable == 1 & MATE_Mappable == 1 & SVTYPE == "BND" & Tier1 == 0 & Tier2 == 0) ~ 1,
  TRUE ~ 0
))

domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Tier4 = case_when(
  (SVTYPE != "BND" & (Umap36_PCT_overlap_START != 1 & Umap36_PCT_overlap_END == 1))  ~ 1,
  (SVTYPE != "BND" & (Umap36_PCT_overlap_START == 1 & Umap36_PCT_overlap_END != 1))  ~ 1,
  (Mappable == 1 & MATE_Mappable == 0 & SVTYPE == "BND") ~ 1,
  (Mappable == 0 & MATE_Mappable == 1 & SVTYPE == "BND") ~ 1,
  TRUE ~ 0
))

domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Tier5 = case_when(
  (Tier1 == 0 & Tier2 == 0 & Tier3 == 0 & Tier4 == 0) ~ 1,
  TRUE ~ 0
))

# Add Tier field
domain1_SVs_fix <- domain1_SVs_fix %>% mutate(Tier = case_when(
  Tier1 == 1 ~ "Tier 1",
  Tier2 == 1 ~ "Tier 2",
  Tier3 == 1 ~ "Tier 3",
  Tier4 == 1 ~ "Tier 4",
  Tier5 == 1 ~ "Tier 5",
  TRUE ~ ""
))



# Sanity checks
table(domain1_SVs_fix$Tier, exclude = NULL)
table(domain1_SVs_fix$Tier, domain1_SVs_fix$SVTYPE)

table(domain1_SVs_fix$Tier1, domain1_SVs_fix$SVTYPE)
table(domain1_SVs_fix$Tier2, domain1_SVs_fix$SVTYPE)
table(domain1_SVs_fix$Tier3, domain1_SVs_fix$SVTYPE)
table(domain1_SVs_fix$Tier4, domain1_SVs_fix$SVTYPE)
table(domain1_SVs_fix$Tier5, domain1_SVs_fix$SVTYPE)

domain1_SVs_fix %>% filter(Tier3 == 1) %>% dplyr::select(SAMPLE_WELL_ID, KEY, MATE_KEY, Mappable, MATE_Mappable, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END, Tier1, Tier2, Tier3)
domain1_SVs_fix %>% filter(Tier4 == 1, SVTYPE != "BND") %>% dplyr::select(SAMPLE_WELL_ID, KEY, Mappable, Umap36_PCT_overlap_START, Umap36_PCT_overlap_END, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)
domain1_SVs_fix %>% filter(Tier4 == 1, SVTYPE == "BND") %>% dplyr::select(SAMPLE_WELL_ID, KEY, MATE_KEY, Mappable, MATE_Mappable, Umap36_PCT_overlap_START, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)
domain1_SVs_fix %>% filter(Tier5 == 1, SVTYPE != "BND") %>% dplyr::select(SAMPLE_WELL_ID, KEY, Mappable, Umap36_PCT_overlap_START, Umap36_PCT_overlap_END, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)
domain1_SVs_fix %>% filter(Tier5 == 1, SVTYPE == "BND") %>% dplyr::select(SAMPLE_WELL_ID, KEY, MATE_KEY, Mappable, MATE_Mappable, Umap36_PCT_overlap_START, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)


# List of Tier 4 non-TR SVs
domain1_SVs_fix %>% filter(Tier == "Tier 4", SVTYPE != "BND") %>% 
  dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, CIEND, SOMATICSCORE, Umap36_PCT_overlap_START, Umap36_PCT_overlap_END, wMasker_PCT_overlap_START, wMasker_PCT_overlap_END, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)

# List of Tier 4 TR SVs
domain1_SVs_fix %>% filter(Tier == "Tier 4", SVTYPE == "BND") %>% 
  dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, CIPOS, SOMATICSCORE, Umap36_PCT_overlap_START, wMasker_PCT_overlap_START, simpleRepeat_PCT_overlap_START, segDups_PCT_overlap_START)


# List of Tier 5 non-TR SVs
domain1_SVs_fix %>% filter(Tier == "Tier 5", SVTYPE != "BND") %>% 
  dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, SVLEN, CIPOS, CIEND, SOMATICSCORE, Umap36_PCT_overlap_START, Umap36_PCT_overlap_END, wMasker_PCT_overlap_START, wMasker_PCT_overlap_END, simpleRepeat_PCT_overlap_START, simpleRepeat_PCT_overlap_END, segDups_PCT_overlap_START, segDups_PCT_overlap_END)


# List of Tier 5 TR SVs
domain1_SVs_fix %>% filter(Tier == "Tier 5", SVTYPE == "BND") %>% 
  dplyr::select(SAMPLE_WELL_ID, IMPRECISE, KEY, MATE_KEY, VAF, recurrent, PR_REF, PR_ALT, SR_REF, SR_ALT, CIPOS, SOMATICSCORE, Umap36_PCT_overlap_START, wMasker_PCT_overlap_START, simpleRepeat_PCT_overlap_START, segDups_PCT_overlap_START)









