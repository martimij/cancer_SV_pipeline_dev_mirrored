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
actionable_SVs <- read.table("/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Data/GO/GENOMONCOLOGY_SOLID_TUMOUR.SV.v1.6.tsv", sep = "\t", header = T)
head(actionable_SVs)  # looks ok
length(unique(actionable_SVs$gene_name))  # 106 genes (!?)
length(unique(actionable_SVs$transcript_ID))  # 103 unique transcripts
length(unique(actionable_SVs$gene_ID))  # 103 unique gene IDs
dim(actionable_SVs) # 1421 
sum(duplicated(actionable_SVs))  # 38 duplicate rows
actionable_SVs[duplicated(actionable_SVs),]
# Deduplicate
actionable_SVs <- actionable_SVs[!duplicated(actionable_SVs),]
dim(actionable_SVs)  # 1383
# Write out
write.table(actionable_SVs, file = "./Data/GO/GENOMONCOLOGY_SOLID_TUMOUR.dedup.SV.v1.6.tsv", quote = F, row.names = F, col.names = T, sep = "\t")


# Remove missing fields
actionable_SVs <- actionable_SVs %>% filter(transcript_ID == "")
GO_transcr <- unique(actionable_SVs$transcript_ID)
GO_transcr <- GO_transcr[GO_transcr != ""]
length(GO_transcr)  # 102



####### Annotate SVs #######

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
### Keeps Canvas and Manta calls
### Adds a list of VCF filters to include (default is only "PASS") - NOTE that some might have missing mates/events (if filtered out)
### Removes non-SOMATIC (Canvas REF)

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
  
  vcf_info <- vcf_info %>% filter(FILTER %in% filters_keep, SOMATIC == "TRUE")
  
  
  
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


# Read the list of transcripts (HPC)
actionable_SVs <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/GENOMONCOLOGY_SOLID_TUMOUR.SV.v1.6.tsv", sep = "\t", header = T)

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



# Run (HPC) --- test
setwd("/home/mmijuskovic/SV_dev/SV_tiering/tiered_results/test")

library(dplyr)
library(VariantAnnotation)
library(ensembldb)

# Read the list of samples with paths to SV VCF
FF_list <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/FF_PCRfree_wClinical_forSVtiering_wPaths.tsv", sep = "\t", header = T)
FF_list$SV_vcf_path <- as.character(FF_list$SV_vcf_path)

# Read the list of transcripts (HPC)
actionable_SVs <- read.table("/home/mmijuskovic/SV_dev/SV_tiering/GENOMONCOLOGY_SOLID_TUMOUR.SV.v1.6.tsv", sep = "\t", header = T)
# Get unique list of transcripts
GO_transcr <- unique(actionable_SVs$transcript_ID)
GO_transcr <- GO_transcr[GO_transcr != ""]

domain1_SVs <- lapply(FF_list$SV_vcf_path[1:50], getSVsInTranscripts_v2, GO_transcr)  # job 1
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
  NUM_OBS = sapply(unique_keys, function(x){ sum(domain1_SVs$KEY == x)}),
  NUM_PATIENTS = sapply(unique_keys, function(x){ length(unique(domain1_SVs %>% filter(KEY == x) %>% .$SAMPLE_WELL_ID)) })
  )

# Add variant frequency info to the main table



### Annotation function adding repeat overlaps
### Works on SV VCF info table
### Overlaps with simpleRepeats and windowMasker









