# Martina Mijuskovic
# SV pipeline development
# Analysis of TMPRSS2-ERG fusions in prostate cancer cohort
# Extraction and analysis of TMPRSS2-ERG and related fusions
# Oct 2017

library(dplyr)
library(ggplot2)
library(data.table)
library(jsonlite)
library(VariantAnnotation)
library(ensembldb)

#### Install GRCh38 Ensembl transcriptome ####

# Get GRCh38 transcriptome from Bioconductor
# source("https://bioconductor.org/biocLite.R")


# Build Ensembl transcript database object (pre-downloaded GTF from here: ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz)
#biocLite("ensembldb")  # error with package "expands"
#library(ensembldb)

# From https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/
txdb_pth <- './Homo_sapiens.GRCh38.84.sqlite'
if (!file.exists(txdb_pth)) {
  txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
}
txdb <- EnsDb(txdb_pth)
txdb  # Preview the metadata


### Old strategy
# biocLite("BSgenome.Hsapiens.NCBI.GRCh38")  # install GRCh38  ----(error with package "expands")
# library(BSgenome.Hsapiens.NCBI.GRCh38)
# library(GenomicFeatures)
# biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)   #---> H. Sapiens transcript database
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene  #---> rename for convenience

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)
# Bigger plot text
bigger <- theme(legend.text=element_text(size=15), legend.title = element_text(size=15), axis.title = element_text(size=15), axis.text = element_text(size=15))
# Tilted x-axis labels
tiltedX <- theme(axis.text.x=element_text(angle=45,hjust=1))



######### Load sample manifest #########

samples <- read.csv("./Data/cancer_samples_prostate_cleaned.csv")


######### Annotate SVs (DEV) ######### 

# Get transcripts involved (regions)
transcripts_all <- read.table("./Data/LIST_OF_CANONICAL_TRANSCRIPTS.tsv", header = T)

# Transcripts of genes involved in fusions
transcripts_fusion <- transcripts_all %>% filter(gene_name %in% c("TMPRSS2", "ERG", "ETV1", "ETV4", "ETV5"))

# Read SV VCF
#vcf <- readVcf("./Data/test/LP3000178-DNA_E10_LP3000200-DNA_B03.somatic.SV.vcf.gz")
#vcf <- readVcf("./Data/test/LP3000169-DNA_C05_LP3000200-DNA_G02.somatic.SV.vcf.gz")


# Get VCF info fields
vcf_info <- as.data.frame(info(vcf))
# Get FILTER field
vcf_info$FILTER <- rowRanges(vcf)$FILTER
# Create Application variable (Canvas or Manta?)
vcf_info$Application <- ""
vcf_info[grepl("Canvas", rownames(vcf_info)),]$Application <- "Canvas"
vcf_info[grepl("Manta", rownames(vcf_info)),]$Application <- "Manta"
# Extract ID
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


### Look into Illumina's SV annotations (listing genes encompassing the whole length of the SV)

# ### Look for specific transcripts in the Illumina VCF annotation --- abandoning this strategy; only breakpoint annotations should be used
# vcf_info$TMPRSS2 <- grepl("ENST00000398585", vcf_info$CSQT)
# vcf_info$ERG <- grepl("ENST00000417133", vcf_info$CSQT)
# vcf_info$ETV1 <- grepl("ENST00000430479", vcf_info$CSQT)  #  ENST00000242066 is on our list of canonical transcripts
# vcf_info$ETV4 <- grepl("ENST00000591713", vcf_info$CSQT)
# vcf_info$ETV5 <- grepl("ENST00000306376", vcf_info$CSQT)
# # Summary
# sum(vcf_info$TMPRSS2)
# sum(vcf_info$ERG)
# sum(vcf_info$ETV1)
# sum(vcf_info$ETV4)
# sum(vcf_info$ETV5)
# # Extracting SVs potentially related to the fusion
# vcf_info_fusion <- vcf_info %>% filter(TMPRSS2 == 1 | ERG == 1 | ETV1 == 1 | ETV4 == 1 | ETV5 == 1)


# Get transcripts' locations from Ensembl transcript database
fusion_transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts_fusion$transcript_ID))
transcripts_fusion$gene_name <- as.character(transcripts_fusion$gene_name)
fusion_transcripts_gr$gene_name <- transcripts_fusion[match(names(fusion_transcripts_gr), transcripts_fusion$transcript_ID),]$gene_name

### Add transcript overlap annotations to SV breakpoints

# Change chr names in the VCF info table
vcf_info$CHR <- sub("chr", "", vcf_info$CHR)

# Overlap with START
# findOverlapPairs(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), fusion_transcripts_gr)
# vcf_info[overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), fusion_transcripts_gr),]
# mergeByOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), fusion_transcripts_gr)

# Annotate start
#vcf_info$TMPRSS2_ERG_related <- overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), fusion_transcripts_gr)

# Flag start
vcf_info$TMPRSS2_ERG_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), fusion_transcripts_gr))

# Flag end
#overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$END, width = 1)), fusion_transcripts_gr)
vcf_info$TMPRSS2_ERG_related_end <- sapply(1:dim(vcf_info)[1], function(x){
                                      if(is.na(vcf_info$END[x])){
                                        return(NA)
                                        } 
                                      else {
                                        if(overlapsAny(GRanges(seqnames=vcf_info$CHR[x], ranges=IRanges(start = vcf_info$END[x], width = 1)), fusion_transcripts_gr)) {
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

# Find overlaps of start and end position with transcripts
overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), fusion_transcripts_gr)
overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$TMPRSS2_ERG_related_end),]$CHR[x], ranges=IRanges(start = vcf_info[!is.na(vcf_info$TMPRSS2_ERG_related_end),]$END[x], width = 100000)), fusion_transcripts_gr)

# Add start position overlaps to table
if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
  vcf_info[queryHits(overlaps_start),]$start_ann <- sapply(seq(length(overlaps_start)), function(x){
                                                          i <- queryHits(overlaps_start)[x]
                                                          j <- subjectHits(overlaps_start)[x]
                                                          vcf_info[i,]$start_ann <- fusion_transcripts_gr@elementMetadata$gene_name[j]
                                                        })
}

# Add end position overlaps to table
if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
  vcf_info[queryHits(overlaps_end),]$end_ann <- sapply(seq(length(overlaps_end)), function(x){
    i <- queryHits(overlaps_end)[x]
    j <- subjectHits(overlaps_end)[x]
    vcf_info[i,]$end_ann <- fusion_transcripts_gr@elementMetadata$gene_name[j]
  })
}



##### Annotation function ##### 

### Function to extract VCF SVs with breakpoints overlapping specified transcripts and add gene name annotation
### Needs Homo_sapiens.GRCh38.84.sqlite in the home directory (see https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/)
### No filtering based on FILTER
### Canvas + Manta kept
### Paired SVs within same EVENT and BND mates kept (even if not overlapping transcripts)

# Use locally (fixed)
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
  txdb_pth <- './Homo_sapiens.GRCh38.84.sqlite'
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
    #extracted_SVs$SAMPLE_WELL_ID <- strsplit(strsplit(vcf_path, split = "Cancer")[[1]][2], split = "_Normal")[[1]][1]  # works for HPC path, not local
    # Use locally
    extracted_SVs$SAMPLE_WELL_ID <- paste0("LP", strsplit(strsplit(vcf_path, ".somatic")[[1]][1], "_LP")[[1]][2])
    
    return(extracted_SVs)
    }
}

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



# Test function
#getSVsInTranscripts("./Data/test/LP3000169-DNA_C05_LP3000200-DNA_G02.somatic.SV.vcf.gz", transcripts_fusion$transcript_ID)





### Run (HPC,  module load R/3.4.0 libxml2/2.9.4)

# Setup
#source("https://bioconductor.org/biocLite.R") # install ensembldb
#biocLite("ensembldb")   # install ensembldb
library(ensembldb)
library(VariantAnnotation)
library(dplyr)

setwd("/home/mmijuskovic/Fusions")

# Read sample manifest, remove excluded
samples <- read.csv("./cancer_samples_prostate_cleaned.csv")
samples <- samples %>% filter(Excluded == 0)  # 59 total
# Read transcripts, extract relevant ones
tr <- read.table("/home/mmijuskovic/Fusions/LIST_OF_CANONICAL_TRANSCRIPTS.tsv", sep = "\t", header = T)
tr_fusion <- tr %>% filter(gene_name %in% c("TMPRSS2", "ERG", "ETV1", "ETV4", "ETV5"))
# Get VCF paths
samples$VCF_path <- as.character(sapply(samples$Path, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))
# Remove sample with missing path
samples[samples$sampleId == "LP2000907-DNA_A02",]$Excluded <- 1
samples <- samples %>% filter(Excluded == 0)  # 58 total

### Run function to extract all SVs overlapping with TMRPSS2-ERG and related transcripts
fusions <- lapply(samples$VCF_path, getSVsInTranscripts, tr_fusion$transcript_ID)
#fusions <- rbind_all(fusions)  
fusions <- bind_rows(fusions) # ok
# Remove list fields
fusions <- fusions %>% dplyr::select(-(CSQR), -(CSQT), -(phyloP), -(EVS), -(clinvar), -(cosmic), -(AA), -(GMAF), -(AF1000G))
# Write Rdata
save.image("fusions.RData")
#write.table(fusions, file = "TMPRSS2_ERG_fusion_SVs_2017-10-23.tsv", sep = "\t", quote = F, row.names = F) # issue with lists


### Run locally

# Get local VCF paths
samples$local_VCF_path <- sapply(1:dim(samples)[1], function(x){
  command <- paste("find", paste0("/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Data/VCFs/", samples$Platekey_normal[x], "_", samples$sampleId[x], ".somatic.SV.vcf.gz"), sep = " ")
  system(command, intern = T)
})

# Extract SVs overlapping fusion-related transcripts ("samples" list is already based on non-excluded samples)
fusions <- lapply(samples$local_VCF_path, getSVsInTranscripts, transcripts_fusion$transcript_ID)
fusions <- rbind_all(fusions)



#### Add gene names to extracted fusions ####

# Reduce fusion table to relevant columns
fusions2 <- fusions %>% dplyr::select(-(LEFT_SVINSSEQ), -(RIGHT_SVINSSEQ), -(AF1000G), -(AA), -(GMAF), -(cosmic), -(clinvar), -(EVS), -(RefMinor), -(phyloP), -(CSQT), -(CSQR))

# Annotate fusion table (VCF info data frame) with all canonical ENSEMBL transcripts
annGenes <- function(vcf_info, transcripts){
  
  require(ensembldb)
  require(dplyr)
  
  txdb_pth <- '/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Homo_sapiens.GRCh38.84.sqlite'
  if (!file.exists(txdb_pth)) {
    #txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
    print("No Homo_sapiens.GRCh38.84.sqlite found. See https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/ for installation")
  }
  txdb <- EnsDb(txdb_pth)
  
  # Get transcript ranges from Ensembl transcript database
  transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts))
  
  # Flag if start overlaps with any transcripts
  vcf_info$all_transcript_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr))
  
  # Flag if end overlaps with any transcripts
  vcf_info$all_transcript_related_end <- sapply(1:dim(vcf_info)[1], function(x){
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
  
  
  
  # Initiate gene annotations
  vcf_info$ann_start <- ""
  vcf_info$ann_end  <- ""
  
  # Find exact overlaps of start and end position with transcripts
  #overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), transcripts_gr)  # testing
  overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr)
  overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$all_transcript_related_end),]$CHR, ranges=IRanges(start = vcf_info[!is.na(vcf_info$all_transcript_related_end),]$END, width = 1)), transcripts_gr)
  
  # Add start position overlaps to table
  if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
    vcf_info[queryHits(overlaps_start),]$ann_start <- sapply(seq(length(overlaps_start)), function(x){
      i <- queryHits(overlaps_start)[x]
      j <- subjectHits(overlaps_start)[x]
      vcf_info[i,]$ann_start <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  # Add end position overlaps to table --- FIXED
  if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
    vcf_info[!is.na(vcf_info$all_transcript_related_end),][queryHits(overlaps_end),]$ann_end <- sapply(seq(length(overlaps_end)), function(x){
      i <- queryHits(overlaps_end)[x]
      j <- subjectHits(overlaps_end)[x]
      vcf_info[!is.na(vcf_info$all_transcript_related_end),][i,]$ann_end <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  
  return(vcf_info)  
}

# Add annotations for all breakpoints
fusions2 <- annGenes(fusions2, transcripts_all$transcript_ID)

# Add gene names
fusions2$gene_start <- ""
fusions2$gene_end <- ""
transcripts_all$gene_name <- as.character(transcripts_all$gene_name)

fusions2[fusions2$ann_start != "",]$gene_start <- transcripts_all[match(fusions2[fusions2$ann_start != "",]$ann_start, transcripts_all$transcript_ID,),]$gene_name
fusions2[fusions2$ann_end != "",]$gene_end <- transcripts_all[match(fusions2[fusions2$ann_end != "",]$ann_end, transcripts_all$transcript_ID,),]$gene_name

# Table of fusions (precise breakpoints)
table(fusions2[fusions2$Application == "Manta",]$gene_start, fusions2[fusions2$Application == "Manta",]$gene_end)


#### TMPRSS2-ERG fusions #### 

# Manta calls
fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta")
# Number of samples with a manta call
length(unique(fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta") %>% pull(SAMPLE_WELL_ID)))  # 12

# Samples with a Canvas call
temp <- fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Canvas") %>% pull(SAMPLE_WELL_ID)
fusions2 %>% filter(SAMPLE_WELL_ID %in% temp) %>%  filter(gene_start == "ERG", gene_end == "TMPRSS2")

# Purity of samples with a call
temp <- fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta") %>% pull(SAMPLE_WELL_ID)
samples %>% filter(sampleId %in% temp) %>% select(sampleId, TUMOUR_PURITY)

# ERG breakpoints
ggplot((fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta")), aes(x = START, y = 0)) +
  geom_point(alpha = 0.4, size = 2) +
  tiltedX +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank())
  
# How many ERG breakpoints fall within the Weier et al ERG hotspot (chr21:38482768-38510697)
fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta", START > 38482768, START < 38510697) # 7/13

# TMPRSS breakpoints
ggplot((fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta")), aes(x = END, y = 0)) +
  geom_point(alpha = 0.4, size = 2) +
  tiltedX +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x = element_blank())

### Samples with colocalized Canvas call

# Canvas call overlapping the transcripts
temp <- fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Canvas") %>% pull(SAMPLE_WELL_ID)
fusions2 %>% filter(SAMPLE_WELL_ID %in% temp) %>%  filter(gene_start == "ERG", gene_end == "TMPRSS2")
# Other colocalized Canvas call
fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta", ColocalizedCanvas == TRUE, !SAMPLE_WELL_ID %in% temp)




#### TMPRSS2 translocations #### 

# All TMPRSS2
fusions2 %>% filter((gene_start == "TMPRSS2" | gene_end == "TMPRSS2"), Application == "Manta")

# Translocation partners
temp <- fusions2 %>% filter((gene_start == "TMPRSS2" | gene_end == "TMPRSS2"), SVTYPE == "BND") %>% pull(MATEID)
fusions2 %>% filter(ID %in% temp)

# Breakpoint location range (all)
range(c(fusions2 %>% filter(gene_start == "TMPRSS2", Application == "Manta") %>% pull(START)), (fusions2 %>% filter(gene_end == "TMPRSS2", Application == "Manta") %>% pull(END)))
# chr21:41464721-41507611
# Breakpoint location range (TMPRSS2-ERG fusions)
range(fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta") %>% pull(END))
# chr21:41485614-41507611
range(fusions2 %>% filter(gene_start == "ERG", gene_end == "TMPRSS2", Application == "Manta") %>% pull(START))
# chr21:38448535-38569131



#### ERG translocations #### 

# All ERG
fusions2 %>% filter((gene_start == "ERG" | gene_end == "ERG"), Application == "Manta")

# ERG translocation partners
temp2 <- fusions2 %>% filter((gene_start == "ERG" | gene_end == "ERG"), SVTYPE == "BND") %>% pull(MATEID)
fusions2 %>% filter(ID %in% temp2) %>% pull(gene_start)


#### Other rearrangements #### 

#### ETV1 rearrangements 

# All ETV1 
fusions2 %>% filter((gene_start == "ETV1" | gene_end == "ETV1"), Application == "Manta")

# ETV1 translocation partners
temp3 <- fusions2 %>% filter((gene_start == "ETV1" | gene_end == "ETV1"), SVTYPE == "BND") %>% pull(MATEID)
fusions2 %>% filter(ID %in% temp3)


#### ETV4 rearrangements 

# All ETV4 
fusions2 %>% filter((gene_start == "ETV4" | gene_end == "ETV4"), Application == "Manta")  # None


#### ETV5 rearrangements 

# All ETV5
fusions2 %>% filter((gene_start == "ETV5" | gene_end == "ETV5"), Application == "Manta")  # None






