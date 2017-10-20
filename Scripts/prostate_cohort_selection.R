# Martina Mijuskovic
# SV pipeline development
# Analysis of TMPRSS2-ERG fusions in prostate cancer cohort
# Selection of the prostate cancer cohort
# Oct 2017

library(dplyr)
library(ggplot2)
library(data.table)
library(jsonlite)

today <- Sys.Date()

#### Functions #### 

getCatalogParticipants <- function(PATIENT_ID, sessionID, studyID="1000000041"){
  require(jsonlite)
  require(dplyr)
  # Create Catalog 1.0 search command line
  #curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/individuals/search?sid=eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJiZXJ0aGEiLCJhdWQiOiJPcGVuQ0dBIHVzZXJzIiwiaWF0IjoxNTA3NjQwNjU5LCJleHAiOjE1MDc2NDI0NTl9.kGXQJE-xRLzL4DfB1dwHwyiUF9gkfUJBmjtkm6SzTYw&study=1000000041&skipCount=false&limit=1"
  #curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/individuals/search?sid=eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJiZXJ0aGEiLCJhdWQiOiJPcGVuQ0dBIHVzZXJzIiwiaWF0IjoxNTA3NjQwNjU5LCJleHAiOjE1MDc2NDI0NTl9.kGXQJE-xRLzL4DfB1dwHwyiUF9gkfUJBmjtkm6SzTYw&study=1000000041&name=217000038&skipCount=false&limit=1"
  
  command <- paste0('curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/individuals/search?sid=', sessionID, '&study=', studyID, '&name=', PATIENT_ID, '&skipCount=false&limit=1"')
  # Read in the json file containing the QC metrics for the specific sample
  bam_json <- fromJSON(system(command, intern = T), flatten = T)
  
  if (bam_json$response$numResults == 0) {
    result <- data.frame(
      PATIENT_ID = PATIENT_ID,
      Platekey_normal = NA,
      Platekey = NA)
    return(result)
  }
  result <- data.frame(PATIENT_ID = bam_json$response$result[[1]]$name, 
                       Platekey_normal = bam_json$response$result[[1]]$annotationSets[[1]]$annotations[[2]][bam_json$response$result[[1]]$annotationSets[[1]]$annotations[[2]]$name == "germlineSampleId",]$value,
                       Platekey = bam_json$response$result[[1]]$annotationSets[[1]]$annotations[[2]][bam_json$response$result[[1]]$annotationSets[[1]]$annotations[[2]]$name == "tumorSampleId",]$value  # Needs fixing in case FF and FFPE exist for this patient
                        )
  return(result)
}


#### Explore the prostate cohort ####

### Read a recent upload report
upload <- fread("./Data/upload_report.2017-09-25.txt", skip = 14, header = T)
### Read today's upload report
today <- Sys.Date()
system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
upload <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")
colnames(upload) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
upload <- upload %>% filter(Type != "rare disease", `Delivery Version` == "V4")

upload_full <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")  # unfiltered version
colnames(upload_full) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))


### Table of cancer samples with tumour types (from Alona)
cancer_samples <- read.csv("./Data/cancer_samples.csv")

# Select prostate cohort
table(cancer_samples$tumorType, exclude = NULL)
cancer_samples_prostate <- cancer_samples %>% filter(tumorType == "Prostate")  # 105
table(cancer_samples_prostate$preservationMethod, exclude = NULL)  # 55 FF, 50 FFPE

# Check for duplicates (FFPE trios?)
sum(duplicated(cancer_samples_prostate$gelId))  # 7 dups
cancer_samples_prostate[cancer_samples_prostate$gelId %in% cancer_samples_prostate[duplicated(cancer_samples_prostate$gelId),]$gelId,] # Some look like FFPE trios, others just duplicated samples

### Already interpreted cases (batches 1-4, fast track samples): collect patients from the Catalog study 
old_batches <- bind_rows(lapply(cancer_samples_prostate$gelId, getCatalogParticipants, sessionID = "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJiZXJ0aGEiLCJhdWQiOiJPcGVuQ0dBIHVzZXJzIiwiaWF0IjoxNTA4MTc0OTE2LCJleHAiOjE1MDgxNzY3MTZ9.LqicAMwwjiiJEuaWLJAdMmrC2K6LuhZT_5bBSZzWLg0"))
summary(old_batches)  # all 105 registered in catalog

### Batch 5 (uninterpreted cases): collect patients from the DDF payload
batch5 <- read.csv("./Data/Batch 5 Reviewed, Filtered, SampleIDs, Payloads.csv")

dim(batch5)  # 326
summary(batch5)

# Find prostate samples
batch5$Prostate <- sapply(1:dim(batch5)[1], function(x){
  grepl("prostate", batch5$Payload[x], ignore.case = T)
})

# Summary
dim(batch5 %>% filter(Prostate == 1))  # 17 samples

# Extract payload jsons from prostate cases
sapply(1:dim(batch5)[1], function(x){
  if(grepl("prostate", batch5$Payload[x], ignore.case = T)){
    json <- batch5$Payload[x]
    write.csv(json, file = paste0("./Data/batch5/", batch5$individualId[x], ".json"), quote = F, row.names = F, col.names = F)
  }
})

# Get all clinical data for the prostate batch5
batch5clin <- lapply(list.files("./Data/batch5/"), function(x){
  json <- fromJSON(readLines(paste0("./Data/batch5/", x))[2], flatten = T)
  data.frame(PATIENT_ID = json$individualId, 
             CENTER_CODE = json$center, 
             TumourType = json$primaryDiagnosisDisease, 
             SAMPLE_TYPE = json$tumourSamples$preparationMethod, 
             GERMLINE_SOURCE = json$germlineSamples$source, 
             Platekey = json$matchedSamples$tumourSampleId,
             Platekey_normal = json$matchedSamples$germlineSampleId
             )
})
batch5clin <- bind_rows(batch5clin)


# Any batch5 (=potential batch5) samples in the cancer samples table?
sum(batch5clin$Platekey %in% cancer_samples_prostate$sampleId)  # 9


#### Prostate cohort selection #### 

# Find prostate cases in the upload report
prostate_ids <- unique(cancer_samples_prostate$sampleId)
length(prostate_ids) # 105
dim(upload %>% filter(Platekey %in% prostate_ids))  # 69 only (36 missing)
dim(upload_full %>% filter(Platekey %in% prostate_ids))  # 98 samples (1 qc-failed, 1 unknown version)
upload_prostate <- upload_full %>% filter(Platekey %in% prostate_ids)
table(upload_prostate$`Delivery Version`, upload_prostate$Status, exclude = NULL)
dim(upload_prostate)  # 98 total
sum(duplicated(upload_prostate$Platekey))  # 29
# Are V2 samples all lifted? - yes, all except 1 qc-failed sample (total = 28 lifted)
dups <- as.character(upload_prostate$Platekey[duplicated(upload_prostate$Platekey)])
table(upload_prostate[upload_prostate$Platekey %in% dups,]$`Delivery Version`, upload_prostate[upload_prostate$Platekey %in% dups,]$Status, exclude = NULL)
# Check one "unknown" version
upload_prostate %>% filter(`Delivery Version` == "unknown")  # LP3000212-DNA_B03
upload_prostate %>% filter(Platekey == "LP3000212-DNA_B03")  # This sample is actually listed as germline in V4
# Check QC failed sample
upload_prostate %>% filter(Status == "qc_failed") # LP3000079-DNA_G04 (V2)
upload_full %>% filter(Platekey == "LP3000079-DNA_G04") # V4 version has passed QC (!)
# Remove duplicated V2 samples that are lifted to V4
rm <- as.character(upload_prostate %>% filter(Platekey %in% dups, `Delivery Version` != "V4") %>% pull(DeliveryID))
upload_prostate <- upload_prostate %>% filter(!DeliveryID %in% rm) # 69 left
table(upload_prostate[upload_prostate$Platekey %in% dups,]$`Delivery Version`, upload_prostate[upload_prostate$Platekey %in% dups,]$Status, exclude = NULL)
table(upload_prostate$`Delivery Version`, upload_prostate$Status)

# Excluding sample LP3000212-DNA_B03 (too hard to figure out what's going on with it)
cancer_samples_prostate <- cancer_samples_prostate %>% filter(sampleId != "LP3000212-DNA_B03") # 104 left
upload_prostate <- upload_prostate %>% filter(Platekey != "LP3000212-DNA_B03") # 68 left 

# Remove samples not in the upload report
dim(cancer_samples_prostate)  # 104
cancer_samples_prostate <- cancer_samples_prostate %>% filter(sampleId %in% upload_prostate$Platekey)
dim(cancer_samples_prostate)  # 68

# Merge upload report with clinical inf
cancer_samples_prostate <- full_join(cancer_samples_prostate, upload_prostate, by = c("sampleId" = "Platekey"))

# Mark FFPE trios
# Flag trios
trios <- sapply(unique(cancer_samples_prostate$gelId), function(x){
  if (sum(c("FFPE", "FF") %in% cancer_samples_prostate[cancer_samples_prostate$gelId == x,]$preservationMethod) == 2) {
    return(1)
  }
  else {
    return(0)
  }
})
sum(trios)  # 2 trios
names(trios) <- unique(cancer_samples_prostate$gelId)
cancer_samples_prostate$Trio <- 0
cancer_samples_prostate[cancer_samples_prostate$gelId %in% names(trios[trios == 1]),]$Trio <- 1
# Sanity check
cancer_samples_prostate %>% filter(Trio == 1)
cancer_samples_prostate$preservationMethod <- as.character(cancer_samples_prostate$preservationMethod)
# Summary
table(cancer_samples_prostate$Trio, cancer_samples_prostate$preservationMethod)

# Remove FFPE sample from a trio
cancer_samples_prostate <- cancer_samples_prostate %>% filter(!(Trio == 1 & preservationMethod == "FFPE"))
table(cancer_samples_prostate$Trio, cancer_samples_prostate$preservationMethod)

### Check tumour contamination for selected 68 samples
# Add germline Platekey
cancer_samples_prostate$Path <- as.character(cancer_samples_prostate$Path)
cancer_samples_prostate$Platekey_normal <- sapply(cancer_samples_prostate$Path, function(x){strsplit(strsplit(x, split = "/")[[1]][6], split = "Normal")[[1]][2]})
# Add upload report info for germline
cancer_samples_prostate$Platekey_normal[!cancer_samples_prostate$Platekey_normal %in% upload$Platekey] # no germline missing from upload report
dim(upload %>% filter(Platekey %in% cancer_samples_prostate$Platekey_normal))  # 66 (all there)
upload_gl <- upload %>% filter(Platekey %in% cancer_samples_prostate$Platekey_normal)
names(upload_gl) <- paste0(names(upload_gl), "_GL")
sum(duplicated(upload_gl$Platekey_GL))
cancer_samples_prostate <- full_join(cancer_samples_prostate, upload_gl, by = c("Platekey_normal" = "Platekey_GL"))

# Read latest QC table (from Alona)
QC <- read.csv("./Data/QC/ready_to_upload.09-24-17.csv")

dim(QC %>% filter(SAMPLE_WELL_ID %in% c(cancer_samples_prostate$sampleId, cancer_samples_prostate$Platekey_normal))) # 132 (ok)

# Add QC info to the manifest
cancer_samples_prostate <- left_join(cancer_samples_prostate, QC, by = c("sampleId" = "SAMPLE_WELL_ID"))

summary(cancer_samples_prostate)

# Exclude contaminated tumours 
cancer_samples_prostate$Excluded <- 0
excl <- as.character(cancer_samples_prostate %>% filter(TUMOUR_CONTAMINATION == "Fail") %>% pull(sampleId))
cancer_samples_prostate[cancer_samples_prostate$sampleId %in% excl,]$Excluded <- 1

table(cancer_samples_prostate$Excluded, cancer_samples_prostate$preservationMethod)

### Check germline contamination for selected samples

table((QC %>% filter(SAMPLE_WELL_ID %in% cancer_samples_prostate$Platekey_normal) %>% pull(GERMLINE_CONTAMINATION))) # All pass


# Write the clean cohort
write.csv(cancer_samples_prostate, file = "./Data/cancer_samples_prostate_cleaned.csv", quote = F, row.names = F)
