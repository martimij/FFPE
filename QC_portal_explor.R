# Look at Alona's QC metrics done on all IIP samples, incl. 26 FFPE trios (samples w germline, fresh frozen and FFPE)

library("data.table")
library("dplyr")
library("ggvis")

setwd("/Users/MartinaMijuskovic/Documents/FFPE")

QC_data <- read.csv("ready_to_upload.batches1-6.v4.csv")

############  Explore QC data ############  

# Number of samples (entries) 
dim(QC_data)[1]   # 1397

# Distribution by center
table(QC_data$CENTER_CODE, exclude = NULL)

# Distribution by sample type
table(QC_data$SAMPLE_TYPE, exclude = NULL)

# Pre-seq QC pass/fail
table(QC_data$PASS_TO_SEQ, exclude = NULL)

# Number of patients involved
length(unique(QC_data$PATIENT_ID))  # 602


############  Clean QC and add flags ############  

QC_data[QC_data$PASS_TO_SEQ == "pass",]$PASS_TO_SEQ <- "Pass"
QC_data$PASS_TO_SEQ <- as.character(QC_data$PASS_TO_SEQ)

# Total passing QC
QC_data %>% filter(PASS_TO_SEQ == "Pass")  # 1091

# FFPE passing QC
QC_data %>% filter(PASS_TO_SEQ == "Pass", SAMPLE_TYPE == "FFPE")  # 192

QC_data$PATIENT_ID <- as.character(QC_data$PATIENT_ID)


############  Add GMC data to the QC table ############  

# Subset of QC table containing only sequenced FFPE samples that pass tumor contamination (192)
QC_FFPE <- QC_data %>% filter(PASS_TO_SEQ == "Pass", SAMPLE_TYPE == "FFPE", TUMOUR_CONTAMINATION == "Pass")
dim(QC_FFPE)  # 105
length(unique(QC_FFPE$PATIENT_ID))  # 99  (WARNING: some patients have multiple FFPE samples sequenced! or duplicates?)

# Clean duplicates
sum(duplicated(QC_FFPE))  # 6 completely duplicate rows (and 6 possible double samples per patient?)
QC_FFPE[QC_FFPE$PATIENT_ID %in% QC_FFPE[duplicated(QC_FFPE),]$PATIENT_ID,]

# Remove row duplicates
QC_FFPE <- QC_FFPE[!duplicated(QC_FFPE),]  # 6 removed

# Samples with duplicated patient IDs (NONE)
#duplicate_IDs <- QC_FFPE[duplicated(QC_FFPE$PATIENT_ID),]$PATIENT_ID
#QC_FFPE %>% filter(PATIENT_ID %in% duplicate_IDs)  # These are probably purposefully duplicate samples (different wells! but same metrics)

# Check for duplicate keys
sum(duplicated(QC_FFPE$PATIENT_ID))  # none
sum(duplicated(QC_FFPE$SAMPLE_LAB_ID))  # none

# # Add a duplicate flag to the second copy
# QC_FFPE$SEQ_DUPLICATE <- 0
# QC_FFPE$SEQ_DUPLICATE[which(duplicated(QC_FFPE$PATIENT_ID))] <- 1


# Read a subset of GMC data from Alona (full file)
gmc <- read.table("collectedsample_2017-01-20_10-42-09.tsv", sep = "\t", fill = T, header = T)

# Clean PATIENT_ID, CENTER_CODE
names(gmc)[3] <- "PATIENT_ID"
names(gmc)[1] <- "CENTER_CODE"
names(gmc)[4] <- "SAMPLE_LAB_ID"

# Subset the gmc table for patient IDs in the QC_FFPE table
gmc_FFPE <- gmc %>% filter(PATIENT_ID %in% QC_FFPE$PATIENT_ID)
dim(gmc_FFPE)  # 329 samples (not all are FFPE)

# Check for duplicated GMC sample IDs (OK)
sum(duplicated(gmc_FFPE$SAMPLE_LAB_ID))   # 0

# Number of matching SAMPLE_IDs
sum(gmc_FFPE$SAMPLE_LAB_ID %in% QC_FFPE$SAMPLE_LAB_ID)   # 51

# Find FFPE sample information in the GMC data: add the FFPE flag to GMC table 
gmc_FFPE$FFPE <- 0
gmc_FFPE[gmc_FFPE$SAMPLE_LAB_ID %in% QC_FFPE$SAMPLE_LAB_ID,]$FFPE <- 1  # FFPE samples that have matching sample IDs between tables
# gmc_FFPE[gmc_FFPE$typeOfFixative != "",]$FFPE <- 1  # assumed FFPE samples with "typeOfFixative" info
# gmc_FFPE[gmc_FFPE$fixationStartDateTime != "",]$FFPE <- 1  # assumed FFPE samples with "fixationStartDateTime" info
#gmc_FFPE[grep("ffpe", gmc_FFPE$sampleID, ignore.case = T),]$FFPE <- 1   # assumed FFPE samples with "ffpe" in the sample ID  (4, already flagged)

# Subset GMC table for samples that are assumed FFPE
gmc_FFPE <- gmc_FFPE %>% filter(FFPE == 1)   # 51 entries total

# Check merging key (Sample ID) for duplicates in both tables again
sum(duplicated(QC_FFPE$SAMPLE_LAB_ID))  # 0
sum(duplicated(gmc_FFPE$SAMPLE_LAB_ID))  # 0
sum(duplicated(QC_FFPE$PATIENT_ID))  # 0
sum(duplicated(gmc_FFPE$PATIENT_ID))  # 0

# Merge QC and GMC tables using SAMPLE_LAB_ID
gmc_FFPE$SAMPLE_LAB_ID <- as.character(gmc_FFPE$SAMPLE_LAB_ID)
QC_FFPE$SAMPLE_LAB_ID <- as.character(QC_FFPE$SAMPLE_LAB_ID)
FFPE_data <- inner_join(QC_FFPE, gmc_FFPE, by = "SAMPLE_LAB_ID")  # 51 observations

# Check that Patient IDs match from both tables
sum(FFPE_data$PATIENT_ID.x != FFPE_data$PATIENT_ID.y)  # 0   (all match!!!)
FFPE_data$CENTER_CODE.x <- as.character(FFPE_data$CENTER_CODE.x)
FFPE_data$CENTER_CODE.y <- as.character(FFPE_data$CENTER_CODE.y)
sum(FFPE_data$CENTER_CODE.x != FFPE_data$CENTER_CODE.y)  # 0  (all match!!!)


############  Summary of 51 matched FFPE samples ############  

# By GMC
as.data.frame(table(FFPE_data$CENTER_CODE.x, exclude = NULL))

# Samples with collection date, fixation start, formalin type and incubation time available
dim(FFPE_data %>% filter(clinicSampleDateTime != ""))  # 33
dim(FFPE_data %>% filter(fixationStartDateTime != ""))  # 28
dim(FFPE_data %>% filter(typeOfFixative != ""))  # 31
dim(FFPE_data %>% filter(processingSchedule != ""))  # 30

############  Summary, plots of OXFORD samples ############  

# Samples with collection date, fixation start, formalin type and incubation time available
dim(FFPE_data %>% filter(CENTER_CODE.x == "RTH", clinicSampleDateTime != ""))  # 30 
dim(FFPE_data %>% filter(CENTER_CODE.x == "RTH", fixationStartDateTime != ""))  # 27
dim(FFPE_data %>% filter(CENTER_CODE.x == "RTH", typeOfFixative != ""))  # 28
dim(FFPE_data %>% filter(CENTER_CODE.x == "RTH", processingSchedule != ""))  # 27

#grep("date", names(gmc), value = T, ignore.case = T)  # "clinicSampleDateTime"  "snapFreexingDateTime"  "fixationStartDateTime" "fixationEndDateTime"
#grep("date", names(gmc), ignore.case = T)  # 6 23 25 26

# Create date of collection variable
FFPE_data$clinicSampleDateTime  <- as.character(FFPE_data$clinicSampleDateTime)
# Remove hours
FFPE_data$SampleCollectionDate <- sapply(1:dim(FFPE_data)[1], function(x){
                strsplit(FFPE_data$clinicSampleDateTime[x], split = " ")[[1]][1]
              })
#sum(is.na(FFPE_data$SampleCollectionDate))

# Convert collection date to Date variable
FFPE_data$SampleCollectionDate <- as.Date(FFPE_data$SampleCollectionDate)

# # Create month of collection variable
# FFPE_data$CollectionMonth <- sapply(1:dim(FFPE_data)[1], function(x){
#   paste(strsplit(FFPE_data$SampleCollectionDate[x], split = "-")[[1]][1], strsplit(FFPE_data$SampleCollectionDate[x], split = "-")[[1]][2], sep = "")
# })
# FFPE_data[FFPE_data$CollectionMonth == "NANA",]$CollectionMonth <- NA
# FFPE_data$CollectionMonth <- as.integer(FFPE_data$CollectionMonth)

# Create "BufferedFormaline" flag
FFPE_data$BufferedFormaline <- NA
FFPE_data[grep("neutral", FFPE_data$typeOfFixative, ignore.case = T),]$BufferedFormaline <- 1
FFPE_data[grep("nonbuffered", FFPE_data$typeOfFixative, ignore.case = T),]$BufferedFormaline <- 0
FFPE_data[grep("saline", FFPE_data$typeOfFixative, ignore.case = T),]$BufferedFormaline <- 0
names(FFPE_data)[74] <- "BufferedFormalin"

# Create "OvernightIncubation" flag
FFPE_data$OvernightIncubation <- NA
FFPE_data[grep("Extended", FFPE_data$processingSchedule, ignore.case = T),]$OvernightIncubation <- 0
FFPE_data[grep("Overnight", FFPE_data$processingSchedule, ignore.case = T),]$OvernightIncubation <- 1

# Oxford subset
FFPE_Ox <- FFPE_data %>% filter(CENTER_CODE.x == "RTH")

# Rename flags
FFPE_Ox[is.na(FFPE_Ox$BufferedFormalin),]$BufferedFormalin <- "Unavailable"
FFPE_Ox[FFPE_Ox$BufferedFormalin == 0,]$BufferedFormalin <- "Non-Buffered"
FFPE_Ox[FFPE_Ox$BufferedFormalin == 1,]$BufferedFormalin <- "Buffered"

FFPE_Ox[is.na(FFPE_Ox$OvernightIncubation),]$OvernightIncubation <- "Unavailable"
FFPE_Ox[FFPE_Ox$OvernightIncubation == 0,]$OvernightIncubation <- "Extended"
FFPE_Ox[FFPE_Ox$OvernightIncubation == 1,]$OvernightIncubation <- "Overnight"

# Plot Oxford samples by time grouped by buffered formalin
FFPE_Ox %>% ggvis(~SampleCollectionDate, ~AT_DROP, fill = ~factor(BufferedFormalin)) %>% layer_points() %>% add_legend("fill", title = "Buffered Formalin") %>% scale_datetime("x", nice = "year")
FFPE_Ox %>% ggvis(~SampleCollectionDate, ~COVERAGE_HOMOGENEITY, fill = ~factor(BufferedFormalin)) %>% layer_points() %>% add_legend("fill", title = "Buffered Formalin") %>% scale_datetime("x", nice = "year")

# Plot Oxford samples by time grouped by overnight incubation
FFPE_Ox %>% ggvis(~SampleCollectionDate, ~AT_DROP, fill = ~factor(OvernightIncubation)) %>% layer_points() %>% add_legend("fill", title = "Overnight Incubation") %>% scale_datetime("x", nice = "year")
FFPE_Ox %>% ggvis(~SampleCollectionDate, ~COVERAGE_HOMOGENEITY, fill = ~factor(OvernightIncubation)) %>% layer_points() %>% add_legend("fill", title = "Overnight Incubation") %>% scale_datetime("x", nice = "year")


