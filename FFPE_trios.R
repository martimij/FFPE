# Martina Mijuskovic
# FFPE project
# Collect QC and clinical data for sequenced FFPE trios

library("data.table")
library("dplyr")
library("ggvis")
library("zoo")
library("ggplot2")

rm(list=ls())

############  Read and clean data (Jan 26 2017) ############  

# Read QC data for all cancer samples
QC_portal <- read.csv("/Users/MartinaMijuskovic/Documents/FFPE/ready_to_upload.01-26-17.csv")

# Summarize and check for dupliactes
dim(QC_portal)  # 1742 entries total
sum(duplicated(QC_portal))   # 0 duplicate rows
sum(duplicated(QC_portal$SAMPLE_WELL_ID))  # 0 duplicates
sum(duplicated(QC_portal$SAMPLE_LAB_ID))  # 81 duplicates

# Investigate SAMPLE_LAB_ID duplicates
QC_portal %>% filter(SAMPLE_LAB_ID %in% QC_portal$SAMPLE_LAB_ID[duplicated(QC_portal$SAMPLE_LAB_ID)])   # Duplicate entries for samples that didn't pass QC

# Any duplicate SAMPLE_LAB_IDs in samples that pass pre-sequencing QC? (NONE)
sum(duplicated(QC_portal %>% filter(SAMPLE_LAB_ID %in% QC_portal$SAMPLE_LAB_ID[duplicated(QC_portal$SAMPLE_LAB_ID)]) %>% filter(PASS_TO_SEQ == "PASS") %>% .$SAMPLE_LAB_ID))  # none


############  Subset trios ############  

# Subset samples that pass QC, pass tumor and germline contamination (remove also those with no data available)
table(QC_portal$PASS_TO_SEQ, exclude = NULL)
QC_portal_passPreSeqQC <- QC_portal %>% filter(PASS_TO_SEQ == "Pass", TUMOUR_CONTAMINATION != "Fail", GERMLINE_CONTAMINATION != "Fail")   # 489 samples (of 1742)
table((QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE")) %>% .$TUMOUR_CONTAMINATION), exclude = NULL)
table((QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE == "GL") %>% .$GERMLINE_CONTAMINATION), exclude = NULL)
# Remove tumor samples (only) with no contamination data available
rmIDs <- QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"), TUMOUR_CONTAMINATION == "") %>% .$SAMPLE_LAB_ID
QC_portal_passPreSeqQC <- QC_portal_passPreSeqQC %>% filter(!SAMPLE_LAB_ID %in% rmIDs)  # 1152 entries

# Flag trios
trios <- sapply(unique(QC_portal_passPreSeqQC$PATIENT_ID), function(x){
  if (sum(c("GL", "FFPE", "FF") %in% QC_portal_passPreSeqQC[QC_portal_passPreSeqQC$PATIENT_ID == x,]$SAMPLE_TYPE) == 3) {
    return(1)
  }
  else {
    return(0)
  }
})
names(trios) <- unique(QC_portal_passPreSeqQC$PATIENT_ID)

QC_portal_passPreSeqQC$Trio <- 0
QC_portal_passPreSeqQC[QC_portal_passPreSeqQC$PATIENT_ID %in% names(trios[trios == 1]),]$Trio <- 1

# Total number of trios that pass QC to sequencing
dim(QC_portal_passPreSeqQC %>% filter(Trio == 1))  # 78
table((QC_portal_passPreSeqQC %>% filter(Trio == 1, SAMPLE_TYPE %in% c("FF", "FFPE")) %>% .$TUMOUR_CONTAMINATION), exclude=NULL)

# Subset data for clean trios
QC_portal_trios <- QC_portal_passPreSeqQC %>% filter(Trio == 1)

# Load trios list from Alona (only now became available!)
Alonas_trios <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/trios.2016-01-26.txt", sep = "\t")

# Check my trios against Alona's
sum(!Alonas_trios[Alonas_trios$V2 == "",]$V1 %in% QC_portal_trios$PATIENT_ID)  # 2 missing  ---- check - do not filter PASS_TO_SEQ "TBDs"

############  Check status of trios in the upload report ############  

