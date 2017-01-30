# Martina Mijuskovic
# FFPE project
# Collect QC and clinical data for sequenced FFPE trios

library("data.table")
library("dplyr")
library("ggvis")
library("zoo")
library("ggplot2")

rm(list=ls())

############  Read and clean data ############  

# Read QC data for all cancer samples ((Jan 26 2017))
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


############  Subset and check trios ############  

# Subset samples that pass QC, pass tumor and germline contamination (remove also those with no data available)
table(QC_portal$PASS_TO_SEQ, exclude = NULL)
# NOTE that I'm not filtering germline for contamination==PASS
QC_portal_passPreSeqQC <- QC_portal %>% filter(PASS_TO_SEQ %in% c("Pass", "TBD"), TUMOUR_CONTAMINATION != "Fail", GERMLINE_CONTAMINATION != "Fail")
table((QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE")) %>% .$TUMOUR_CONTAMINATION), exclude = NULL)
table((QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE == "GL") %>% .$GERMLINE_CONTAMINATION), exclude = NULL)
# Remove tumor samples (only) with no contamination data available
rmIDs <- QC_portal_passPreSeqQC %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"), TUMOUR_CONTAMINATION == "") %>% .$SAMPLE_WELL_ID
QC_portal_passPreSeqQC <- QC_portal_passPreSeqQC %>% filter(!SAMPLE_WELL_ID %in% rmIDs)  # 1170 entries

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
dim(QC_portal_passPreSeqQC %>% filter(Trio == 1))  # 86
table((QC_portal_passPreSeqQC %>% filter(Trio == 1, SAMPLE_TYPE %in% c("FF", "FFPE")) %>% .$TUMOUR_CONTAMINATION), exclude=NULL)

# Subset data for clean trios
QC_portal_trios <- QC_portal_passPreSeqQC %>% filter(Trio == 1)
table(table(QC_portal_trios$PATIENT_ID))  # NOTE that there are 3 samples for 28 patients but 2 patients have 4 samples
table(QC_portal_trios$SAMPLE_TYPE)  # There are 2 samples with extra germline sample
table((QC_portal_trios %>% filter(SAMPLE_TYPE =="GL") %>% .$GERMLINE_CONTAMINATION), exclude = NULL)
# Two patients with 2 GL entries, one entry in each doesn't Pass GERMLINE_CONTAMINATION
QC_portal_trios[QC_portal_trios$PATIENT_ID %in% (QC_portal_trios %>% filter(SAMPLE_TYPE == "GL", GERMLINE_CONTAMINATION != "Pass") %>% .$PATIENT_ID),]
# Remove extra GL entries
rmGLs <- QC_portal_trios %>% filter(SAMPLE_TYPE == "GL", GERMLINE_CONTAMINATION != "Pass") %>% .$SAMPLE_WELL_ID
QC_portal_trios <- QC_portal_trios %>% filter(!SAMPLE_WELL_ID %in% rmGLs)

### IMPORTANT: Remove any samples from WEST MIDLANDS, "RRK" (discontinued - unreliable)
QC_portal_trios <- QC_portal_trios %>% filter(CENTER_CODE != "RRK")


# Load trios list from Alona (only now became available!)
Alonas_trios <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/trios.2016-01-26.txt", sep = "\t", colClasses = rep("character", 2))
names(Alonas_trios) <- c("PATIENT_ID", "STATUS")
Alonas_trios$ToUse <- 0
Alonas_trios[Alonas_trios$STATUS == "",]$ToUse <- 1
table(Alonas_trios$ToUse)  # 26 trios

# Check my trios against Alona's
trioIDs <- Alonas_trios %>% filter(ToUse == 1) %>% .$PATIENT_ID
sum(trioIDs %in% QC_portal_trios$PATIENT_ID)  # 26 (all there)



############   Add tumor type  ############  
tumor_types <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/cancer_disease_information_tum_2017-01-23_20-37-56.tsv", sep = "\t", header = T)
QC_portal_trios$TumorType <- tumor_types[match(QC_portal_trios$PATIENT_ID, tumor_types$participant_identifiers_id),]$disease_type_id
QC_portal_trios$TumorType <- as.character(QC_portal_trios$TumorType)

# Fix the unknown tumor type (from missing Oxford data)
missingTTid <- unique(QC_portal_trios %>% filter(TumorType == "Unknown") %>% .$PATIENT_ID)
missingOx <- read.csv("/Users/MartinaMijuskovic/Documents/FFPE/Oxford_samples_with_missing_data MC.csv", header = T)
QC_portal_trios[QC_portal_trios$PATIENT_ID == missingTTid,]$TumorType <- as.character(missingOx[missingOx$PATIENT_ID == missingTTid,]$X)



############  Add BAM paths ############

# Read the current upload report, restrict to cancer, V4 and qc_passed
today <- Sys.Date()
system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
upload <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")
colnames(upload) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
upload <- upload %>% filter(`Delivery Version` == "V4", Status == "qc_passed", Type %in% c("cancer germline", "cancer tumour"))
upload$Path <- as.character(upload$Path)

# Add BAM paths to the trios table
QC_portal_trios$BamPath <- upload[match(QC_portal_trios$SAMPLE_WELL_ID, upload$Platekey),]$Path


# Write the trios QC portal data
write.csv(QC_portal_trios, file = "QC_portal_trios.csv", quote = F, row.names = F)





############  Explore trios data ############

# Distribution by GMC
as.data.frame(table(as.character(QC_portal_trios$CENTER_CODE), exclude = NULL))

# Distribution by library and sample type
as.data.frame(table(as.character(QC_portal_trios$LIBRARY_TYPE), QC_portal_trios$SAMPLE_TYPE), exclude = NULL)

# Plots with QC metrics, grouped by sample type and colored by GMC
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=AT_DROP, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "A/T Dropout")
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=GC_DROP, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "G/C Dropout")
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=COVERAGE_HOMOGENEITY, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "Unevenness of Coverage")
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=CHIMERIC_PER, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "Chimeric Reads")
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=DEAMINATION_MISMATCHES_PER, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "Deamination Missmatches")
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=AV_FRAGMENT_SIZE_BP, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "Average Fragment Size")
ggplot(QC_portal_trios, aes(x=SAMPLE_TYPE, y=MAPPING_RATE_PER, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Sample type", y = "Mapping Rate")

# Explore FFPE trio data by tumor type
as.data.frame(table(as.character(QC_portal_trios[QC_portal_trios$SAMPLE_TYPE=="GL",]$CENTER_CODE), as.character(QC_portal_trios[QC_portal_trios$SAMPLE_TYPE=="GL",]$TumorType)), exclude = NULL)
ggplot((QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE")), aes(x=TumorType, y=AT_DROP, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Tumor type", y = "A/T Dropout")
ggplot((QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE")), aes(x=TumorType, y=COVERAGE_HOMOGENEITY, colour = CENTER_CODE)) + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank()) + labs(x = "Tumor type", y = "Unevenness of Coverage")








