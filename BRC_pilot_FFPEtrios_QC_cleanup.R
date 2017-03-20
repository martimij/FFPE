### Martina Mijuskovic
# FFPE project
# Compile and clean QC and clinical data for additional BRC pilot FFPE trios

library(dplyr)
library(ggvis)
library(ggplot2)
library(VariantAnnotation)
library(data.table)

rm(list=ls())

### For ggplot
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)


############  Clean BRC pilot data ############  

# Read list of pilot samples (incl. BRC FFPE trios)
pilot <- read.csv("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/cancer_pilot_samples_David.csv")

# Read the two lists with tumor types
pilot_extra1 <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/pilot_Matt_DB_data.txt", header = T)
pilot_extra2 <- read.table("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/Matt_DB_data.txt", header = T, sep = "\t")

# Explore
dim(pilot)  # 1014
table(pilot$study, exclude = NULL)
table(pilot$type, exclude = NULL)
table(pilot$type, exclude = NULL)
table(pilot$study, pilot$type)

# Reduce the list to BRC samples
pilot <- pilot %>% filter(study == "brc")
dim(pilot)  # 164

# Explore further
table(pilot$sex_problem, exclude = NULL)
table(pilot$modifier, pilot$type, exclude = NULL)

# Reduce
pilot <- pilot %>% dplyr::select(manifest_id, study, GelId, type, platekey, base_dir, modifier, assigned_diagnosis)

# Change names to fit 26 previous trios analysis
names(pilot)[3:5] <- c("PATIENT_ID", "SAMPLE_TYPE", "SAMPLE_WELL_ID")

# Flag trios
pilot$SAMPLE_TYPE <- as.character(pilot$SAMPLE_TYPE)
trios <- sapply(unique(pilot$PATIENT_ID), function(x){
  if (sum(c("GL", "FFPE", "FF") %in% pilot[pilot$PATIENT_ID == x,]$SAMPLE_TYPE) == 3) {
    return(1)
  }
  else {
    return(0)
  }
})
names(trios) <- unique(pilot$PATIENT_ID)

pilot$Trio <- 0
pilot[as.character(pilot$PATIENT_ID) %in% names(trios[trios == 1]),]$Trio <- 1

# Number of trios
sum(trios)  # 41
sum(pilot$Trio)/3  # 41.33333  (one FFPE sample extra!)
table(pilot$SAMPLE_TYPE, pilot$Trio, exclude = NULL)

# Find the origin of the extra sample
table(table(pilot$PATIENT_ID))  # 20 patients with 2 samples, 40 with trios, 1 with 4 samples
table(pilot$PATIENT_ID)  # 200000960 has 4 samples
sum(duplicated(pilot$SAMPLE_WELL_ID))  # 0
pilot %>% filter(PATIENT_ID == "200000960")  # this patient has 2 FFPE samples listed 
# (one may have failed but I have no info; older one?; LP2000830-DNA_E01 FFPE sample has only 20% purity)

# NOTE that one sample is denoted "EXPERIMENTAL" - best to remove; also missing diagnosis (PATIENT_ID==200000316)
pilot %>% filter(modifier == "EXPERIMENTAL")

# Explore extra data
table(pilot_extra1$CENTER, pilot_extra1$TUMOR, exclude = NULL)
table(pilot_extra1$CENTER, pilot_extra1$TYPE, exclude = NULL)
table(pilot_extra2$PILOT, pilot_extra2$TUMOR, exclude = NULL)
table(pilot_extra2$PILOT, pilot_extra2$type, exclude = NULL)

# Check if all trio samples have tumor type info
names(pilot_extra1)[1] <- "SAMPLE_WELL_ID"
names(pilot_extra2)[1] <- "SAMPLE_WELL_ID"
names(pilot_extra1)[5] <- "SAMPLE_TYPE"
names(pilot_extra2)[5] <- "SAMPLE_TYPE"
trio_ff_ids <- as.character(pilot %>% filter(Trio == 1, SAMPLE_TYPE == "FF") %>% .$SAMPLE_WELL_ID)
trio_ffpe_ids <- as.character(pilot %>% filter(Trio == 1, SAMPLE_TYPE == "FFPE") %>% .$SAMPLE_WELL_ID)
sum(!(trio_ff_ids %in% pilot_extra1[(!is.na(pilot_extra1$TUMOR)),]$SAMPLE_WELL_ID))  # 7 missing tumor info
sum(!(trio_ff_ids %in% pilot_extra2[(!is.na(pilot_extra2$TUMOR)),]$SAMPLE_WELL_ID))  # 5 missing tumor info
sum(!(trio_ffpe_ids %in% pilot_extra1[(!is.na(pilot_extra1$TUMOR)),]$SAMPLE_WELL_ID))  # 8 missing tumor info
sum(!(trio_ffpe_ids %in% pilot_extra2[(!is.na(pilot_extra2$TUMOR)),]$SAMPLE_WELL_ID))  # 6 missing tumor info

# Reduce extra tables to samples in the BRC pilot
pilot_extra1 <- pilot_extra1 %>% filter(SAMPLE_WELL_ID %in% pilot$SAMPLE_WELL_ID)  # this seems to be an incomplete subset of pilot_extra2
pilot_extra2 <- pilot_extra2 %>% filter(SAMPLE_WELL_ID %in% pilot$SAMPLE_WELL_ID)

# Check
sum(duplicated(pilot_extra2$SAMPLE_WELL_ID))  #0

# Check extra data for FF and FFPE trios
pilot_extra1 %>% filter(SAMPLE_WELL_ID %in% c(trio_ff_ids, trio_ffpe_ids))  # some marked "EXPT" (experimental?) under "PILOT"
pilot_extra2 %>% filter(SAMPLE_WELL_ID %in% c(trio_ff_ids, trio_ffpe_ids))  # some marked "EXPT" under "PILOT"


# Add tumor info from the 2nd file, clean up
pilot <- left_join(pilot, (pilot_extra2 %>% dplyr::select(SAMPLE_WELL_ID, TUMOR, SAMPLE_TYPE, CENTER, PILOT)))
table(pilot$Trio, pilot$TUMOR)
pilot[is.na(pilot$TUMOR),]$TUMOR <- "UNKNOWN"
pilot[pilot$TUMOR == "UNKOWN",]$TUMOR <- "UNKNOWN"
pilot$TUMOR <- as.character(pilot$TUMOR)

# PATIENT_IDs with missing tumor info
unique(pilot %>% filter(Trio == 1, TUMOR == "UNKNOWN") %>% .$PATIENT_ID)  # 200000926 200000929 200000953 200000954 200000960 200001299 200001310
pilot %>% filter(PATIENT_ID %in% (unique(pilot %>% filter(Trio == 1, TUMOR == "UNKNOWN") %>% .$PATIENT_ID)))

# Manually add tumor type using the dianogostic codes (http://www.icd10data.com/)
pilot[pilot$assigned_diagnosis %in% c("C50.9", "C50.3", "C50.8"),]$TUMOR <- "BREAST"
pilot[pilot$assigned_diagnosis == "C34.1",]$TUMOR <- "LUNG"
pilot[pilot$assigned_diagnosis %in% c("C18.7", "C19X", "C18.0"),]$TUMOR <- "COLORECTAL"

# Check tumor types
table(pilot$Trio, pilot$TUMOR, exclude = NULL) # No UNKNOWN tumor types in Trios left

# Subset to keep only trios
pilot <- pilot %>% filter(Trio == 1)  # 124

# Remove all "experimental" samples and one that fails contamination check (200000336)
rm_ids <- as.character(c(unique(pilot %>% filter(PILOT == "EXPT") %>% .$PATIENT_ID), unique(pilot %>% filter(modifier == "EXPERIMENTAL") %>% .$PATIENT_ID), "200000336"))
pilot <- pilot %>% filter(!PATIENT_ID %in% rm_ids)


############  Clean BRC pilot QC data ############ 

# Read QC data for BRC pilot samples
QC_BRC <- read.csv("/Users/MartinaMijuskovic/Documents/FFPE/Pilot data/BRC.QC.csv")

# Subset for clean BRC trios (36 total)
QC_BRC$PATIENT_ID <- as.character(QC_BRC$PATIENT_ID)
QC_BRC <- QC_BRC %>% filter(PATIENT_ID %in% pilot$PATIENT_ID)
dim(QC_BRC)  # 108

# Clean names
names(QC_BRC)[2:3] <- c("SAMPLE_WELL_ID", "SAMPLE_TYPE")

# Merge pilot data with QC data
pilot$PATIENT_ID <- as.character(pilot$PATIENT_ID)
QC_BRC <- inner_join(pilot, QC_BRC)

# Change names to merge with other trios properly
names(QC_BRC)[6] <- "BamPath"
names(QC_BRC)[10] <- "TumorType"

############   Get the correct (v4) BamPath for BRC ############  

# Read the current upload report, restrict to cancer, V4 and qc_passed
today <- Sys.Date()
system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
upload <- read.table(paste0("upload_report.", today, ".txt"), sep = "\t")
colnames(upload) <- as.character(fread(paste0("upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
#upload <- upload %>% filter(`Delivery Version` == "V4", Status == "qc_passed", Type %in% c("cancer germline", "cancer tumour"))
# Samples with non-pass status in bertha (upload report)
#missing <- QC_BRC %>% filter(is.na(BamPath)) %>% .$SAMPLE_WELL_ID
upload %>% filter(Platekey %in% missing, Status != "qc_passed") %>% dplyr::select(Platekey, Status)

upload <- upload %>% filter(`Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"))  # not all have "qc_passed" status, removing it
upload$Path <- as.character(upload$Path)

# # Correct the BamPath (this is the old path, to v2)
# QC_BRC$BamPath <- paste0("/genomes", QC_BRC$BamPath)

# Add new (v4) BAM paths to the BRC table
QC_BRC$BamPath <- upload[match(QC_BRC$SAMPLE_WELL_ID, upload$Platekey),]$Path

# List BRC trio samples with non-pass status
upload %>% filter(Platekey %in% QC_BRC$SAMPLE_WELL_ID, Status != "qc_passed") %>% dplyr::select(Platekey, Type, DeliveryID, `Delivery Date`, `Delivery Version`, `BAM Date`, Status)




############  Add BRC data to initial 26 trios ############ 

# Read QC data of initial 26 trios
QC_portal_trios <- read.csv("QC_portal_trios_final.csv")

# Change tumor type entries to all CAPS
QC_portal_trios$TumorType <- toupper(QC_portal_trios$TumorType)

# Merge BRC with initial 26 trios
QC_portal_trios$PATIENT_ID <- as.character(QC_portal_trios$PATIENT_ID)
QC_portal_trios <- bind_rows(QC_portal_trios, QC_BRC)

# Write out the full FFPE trio QC data table
write.csv(QC_portal_trios, file = "QC_portal_62_trios.csv", quote = F, row.names = F)



