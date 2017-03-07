### Martina Mijuskovic
# FFPE project
# Compile and clean QC and clinical data for additional BRC pilot FFPE trios

library(dplyr)
library(ggvis)
library(ggplot2)
library(VariantAnnotation)

rm(list=ls())

### For ggplot
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)


############  Read and clean data ############  

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
pilot <- pilot %>% select(manifest_id, study, GelId, type, platekey, base_dir, modifier, assigned_diagnosis)

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
pilot[pilot$PATIENT_ID %in% names(trios[trios == 1]),]$Trio <- 1

# Number of trios
sum(trios)  # 41
sum(pilot$Trio)/3  # 41.33333  (one sample extra?)

# Find the origin of the extra sample
sum(duplicated(pilot$SAMPLE_WELL_ID))  # 0
table(table(pilot$PATIENT_ID))  # 20 patients with 2 samples, 40 with trios, 1 with 4 samples
table(pilot$PATIENT_ID)  # 200000960 has 4 samples
pilot %>% filter(PATIENT_ID == "200000960")  # this patient has 2 FFPE samples listed (one may have failed but I have no info; older one?)

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


# Add tumor info from the 2nd file
pilot <- left_join(pilot, (pilot_extra2 %>% select(SAMPLE_WELL_ID, TUMOR, SAMPLE_TYPE, CENTER, PILOT)))
table(pilot$Trio, pilot$TUMOR)
pilot[is.na(pilot$TUMOR),]$TUMOR <- "UNKNOWN"
pilot[pilot$TUMOR == "UNKOWN",]$TUMOR <- "UNKNOWN"
pilot$TUMOR <- as.character(pilot$TUMOR)

# PATIENT_IDs with missing tumor info
unique(pilot %>% filter(Trio == 1, TUMOR == "UNKNOWN") %>% .$PATIENT_ID)  # 200000926 200000929 200000953 200000954 200000960 200001299 200001310
pilot %>% filter(PATIENT_ID %in% (unique(pilot %>% filter(Trio == 1, TUMOR == "UNKNOWN") %>% .$PATIENT_ID)))





