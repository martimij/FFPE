# Martina Mijuskovic
# FFPE project
# Collect QC and clinical data for sequenced FFPE trios

library("data.table")
library("dplyr")
library("ggvis")
library("zoo")
library("ggplot2")

############  Read data (Jan 26 2017) ############  



############  Subset trios ############  

# Flag trios (ALL)
trios <- sapply(unique(QC_data$PATIENT_ID), function(x){
  if (sum(c("GL", "FFPE", "FF") %in% QC_data[QC_data$PATIENT_ID == x,]$SAMPLE_TYPE) == 3) {
    return(1)
  }
  else {
    return(0)
  }
})

QC_data$Trio <- 0
QC_data[QC_data$PATIENT_ID %in% names(trios[trios == 1]),]$Trio <- 1

# Total number of trios in the database
dim(QC_data %>% filter(Trio == 1))  # 325

# Flag trios where all three samples pass QC before sequencing
QC_data_PassToSeq <- QC_data %>% filter(PASS_TO_SEQ == "Pass")   # 1091 samples (of 1397) pass QC for sequencing

triosQC <- sapply(unique(QC_data_PassToSeq$PATIENT_ID), function(x){
  if (sum(c("GL", "FFPE", "FF") %in% QC_data_PassToSeq[QC_data_PassToSeq$PATIENT_ID == x,]$SAMPLE_TYPE) == 3) {
    return(1)
  }
  else {
    return(0)
  }
})

QC_data_PassToSeq$Trio <- 0
QC_data_PassToSeq[QC_data_PassToSeq$PATIENT_ID %in% names(triosQC[triosQC == 1]),]$Trio <- 1

# Total number of trios that pass pre-sequencing QC
QC_trios <- QC_data_PassToSeq %>% filter(Trio == 1)  # 153 (51 patients total) --- NO, carefull, there are DUPLICATES

############  Check status of trios in the upload report ############  

