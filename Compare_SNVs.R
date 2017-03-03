# Martina Mijuskovic
# FFPE project
# Compares SNV calls between FFPE and FF samples from FFPE trios

# NOTE: Installing R packages on HPC: use lib = "~/R/x86_64-pc-linux-gnu-library/3.3"

library(dplyr)
library(reshape)
library(ggplot2)
library(scales)
library(R.utils)
library(jsonlite)

# Working directory on the HPC
setwd("/home/mmijuskovic/FFPE/SV_trio_comparison")

# Today's date
today <- Sys.Date()

### For plots
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)


############# Find Domain1/2 variants, subset VCFs ############# 

# Load the manifest (HPC)
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/CNV_trio_comparison/QC_portal_trios.csv")
# QC_portal_trios <- read.csv("QC_portal_trios_final.csv") # local

# Subset for FF and FFPE samples
QC_portal_trios <- QC_portal_trios %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"))

# Get VCF paths (HPC)   ------ FIX for SNVs!!!!
QC_portal_trios$VCF_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))

# Get json paths (HPC)
QC_portal_trios$json_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", sub("/genomes", "/genomes/analysis", x), "-iname *tiering.json", sep = " ")
  system(command, intern = T)
}))


### Function that extracts Domain 1/2 variant positions from json files, 
# writes them into interval_list files for GATK,
# calls GATK SelectVariants tool to subset VCFs (write to my home on HPC),
# calls GATK GenotypeConcordance tool to find FFPE to FF concordance

compareSNV <- function(patientID){
  
  # Get json paths
  ff_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$json_path
  ffpe_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$json_path
  
  # Read json files: get chromosome, position, VAF
  ff <- fromJSON("LP2000907-DNA_F02_tiering.json")$reportedVariantCancer %>% select(chromosome, position, VAF)
  # Get tier and variant classification ---- continue here
  str(fromJSON("LP2000907-DNA_F02_tiering.json")$reportedVariantCancer$reportEvents[1])

  # Reducte to tier1 and tier2
  
  # Write the interval_list file
  
  # Call GATK to subset VCFs
  
  # Call GATK to calculate concordance


}





