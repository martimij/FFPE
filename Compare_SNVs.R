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
setwd("/home/mmijuskovic/FFPE/SNV_trio_comparison")

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

# Get VCF paths (HPC)   ------ original Illumina SNVs VCF! --- some folders contain two, one ending in "PASS.duprem.atomic.left.split.somatic.vcf.gz"
paths <- unlist(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.somatic.vcf.gz", sep = " ")
  system(command, intern = T)
}))
paths <- paths[!grepl("PASS.duprem.atomic.left.split.somatic.vcf.gz", paths)]
QC_portal_trios$SNV_VCF_path <- as.character(paths)


# Get json paths (HPC)
QC_portal_trios$json_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", sub("/genomes", "/genomes/analysis", x), "-iname *tiering.json", sep = " ")
  system(command, intern = T)
}))


### Function that extracts variant information from json files and calculates concordance between FF and FFPE
compareSNV <- function(patientID){
  
  # Get json paths
  ff_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$json_path
  ffpe_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$json_path
  
  ### Read json files: get chromosome, position, VAF, tier, etc
  # FF
  #ff_json <- fromJSON("LP2000907-DNA_F02_tiering.json", flatten = T) 
  ff_json <- fromJSON(ff_path, flatten = T)
  ff <- ff_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF)
  ff$tier <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ff$class <- sapply(1:dim(ff)[1], function(x){
    ff_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ff) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ff$KEY <- sapply(1:dim(ff)[1], function(x){
    paste(ff$chr[x], ff$pos[x], ff$ref[x], ff$alt[x], sep = "_")
  })

  # FFPE
  #ffpe_json <- fromJSON("LP2000906-DNA_H01_tiering.json", flatten = T)
  ffpe_json <- fromJSON(ffpe_path, flatten = T)
  ffpe <- ffpe_json %>% dplyr::select(
    reportedVariantCancer.chromosome, reportedVariantCancer.position, reportedVariantCancer.reference, reportedVariantCancer.alternate, 
    reportedVariantCancer.VAF)
  ffpe$tier <- sapply(1:dim(ffpe)[1], function(x){
    ffpe_json$reportedVariantCancer.reportEvents[[x]]$tier})
  ffpe$class <- sapply(1:dim(ffpe)[1], function(x){
    ffpe_json$reportedVariantCancer.reportEvents[[x]]$soNames})
  names(ffpe) <- c("chr", "pos", "ref", "alt", "VAF", "tier", "class")
  ffpe$KEY <- sapply(1:dim(ffpe)[1], function(x){
    paste(ffpe$chr[x], ffpe$pos[x], ffpe$ref[x], ffpe$alt[x], sep = "_")
  }) 
 
  # # Sanity check (!!! found duplicates - traced the error to the tiering pipeline)
  # sum(duplicated(ff))  #33
  # sum(duplicated(ff$KEY))  #37
  # sum(duplicated(ffpe))   #38
  # sum(duplicated(ffpe$KEY))  #41
 
  # Deduplicate
  ff <- ff[!duplicated(ff),]
  ffpe <- ffpe[!duplicated(ffpe),]
  
  # # Check duplicate keys
  # sum(duplicated(ff$KEY))
  # sum(duplicated(ffpe$KEY))
  # ff %>% filter(KEY %in% (ff[duplicated(ff$KEY),]$KEY))
  
  # Deduplicate variants with same keys
  ff <- ff[(!duplicated(ff$KEY)),]
  ffpe <- ffpe[(!duplicated(ffpe$KEY)),]
  
  # Summary table with concordance
  result <- data.frame(FF_TOTAL = dim(ff)[1], 
                       FFPE_TOTAL = dim(ffpe)[1], 
                       OVERLAP = sum(ff$KEY %in% ffpe$KEY), 
                       FF_UNIQ = ((dim(ff)[1]) - sum(ff$KEY %in% ffpe$KEY)), 
                       FFPE_UNIQ = ((dim(ffpe)[1]) - sum(ff$KEY %in% ffpe$KEY)), 
                       RECALL = ((sum(ff$KEY %in% ffpe$KEY))/(dim(ff)[1])), 
                       PRECISION = (sum(ff$KEY %in% ffpe$KEY))/(dim(ffpe)[1]))
 #write.table(result, file = paste0(patientID, "_SNV_concord", ".tsv"), row.names = F, quote = F, sep = "\t")
 return(result)
}

# Get all patient IDs
patientIDs <- as.character(unique(QC_portal_trios$PATIENT_ID))

# Remove 2 patient IDs for which I'm still missing FFPE json files
rmIDs <- QC_portal_trios %>% filter(SAMPLE_WELL_ID %in% c("LP3000079-DNA_G02", "LP3000074-DNA_C07")) %>% .$PATIENT_ID
rmIDs <- c(rmIDs, "218000014")
patientIDs <- patientIDs[!patientIDs %in% rmIDs]

# Run CNV comparison for each patient ID
SNV_summary <- lapply(patientIDs, compareSNV)  # errors :/ " no applicable method for 'select_' applied to an object of class "list""
SNV_summary <- lapply( c("212000006", "212000008", "212000009", "212000014", "212000015", "217000011", "217000028", "217000030", "217000052", "217000073", "218000002", "218000013", "218000014", "218000016", "218000017", "218000018", "218000021", "218000022", "218000024", "218000027", "218000032", "218000033", "218000036", "218000039"), compareSNV)
# there must be an error or malformed json somewhere (218000014 ID is problematic, malformed FFPE json file:
# "/genomes/analysis/by_date/2017-01-16/CANCP41020/CancerLP3000162-DNA_E01_NormalLP3000067-DNA_B09/SomaticVariations/1486209216/LP3000162-DNA_E01_tiering.json")

# Running with only 23 trios  ---- DOESN'T work, even more empty JSONs
SNV_summary <- lapply(patientIDs, compareSNV) 
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs

# Write out the resulting table
write.csv(SNV_summary, file = paste0("SNV_summary_", today, ".csv"), row.names = F, quote = F)



