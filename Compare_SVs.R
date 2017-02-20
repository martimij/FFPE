# Martina Mijuskovic
# FFPE project
# Reads the FFPE trio manifest, finds corresponding VCF files and compares Manta SVs between FF and FFPE using bedtools (on HPC)

library(dplyr)
library(VariantAnnotation)

# Working directory on the HPC
setwd("/home/mmijuskovic/FFPE/SV_trio_comparison")

today <- Sys.Date()

# Load the manifest (HPC)
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/CNV_trio_comparison/QC_portal_trios.csv")

# Subset for FF and FFPE samples
QC_portal_trios <- QC_portal_trios %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"))

# Get VCF paths (HPC)
QC_portal_trios$VCF_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))


### Function that compares SVs from FFPE trios
# Extracts only precise SVs that pass Manta filters, from regular chromosomes, with at least 3 supporting reads (either 3 PR or 3 SR minimum)
# Filters out SVs overlapping windowMasker, simple repeats or segmental duplications (either end of SV)
# Writes out the filtered SV table
compareSV <- function(patientID){
  
  # Get FF and FFPE VCF paths for specified patient ID and read into a Large CollapsedVCF object (VariantAnnotation package)
  ff_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$VCF_path
  ffpe_path <- QC_portal_trios %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$VCF_path
  ff_vcf <- readVcf(ff_path)
  ffpe_vcf <- readVcf(ffpe_path)
  
  ### FF info
  # Extract info fields
  SVinfo_ff <- as.data.frame(info(ff_vcf))
  # Get FILTER field
  SVinfo_ff$FILTER <- rowRanges(ff_vcf)$FILTER
  # Create Application variable (Canvas or Manta?)
  SVinfo_ff$Application <- ""
  SVinfo_ff[grepl("Canvas", rownames(SVinfo_ff)),]$Application <- "Canvas"
  SVinfo_ff[grepl("Manta", rownames(SVinfo_ff)),]$Application <- "Manta"
  # Extract ID
  SVinfo_ff$ID <- rownames(SVinfo_ff)
  # Remove filtered entries, keep Canvas only
  Canvas_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Canvas")
  # Extract chr, start, end
  Canvas_ff$Chr <- sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][3]})
  Canvas_ff$Start <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][4]})
  Canvas_ff$End <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][5]})
  Canvas_ff$Type <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][2]})
  # Keep LOSS/GAIN only
  Canvas_ff <- Canvas_ff %>% filter(Type %in% c("LOSS", "GAIN"))
  
  ### FFPE info
  # Extract info fields
  SVinfo_ffpe <- as.data.frame(info(ffpe_vcf))
  # Add filter field
  SVinfo_ffpe$FILTER <- rowRanges(ffpe_vcf)$FILTER
  # Create Application variable (Canvas or Manta?)
  SVinfo_ffpe$Application <- ""
  SVinfo_ffpe[grepl("Canvas", rownames(SVinfo_ffpe)),]$Application <- "Canvas"
  SVinfo_ffpe[grepl("Manta", rownames(SVinfo_ffpe)),]$Application <- "Manta"
  # Extract ID
  SVinfo_ffpe$ID <- rownames(SVinfo_ffpe)
  # Remove filtered entries, keep Canvas only, keep LOSS,GAIN only
  Canvas_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Canvas")
  # Extract chr, start, end
  Canvas_ffpe$Chr <- sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][3]})
  Canvas_ffpe$Start <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][4]})
  Canvas_ffpe$End <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][5]})
  Canvas_ffpe$Type <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][2]})
  # Keep LOSS/GAIN only
  Canvas_ffpe <- Canvas_ffpe %>% filter(Type %in% c("LOSS", "GAIN"))  # 7
  
  
  ### Call bedtools to compare FF and FFPE
  
  # Add new "Type2" ("+" is gain, "-" is loss, to separately compare them via bedtools)
  Canvas_ff$Type2 <- NA
  Canvas_ff$Type2 <- sapply(1:dim(Canvas_ff)[1], function(x){
    if (Canvas_ff$Type[x] == "LOSS") {Canvas_ff$Type2[x] <- "-"}
    else if (Canvas_ff$Type[x] == "GAIN") {Canvas_ff$Type2[x] <- "+"}
  })
  
  Canvas_ffpe$Type2 <- NA
  Canvas_ffpe$Type2 <- sapply(1:dim(Canvas_ffpe)[1], function(x){
    if (Canvas_ffpe$Type[x] == "LOSS") {Canvas_ffpe$Type2[x] <- "-"}
    else if (Canvas_ffpe$Type[x] == "GAIN") {Canvas_ffpe$Type2[x] <- "+"}
  })
  # Make bed files to find number of overlapping and non-overlapping bases (with bedtools)
  ff_bed <- Canvas_ff %>% dplyr::select(Chr, Start, End, ID)
  ff_bed$Score <- ""
  ff_bed <- cbind(ff_bed, (Canvas_ff %>% dplyr::select(Type2)))
  ffpe_bed <- Canvas_ffpe %>% dplyr::select(Chr, Start, End, ID)
  ffpe_bed$Score <- ""
  ffpe_bed <- cbind(ffpe_bed, (Canvas_ffpe %>% dplyr::select(Type2)))
  # Write bed files indicating gain as "+" strand and loss as "-" strand (to enable comparing them at once in bedtools)
  write.table(ff_bed, file = paste0(patientID, "_CNV_ff.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
  write.table(ffpe_bed, file = paste0(patientID, "_CNV_ffpe.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
  
  # Call bedtools to get FFPE overlap with FF
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -s -a", paste0(patientID, "_CNV_ff.bed"), "-b", paste0(patientID, "_CNV_ffpe.bed"), ">", paste0(patientID, "_CNV_ff_overlap.bed")), intern = T)
  ff_overlap <- read.table(paste0(patientID, "_CNV_ff_overlap.bed"), sep = "\t")
  names(ff_overlap) <- c(names(ff_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  
  # Call bedtools to get FF overlap with FFPE
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -s -a", paste0(patientID, "_CNV_ffpe.bed"), "-b", paste0(patientID, "_CNV_ff.bed"), ">", paste0(patientID, "_CNV_ffpe_overlap.bed")), intern = T)
  ffpe_overlap <- read.table(paste0(patientID, "_CNV_ffpe_overlap.bed"), sep = "\t")
  names(ffpe_overlap) <- c(names(ffpe_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  
  # Summary table
  result <- data.frame(PERCENT_RECALL_FF = sum(ff_overlap$BPoverlap) / sum(ff_overlap$BPTotal), PERCENT_RECALL_FFPE = sum(ffpe_overlap$BPoverlap) / sum(ffpe_overlap$BPTotal), BP_OVERLAP = sum(ff_overlap$BPoverlap), BP_FF_ONLY = (sum(ff_overlap$BPTotal) - sum(ff_overlap$BPoverlap)), BP_FFPE_ONLY = (sum(ffpe_overlap$BPTotal) - sum(ffpe_overlap$BPoverlap)))
  return(result)
}

# Get all patient IDs
patientIDs <- unique(QC_portal_trios$PATIENT_ID)

# Run CNV comparison for each patient ID
CNV_summary <- lapply(patientIDs, compareCNV)
CNV_summary <- bind_rows(CNV_summary)
CNV_summary$PATIENT_ID <- patientIDs

# Write out the resulting table
write.csv(CNV_summary, file = paste0("CNV_summary_", today, ".csv"), row.names = F, quote = F)

# Read the result table (local copy)
#CNV_summary <- read.csv(paste0("CNV_summary_", today, ".csv"))



### Put all data together (QC, CNV)

# Add QC details to the CNV summary table
QC_table <- QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE") %>% select(PATIENT_ID, CENTER_CODE, TumorType, SAMPLE_WELL_ID, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER)
names(QC_table)[4:11] <- paste0("FFPE_", names(QC_table)[4:11])
QC_table <- left_join(QC_table, (QC_portal_trios %>% filter(SAMPLE_TYPE == "FF") %>% select(PATIENT_ID, SAMPLE_WELL_ID, LIBRARY_TYPE, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER)), by = "PATIENT_ID")
names(QC_table)[12:20] <- paste0("FF_", names(QC_table)[12:20])
CNV_summary <- left_join(CNV_summary, QC_table, by = "PATIENT_ID")

# Calculate normalized overlap
CNV_summary$BP_OVERLAP <- as.numeric(CNV_summary$BP_OVERLAP)
CNV_summary$BP_FF_ONLY <- as.numeric(CNV_summary$BP_FF_ONLY)
CNV_summary$BP_FFPE_ONLY <- as.numeric(CNV_summary$BP_FFPE_ONLY)
CNV_summary$TOTAL_BP <- CNV_summary$BP_OVERLAP + CNV_summary$BP_FF_ONLY + CNV_summary$BP_FFPE_ONLY

# Write the full table
write.csv(CNV_summary, file = paste0("Full_CNV_summary_", today, ".csv"), row.names = F, quote = F)
