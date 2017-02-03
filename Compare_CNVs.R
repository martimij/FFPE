# Martina Mijuskovic
# FFPE project
# Reads the FFPE trio manifest, finds corresponding VCF files and compares Canvas CNVs between FF and FFPE using bedtools (on HPC)

library(dplyr)
library(VariantAnnotation)

setwd("/home/mmijuskovic/FFPE")

# Load the manifest
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/QC_portal_trios.csv")

# Get VCF paths
QC_portal_trios$VCF_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))

<<<<<<< HEAD
# # # Test get VCF paths
# dummy_bam_path <- "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/"
# QC_portal_trios2 <- QC_portal_trios %>% filter(SAMPLE_WELL_ID %in% c("LP2000906-DNA_A02", "LP2000907-DNA_B01", "LP3000070-DNA_G01", "LP3000079-DNA_B02"))
# QC_portal_trios2$VCF_path <- sapply(rep(dummy_bam_path, 4), function(x){
#   command <- paste("find", x, "-iname *.SV.vcf.gz", sep = " ")
#   system(command, intern = T)
# })
# QC_portal_trios2$VCF_path <- QC_portal_trios2$VCF_path[1:4]

# Function that compares CNVs from FFPE trios (VCFs copied locally)
=======
# # Test get VCF paths
dummy_bam_path <- "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/"
QC_portal_trios2 <- QC_portal_trios %>% filter(SAMPLE_WELL_ID %in% c("LP2000906-DNA_A02", "LP2000907-DNA_B01", "LP3000070-DNA_G01", "LP3000079-DNA_B02"))
QC_portal_trios2$VCF_path <- sapply(rep(dummy_bam_path, 4), function(x){
  command <- paste("find", x, "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
})
QC_portal_trios2$VCF_path <- QC_portal_trios2$VCF_path[1:4]

# Function that compares CNVs from FFPE trios (VCFs copied locally)

>>>>>>> 775483df6c951699042bffdf0dd2a6fac8950712
compareCNV <- function(patientID){
  
  # Get FF and FFPE VCF paths for specified patient ID and read into a Large CollapsedVCF object (VariantAnnotation package)
  ff_path <- QC_portal_trios2 %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$VCF_path
  ffpe_path <- QC_portal_trios2 %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$VCF_path
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
  
  # Add new type coding ("+" is gain, "-" is loss)
  Canvas_ff$Type2 <- NA
  Canvas_ff$Type2 <- sapply(1:dim(Canvas_ff)[1], function(x){
    if (Canvas_ff$Type[x] == "LOSS") {Canvas_ff$Type2[x] <- "-"}
    else if (Canvas_ff$Type[x] == "GAIN") {Canvas_ff$Type2[x] <- "+"}
  })
  # Add new type coding ("+" is gain, "-" is loss)
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
  names(ffpe_overlap) <- c(names(ff_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  
  # Summary table
  result <- data.frame(PERCENT_RECALL_FF = sum(ff_overlap$BPoverlap) / sum(ff_overlap$BPTotal), PERCENT_RECALL_FFPE = sum(ffpe_overlap$BPoverlap) / sum(ffpe_overlap$BPTotal), BP_OVERLAP = sum(ff_overlap$BPoverlap), BP_FF_ONLY = (sum(ff_overlap$BPTotal) - sum(ff_overlap$BPoverlap)), BP_FFPE_ONLY = (sum(ffpe_overlap$BPTotal) - sum(ffpe_overlap$BPoverlap)))
  return(result)
}

# Get all patient IDs
<<<<<<< HEAD
patientIDs <- unique(QC_portal_trios$PATIENT_ID)
=======
patientIDs <- unique(QC_portal_trios2$PATIENT_ID)
>>>>>>> 775483df6c951699042bffdf0dd2a6fac8950712

# Run CNV comparison for each patient ID
CNV_summary <- lapply(patientIDs, compareCNV)
names(CNV_summary) <- patientIDs
<<<<<<< HEAD
# Convert to data.frame --- continue here

# Write the output table


=======
>>>>>>> 775483df6c951699042bffdf0dd2a6fac8950712

# ### Test function locally
# 
# compareCNV <- function(patientID){
#   
#   # Get FF and FFPE VCF paths for specified patient ID and read into a Large CollapsedVCF object (VariantAnnotation package)
#   ff_path <- QC_portal_trios2 %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FF") %>% .$VCF_path
#   ffpe_path <- QC_portal_trios2 %>% filter(PATIENT_ID == patientID, SAMPLE_TYPE == "FFPE") %>% .$VCF_path
#   ff_vcf <- readVcf(ff_path)
#   ffpe_vcf <- readVcf(ffpe_path)
#   
#   ### FF info
#   # Extract info fields
#   SVinfo_ff <- as.data.frame(info(ff_vcf))
#   # Get FILTER field
#   SVinfo_ff$FILTER <- rowRanges(ff_vcf)$FILTER
#   # Create Application variable (Canvas or Manta?)
#   SVinfo_ff$Application <- ""
#   SVinfo_ff[grepl("Canvas", rownames(SVinfo_ff)),]$Application <- "Canvas"
#   SVinfo_ff[grepl("Manta", rownames(SVinfo_ff)),]$Application <- "Manta"
#   # Extract ID
#   SVinfo_ff$ID <- rownames(SVinfo_ff)
#   # Remove filtered entries, keep Canvas only
#   Canvas_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Canvas")
#   # Extract chr, start, end
#   Canvas_ff$Chr <- sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][3]})
#   Canvas_ff$Start <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][4]})
#   Canvas_ff$End <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][5]})
#   Canvas_ff$Type <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][2]})
#   # Keep LOSS/GAIN only
#   Canvas_ff <- Canvas_ff %>% filter(Type %in% c("LOSS", "GAIN"))
#   
#   ### FFPE info
#   # Extract info fields
#   SVinfo_ffpe <- as.data.frame(info(ffpe_vcf))
#   # Add filter field
#   SVinfo_ffpe$FILTER <- rowRanges(ffpe_vcf)$FILTER
#   # Create Application variable (Canvas or Manta?)
#   SVinfo_ffpe$Application <- ""
#   SVinfo_ffpe[grepl("Canvas", rownames(SVinfo_ffpe)),]$Application <- "Canvas"
#   SVinfo_ffpe[grepl("Manta", rownames(SVinfo_ffpe)),]$Application <- "Manta"
#   # Extract ID
#   SVinfo_ffpe$ID <- rownames(SVinfo_ffpe)
#   # Remove filtered entries, keep Canvas only, keep LOSS,GAIN only
#   Canvas_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Canvas")
#   # Extract chr, start, end
#   Canvas_ffpe$Chr <- sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][3]})
#   Canvas_ffpe$Start <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][4]})
#   Canvas_ffpe$End <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][5]})
#   Canvas_ffpe$Type <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][2]})
#   # Keep LOSS/GAIN only
#   Canvas_ffpe <- Canvas_ffpe %>% filter(Type %in% c("LOSS", "GAIN"))  # 7
#   
#   
#   ### Call bedtools to compare FF and FFPE
#   
#   # Add new type coding ("+" is gain, "-" is loss)
#   Canvas_ff$Type2 <- NA
#   Canvas_ff$Type2 <- sapply(1:dim(Canvas_ff)[1], function(x){
#     if (Canvas_ff$Type[x] == "LOSS") {Canvas_ff$Type2[x] <- "-"}
#     else if (Canvas_ff$Type[x] == "GAIN") {Canvas_ff$Type2[x] <- "+"}
#   })
#   # Add new type coding ("+" is gain, "-" is loss)
#   Canvas_ffpe$Type2 <- NA
#   Canvas_ffpe$Type2 <- sapply(1:dim(Canvas_ffpe)[1], function(x){
#     if (Canvas_ffpe$Type[x] == "LOSS") {Canvas_ffpe$Type2[x] <- "-"}
#     else if (Canvas_ffpe$Type[x] == "GAIN") {Canvas_ffpe$Type2[x] <- "+"}
#   })
#   # Make bed files to find number of overlapping and non-overlapping bases (with bedtools)
#   ff_bed <- Canvas_ff %>% dplyr::select(Chr, Start, End, ID)
#   ff_bed$Score <- ""
#   ff_bed <- cbind(ff_bed, (Canvas_ff %>% dplyr::select(Type2)))
#   ffpe_bed <- Canvas_ffpe %>% dplyr::select(Chr, Start, End, ID)
#   ffpe_bed$Score <- ""
#   ffpe_bed <- cbind(ffpe_bed, (Canvas_ffpe %>% dplyr::select(Type2)))
#   # Write bed files indicating gain as "+" strand and loss as "-" strand (to enable comparing them at once in bedtools)
#   write.table(ff_bed, file = paste0(patientID, "_CNV_ff.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
#   write.table(ffpe_bed, file = paste0(patientID, "_CNV_ffpe.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
#   
#   # Call bedtools to get FFPE overlap with FF
#   system(paste("bedtools coverage -s -a", paste0(patientID, "_CNV_ff.bed"), "-b", paste0(patientID, "_CNV_ffpe.bed"), ">", paste0(patientID, "_CNV_ff_overlap.bed")), intern = T)
#   ff_overlap <- read.table(paste0(patientID, "_CNV_ff_overlap.bed"), sep = "\t")
#   names(ff_overlap) <- c(names(ff_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
#   
#   # Call bedtools to get FF overlap with FFPE
#   system(paste("bedtools coverage -s -a", paste0(patientID, "_CNV_ffpe.bed"), "-b", paste0(patientID, "_CNV_ff.bed"), ">", paste0(patientID, "_CNV_ffpe_overlap.bed")), intern = T)
#   ffpe_overlap <- read.table(paste0(patientID, "_CNV_ffpe_overlap.bed"), sep = "\t")
#   names(ffpe_overlap) <- c(names(ff_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
#   
#   # Summary table
#   result <- data.frame(PATIENT_ID = patientID, PERCENT_RECALL_FF = sum(ff_overlap$BPoverlap) / sum(ff_overlap$BPTotal), PERCENT_RECALL_FFPE = sum(ffpe_overlap$BPoverlap) / sum(ffpe_overlap$BPTotal), BP_OVERLAP = sum(ff_overlap$BPoverlap), BP_FF_ONLY = (sum(ff_overlap$BPTotal) - sum(ff_overlap$BPoverlap)), BP_FFPE_ONLY = (sum(ffpe_overlap$BPTotal) - sum(ffpe_overlap$BPoverlap)))
#   return(result)
# }
# 
