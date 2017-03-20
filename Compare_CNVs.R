# Martina Mijuskovic
# FFPE project
# Reads the FFPE trio manifest, finds corresponding VCF files and compares Canvas CNVs between FF and FFPE using bedtools (on HPC)


### Install R packages if not present
# install.packages("dplyr", lib = "~/R/x86_64-pc-linux-gnu-library/3.3")
# NOTE that to install XML and Rcurl, libcurl4-openssl-dev and libxml2-dev are needed on the system
# See: http://stackoverflow.com/questions/10965755/genomicfeatures-package-installation-trouble
# install.packages(c("XML", "RCurl"), lib = "~/R/x86_64-pc-linux-gnu-library/3.3")
# source("https://bioconductor.org/biocLite.R") # Install in the local library "~/R/x86_64-pc-linux-gnu-library/3.3"
# biocLite("VariantAnnotation", lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.3", lib = "~/R/x86_64-pc-linux-gnu-library/3.3")
# install.packages("VennDiagram", lib = "~/R/x86_64-pc-linux-gnu-library/3.3")

library(dplyr)
library(VariantAnnotation)
library(ggplot2)
library(reshape)
library(scales)

# Working directory on the HPC
setwd("/home/mmijuskovic/FFPE/CNV_trio_comparison")

today <- Sys.Date()

# Load the manifest (HPC)
#QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/CNV_trio_comparison/QC_portal_trios.csv") # This is the file with initial 26 trios
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/QC_portal_62_trios.csv")  # This is the file with final 62 trios (BRC pilot FFPE trios added)
  
# Subset for FF and FFPE samples
QC_portal_trios <- QC_portal_trios %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"))

# Get VCF paths (HPC)
QC_portal_trios$VCF_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))


# Function that compares CNVs from FFPE trios
compareCNV <- function(patientID){
  
  # Get FF and FFPE VCF paths for specified patient ID and read into a Large CollapsedVCF object (VariantAnnotation package)
  # NOTE that some SV.VCFs are original ones from Illumina (somatic, but contain also NORMAL sample) and some are processed by BERTHA, removing the NORMAL sample)
  # This does not affect this comparison though, as it uses only the variant INFO, not genotypes
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
  
  # Check that both ff and ffpe have calls; if not create empty result table and move to the next sample
  # I'm putting 1 bp call in respective sample for a placeholder (and % calculations later on)
  QC_portal_trios$PATIENT_ID <- as.character(QC_portal_trios$PATIENT_ID)
  if (dim(Canvas_ff)[1] == 0) {
    print(paste0("No ff calls in ", patientID))
    result <- data.frame(PERCENT_RECALL_FF = NA, PERCENT_RECALL_FFPE = 0, BP_OVERLAP = 0, BP_FF_ONLY = 0, BP_FFPE_ONLY = 1)
    return(result)
  }
  if (dim(Canvas_ffpe)[1] == 0) {
    print(paste0("No ffpe calls in ", patientID))
    result <- data.frame(PERCENT_RECALL_FF = 0, PERCENT_RECALL_FFPE = NA, BP_OVERLAP = 0, BP_FF_ONLY = 1, BP_FFPE_ONLY = 0)
    return(result)
  }
  
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
CNV_summary$PATIENT_ID <- as.character(CNV_summary$PATIENT_ID)
QC_table <- QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE") %>% dplyr::select(PATIENT_ID, CENTER_CODE, TumorType, SAMPLE_WELL_ID, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER)
names(QC_table)[4:11] <- paste0("FFPE_", names(QC_table)[4:11])
QC_table <- left_join(QC_table, (QC_portal_trios %>% filter(SAMPLE_TYPE == "FF") %>% dplyr::select(PATIENT_ID, SAMPLE_WELL_ID, LIBRARY_TYPE, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER)), by = "PATIENT_ID")
names(QC_table)[12:20] <- paste0("FF_", names(QC_table)[12:20])
CNV_summary <- left_join(CNV_summary, QC_table, by = "PATIENT_ID")

# Calculate normalized overlap
CNV_summary$BP_OVERLAP <- as.numeric(CNV_summary$BP_OVERLAP)
CNV_summary$BP_FF_ONLY <- as.numeric(CNV_summary$BP_FF_ONLY)
CNV_summary$BP_FFPE_ONLY <- as.numeric(CNV_summary$BP_FFPE_ONLY)
CNV_summary$TOTAL_BP <- CNV_summary$BP_OVERLAP + CNV_summary$BP_FF_ONLY + CNV_summary$BP_FFPE_ONLY

# Write the full table
write.csv(CNV_summary, file = paste0("Full_CNV_summary_", today, ".csv"), row.names = F, quote = F)



### Create plots

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)

# Barplots of overlapping and unique bp (normalized) for all 26 trios
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
CNV_summary_m <- as.data.frame(t(CNV_summary %>% dplyr::select(PATIENT_ID, BP_OVERLAP, BP_FF_ONLY, BP_FFPE_ONLY)))
#names(CNV_summary_m) <- CNV_summary_m[1,]
names(CNV_summary_m) <- as.matrix(CNV_summary_m)[1,]
CNV_summary_m <- CNV_summary_m[2:4,]
CNV_summary_m <- melt(cbind(CNV_summary_m, ind = rownames(CNV_summary_m)), id.vars = c('ind')) # attributes are not identical across measure variables; they will be dropped
CNV_summary_m$value <- as.numeric(CNV_summary_m$value)
# Plot (needs package "scales")
ggplot(CNV_summary_m,aes(x = variable, y = value, fill = ind)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + blank


# Recall of FF vs recall of FFPE
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = PERCENT_RECALL_FFPE)) + geom_jitter() + geom_smooth(method = "lm", se = T) + labs(x = "Percent Recall of FF", y = "Percent Recall of FFPE") 
# Spearman correlation coefficient (remove NAs!)
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$PERCENT_RECALL_FFPE, method = "spearman")  #  0.5186325

# RECALL of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, colored by center
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$FFPE_AT_DROP, method = "spearman")  # 0.1863566
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.05025641
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # 0.1108261
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # -0.1432479
# Recall of FF by AT dropout of FFPE, colored by tumor type
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_AT_DROP, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
# Recall by tumor purity
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$FFPE_TUMOUR_PURITY, method = "spearman")  # 0.2920236
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FF_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FF Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(CNV_summary$PERCENT_RECALL_FF, CNV_summary$FF_TUMOUR_PURITY, method = "spearman")  # 0.373867
# # Recall by difference in tumor purity between FF and FFPE
# CNV_summary$PurityDiff <- CNV_summary$FFPE_TUMOUR_PURITY-CNV_summary$FF_TUMOUR_PURITY
# ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = PurityDiff, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "Difference In Tumour Purity")  + theme(legend.title=element_blank())
# ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_AT_DROP, col = factor(PurityDiff > 15))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout")  + theme(legend.title=element_blank())
# Recall by high vs low FF and FFPE purity (both > 50)
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_AT_DROP, col = factor((FFPE_TUMOUR_PURITY > 50) & (FF_TUMOUR_PURITY > 50)))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
# Divide samples into good and bad FF recall (>50%), see which variables make a difference
ggplot(data=CNV_summary, aes(x=(PERCENT_RECALL_FF>0.5), y=FFPE_AT_DROP, fill=(PERCENT_RECALL_FF>0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High FF CNV Recall (>50%)", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
ggplot(data=CNV_summary, aes(x=(PERCENT_RECALL_FF>0.5), y=FFPE_COVERAGE_HOMOGENEITY, fill=(PERCENT_RECALL_FF>0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High FF CNV Recall (>50%)", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
ggplot(data=CNV_summary, aes(x=(PERCENT_RECALL_FF>0.5), y=FFPE_TUMOUR_PURITY, fill=(PERCENT_RECALL_FF>0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High FF CNV Recall (>50%)", y = "FFPE Tumour Purity") + theme(legend.title=element_blank())
ggplot(data=CNV_summary, aes(x=(PERCENT_RECALL_FF>0.5), y=FF_TUMOUR_PURITY, fill=(PERCENT_RECALL_FF>0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High FF CNV Recall (>50%)", y = "FF Tumour Purity") + theme(legend.title=element_blank())
# Recall by % chimeric fragments
ggplot(data=CNV_summary, aes(x=(PERCENT_RECALL_FF>0.5), y=FFPE_CHIMERIC_PER, fill=(PERCENT_RECALL_FF>0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High FF CNV Recall (>50%)", y = "FFPE Chimeric Frag") + theme(legend.title=element_blank())
# Recall by tumor type (not enough data...)
ggplot(data=CNV_summary, aes(x=TumorType, y=PERCENT_RECALL_FF)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "Tumor Type", y = "Percent CNV Recall") + theme(legend.title=element_blank())
# Recall by GMC, color by tumor type
ggplot(data=CNV_summary, aes(x=CENTER_CODE, y=PERCENT_RECALL_FF)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "GMC", y = "Percent CNV Recall") + theme(legend.title=element_blank())
# Recall by FF library type
ggplot(data=CNV_summary, aes(x=FF_LIBRARY_TYPE, y=PERCENT_RECALL_FF)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "FF Library Type", y = "Percent CNV Recall") + theme(legend.title=element_blank())
# Recall by fragment size
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FF, y = FFPE_AV_FRAGMENT_SIZE_BP, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Fragment Size") + theme(legend.title=element_blank())



# PRECISION by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, colored by center
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FFPE, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "FFPE Precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FFPE, CNV_summary$FFPE_AT_DROP, method = "spearman")  # 0.6209609
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FFPE, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "FFPE Precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FFPE, CNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.2841026
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FFPE, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "FFPE Precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FFPE, CNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # 0.3519754
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FFPE, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "FFPE Precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(CNV_summary$PERCENT_RECALL_FFPE, CNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # -0.2622222
# Recall by tumor purity
ggplot(CNV_summary, aes(x = PERCENT_RECALL_FFPE, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "FFPE Precision", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(CNV_summary$PERCENT_RECALL_FFPE, CNV_summary$FFPE_TUMOUR_PURITY, method = "spearman")  # 0.516262
# PRECISION by GMC, color by tumor type
ggplot(data=CNV_summary, aes(x=CENTER_CODE, y=PERCENT_RECALL_FFPE)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "GMC", y = "FFPE Precision") + theme(legend.title=element_blank())
# PRECISION by FF library type
ggplot(data=CNV_summary, aes(x=FF_LIBRARY_TYPE, y=PERCENT_RECALL_FFPE)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "FF Library Type", y = "FFPE Precision") + theme(legend.title=element_blank())


### FF vs FFPE tumor purity
ggplot(data=CNV_summary, aes(x=FFPE_TUMOUR_PURITY, y=FF_TUMOUR_PURITY)) + geom_jitter(aes(col = TumorType)) + labs(x = "FFPE Tumor Purity", y = "FF Tumor Purity") + theme(legend.title=element_blank()) + regr_line
cor(CNV_summary$FFPE_TUMOUR_PURITY, CNV_summary$FF_TUMOUR_PURITY, method = "spearman")  #  0.5927387








### Test function locally

# patientIDs <- c("212000009", "217000028")
# 
# # # Test get VCF paths
# dummy_bam_path <- "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/"
# QC_portal_trios2 <- QC_portal_trios %>% filter(SAMPLE_WELL_ID %in% c("LP2000906-DNA_A02", "LP2000907-DNA_B01", "LP3000070-DNA_G01", "LP3000079-DNA_B02"))
# QC_portal_trios2$VCF_path <- sapply(rep(dummy_bam_path, 4), function(x){
#   command <- paste("find", x, "-iname *.SV.vcf.gz", sep = " ")
#   system(command, intern = T)
# })
# QC_portal_trios2$VCF_path <- QC_portal_trios2$VCF_path[1:4]
# 
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
#   result <- data.frame(PERCENT_RECALL_FF = sum(ff_overlap$BPoverlap) / sum(ff_overlap$BPTotal), PERCENT_RECALL_FFPE = sum(ffpe_overlap$BPoverlap) / sum(ffpe_overlap$BPTotal), BP_OVERLAP = sum(ff_overlap$BPoverlap), BP_FF_ONLY = (sum(ff_overlap$BPTotal) - sum(ff_overlap$BPoverlap)), BP_FFPE_ONLY = (sum(ffpe_overlap$BPTotal) - sum(ffpe_overlap$BPoverlap)))
#   return(result)
# }
# 
# # Run CNV comparison for each patient ID
# CNV_summary <- lapply(patientIDs, compareCNV)
# CNV_summary <- bind_rows(CNV_summary)
# CNV_summary$PATIENT_ID <- patientIDs
# 
# # Write out the resulting table
# write.csv(CNV_summary, col.names = T, row.names = F, quote = F)
