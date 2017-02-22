# Martina Mijuskovic
# FFPE project
# Collect QC and clinical data for sequenced FFPE trios
# Develop the comparison method for Canvas CNV output between FF and FFPE samples

library("data.table")
library("dplyr")
library("ggvis")
library("ggplot2")

# source("https://bioconductor.org/biocLite.R")
# biocLite("VariantAnnotation")
library(VariantAnnotation)

rm(list=ls())

### For ggplot
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)


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

# Manually add the BAM path for one sample (FFPE from 217000030, LP3000079-DNA_H01) - passes QC in v2 but not in v4, we are keeping it
#upload %>% filter(Platekey == "LP3000079-DNA_H01") # Chose v4 path from unfiltered upload report
QC_portal_trios[QC_portal_trios$SAMPLE_WELL_ID == "LP3000079-DNA_H01",]$BamPath <- "/genomes/by_date/2016-12-13/CANCP40747/CancerLP3000079-DNA_H01_NormalLP3000067-DNA_H08"

# Write the trios QC portal data
write.csv(QC_portal_trios, file = "QC_portal_trios.csv", quote = F, row.names = F)


############  Add PCR dup percentage (HPC) ############

# Get Metrics file paths (HPC) - NOTE that there is none for the germline
QC_portal_trios$Metrics_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/*Metrics.csv"), sep = " ")
  system(command, intern = T)
}))

# Get 2 missing paths manually
# /genomes/by_date/2017-01-16/CANCP41020/CancerLP3000162-DNA_E01_NormalLP3000067-DNA_B09/LP3000067-DNA_B09_LP3000162-DNA_E01.summary.csv (LP3000162-DNA_E01)
# /genomes/by_date/2017-01-16/CANCP41019/CancerLP3000162-DNA_D01_NormalLP3000067-DNA_D04/LP3000067-DNA_D04_LP3000162-DNA_D01.summary.csv (LP3000162-DNA_D01)
QC_portal_trios[QC_portal_trios$SAMPLE_WELL_ID == "LP3000162-DNA_E01",]$Metrics_path <- "/genomes/by_date/2017-01-16/CANCP41020/CancerLP3000162-DNA_E01_NormalLP3000067-DNA_B09/LP3000067-DNA_B09_LP3000162-DNA_E01.summary.csv"
QC_portal_trios[QC_portal_trios$SAMPLE_WELL_ID == "LP3000162-DNA_D01",]$Metrics_path <- "/genomes/by_date/2017-01-16/CANCP41019/CancerLP3000162-DNA_D01_NormalLP3000067-DNA_D04/LP3000067-DNA_D04_LP3000162-DNA_D01.summary.csv"


# Read the csv file with the metrics and add PCR % info to the QC table
QC_portal_trios$PCR_DUPL <- NA
QC_portal_trios[QC_portal_trios$SAMPLE_TYPE != "GL",]$PCR_DUPL <- sapply(1:dim(QC_portal_trios[QC_portal_trios$SAMPLE_TYPE != "GL",])[1], function(x){
  as.character(read.csv(QC_portal_trios[QC_portal_trios$SAMPLE_TYPE != "GL",]$Metrics_path[x], header = F, skip = 13, nrows=1)$V2)
})
# Correct for those with shifted rows
QC_portal_trios$PCR_DUPL2 <- NA
QC_portal_trios[QC_portal_trios$SAMPLE_TYPE != "GL",]$PCR_DUPL2 <- sapply(1:dim(QC_portal_trios[QC_portal_trios$SAMPLE_TYPE != "GL",])[1], function(x){
  as.character(read.csv(QC_portal_trios[QC_portal_trios$SAMPLE_TYPE != "GL",]$Metrics_path[x], header = F, skip = 15, nrows=1)$V2)
})
# Clean up PCR_DUPL
QC_portal_trios[!grepl("%", QC_portal_trios$PCR_DUPL),]$PCR_DUPL <- QC_portal_trios[!grepl("%", QC_portal_trios$PCR_DUPL),]$PCR_DUPL2
x <- dim(QC_portal_trios)[2]
QC_portal_trios <- QC_portal_trios[-x]

# Write the trios QC portal data
write.csv(QC_portal_trios, file = "QC_portal_trios_final.csv", quote = F, row.names = F)



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


############  Explore SVs from VCF ############  

# Load trios QC data
QC_portal_trios <- read.csv("QC_portal_trios.csv")

ff_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06_LP3000070-DNA_G01.somatic.SV.vcf.gz")
ffpe_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06_LP3000079-DNA_B02.somatic.SV.vcf.gz")

# SV types by filter
t(table(ff_vcf@info@listData$SVTYPE, ff_vcf@fixed@listData$FILTER, exclude = NULL))
t(table(ffpe_vcf@info@listData$SVTYPE, ffpe_vcf@fixed@listData$FILTER, exclude = NULL))



############  Compare CANVAS CNVs ############  

### FF
# Extract INFO table (all SVs)
SVinfo_ff <- as.data.frame(info(ff_vcf))
# Add filter field
SVinfo_ff$FILTER <- rowRanges(ff_vcf)$FILTER
# Create Application variable (Canvas or Manta?)
SVinfo_ff$Application <- ""
SVinfo_ff[grepl("Canvas", rownames(SVinfo_ff)),]$Application <- "Canvas"
SVinfo_ff[grepl("Manta", rownames(SVinfo_ff)),]$Application <- "Manta"
# Extract ID
SVinfo_ff$ID <- rownames(SVinfo_ff)

# table((SVinfo_ff %>% filter(FILTER == "PASS") %>% .$SOMATIC), (SVinfo_ff %>% filter(FILTER == "PASS") %>% .$Application))
######## WARNING !!!
# Do not use the "SOMATIC" field to filter CANVAS variants - they are all set to FALSE 
# as this field does not exist in the Canvas VCF entries. 
# Manta entries are all set to TRUE so the field is actually useless.

# Remove filtered entries, keep Canvas only
Canvas_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Canvas")  # 59 

# Extract chr, start, end
Canvas_ff$Chr <- sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][3]})
Canvas_ff$Start <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][4]})
Canvas_ff$End <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][5]})
Canvas_ff$Type <-  sapply(1:dim(Canvas_ff)[1], function(x){strsplit(Canvas_ff$ID[x], split = ":")[[1]][2]})

# Distribution of types
table(Canvas_ff$Type)  # Only 5 loss, and 54 REF (What is "REF"???? Non-somatic???)

# Keep LOSS/GAIN only
Canvas_ff <- Canvas_ff %>% filter(Type %in% c("LOSS", "GAIN"))  # 5


### FFPE
# Extract INFO table (all SVs)
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
Canvas_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Canvas")  # 54

# Extract chr, start, end
Canvas_ffpe$Chr <- sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][3]})
Canvas_ffpe$Start <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][4]})
Canvas_ffpe$End <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][5]})
Canvas_ffpe$Type <-  sapply(1:dim(Canvas_ffpe)[1], function(x){strsplit(Canvas_ffpe$ID[x], split = ":")[[1]][2]})

# Distribution of types
table(Canvas_ffpe$Type)  # Only 7 loss, and 47 REF (What is "REF"???? Non-somatic??? Or Ns?)

# Keep LOSS/GAIN only
Canvas_ffpe <- Canvas_ffpe %>% filter(Type %in% c("LOSS", "GAIN"))  # 7


### Compare ff and ffpe

# List SVs
Canvas_ff %>% dplyr::select(ID, Chr, Start, End, Type)
Canvas_ffpe %>% dplyr::select(ID, Chr, Start, End, Type)

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
write.table(ff_bed, file = "ff.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(ffpe_bed, file = "ffpe.bed", quote = F, row.names = F, col.names = F, sep = "\t")

# Call bedtools to get FFPE overlap with FF
system('bedtools coverage -s -a ff.bed -b ffpe.bed > ff_overlap.bed', intern = T)
ff_overlap <- read.table("ff_overlap.bed", sep = "\t")
names(ff_overlap) <- c(names(ff_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")

# Call bedtools to get FF overlap with FFPE
system('bedtools coverage -s -a ffpe.bed -b ff.bed > ffpe_overlap.bed', intern = T)
ffpe_overlap <- read.table("ffpe_overlap.bed", sep = "\t")
names(ffpe_overlap) <- c(names(ff_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")

# Summary table
result <- data.frame(PATIENT_ID = "217000028", PERCENT_RECALL_FF = sum(ff_overlap$BPoverlap) / sum(ff_overlap$BPTotal), PERCENT_RECALL_FFPE = sum(ffpe_overlap$BPoverlap) / sum(ffpe_overlap$BPTotal), BP_OVERLAP = sum(ff_overlap$BPoverlap), BP_FF_ONLY = (sum(ff_overlap$BPTotal) - sum(ff_overlap$BPoverlap)), BP_FFPE_ONLY = (sum(ffpe_overlap$BPTotal) - sum(ffpe_overlap$BPoverlap)))


# Write BAM paths into a file so I can copy the VCFs locally (not used)
#write.table(QC_portal_trios$BamPath, file = "FFPEtrio_bamPaths.txt", quote = F, row.names = F, col.names = F)





############  Explore and clean MANTA SVs ############

# Load QC_portal_trios with PCR_DUPL added
QC_portal_trios <- read.csv("QC_portal_trios_final.csv", header = T)

# Load original test VCFs
ff_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06_LP3000070-DNA_G01.somatic.SV.vcf.gz")
ffpe_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06_LP3000079-DNA_B02.somatic.SV.vcf.gz")

### FF
# Extract INFO table (all SVs)
SVinfo_ff <- as.data.frame(info(ff_vcf))
# Add filter field
SVinfo_ff$FILTER <- rowRanges(ff_vcf)$FILTER
# Create Application variable (Canvas or Manta?)
SVinfo_ff$Application <- ""
SVinfo_ff[grepl("Canvas", rownames(SVinfo_ff)),]$Application <- "Canvas"
SVinfo_ff[grepl("Manta", rownames(SVinfo_ff)),]$Application <- "Manta"
# Extract ID
SVinfo_ff$ID <- rownames(SVinfo_ff)
# Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
SVinfo_ff$START <- as.data.frame(ranges(ff_vcf))$start
SVinfo_ff$CHR <- as.character(seqnames(ff_vcf))

# Examine some variants, look at false positives
#SVinfo_ff %>% filter(FILTER == "PASS", IMPRECISE == "FALSE") %>% select(CHR, START, END, SVLEN, SVTYPE, IMPRECISE, FILTER)

# Add number of supporting paired and reads
#ff_vcf@assays$data$PR
#geno(ff_vcf)[[1]][1:5]
SVinfo_ff$PR_REF <- sapply(1:dim(SVinfo_ff)[1], function(x){
  geno(ff_vcf)[[1]][x][[1]][1]
})
SVinfo_ff$PR_ALT <- sapply(1:dim(SVinfo_ff)[1], function(x){
  geno(ff_vcf)[[1]][x][[1]][2]
})

# Add number of supporting split reads
#ff_vcf@assays$data$SR
#geno(ff_vcf)[[2]][1:5]
SVinfo_ff$SR_REF <- sapply(1:dim(SVinfo_ff)[1], function(x){
  geno(ff_vcf)[[2]][x][[1]][1]
})
SVinfo_ff$SR_ALT <- sapply(1:dim(SVinfo_ff)[1], function(x){
  geno(ff_vcf)[[2]][x][[1]][2]
})

# Add PR/SR flag (1=evidence based on PR/SR exists for somatic variant)
SVinfo_ff$PR_EV <- as.numeric(SVinfo_ff$PR_ALT > 0)
SVinfo_ff$SR_EV <- as.numeric(SVinfo_ff$SR_ALT > 0)  

# Remove filtered entries, keep Manta only (WARNING: some "good" >10kb SVs might be filtered out too, filter "MGE10kb", 85 in total)
Manta_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Manta")  # 207 

# Look at some filtered out >10kb SVs (looked at ~20 manually, none look real, most overlap with data base of known structural variants)
#SVinfo_ff %>% filter(Application == "Manta", FILTER == "MGE10kb") %>% select(CHR, START, END, IMPRECISE, SVLEN, PR_ALT, SR_ALT)

# Remove unknown CHR from the table and add chrY
#normal_chr <- c(grep("chrUn", grep("chr", levels(factor(Manta_ff$CHR)), value = T), value = T, invert = T), "chrY")  # missing chr, not clean
normal_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
Manta_ff <- Manta_ff %>% filter(CHR %in% normal_chr)  # 200

# Compare precise vs non-precise, look at range of supporting reads
table(Manta_ff$IMPRECISE, Manta_ff$SVTYPE)
sum(Manta_ff$IMPRECISE)/dim(Manta_ff)[1]  # 0.11 (percentage of imprecise)

# Filter out imprecise SVs for this purpose
Manta_ff <- Manta_ff %>% filter(IMPRECISE == FALSE)  # 178

# Look at range of supporting reads
table( Manta_ff$PR_ALT <3, Manta_ff$SVTYPE)
table( Manta_ff$SR_ALT <3, Manta_ff$SVTYPE)
table( ((Manta_ff$PR_ALT <3) & (Manta_ff$SR_ALT <3)), Manta_ff$SVTYPE)  # 2 BND with < 3 both SR and PR

# Filter out SVs with <3 supporting PR AND SR reads
Manta_ff <- Manta_ff[!((Manta_ff$PR_ALT <3) & (Manta_ff$SR_ALT <3)),]  # 176






### FFPE
# Extract INFO table (all SVs)
SVinfo_ffpe <- as.data.frame(info(ffpe_vcf))
# Add filter field
SVinfo_ffpe$FILTER <- rowRanges(ffpe_vcf)$FILTER
# Create Application variable (Canvas or Manta?)
SVinfo_ffpe$Application <- ""
SVinfo_ffpe[grepl("Canvas", rownames(SVinfo_ffpe)),]$Application <- "Canvas"
SVinfo_ffpe[grepl("Manta", rownames(SVinfo_ffpe)),]$Application <- "Manta"

# Extract ID
SVinfo_ffpe$ID <- rownames(SVinfo_ffpe)

# Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
SVinfo_ffpe$START <- as.data.frame(ranges(ffpe_vcf))$start
SVinfo_ffpe$CHR <- as.character(seqnames(ffpe_vcf))

# Add number of supporting paired and reads
SVinfo_ffpe$PR_REF <- sapply(1:dim(SVinfo_ffpe)[1], function(x){
  geno(ffpe_vcf)[[1]][x][[1]][1]
})
SVinfo_ffpe$PR_ALT <- sapply(1:dim(SVinfo_ffpe)[1], function(x){
  geno(ffpe_vcf)[[1]][x][[1]][2]
})

# Add number of supporting split reads
SVinfo_ffpe$SR_REF <- sapply(1:dim(SVinfo_ffpe)[1], function(x){
  geno(ffpe_vcf)[[2]][x][[1]][1]
})
SVinfo_ffpe$SR_ALT <- sapply(1:dim(SVinfo_ffpe)[1], function(x){
  geno(ffpe_vcf)[[2]][x][[1]][2]
})

# Add PR/SR flag (1=evidence based on PR/SR exists for somatic variant)
SVinfo_ffpe$PR_EV <- as.numeric(SVinfo_ffpe$PR_ALT > 0)
SVinfo_ffpe$SR_EV <- as.numeric(SVinfo_ffpe$SR_ALT > 0)  

# Remove filtered entries, keep Manta only
Manta_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Manta")  # 1239

# Remove unknown CHR from the table and add chrY
Manta_ffpe <- Manta_ffpe %>% filter(CHR %in% normal_chr)  # 1221

# Compare precise vs non-precise, look at range of supporting reads
table(Manta_ffpe$IMPRECISE, Manta_ffpe$SVTYPE)
sum(Manta_ffpe$IMPRECISE)/dim(Manta_ffpe)[1]  # 0.02457002 (percentage of imprecise)

# Filter out imprecise SVs for this purpose
Manta_ffpe <- Manta_ffpe %>% filter(IMPRECISE == FALSE)  # 1191

# Look at range of supporting reads
table( Manta_ffpe$PR_ALT <3, Manta_ffpe$SVTYPE)
table( Manta_ffpe$SR_ALT <3, Manta_ffpe$SVTYPE)
table( ((Manta_ffpe$PR_ALT <3) & (Manta_ffpe$SR_ALT <3)), Manta_ffpe$SVTYPE)  # 2 BND with < 3 both SR and PR

# Filter out SVs with <3 supporting PR AND SR reads
Manta_ffpe <- Manta_ffpe[!((Manta_ffpe$PR_ALT <3) & (Manta_ffpe$SR_ALT <3)),]  # 806



############  Plots comparing FF and FFPE read support ############  

# Plot supporting reads per SV type (FF)
ggplot(Manta_ff, aes(x=factor(SVTYPE), y=PR_ALT)) + geom_boxplot() + geom_jitter(aes(colour = factor(SVTYPE))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads (FF)")
ggplot(Manta_ff, aes(x=factor(SVTYPE), y=SR_ALT)) + geom_boxplot() + geom_jitter(aes(colour = factor(SVTYPE))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads (FF)")
ggplot(Manta_ff, aes(x=PR_ALT, y=SR_ALT, colour=SVTYPE)) + geom_jitter() + blank + regr_line + ggtitle("Split vs Paired Reads (FF)")
cor(Manta_ff$PR_ALT, Manta_ff$SR_ALT, method = "pearson")  # 0.5392241

# Plot supporting reads per SV type (FFPE)
ggplot(Manta_ffpe, aes(x=factor(SVTYPE), y=PR_ALT)) + geom_boxplot() + geom_jitter(aes(colour = factor(SVTYPE))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads (FFPE)")
ggplot(Manta_ffpe, aes(x=factor(SVTYPE), y=SR_ALT)) + geom_boxplot() + geom_jitter(aes(colour = factor(SVTYPE))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads (FFPE)")
ggplot(Manta_ffpe, aes(x=PR_ALT, y=SR_ALT, colour=SVTYPE)) + geom_jitter() + blank + regr_line + ggtitle("Split vs Paired Reads (FFPE)")
cor(Manta_ffpe$PR_ALT, Manta_ffpe$SR_ALT, method = "pearson")  # 0.3317059

# FF vs FFPE
Manta_ff$FF <- "FF"
Manta_ffpe$FF <- "FFPE"
ff_ffpe_merged <- rbind(Manta_ff, Manta_ffpe)
# PR
ggplot(ff_ffpe_merged, aes(x=factor(SVTYPE), y=PR_ALT)) + geom_jitter(aes(colour = factor(FF)), alpha = 0.6) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads")
ggplot(ff_ffpe_merged, aes(x=factor(SVTYPE), y=PR_ALT)) + geom_boxplot(aes(colour = factor(FF))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads")
# SR
ggplot(ff_ffpe_merged, aes(x=factor(SVTYPE), y=SR_ALT)) + geom_jitter(aes(colour = factor(FF)), alpha = 0.6) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads")
ggplot(ff_ffpe_merged, aes(x=factor(SVTYPE), y=SR_ALT)) + geom_boxplot(aes(colour = factor(FF))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads")

# Only concordant reads
ggplot(ff_ffpe_merged[ff_ffpe_merged$CONCORDANT ==1,], aes(x=factor(SVTYPE), y=PR_ALT)) + geom_jitter(aes(colour = factor(FF)), alpha = 0.6) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads")
ggplot(ff_ffpe_merged[ff_ffpe_merged$CONCORDANT ==1,], aes(x=factor(SVTYPE), y=PR_ALT)) + geom_boxplot(aes(colour = factor(FF))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads")
ggplot(ff_ffpe_merged[ff_ffpe_merged$CONCORDANT ==1,], aes(x=factor(SVTYPE), y=SR_ALT)) + geom_jitter(aes(colour = factor(FF)), alpha = 0.6) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads")
ggplot(ff_ffpe_merged[ff_ffpe_merged$CONCORDANT ==1,], aes(x=factor(SVTYPE), y=SR_ALT)) + geom_boxplot(aes(colour = factor(FF))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads")


############  Compare non-filtered MANTA SVs ############

### This is the direct comparison (exact breakpoint)

# Create CHR_START_END_TYPE key and look for dups
ff_ffpe_merged$KEY <- sapply(1:dim(ff_ffpe_merged)[1], function(x){
  paste(ff_ffpe_merged$CHR[x], ff_ffpe_merged$START[x], ff_ffpe_merged$END[x], ff_ffpe_merged$SVTYPE[x], sep = "-")
})
sum(duplicated(ff_ffpe_merged$KEY))  # 54 exact same events

# Check that dups are among different samples
sum(duplicated(ff_ffpe_merged[ff_ffpe_merged$FF == "FF",]$KEY))  # 0
sum(duplicated(ff_ffpe_merged[ff_ffpe_merged$FF == "FFPE",]$KEY))  # 0

# Add CONCORDANT flag to FF and FFPE calls
concordant_keys <- ff_ffpe_merged[duplicated(ff_ffpe_merged$KEY),]$KEY
ff_ffpe_merged$CONCORDANT <- 0
ff_ffpe_merged[ff_ffpe_merged$KEY %in% concordant_keys,]$CONCORDANT <- 1

# Concordant calls by SV type
table(ff_ffpe_merged$CONCORDANT, ff_ffpe_merged$SVTYPE)[2,]/2

# Total SVs in FF and FFPE
table(ff_ffpe_merged$FF, ff_ffpe_merged$SVTYPE)

# FFPE recall and precision 
(table(ff_ffpe_merged$CONCORDANT, ff_ffpe_merged$SVTYPE)[2,]/2) / table(ff_ffpe_merged$FF, ff_ffpe_merged$SVTYPE)[1,] # Recall
(table(ff_ffpe_merged$CONCORDANT, ff_ffpe_merged$SVTYPE)[2,]/2) / table(ff_ffpe_merged$FF, ff_ffpe_merged$SVTYPE)[2,] # Precision



############  Compare filtered MANTA SVs ############

### Prepare filtering (on HPC)

# Convert WindowMasker to BED format
wMasker <- read.table("/home/mmijuskovic/FFPE/windowmaskerSdust.hg38.txt", header = F)
wMasker <- wMasker %>% select(-(V1))
write.table(wMasker, file = "windowmaskerSdust.hg38.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Simple repeats
repeats <- read.table("/home/mmijuskovic/FFPE/simpleRepeat.hg38.txt", header = F)
repeats <- repeats %>% select(-(V1))
repeats <- repeats %>% select(1:4)
write.table(repeats, file = "/home/mmijuskovic/FFPE/simpleRepeat.hg38.bed", row.names = F, col.names = F, quote = F, sep = "\t")

# Segmental duplications
seg_dups <- read.table("/home/mmijuskovic/FFPE/genomicSuperDups.hg38.txt", header = F)
seg_dups <- seg_dups %>% select(-(V1))
seg_dups <- seg_dups %>% select(1:4)
write.table(seg_dups, file = "/home/mmijuskovic/FFPE/genomicSuperDups.hg38.bed", row.names = F, col.names = F,  quote = F, sep = "\t")



# Add new type coding ("+" is FF, "-" is FFPE)
ff_ffpe_merged$Type2 <- NA
ff_ffpe_merged$Type2 <- sapply(1:dim(ff_ffpe_merged)[1], function(x){
  if (ff_ffpe_merged$FF[x] == "FF") {ff_ffpe_merged$Type2[x] <- "+"}
  else if (ff_ffpe_merged$FF[x] == "FFPE") {ff_ffpe_merged$Type2[x] <- "-"}
})

# Create a bed file with SVs, start and end separately, remove NAs, code FF and FFPE as strand ("+" for FF, "-" for FFPE)
# Note that START has to be adjusted to 0-based (-1) and end is not included (stays same)
# 170216 amended to add 150 bp on either side of the breakpoint
start_bed <- cbind((ff_ffpe_merged %>% dplyr::select(CHR, START)), (ff_ffpe_merged %>% dplyr::select(START, KEY, Type2)))
names(start_bed)[3] <- "end"
start_bed$Score <- ""
start_bed <- start_bed %>% select(CHR, START, end, KEY, Score, Type2)
start_bed <- start_bed %>% filter(!is.na(START))
#start_bed$START <- start_bed$START-1
start_bed$START <- start_bed$START - 151
start_bed$end <- start_bed$end + 150
write.table(start_bed, file = "sv_start.bed", quote = F, row.names = F, col.names = F, sep = "\t")

end_bed <- cbind((ff_ffpe_merged %>% dplyr::select(CHR, END)), (ff_ffpe_merged %>% dplyr::select(END, KEY, Type2)))
names(end_bed)[2] <- "start"
end_bed$Score <- ""
end_bed <- end_bed %>% select(CHR, start, END, KEY, Score, Type2)
end_bed <- end_bed %>% filter(!is.na(END))
#end_bed$start <- end_bed$start-1
end_bed$start <- end_bed$start - 151
end_bed$END <- end_bed$END + 150
write.table(end_bed, file = "sv_end.bed", quote = F, row.names = F, col.names = F, sep = "\t")


### Call bedtools to find overlaps with WindowMasker 
# start
system('bedtools coverage -a sv_start.bed -b windowmaskerSdust.hg38.bed > sv_wMasker_start_overlap.bed', intern = T)
sv_wMasker_start <- read.table("sv_wMasker_start_overlap.bed", sep = "\t")
names(sv_wMasker_start) <- c(names(start_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system('bedtools coverage -a sv_end.bed -b windowmaskerSdust.hg38.bed > sv_wMasker_end_overlap.bed', intern = T)
sv_wMasker_end <- read.table("sv_wMasker_end_overlap.bed", sep = "\t")
names(sv_wMasker_end) <- c(names(end_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")

# Collect KEYs with wMasker overlap
wMasker_keys <- unique(c(as.character(sv_wMasker_start %>% filter(NumOverlap != 0) %>% .$KEY), as.character(sv_wMasker_end %>% filter(NumOverlap != 0) %>% .$KEY)))
length(wMasker_keys) # 916 of 982 filtered out

# Filter out SVs where START or END overlaps with WindowMasker
ff_ffpe_merged$wMasker_filtered <- 0
ff_ffpe_merged[(ff_ffpe_merged$KEY %in% wMasker_keys),]$wMasker_filtered <- 1


### Call bedtools to find overlaps with simple repeats

# start
system('bedtools coverage -a sv_start.bed -b simpleRepeat.hg38.bed > sv_repeats_start_overlap.bed', intern = T)
sv_repeats_start <- read.table("sv_repeats_start_overlap.bed", sep = "\t")
names(sv_repeats_start) <- c(names(start_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system('bedtools coverage -a sv_end.bed -b simpleRepeat.hg38.bed > sv_repeats_end_overlap.bed', intern = T)
sv_repeats_end <- read.table("sv_repeats_end_overlap.bed", sep = "\t")
names(sv_repeats_end) <- c(names(end_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")

# Collect KEYs with overlap
repeats_keys <- unique(c(as.character(sv_repeats_start %>% filter(NumOverlap != 0) %>% .$KEY), as.character(sv_repeats_end %>% filter(NumOverlap != 0) %>% .$KEY)))
length(repeats_keys) # 186 of 982 filtered out

# Filter out SVs where START or END overlaps with repeats
ff_ffpe_merged$repeats_filtered <- 0
ff_ffpe_merged[(ff_ffpe_merged$KEY %in% repeats_keys),]$repeats_filtered <- 1



### Call bedtools to find overlaps with segmental duplications

# start
system('bedtools coverage -a sv_start.bed -b genomicSuperDups.hg38.bed > sv_segdups_start_overlap.bed', intern = T)
sv_segdups_start <- read.table("sv_segdups_start_overlap.bed", sep = "\t")
names(sv_segdups_start) <- c(names(start_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# end
system('bedtools coverage -a sv_end.bed -b genomicSuperDups.hg38.bed > sv_segdups_end_overlap.bed', intern = T)
sv_segdups_end <- read.table("sv_segdups_end_overlap.bed", sep = "\t")
names(sv_segdups_end) <- c(names(end_bed), "NumOverlap", "BPoverlap", "BPTotal", "PCT")

# Collect KEYs with overlap
segdups_keys <- unique(c(as.character(sv_segdups_start %>% filter(NumOverlap != 0) %>% .$KEY), as.character(sv_segdups_end %>% filter(NumOverlap != 0) %>% .$KEY)))
length(segdups_keys) # 45 of 982 filtered out

# Filter out SVs where START or END overlaps with repeats
ff_ffpe_merged$segdups_filtered <- 0
ff_ffpe_merged[(ff_ffpe_merged$KEY %in% segdups_keys),]$segdups_filtered <- 1 



### Apply all filtering

# Check the distribution of filters
table(ff_ffpe_merged$wMasker_filtered, ff_ffpe_merged$repeats_filtered, ff_ffpe_merged$segdups_filtered)

# Add a general FILTERED field
ff_ffpe_merged$FILTERED <- 0
filtered_id <- ff_ffpe_merged %>% filter((wMasker_filtered == 1) | (repeats_filtered == 1) | (segdups_filtered == 1)) %>% .$ID
ff_ffpe_merged[ff_ffpe_merged$ID %in% filtered_id,]$FILTERED <- 1

# Flag breakends whose mates are filtered out
ff_ffpe_merged$MATEID <- as.character(ff_ffpe_merged$MATEID)
ff_ffpe_merged$BND_MATE_FILTERED <- NA
ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",]$BND_MATE_FILTERED <- sapply(1:dim(ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",])[1], function(x){
  if (ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",]$MATEID[x] %in% filtered_id) { return(1)}
  else {return(0)}
})

# Remove filtered reads
ff_ffpe_merged_fil <- ff_ffpe_merged %>% filter(wMasker_filtered == 0, repeats_filtered == 0, segdups_filtered == 0) %>% filter(is.na(BND_MATE_FILTERED) | BND_MATE_FILTERED == 0)




### Summary SV comparison

# Overview of leftover high quality SV candidates
table(ff_ffpe_merged_fil$FF, ff_ffpe_merged_fil$SVTYPE)

### Comparison of filtered SV candidates
# Concordant calls by SV type (NONE left for this sample)
table(ff_ffpe_merged_fil$CONCORDANT, ff_ffpe_merged_fil$SVTYPE)[2,]/2

# FFPE recall and precision (all zero for this sample, no concordant SVs)
(table(ff_ffpe_merged_fil$CONCORDANT, ff_ffpe_merged_fil$SVTYPE)[2,]/2) / table(ff_ffpe_merged_fil$FF, ff_ffpe_merged_fil$SVTYPE)[1,] # Recall
(table(ff_ffpe_merged_fil$CONCORDANT, ff_ffpe_merged_fil$SVTYPE)[2,]/2) / table(ff_ffpe_merged_fil$FF, ff_ffpe_merged_fil$SVTYPE)[2,] # Precision

# Check concordant filtered calls
ff_ffpe_merged_fil %>% filter(CONCORDANT == 1) %>% select(KEY, FF, PR_ALT, SR_ALT, BND_DEPTH, MATE_BND_DEPTH, SVLEN)

# Check all filtered SVs
ff_ffpe_merged_fil %>% select(KEY, FF, PR_ALT, SR_ALT, BND_DEPTH, MATE_BND_DEPTH, SVLEN)







################# Check germline VCF ################# 

# Read germline VCF
gl_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06.SV.vcf.gz")

# Extract INFO table (all SVs)
SVinfo_gl <- as.data.frame(info(gl_vcf))
# Add filter field
SVinfo_gl$FILTER <- rowRanges(gl_vcf)$FILTER
# Create Application variable (Canvas or Manta?)
SVinfo_gl$Application <- ""
SVinfo_gl[grepl("Canvas", rownames(SVinfo_gl)),]$Application <- "Canvas"
SVinfo_gl[grepl("Manta", rownames(SVinfo_gl)),]$Application <- "Manta"
# Extract ID
SVinfo_gl$ID <- rownames(SVinfo_gl)

# Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
SVinfo_gl$START <- as.data.frame(ranges(gl_vcf))$start
SVinfo_gl$CHR <- as.character(seqnames(gl_vcf))

# Add number of supporting paired and reads
# gl_vcf@assays$data$PR
# geno(gl_vcf)[[5]][1:5]  # PR
# geno(gl_vcf)[[6]][1:5]  # SR
SVinfo_gl$PR_REF <- sapply(1:dim(SVinfo_gl)[1], function(x){
  geno(gl_vcf)[[5]][x][[1]][1]
})
SVinfo_gl$PR_ALT <- sapply(1:dim(SVinfo_gl)[1], function(x){
  geno(gl_vcf)[[5]][x][[1]][2]
})

# Add number of supporting split reads
SVinfo_gl$SR_REF <- sapply(1:dim(SVinfo_gl)[1], function(x){
  geno(gl_vcf)[[6]][x][[1]][1]
})
SVinfo_gl$SR_ALT <- sapply(1:dim(SVinfo_gl)[1], function(x){
  geno(gl_vcf)[[6]][x][[1]][2]
})

# Add PR/SR flag (1=evidence based on PR/SR exists for somatic variant)
SVinfo_gl$PR_EV <- as.numeric(SVinfo_gl$PR_ALT > 0)
SVinfo_gl$SR_EV <- as.numeric(SVinfo_gl$SR_ALT > 0)  

# Keep Manta only (no filtering!)
Manta_gl <- SVinfo_gl %>% filter(FILTER == "PASS", Application == "Manta")  #  11143 

# Create KEY to look at precise events 
# Create CHR_START_END_TYPE key and look for dups
Manta_gl$KEY <- sapply(1:dim(Manta_gl)[1], function(x){
  paste(Manta_gl$CHR[x], Manta_gl$START[x], Manta_gl$END[x], Manta_gl$SVTYPE[x], sep = "-")
})
sum(duplicated(Manta_gl$KEY))  # 23 exact same events (WHERE do these come from?!)

# Look for SVs found in FF and FFPE
concordant_keys2 <- ff_ffpe_merged_fil %>% filter(CONCORDANT == 1) %>% .$KEY
Manta_gl %>% filter(KEY %in% concordant_keys2)  # none there
# Look for Chr11 breakends
Manta_gl %>% filter(CHR == "chr11", SVTYPE == "BND")



# # Remove unknown CHR from the table and add chrY
# Manta_gl <- Manta_gl %>% filter(CHR %in% normal_chr)  # 
# 
# # Compare precise vs non-precise, look at range of supporting reads
# table(Manta_gl$IMPRECISE, Manta_gl$SVTYPE)
# sum(Manta_gl$IMPRECISE)/dim(Manta_gl)[1]  # 0.11 (percentage of imprecise)
# 
# # Filter out imprecise SVs for this purpose
# Manta_gl <- Manta_gl %>% filter(IMPRECISE == FALSE)  # 
# 
# # Look at range of supporting reads
# table( Manta_gl$PR_ALT <3, Manta_gl$SVTYPE)
# table( Manta_gl$SR_ALT <3, Manta_gl$SVTYPE)
# table( ((Manta_gl$PR_ALT <3) & (Manta_gl$SR_ALT <3)), Manta_gl$SVTYPE)  # 2 BND with < 3 both SR and PR
# 
# # Filter out SVs with <3 supporting PR AND SR reads
# Manta_gl <- Manta_gl[!((Manta_gl$PR_ALT <3) & (Manta_gl$SR_ALT <3)),]  # 


# Function to extract SV ID given the SV KEY
keyToID <- function(x){
  key <- x
  ff_ffpe_merged %>% filter(KEY == key) %>% .$ID
}



