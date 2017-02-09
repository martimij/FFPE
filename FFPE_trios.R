# Martina Mijuskovic
# FFPE project
# Collect QC and clinical data for sequenced FFPE trios
# Develop the comparison method for Canvas CNV output between FF and FFPE samples

library("data.table")
library("dplyr")
library("ggvis")
library("ggplot2")

source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
library(VariantAnnotation)

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

# Manually add the BAM path for one sample (FFPE from 217000030, LP3000079-DNA_H01) - passes QC in v2 but not in v4, we are keeping it
#upload %>% filter(Platekey == "LP3000079-DNA_H01") # Chose v4 path from unfiltered upload report
QC_portal_trios[QC_portal_trios$SAMPLE_WELL_ID == "LP3000079-DNA_H01",]$BamPath <- "/genomes/by_date/2016-12-13/CANCP40747/CancerLP3000079-DNA_H01_NormalLP3000067-DNA_H08"

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


############  Explore SVs from VCF ############  

# Load trios QC data
QC_portal_trios <- read.csv("QC_portal_trios.csv")

ff_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06_LP3000070-DNA_G01.somatic.SV.vcf.gz")
ffpe_vcf <- readVcf(file = "/Users/MartinaMijuskovic/FFPE/Trio_VCFs/LP3000067-DNA_E06_LP3000079-DNA_B02.somatic.SV.vcf.gz")

# Filter data
as.data.frame(table(ff_vcf@fixed@listData$FILTER))
as.data.frame(table(ffpe_vcf@fixed@listData$FILTER))

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





############  Explore MANTA SVs ############

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

# Remove filtered entries, keep Manta only
Manta_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Manta")  # 207 

# Add chr, start, end positions
as.data.frame(ranges(ff_vcf))  # gets start, end (by width), IDs (names)
as.data.frame(seqnames(ff_vcf))  # gets chr name










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

# Remove filtered entries, keep Manta only
Manta_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Manta")  # 1239



### Plots to explore Manta SVs

# Variant qualities by SV type
Manta_ff %>% group_by(SVTYPE) %>% ggvis(~SOMATICSCORE, fill = ~SVTYPE) %>% layer_densities()


