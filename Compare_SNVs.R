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


############# Calculate concordance ############# 

# Load the manifest (HPC)
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/CNV_trio_comparison/QC_portal_trios.csv")  # original 26 trios
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
system("touch /home/mmijuskovic/FFPE/FFPE_trios_timestamp -d 2017-03-29") # timestamp to help select only new tiering json files (produced after March 29th 2017)
QC_portal_trios$json_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", sub("/genomes", "/genomes/analysis", x), "-newer /home/mmijuskovic/FFPE/FFPE_trios_timestamp -iname *tiering.json", sep = " ")
  system(command, intern = T)
}))


### Function that extracts variant information from json files and calculates concordance between FF and FFPE
compareSNV <- function(patientID, var_freq = NULL){
  
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
  # sum(duplicated(ff))  #33  # after OpenCGA fix -> 1
  # sum(duplicated(ff$KEY))  #37  # after OpenCGA fix -> 1
  # # Checking the leftover duplicate
  # dup_key <- ff[duplicated(ff$KEY),]$KEY
  # ff %>% filter(KEY == dup_key)

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
  
  # Subset the variants by frequency if frequency argument is provided
  if (!is.null(var_freq)) {
    ff <- ff %>% filter(VAF > var_freq)
    ffpe <- ffpe %>% filter(VAF > var_freq)
  }
  
  # Summary table with concordance
  result <- data.frame(FF_TOTAL = dim(ff)[1],
                       FFPE_TOTAL = dim(ffpe)[1],
                       OVERLAP = sum(ff$KEY %in% ffpe$KEY), # bug
                       FF_UNIQ = ((dim(ff)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       FFPE_UNIQ = ((dim(ffpe)[1]) - sum(ff$KEY %in% ffpe$KEY)),
                       RECALL = ((sum(ff$KEY %in% ffpe$KEY))/(dim(ff)[1])),
                       PRECISION = (sum(ff$KEY %in% ffpe$KEY))/(dim(ffpe)[1]))

 #write.table(result, file = paste0(patientID, "_SNV_concord", ".tsv"), row.names = F, quote = F, sep = "\t")
 return(result)
}

# Get all patient IDs
patientIDs <- as.character(unique(QC_portal_trios$PATIENT_ID))

# # Remove 2 patient IDs for which I'm still missing FFPE json files
# rmIDs <- QC_portal_trios %>% filter(SAMPLE_WELL_ID %in% c("LP3000079-DNA_G02", "LP3000074-DNA_C07")) %>% .$PATIENT_ID
# rmIDs <- c(rmIDs, "218000014")
# patientIDs <- patientIDs[!patientIDs %in% rmIDs]

# Run CNV comparison for each patient ID (all allele frequencies)
SNV_summary <- lapply(patientIDs, compareSNV)  # errors :/ " no applicable method for 'select_' applied to an object of class "list""
#SNV_summary <- lapply( c("212000006", "212000008", "212000009", "212000014", "212000015", "217000011", "217000028", "217000030", "217000052", "217000073", "218000002", "218000013", "218000014", "218000016", "218000017", "218000018", "218000021", "218000022", "218000024", "218000027", "218000032", "218000033", "218000036", "218000039"), compareSNV)
# there must be an error or malformed json somewhere (218000014 ID is problematic, malformed FFPE json file:
# "/genomes/analysis/by_date/2017-01-16/CANCP41020/CancerLP3000162-DNA_E01_NormalLP3000067-DNA_B09/SomaticVariations/1486209216/LP3000162-DNA_E01_tiering.json")
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs

# Write out the resulting table
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_allFreq", today, ".csv"), row.names = F, quote = F)

# Re-run using variants with >5% frequency only
SNV_summary <- lapply(patientIDs, compareSNV, 0.05)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_5pctFreq", today, ".csv"), row.names = F, quote = F)

# Re-run using variants with >20% frequency only
SNV_summary <- lapply(patientIDs, compareSNV, 0.2)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_20pctFreq", today, ".csv"), row.names = F, quote = F)

# Re-run with other frequencies (to get 0, 5, 10, 15, 20, 25, and 30)

# 10
SNV_summary <- lapply(patientIDs, compareSNV, 0.1)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_10pctFreq", today, ".csv"), row.names = F, quote = F)

# 15
SNV_summary <- lapply(patientIDs, compareSNV, 0.15)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_15pctFreq", today, ".csv"), row.names = F, quote = F)

# 25
SNV_summary <- lapply(patientIDs, compareSNV, 0.25)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_25pctFreq", today, ".csv"), row.names = F, quote = F)


# 30
SNV_summary <- lapply(patientIDs, compareSNV, 0.3)
SNV_summary <- bind_rows(SNV_summary)
SNV_summary$PATIENT_ID <- patientIDs
write.csv(SNV_summary, file = paste0("SNV_summary_26trios_30pctFreq", today, ".csv"), row.names = F, quote = F)




################ SNV summary plots (original 26 trios) ################ 

# Load the data (>5% variant freq only) - local
SNV_summary <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_5pctFreq2017-04-19.csv")

# Add QC details to the SNV summary table
SNV_summary$PATIENT_ID <- as.character(SNV_summary$PATIENT_ID)
QC_table <- QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE") %>% dplyr::select(PATIENT_ID, CENTER_CODE, TumorType, SAMPLE_WELL_ID, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, DEAMINATION_MISMATCHES_PER)
names(QC_table)[4:12] <- paste0("FFPE_", names(QC_table)[4:12])
QC_table <- left_join(QC_table, (QC_portal_trios %>% filter(SAMPLE_TYPE == "FF") %>% dplyr::select(PATIENT_ID, SAMPLE_WELL_ID, LIBRARY_TYPE, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, DEAMINATION_MISMATCHES_PER)), by = "PATIENT_ID")
names(QC_table)[13:22] <- paste0("FF_", names(QC_table)[13:22])
QC_table$PATIENT_ID <- as.character(QC_table$PATIENT_ID)
SNV_summary <- left_join(SNV_summary, QC_table, by = "PATIENT_ID")

# Load the Domain1, Domain2 gene lists
domain1 <- read.table("GENOMONCOLOGY_SOLID_TUMOUR.v1.4.tsv", sep = "\t", header = T)
domain1 <- domain1[,1:4]
domain1 <- unique(domain1)
length(unique(domain1$gene_name))  # 72 genes
domain2 <-  read.table("CANCER_CENSUS_GENES.v1.4.tsv", sep = "\t", header = T)
length(unique(domain2$gene_name))  # 544 genes

### Summary plots

# Distribution of RECALL and PRECISION
hist(SNV_summary$RECALL)
hist(SNV_summary$PRECISION)

# Overlap of Domain 1-3 SNVs >5% frequency per patient
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SNV_summary_m <- as.data.frame(t(SNV_summary %>% arrange(PRECISION) %>% dplyr::select(PATIENT_ID, OVERLAP, FF_UNIQ, FFPE_UNIQ)))
names(SNV_summary_m) <- as.matrix(SNV_summary_m[1,])
SNV_summary_m <- SNV_summary_m[2:4,]
SNV_summary_m <- melt(cbind(SNV_summary_m, ind = rownames(SNV_summary_m)), id.vars = c('ind'))
SNV_summary_m$value <- as.numeric(SNV_summary_m$value)
# Overap plot (needs package "scales")
ggplot(SNV_summary_m, aes(x = variable, y = value, fill = ind)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + blank


# RECALL of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary, aes(x = RECALL, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_AT_DROP, method = "spearman")  # 0.01094204

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.06393162

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # -0.08038311

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.07213675

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # -0.1070085

ggplot(SNV_summary, aes(x = RECALL, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SNV_summary$RECALL, SNV_summary$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.3930165


# PRECISION of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_AT_DROP, method = "spearman")  # -0.3135579

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.3996581

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # -0.04241492

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.2205128

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary$PRECISION, SNV_summary$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # 0.04820513

ggplot(SNV_summary, aes(x = PRECISION, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Precision", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SNV_summary$PRECISION, SNV_summary$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.4334136


######## Cluster samples 

summary(SNV_summary$PRECISION)
ggplot(SNV_summary, aes(x = CENTER_CODE, y = PRECISION, col = factor(CENTER_CODE))) + geom_jitter() + geom_boxplot(alpha = 0) + theme(legend.title=element_blank())
summary(SNV_summary$RECALL)
ggplot(SNV_summary, aes(x = CENTER_CODE, y = RECALL, col = factor(CENTER_CODE))) + geom_jitter() + geom_boxplot(alpha = 0) + theme(legend.title=element_blank())


# Recall of FF vs recall of FFPE (recall vs precision)
#ggplot(SNV_summary, aes(x = RECALL, y = PRECISION, col = factor(CENTER_CODE))) + geom_jitter() + geom_smooth(method = "lm", se = F) + labs(x = "Percent Recall of FF", y = "Percent Recall of FFPE") 
ggplot(SNV_summary, aes(x = RECALL, y = PRECISION, col = factor(CENTER_CODE))) + geom_jitter() + labs(x = "Recall", y = "Precision")  + theme(legend.title=element_blank())
cor(SNV_summary$RECALL, SNV_summary$PRECISION, method = "spearman")  # 0.6690598

summary(SNV_summary$FFPE_TUMOUR_PURITY)
summary(SNV_summary$FF_TUMOUR_PURITY)



################ Overlap plots with different VAFs ################ 

# Load data, add to QC data
### 20% VAF
SNV_summary <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_20pctFreq2017-04-20.csv")

# Add QC details to the SNV summary table
SNV_summary$PATIENT_ID <- as.character(SNV_summary$PATIENT_ID)
QC_table <- QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE") %>% dplyr::select(PATIENT_ID, CENTER_CODE, TumorType, SAMPLE_WELL_ID, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, DEAMINATION_MISMATCHES_PER)
names(QC_table)[4:12] <- paste0("FFPE_", names(QC_table)[4:12])
QC_table <- left_join(QC_table, (QC_portal_trios %>% filter(SAMPLE_TYPE == "FF") %>% dplyr::select(PATIENT_ID, SAMPLE_WELL_ID, LIBRARY_TYPE, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, DEAMINATION_MISMATCHES_PER)), by = "PATIENT_ID")
names(QC_table)[13:22] <- paste0("FF_", names(QC_table)[13:22])
QC_table$PATIENT_ID <- as.character(QC_table$PATIENT_ID)
SNV_summary <- left_join(SNV_summary, QC_table, by = "PATIENT_ID")

# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SNV_summary_m <- as.data.frame(t(SNV_summary %>% arrange(PRECISION) %>% dplyr::select(PATIENT_ID, OVERLAP, FF_UNIQ, FFPE_UNIQ)))
names(SNV_summary_m) <- as.matrix(SNV_summary_m[1,])
SNV_summary_m <- SNV_summary_m[2:4,]
SNV_summary_m <- melt(cbind(SNV_summary_m, ind = rownames(SNV_summary_m)), id.vars = c('ind'))
SNV_summary_m$value <- as.numeric(SNV_summary_m$value)
# Overap plot (needs package "scales")
ggplot(SNV_summary_m, aes(x = variable, y = value, fill = ind)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + blank


######### Correlation of overlap with tumor purity, >20% VAF ######### 

# Outliers in the overlap plot
SNV_summary %>% filter(PATIENT_ID %in% c("217000052", "217000011")) %>% dplyr::select(PATIENT_ID, TumorType, RECALL, PRECISION, FFPE_TUMOUR_PURITY, FF_TUMOUR_PURITY, FFPE_GC_DROP, FFPE_COVERAGE_HOMOGENEITY, FFPE_DEAMINATION_MISMATCHES_PER, FF_LIBRARY_TYPE)

# Samples with FF and FFPE tumor purity > 50%
dim(SNV_summary %>% filter(FFPE_TUMOUR_PURITY > 50, FF_TUMOUR_PURITY > 50))  # 8 (FF and FFPE)
dim(SNV_summary %>% filter(FFPE_TUMOUR_PURITY > 50))  # 10 (FFPE only)

# Distribution of tumor purity
ggplot(QC_portal_trios, aes(x = SAMPLE_TYPE, y = TUMOUR_PURITY, col = factor(SAMPLE_TYPE))) + geom_jitter() + geom_boxplot(alpha = 0) + theme(legend.title=element_blank())



#########  Only samples with FFPE tumor purity >50% ######### 


# RECALL by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = RECALL, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Recall", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$RECALL, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_AT_DROP, method = "spearman")  #  0.05454545

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = RECALL, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Recall", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$RECALL, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.1757576

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = RECALL, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Recall", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$RECALL, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_CHIMERIC_PER, method = "spearman")  # -0.1515152

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = RECALL, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Recall", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$RECALL, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.4424242

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = RECALL, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Recall", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$RECALL, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # 0.4181818


# PRECISION by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, etc
ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = PRECISION, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$PRECISION, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_AT_DROP, method = "spearman")  #  -0.5757576

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = PRECISION, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$PRECISION, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.5636364

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = PRECISION, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$PRECISION, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_CHIMERIC_PER, method = "spearman")  # -0.4424242

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = PRECISION, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$PRECISION, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.1636364

ggplot(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,], aes(x = PRECISION, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$PRECISION, SNV_summary[SNV_summary$FFPE_TUMOUR_PURITY > 50,]$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # -0.3333333



########## Overlap plots by allele frequency ########## 

# Read SNV overlap data for difference allele frequency thresholds
SNV_summary_0 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_allFreq2017-04-19.csv")
SNV_summary_5 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_5pctFreq2017-04-19.csv")
SNV_summary_10 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_10pctFreq2017-04-24.csv")
SNV_summary_15 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_15pctFreq2017-04-24.csv")
SNV_summary_20 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_20pctFreq2017-04-20.csv")
SNV_summary_25 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_25pctFreq2017-04-24.csv")
SNV_summary_30 <- read.csv("/Users/MartinaMijuskovic/FFPE/SNV_summary_26trios_30pctFreq2017-04-24.csv")

# Create a table with overlap % for each freq, each patient
SNV_summary_0$OVERLAP_PCT <- SNV_summary_0$OVERLAP/(SNV_summary_0$FF_UNIQ + SNV_summary_0$FFPE_UNIQ + SNV_summary_0$OVERLAP)
SNV_summary_5$OVERLAP_PCT <- SNV_summary_5$OVERLAP/(SNV_summary_5$FF_UNIQ + SNV_summary_5$FFPE_UNIQ + SNV_summary_5$OVERLAP)
SNV_summary_10$OVERLAP_PCT <- SNV_summary_10$OVERLAP/(SNV_summary_10$FF_UNIQ + SNV_summary_10$FFPE_UNIQ + SNV_summary_10$OVERLAP)
SNV_summary_15$OVERLAP_PCT <- SNV_summary_15$OVERLAP/(SNV_summary_15$FF_UNIQ + SNV_summary_15$FFPE_UNIQ + SNV_summary_15$OVERLAP)
SNV_summary_20$OVERLAP_PCT <- SNV_summary_20$OVERLAP/(SNV_summary_20$FF_UNIQ + SNV_summary_20$FFPE_UNIQ + SNV_summary_20$OVERLAP)
SNV_summary_25$OVERLAP_PCT <- SNV_summary_25$OVERLAP/(SNV_summary_25$FF_UNIQ + SNV_summary_25$FFPE_UNIQ + SNV_summary_25$OVERLAP)
SNV_summary_30$OVERLAP_PCT <- SNV_summary_30$OVERLAP/(SNV_summary_30$FF_UNIQ + SNV_summary_30$FFPE_UNIQ + SNV_summary_30$OVERLAP)
SNV_summary_0$VAF <- "0"
SNV_summary_5$VAF <- "0.05"
SNV_summary_10$VAF <- "0.10"
SNV_summary_15$VAF <- "0.15"
SNV_summary_20$VAF <- "0.20"
SNV_summary_25$VAF <- "0.25"
SNV_summary_30$VAF <- "0.30"
SNV_summary_all <- rbind(SNV_summary_0, SNV_summary_5, SNV_summary_10, SNV_summary_15, SNV_summary_20, SNV_summary_25, SNV_summary_30)

# Plot all data
ggplot(SNV_summary_all, aes(x = VAF, y = OVERLAP_PCT)) + geom_boxplot() + labs(x = "VAF", y = "Overlap %") + theme(legend.title=element_blank())

# Plot stratified by tumor type and purity
SNV_summary_all$TumorType <- QC_table[match(SNV_summary_all$PATIENT_ID, QC_table$PATIENT_ID),]$TumorType
SNV_summary_all$FF_TUMOUR_PURITY <- QC_table[match(SNV_summary_all$PATIENT_ID, QC_table$PATIENT_ID),]$FF_TUMOUR_PURITY
SNV_summary_all$FFPE_TUMOUR_PURITY <- QC_table[match(SNV_summary_all$PATIENT_ID, QC_table$PATIENT_ID),]$FFPE_TUMOUR_PURITY

# Boxplots by tumor type
ggplot(SNV_summary_all, aes(x = VAF, y = OVERLAP_PCT, col = factor(TumorType))) + geom_boxplot(aes(col = factor(TumorType))) + labs(x = "VAF", y = "Overlap %") + theme(legend.title=element_blank())

# Boxplots by tumor purity
ggplot(SNV_summary_all, aes(x = VAF, y = OVERLAP_PCT, col = factor(FFPE_TUMOUR_PURITY > 50))) + geom_boxplot(aes(col = factor(FFPE_TUMOUR_PURITY > 50))) + labs(x = "VAF", y = "Overlap %") + theme(legend.title=element_blank())
ggplot(SNV_summary_all, aes(x = VAF, y = OVERLAP_PCT, col = factor(FF_TUMOUR_PURITY > 50))) + geom_boxplot(aes(col = factor(FF_TUMOUR_PURITY > 50))) + labs(x = "VAF", y = "Overlap %") + theme(legend.title=element_blank())

# Plots by sample
ggplot(SNV_summary_all, aes(x = VAF, y = OVERLAP_PCT, group = PATIENT_ID)) + geom_line(aes(col = FF_TUMOUR_PURITY), size = 1) + labs(x = "VAF", y = "Overlap %")

# List of samples by purity and tumor type
SNV_summary_all %>% filter(VAF == "0") %>% select(PATIENT_ID, TumorType, FF_TUMOUR_PURITY, FFPE_TUMOUR_PURITY)

# Summary of overlap by different VAF cutoffs
SNV_summary_all %>% group_by(VAF) %>% summarise(MEAN_OVERLAP_PCT = mean(OVERLAP_PCT), MIN_OVERLAP_PCT = min(OVERLAP_PCT), MAX_OVERLAP_PCT = max(OVERLAP_PCT))


# Plot FF unique # by VAF
ggplot(SNV_summary_all, aes(x = VAF, y = FF_UNIQ, group = PATIENT_ID)) + geom_line(aes(col = FF_TUMOUR_PURITY), size = 1) + labs(x = "VAF", y = "FF Unique Variants")
# Plot FFPE unique # by VAF
ggplot(SNV_summary_all, aes(x = VAF, y = FFPE_UNIQ, group = PATIENT_ID)) + geom_line(aes(col = FF_TUMOUR_PURITY), size = 1) + labs(x = "VAF", y = "FFPE Unique Variants")


######### Correlation of max overlap with QC metrics ######### 

# I will use the maximum overlap % (calculated on per-sample basis by using different VAF cutoffs) to calculate correlations with QC metrics
# The idea is that, since overlap % also depends on tumor purity, clonality and VAF, maximum overlap will normalize some of these effects
SNV_summary_all$PATIENT_ID <- as.character(SNV_summary_all$PATIENT_ID)
SNV_summary_maxO <- as.data.frame(SNV_summary_all %>% group_by(PATIENT_ID) %>% summarise(MAX_OVERLAP = max(OVERLAP_PCT)))

# Add maximum overlap % to QC data
SNV_summary <- full_join(SNV_summary, SNV_summary_maxO, by = "PATIENT_ID")

# Plot maximum overlap % correlation with QC metrics
ggplot(SNV_summary, aes(x = MAX_OVERLAP, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum FF-FFPE overlap", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_OVERLAP, SNV_summary$FFPE_AT_DROP, method = "spearman")  # 0.128569

ggplot(SNV_summary, aes(x = MAX_OVERLAP, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum FF-FFPE overlap", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_OVERLAP, SNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.1117949

ggplot(SNV_summary, aes(x = MAX_OVERLAP, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum FF-FFPE overlap", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_OVERLAP, SNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  # -0.1641868

ggplot(SNV_summary, aes(x = MAX_OVERLAP, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum FF-FFPE overlap", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_OVERLAP, SNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # -0.1405128

ggplot(SNV_summary, aes(x = MAX_OVERLAP, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum FF-FFPE overlap", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_OVERLAP, SNV_summary$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # -0.117265

# Plot maximum precision correlation with QC metrics
SNV_summary_maxP <- as.data.frame(SNV_summary_all %>% group_by(PATIENT_ID) %>% summarise(MAX_PRECISION = max(PRECISION)))
SNV_summary <- full_join(SNV_summary, SNV_summary_maxP, by = "PATIENT_ID")

ggplot(SNV_summary, aes(x = MAX_PRECISION, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_PRECISION, SNV_summary$FFPE_AT_DROP, method = "spearman")  # 0.01060192

ggplot(SNV_summary, aes(x = MAX_PRECISION, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_PRECISION, SNV_summary$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.006496837

ggplot(SNV_summary, aes(x = MAX_PRECISION, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_PRECISION, SNV_summary$FFPE_CHIMERIC_PER, method = "spearman")  #-0.08997606

ggplot(SNV_summary, aes(x = MAX_PRECISION, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_PRECISION, SNV_summary$FFPE_MAPPING_RATE_PER, method = "spearman")  # -0.09813644

ggplot(SNV_summary, aes(x = MAX_PRECISION, y = FFPE_DEAMINATION_MISMATCHES_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Maximum precision", y = "FFPE Deamination Mismatches") + theme(legend.title=element_blank())
cor(SNV_summary$MAX_PRECISION, SNV_summary$FFPE_DEAMINATION_MISMATCHES_PER, method = "spearman")  # -0.1258335




######### Purity, clonality simulation ######### 

# SIM 1: 50% purity both

simSNV <- function(VAF_cutoff, PUR_FF = 50, PUR_FFPE = 50){
  ff <- data.frame(KEY = paste0("VAR", 1:2000), VAF = rnorm(2000, (PUR_FF/200), 0.02) )
  ffpe <- data.frame(KEY = paste0("VAR", 1:2000), VAF = rnorm(2000, (PUR_FFPE/200), 0.02) )
  ff <- ff %>% filter(VAF > VAF_cutoff)
  ffpe <- ffpe %>% filter(VAF > VAF_cutoff)
  overlap <- sum(ff$KEY %in% ffpe$KEY)
  result <- data.frame(FF_TOTAL = dim(ff)[1], 
                       FFPE_TOTAL = dim(ffpe)[1], 
                       OVERLAP = overlap, 
                       FF_UNIQ = ((dim(ff)[1]) - overlap), 
                       FFPE_UNIQ = ((dim(ffpe)[1]) - overlap), 
                       RECALL = overlap/dim(ff)[1], 
                       PRECISION = overlap/dim(ffpe)[1])
  return(result)
}


set.seed(100)
sim <- as.data.frame(sapply(seq(0, 0.3, 0.02), simSNV))
sim <- as.data.frame(t(sim))
sim$VAF <- as.character(seq(0, 0.3, 0.02))
sim$ID <- 1









