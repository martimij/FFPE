# Martina Mijuskovic
# FFPE project
# Reads the FFPE trio manifest, finds corresponding VCF files and compares Manta SVs between FF and FFPE using bedtools (on HPC)
# Installing R packages: use lib = "~/R/x86_64-pc-linux-gnu-library/3.3"

library(dplyr)
library(VariantAnnotation)
library(reshape)
library(ggplot2)
library(scales)

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
  # Add START, CHR (correct END already exists in the info table, note that for insertions - start/end are the same)
  SVinfo_ff$START <- as.data.frame(ranges(ff_vcf))$start
  SVinfo_ff$CHR <- as.character(seqnames(ff_vcf))
  
  # Add number of supporting paired and reads
  SVinfo_ff$PR_REF <- sapply(1:dim(SVinfo_ff)[1], function(x){
    geno(ff_vcf)[[1]][x][[1]][1]
  })
  SVinfo_ff$PR_ALT <- sapply(1:dim(SVinfo_ff)[1], function(x){
    geno(ff_vcf)[[1]][x][[1]][2]
  })
  
  # Add number of supporting split reads
  SVinfo_ff$SR_REF <- sapply(1:dim(SVinfo_ff)[1], function(x){
    geno(ff_vcf)[[2]][x][[1]][1]
  })
  SVinfo_ff$SR_ALT <- sapply(1:dim(SVinfo_ff)[1], function(x){
    geno(ff_vcf)[[2]][x][[1]][2]
  })
  
  # Add PR/SR flag (1=evidence based on PR/SR exists for somatic variant)
  SVinfo_ff$PR_EV <- as.numeric(SVinfo_ff$PR_ALT > 0)
  SVinfo_ff$SR_EV <- as.numeric(SVinfo_ff$SR_ALT > 0)  
  
  # Remove filtered entries, keep Manta only (WARNING: some "good" >10kb SVs might be filtered out too, filter "MGE10kb")
  Manta_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Manta")
  
  # Remove unknown CHR from the table and add chrY
  normal_chr <- c(paste0("chr", 1:22), "chrX", "chrY")  
  Manta_ff <- Manta_ff %>% filter(CHR %in% normal_chr)
  
  # Filter out imprecise SVs for this purpose
  Manta_ff <- Manta_ff %>% filter(IMPRECISE == FALSE)
  
  # Filter out SVs with <3 supporting PR AND SR reads
  Manta_ff <- Manta_ff[!((Manta_ff$PR_ALT <3) & (Manta_ff$SR_ALT <3)),]
  
  
  
  ### FFPE info
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
  Manta_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Manta")
  
  # Remove unknown CHR from the table and add chrY
  Manta_ffpe <- Manta_ffpe %>% filter(CHR %in% normal_chr)
  
  # Filter out imprecise SVs for this purpose
  Manta_ffpe <- Manta_ffpe %>% filter(IMPRECISE == FALSE)
  
  # Filter out SVs with <3 supporting PR AND SR reads
  Manta_ffpe <- Manta_ffpe[!((Manta_ffpe$PR_ALT <3) & (Manta_ffpe$SR_ALT <3)),]  # 806
  
  # Merge FF and FFPE into one table, keep info on SV source (FF or FFPE?)
  Manta_ff$FF <- "FF"
  Manta_ffpe$FF <- "FFPE"
  ff_ffpe_merged <- rbind(Manta_ff, Manta_ffpe)
  # Add new type coding ("+" is FF, "-" is FFPE), to keep info in bed files downstream
  ff_ffpe_merged$Type2 <- NA
  ff_ffpe_merged$Type2 <- sapply(1:dim(ff_ffpe_merged)[1], function(x){
    if (ff_ffpe_merged$FF[x] == "FF") {ff_ffpe_merged$Type2[x] <- "+"}
    else if (ff_ffpe_merged$FF[x] == "FFPE") {ff_ffpe_merged$Type2[x] <- "-"}
  })
  # Create CHR_START_END_TYPE key
  ff_ffpe_merged$KEY <- sapply(1:dim(ff_ffpe_merged)[1], function(x){
    paste(ff_ffpe_merged$CHR[x], ff_ffpe_merged$START[x], ff_ffpe_merged$END[x], ff_ffpe_merged$SVTYPE[x], sep = "-")
  })
  
  ### Filter out likely false positives and keep high quality SV candidates only (for comparison)
  ### (Remove SVs where one or both breaksites +/- 150 bp overlap windowMasker, simple repeats and/or segmental duplications)
  
  # Create a bed file with SVs, start and end separately, remove NAs, code FF and FFPE as strand ("+" for FF, "-" for FFPE)
  # Note that for BED format START has to be adjusted to 0-based (-1 bp) and END position is not included (stays same)
  # Add 150 bp on either side of the breakpoint
  start_bed <- cbind((ff_ffpe_merged %>% dplyr::select(CHR, START)), (ff_ffpe_merged %>% dplyr::select(START, KEY, Type2)))
  names(start_bed)[3] <- "end"
  start_bed$Score <- ""
  start_bed <- start_bed %>% dplyr::select(CHR, START, end, KEY, Score, Type2)
  start_bed <- start_bed %>% filter(!is.na(START))
  # Adjust window around breaksite
  start_bed$START <- start_bed$START - 151
  start_bed$end <- start_bed$end + 150
  write.table(start_bed, file = paste0(patientID, "_sv_start.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
  
  end_bed <- cbind((ff_ffpe_merged %>% dplyr::select(CHR, END)), (ff_ffpe_merged %>% dplyr::select(END, KEY, Type2)))
  names(end_bed)[2] <- "start"
  end_bed$Score <- ""
  end_bed <- end_bed %>% dplyr::select(CHR, start, END, KEY, Score, Type2)
  end_bed <- end_bed %>% filter(!is.na(END))
  # Adjust window around breaksite
  end_bed$start <- end_bed$start - 151
  end_bed$END <- end_bed$END + 150
  write.table(end_bed, file = paste0(patientID, "_sv_end.bed"), quote = F, row.names = F, col.names = F, sep = "\t")
  
  
  #############
  ### Call bedtools to find overlaps with WindowMasker 
  # start
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -a", paste0(patientID, "_sv_start.bed"), "-b /home/mmijuskovic/FFPE/windowmaskerSdust.hg38.bed >", paste0(patientID, "_sv_wMasker_start_overlap.bed")), intern = T)
  sv_wMasker_start <- read.table(paste0(patientID, "_sv_wMasker_start_overlap.bed"), sep = "\t")
  names(sv_wMasker_start) <- c("CHR", "START", "END", "KEY", "Score", "Type2", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  # end
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -a", paste0(patientID, "_sv_end.bed"), "-b /home/mmijuskovic/FFPE/windowmaskerSdust.hg38.bed >", paste0(patientID, "_sv_wMasker_end_overlap.bed")), intern = T)
  sv_wMasker_end <- read.table(paste0(patientID, "_sv_wMasker_end_overlap.bed"), sep = "\t")
  names(sv_wMasker_end) <- c("CHR", "START", "END", "KEY", "Score", "Type2", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  
  # FLAG SVs where START or END overlaps with WindowMasker
  wMasker_keys <- unique(c(as.character(sv_wMasker_start %>% filter(NumOverlap != 0) %>% .$KEY), as.character(sv_wMasker_end %>% filter(NumOverlap != 0) %>% .$KEY)))
  ff_ffpe_merged$wMasker_filtered <- 0
  ff_ffpe_merged[(ff_ffpe_merged$KEY %in% wMasker_keys),]$wMasker_filtered <- 1
  
  
  ### Call bedtools to find overlaps with simple repeats
  
  # start
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -a", paste0(patientID, "_sv_start.bed"), "-b /home/mmijuskovic/FFPE/simpleRepeat.hg38.bed >", paste0(patientID, "_sv_repeats_start_overlap.bed")), intern = T)
  sv_repeats_start <- read.table(paste0(patientID, "_sv_repeats_start_overlap.bed"), sep = "\t")
  names(sv_repeats_start) <- c("CHR", "START", "END", "KEY", "Score", "Type2", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  # end
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -a", paste0(patientID, "_sv_end.bed"), "-b /home/mmijuskovic/FFPE/simpleRepeat.hg38.bed >", paste0(patientID, "_sv_repeats_end_overlap.bed")), intern = T)
  sv_repeats_end <- read.table(paste0(patientID, "_sv_repeats_end_overlap.bed"), sep = "\t")
  names(sv_repeats_end) <- c("CHR", "START", "END", "KEY", "Score", "Type2", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  
  # FLAG SVs where START or END overlaps with repeats
  repeats_keys <- unique(c(as.character(sv_repeats_start %>% filter(NumOverlap != 0) %>% .$KEY), as.character(sv_repeats_end %>% filter(NumOverlap != 0) %>% .$KEY)))
  ff_ffpe_merged$repeats_filtered <- 0
  ff_ffpe_merged[(ff_ffpe_merged$KEY %in% repeats_keys),]$repeats_filtered <- 1
  
  
  
  ### Call bedtools to find overlaps with segmental duplications
  
  # start
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -a", paste0(patientID, "_sv_start.bed"), "-b /home/mmijuskovic/FFPE/genomicSuperDups.hg38.bed >", paste0(patientID, "_sv_segdups_start_overlap.bed")), intern = T)
  sv_segdups_start <- read.table(paste0(patientID, "_sv_segdups_start_overlap.bed"), sep = "\t")
  names(sv_segdups_start) <- c("CHR", "START", "END", "KEY", "Score", "Type2", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  # end
  system(paste("/home/mmijuskovic/bedtools2/bin/bedtools coverage -a", paste0(patientID, "_sv_end.bed"), "-b /home/mmijuskovic/FFPE/genomicSuperDups.hg38.bed >", paste0(patientID, "_sv_segdups_end_overlap.bed")), intern = T)
  sv_segdups_end <- read.table(paste0(patientID, "_sv_segdups_end_overlap.bed"), sep = "\t")
  names(sv_segdups_end) <- c("CHR", "START", "END", "KEY", "Score", "Type2", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
  
  # FLAG SVs where START or END overlaps with repeats
  segdups_keys <- unique(c(as.character(sv_segdups_start %>% filter(NumOverlap != 0) %>% .$KEY), as.character(sv_segdups_end %>% filter(NumOverlap != 0) %>% .$KEY)))
  ff_ffpe_merged$segdups_filtered <- 0
  ff_ffpe_merged[(ff_ffpe_merged$KEY %in% segdups_keys),]$segdups_filtered <- 1 

  
  
  ### FILTER low confidence SVs
  
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
  filtered_id <- unique(c(filtered_id, (ff_ffpe_merged %>% filter(BND_MATE_FILTERED == 1) %>% .$ID)))
  ff_ffpe_merged[ff_ffpe_merged$ID %in% filtered_id,]$FILTERED <- 1
  
  # Add CONCORDANT flag to FF and FFPE calls
  concordant_keys <- ff_ffpe_merged[duplicated(ff_ffpe_merged$KEY),]$KEY
  ff_ffpe_merged$CONCORDANT <- 0
  ff_ffpe_merged[ff_ffpe_merged$KEY %in% concordant_keys,]$CONCORDANT <- 1
  
  # Write out the filtered SVs
  ff_ffpe_merged <- ff_ffpe_merged %>% dplyr::select(-(CSQT))  # Get rid of transcripts to save space & avoid parsing problems
  ff_ffpe_merged$PATIENT_ID <- patientID
  write.table(ff_ffpe_merged, file = paste0(patientID, "_SV_filtered", ".tsv"), row.names = F, quote = F, sep = "\t")

}  

# Get all patient IDs
patientIDs <- unique(QC_portal_trios$PATIENT_ID)

# Run SV function for each patient ID
sapply(patientIDs, compareSV)

# Get all data together
result <- data.frame()
result <- lapply(patientIDs, function(x){
  table <- read.table(paste0(x, "_SV_filtered", ".tsv"), sep = "\t", header = T)
  result <- rbind(result, table)
})

# Merge into one data frame
result <- merge_recurse(result)

# Write out the result
write.table(result, file = paste0(today, "_SV_all", ".tsv"), row.names = F, quote = F, sep = "\t")

# Read in the result
result <- read.table("2017-02-22_SV_all.tsv", sep = "\t", header = T)

# Overview of filtered concordant SVs (NOTE that BNDs are listed twice - each side as a separate BND, need to adjust for that)
table(result[result$FILTERED ==0,]$PATIENT_ID, result[result$FILTERED ==0,]$CONCORDANT)
result %>% filter(FILTERED == 0) %>% dplyr::select(PATIENT_ID, KEY, SVTYPE, SVLEN, FF, ID, MATEID, PR_ALT, SR_ALT, CONCORDANT)

### Look at the recurrent DEL (MIR7155 gene)
result %>% filter(FILTERED == 0, KEY == "chr11-64341843-64341923-DEL") %>% dplyr::select(PATIENT_ID, KEY, SVTYPE, SVLEN, FF, ID, MATEID, PR_ALT, SR_ALT, CONCORDANT)
# By PATIENT_ID
length(unique(result %>% filter(FILTERED == 0, KEY == "chr11-64341843-64341923-DEL") %>% .$PATIENT_ID))
# By FF/FFPE (non-condordant)
table(result %>% filter(FILTERED == 0, CONCORDANT == 0, KEY == "chr11-64341843-64341923-DEL") %>% .$FF)

# For each SV, check how many unique PATIENT IDs is IN
dim(result %>% filter(FILTERED == 0))  # 139 total fitered SVs
unique_keys <- as.character(unique(result %>% filter(FILTERED == 0) %>% .$KEY))  # 62 unique KEYs
result$KEY <- as.character(result$KEY)
recurrent_SVs <- data.frame(
  NUM_OBS = sapply(unique_keys, function(x){ sum(result$KEY == x)}), 
  NUM_PATIENTS = sapply(unique_keys, function(x){ length(unique(result %>% filter(KEY == x) %>% .$PATIENT_ID)) })
  )
recurrent_SVs$KEY <- rownames(recurrent_SVs)
rownames(recurrent_SVs) <- NULL
table(recurrent_SVs$NUM_PATIENTS)  # 3 observed in >1 patient (3, 16 and 22 patients)
recurr_KEYs <- recurrent_SVs %>% filter(NUM_PATIENTS >1) %>% .$KEY
recurrent_SVs %>% filter(KEY %in% recurr_KEYs)
result %>% filter(KEY %in% recurr_KEYs) %>% dplyr::select(KEY, SVLEN, SR_REF, SR_ALT, CONCORDANT, PATIENT_ID)

### Remove 3 recurrent SVs and plot recall/precision per patient
result_filt <- result %>% filter(FILTERED == 0, !(KEY %in% recurr_KEYs))  # 90 left

# Table of concordance vs SVTYPE (note that concordant SV number is half of the number listed)
table(result_filt$CONCORDANT, result_filt$SVTYPE)
# Discordant SVs by sample type
table(result_filt[result_filt$CONCORDANT == 0,]$FF, result_filt[result_filt$CONCORDANT ==0,]$SVTYPE)

# Concordance by patient and SV type
table(result_filt$CONCORDANT, result_filt$SVTYPE, result_filt$PATIENT_ID)




### Put all data together (QC, SV)

# Keep only one of BND mates (remove all with ID ending with "0")
# NOTE that one SV has a mate missing from the original results file, possibly missing in the VCF --- CHECK THIS (MantaBND:341610:0:1:3:1:0:0)
result_filt$ID <- as.character(result_filt$ID)
result_filt$MATE <- sapply(1:dim(result_filt)[1], function(x){
  strsplit(result_filt$ID[x], split = "[[:punct:]]")[[1]][length(strsplit(result_filt$ID[x], split = "[[:punct:]]")[[1]])]
})
result_filt <- result_filt %>% filter(!(SVTYPE == "BND" & MATE == "0"))  # 11 BND mates removed

# Summarize SVs (SV count total, count FF only, count FFPE only, count overlap, recall, precision - per PATIENT_ID)
SV_summary <- data.frame(
  PATIENT_ID=(result_filt %>% group_by(PATIENT_ID) %>% summarise(n()))$PATIENT_ID,
  TOTAL_SV=(result_filt %>% group_by(PATIENT_ID) %>% summarise(n()))$`n()`,
  SV_OVERLAP=table(result_filt$PATIENT_ID, result_filt$CONCORDANT)[,2]/2,  # Fix, IDs skipped if no entries
  SV_FF_ONLY=table(result_filt[result_filt$CONCORDANT == 0,]$PATIENT_ID, result_filt[result_filt$CONCORDANT == 0,]$FF)[,1]
  SV_FFPE_ONLY=table(result_filt[result_filt$CONCORDANT == 0,]$PATIENT_ID, result_filt[result_filt$CONCORDANT == 0,]$FF)[,2]
  )

  
######### ------ continue here 
  
# Add QC details to the SV summary table
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



### Create plots

# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)

# Barplots of overlapping and unique SVs (normalized) for all 26 trios
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
CNV_summary_m <- as.data.frame(t(CNV_summary %>% select(PATIENT_ID, BP_OVERLAP, BP_FF_ONLY, BP_FFPE_ONLY)))
names(CNV_summary_m) <- CNV_summary_m[1,]
CNV_summary_m <- CNV_summary_m[2:4,]
CNV_summary_m <- melt(cbind(CNV_summary_m, ind = rownames(CNV_summary_m)), id.vars = c('ind'))
# Plot (needs package "scales")
ggplot(CNV_summary_m,aes(x = variable, y = value, fill = ind)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + blank







