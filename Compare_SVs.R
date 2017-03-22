# Martina Mijuskovic
# FFPE project
# Reads the FFPE trio manifest, finds corresponding VCF files and compares Manta SVs between FF and FFPE using bedtools (on HPC)
# Various SV summaries and plots


# NOTE: Installing R packages on HPC: use lib = "~/R/x86_64-pc-linux-gnu-library/3.3"

library(dplyr)
library(VariantAnnotation)
library(reshape)
library(ggplot2)
library(scales)
library(R.utils)

# Working directory on the HPC
setwd("/home/mmijuskovic/FFPE/SV_trio_comparison")

### For plots
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)


######### Compile and filter Manta SVs ######### 

today <- Sys.Date()

# Load the manifest (HPC)
#QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/CNV_trio_comparison/QC_portal_trios.csv")  # old
#QC_portal_trios <- read.csv("/Users/MartinaMijuskovic/FFPE/QC_portal_62_trios.csv")  # local
QC_portal_trios <- read.csv("/home/mmijuskovic/FFPE/QC_portal_62_trios.csv")

# Subset for FF and FFPE samples
QC_portal_trios <- QC_portal_trios %>% filter(SAMPLE_TYPE %in% c("FF", "FFPE"))

# Get VCF paths (HPC)
QC_portal_trios$VCF_path <- as.character(sapply(QC_portal_trios$BamPath, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.SV.vcf.gz", sep = " ")
  system(command, intern = T)
}))


# Function that compares SVs from FFPE trios
# Extracts only precise SVs that pass Manta filters, from regular chromosomes, with at least 3 supporting reads (either 3 PR or 3 SR minimum)
# Filters out SVs overlapping windowMasker, simple repeats or segmental duplications (either end of SV)
# Writes out the filtered SV table
# 20170321 Modified to include SVs >10kb
# 20170321 Fixed to allow for no BND
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
  # 20170321 Keep PASS and MGE10kb
  Manta_ff <- SVinfo_ff %>% filter(FILTER %in% c("PASS","MGE10kb"), Application == "Manta")
  
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
  # 20170321 Keep PASS and MGE10kb
  Manta_ffpe <- SVinfo_ffpe %>% filter(FILTER %in% c("PASS","MGE10kb"), Application == "Manta")
  
  # Remove unknown CHR from the table and add chrY
  Manta_ffpe <- Manta_ffpe %>% filter(CHR %in% normal_chr)
  
  # Filter out imprecise SVs for this purpose
  Manta_ffpe <- Manta_ffpe %>% filter(IMPRECISE == FALSE)
  
  # Filter out SVs with <3 supporting PR AND SR reads
  Manta_ffpe <- Manta_ffpe[!((Manta_ffpe$PR_ALT <3) & (Manta_ffpe$SR_ALT <3)),] 
  
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
  # Check if there are any BND
  if ( dim(ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",])[1] != 0 ) {
    ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",]$BND_MATE_FILTERED <- sapply(1:dim(ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",])[1], function(x){
      if (ff_ffpe_merged[ff_ffpe_merged$SVTYPE == "BND",]$MATEID[x] %in% filtered_id) { return(1)}
      else {return(0)}
    })
  }
  filtered_id <- unique(c(filtered_id, (ff_ffpe_merged %>% filter(BND_MATE_FILTERED == 1) %>% .$ID)))
  ff_ffpe_merged[ff_ffpe_merged$ID %in% filtered_id,]$FILTERED <- 1
  
  # Add CONCORDANT flag to FF and FFPE calls
  concordant_keys <- ff_ffpe_merged[duplicated(ff_ffpe_merged$KEY),]$KEY
  ff_ffpe_merged$CONCORDANT <- 0
  # Add a check for any concordant calls
  if ( length(concordant_keys) != 0 ) {
    ff_ffpe_merged[ff_ffpe_merged$KEY %in% concordant_keys,]$CONCORDANT <- 1
  }
  
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
write.table(result, file = paste0(today, "_62_FFPEtrios_SV_all", ".tsv"), row.names = F, quote = F, sep = "\t")


######### Check and clean SVs ######### 

# Read in the result
result <- read.table("2017-03-22_62_FFPEtrios_SV_all.tsv", sep = "\t", header = T)
result$PATIENT_ID <- as.character(result$PATIENT_ID)
# Sanity check & summary
dim(result)  # 28064
length(unique(result$PATIENT_ID))  # 62
sum(duplicated(result$KEY))   # 7926
table(result$CONCORDANT, result$SVTYPE)
table(result$FF, result$SVTYPE)
table(result$FILTER, result$SVTYPE)
table(result$ColocalizedCanvas, result$FILTER)

# Read in the QC data
QC_portal_trios <- read.csv("/Users/MartinaMijuskovic/FFPE/QC_portal_62_trios.csv")

# Overview of filtered concordant SVs (NOTE that BNDs are listed twice - each side as a separate BND, need to adjust for that)
table(result[result$FILTERED ==0,]$PATIENT_ID, result[result$FILTERED ==0,]$CONCORDANT) # concordant events counted twice
result %>% filter(FILTERED == 0) %>% dplyr::select(PATIENT_ID, KEY, SVTYPE, SVLEN, FF, ID, MATEID, PR_ALT, SR_ALT, CONCORDANT)

########### Recurrent SVs ########### 

# ### Look at the recurrent DEL (MIR7155 gene)
# result %>% filter(FILTERED == 0, KEY == "chr11-64341843-64341923-DEL") %>% dplyr::select(PATIENT_ID, KEY, SVTYPE, SVLEN, FF, ID, MATEID, PR_ALT, SR_ALT, CONCORDANT)
# # By PATIENT_ID
# length(unique(result %>% filter(FILTERED == 0, KEY == "chr11-64341843-64341923-DEL") %>% .$PATIENT_ID))
# # By FF/FFPE (non-condordant)
# table(result %>% filter(FILTERED == 0, CONCORDANT == 0, KEY == "chr11-64341843-64341923-DEL") %>% .$FF)


### UNFILTERED SVs

# For each SV, check how many unique PATIENT IDs is IN
dim(result)  # 28064 total unfitered SVs (62 trios), translocation BNDs listed separately
unique_keys <- as.character(unique(result %>% filter(FILTERED == 0) %>% .$KEY))  # 133 unique KEYs (62 trios)
result$KEY <- as.character(result$KEY)
recurrent_SVs <- data.frame(
  NUM_OBS = sapply(unique_keys, function(x){ sum(result$KEY == x)}), 
  NUM_PATIENTS = sapply(unique_keys, function(x){ length(unique(result %>% filter(KEY == x) %>% .$PATIENT_ID)) })
)
recurrent_SVs$KEY <- rownames(recurrent_SVs)
rownames(recurrent_SVs) <- NULL
table(recurrent_SVs$NUM_PATIENTS)  # 4 observed in >1 patient (62 trios)
recurr_KEYs <- recurrent_SVs %>% filter(NUM_PATIENTS >1) %>% .$KEY
recurrent_SVs %>% filter(KEY %in% recurr_KEYs)
result %>% filter(KEY %in% recurr_KEYs) %>% dplyr::select(KEY, SVLEN, SR_REF, SR_ALT, CONCORDANT, PATIENT_ID)
# Concordance of recurrent SVs
table((result %>% filter(KEY %in% recurr_KEYs) %>% .$CONCORDANT))

# ### FILTERED SVs
# 
# # For each SV, check how many unique PATIENT IDs is IN
# dim(result %>% filter(FILTERED == 0))  # 326 total fitered SVs (62 trios)
# unique_keys <- as.character(unique(result %>% filter(FILTERED == 0) %>% .$KEY))  # 133 unique KEYs (62 trios)
# result$KEY <- as.character(result$KEY)
# recurrent_SVs <- data.frame(
#   NUM_OBS = sapply(unique_keys, function(x){ sum(result$KEY == x)}), 
#   NUM_PATIENTS = sapply(unique_keys, function(x){ length(unique(result %>% filter(KEY == x) %>% .$PATIENT_ID)) })
#   )
# recurrent_SVs$KEY <- rownames(recurrent_SVs)
# rownames(recurrent_SVs) <- NULL
# table(recurrent_SVs$NUM_PATIENTS)
# recurr_KEYs <- recurrent_SVs %>% filter(NUM_PATIENTS >1) %>% .$KEY
# recurrent_SVs %>% filter(KEY %in% recurr_KEYs)
# result %>% filter(KEY %in% recurr_KEYs) %>% dplyr::select(KEY, SVLEN, SR_REF, SR_ALT, CONCORDANT, PATIENT_ID)


#########  Summarize SVs ######### 

# Add MATE_KEY to the main result table
result$KEY <- as.character(result$KEY)
result$ID <- as.character(result$ID)
result$MATEID <- as.character(result$MATEID)
result$MATE_KEY <- NA
result[result$SVTYPE == "BND",]$MATE_KEY <- sapply(1:dim(result[result$SVTYPE == "BND",])[1], function(x){
  result[result$ID == (result[result$SVTYPE == "BND",]$MATEID[x]),]$KEY
})


# Extract the common ID part for both BND mates
result$ID <- as.character(result$ID)
result$MATEID <- as.character(result$MATEID)
result$MAIN_ID  <- sapply(1:dim(result)[1], function(x){
  stop <- nchar(result$ID[x])-1
  substr(result$ID[x], 1, stop)
})

# Check that there are always 2 BND for each MAIN_ID
result %>% group_by(MAIN_ID) %>% summarise(n())  # some BND are missing the mate (mates on uknown CHR, filtered out previously)
# NOTE that some MAIN_IDs are duplicated across samples (different SVs in different samples); this is not the best way to filter BNDs


### Summaries by filter

# Create a new filter variable
result$FILTER2 <- sapply(1:dim(result)[1], function(x) {
  if (result$wMasker_filtered[x] == 1) {"wMasker"} 
  else if (result$repeats_filtered[x] == 1) {"simple_rep"} 
  else if (result$segdups_filtered[x] == 1) {"seg_dupl"} 
  else if ((!is.na(result$BND_MATE_FILTERED[x])) & result$BND_MATE_FILTERED[x] == 1) {"mate_filtered"}
  else {"PASS"}
})
# Overview of all SVs by FILTER  
table(result$FILTER2, result$SVTYPE)
# Overview of FF/FFPE SVs separately by FILTER
table(result[result$FF == "FF",]$FILTER2, result[result$FF == "FF",]$SVTYPE)
table(result[result$FF == "FFPE",]$FILTER2, result[result$FF == "FFPE",]$SVTYPE)

# > sum(is.na(result[result$FF == "FFPE",]$MATE_KEY))
# [1] 8615
# > dim(result[result$FF == "FFPE",])
# [1] 14560    56

##################### ------continue HERE



# # Remove the BND with a filtered mate (discovered above) -- checked VCF; mate maps to unkn chr
# result <- result %>% filter(ID != "MantaBND:341610:0:1:3:1:0:1")

#### Remove BNDs with no mate (mate filtered previously!!!)  ------ continue here



############ Filter SVs ############ 

### Remove 4 recurrent SVs and plot recall/precision per patient
result_filt <- result %>% filter(FILTERED == 0, !(KEY %in% recurr_KEYs))  # 260 left

# Keep only one BND mate
result_filt <- result_filt[!duplicated(result_filt$MAIN_ID),]

# Add TUMOR TYPE
result_filt$TumorType <- QC_portal_trios[match(result_filt$PATIENT_ID, QC_portal_trios$PATIENT_ID),]$TumorType








# List unique SVs to be reviewed manually (51 total)
result_filt[(!duplicated(result_filt$KEY)),] %>% dplyr::select(PATIENT_ID, TumorType, FF, KEY, MATE_KEY, ID, SVLEN, CIGAR, HOMLEN, LEFT_SVINSSEQ, RIGHT_SVINSSEQ, SOMATICSCORE, PR_ALT, SR_ALT, CONCORDANT)
#write.table((result_filt[(!duplicated(result_filt$KEY)),] %>% dplyr::select(PATIENT_ID, TumorType, FF, KEY, MATE_KEY, ID, SVLEN, CIGAR, HOMLEN, LEFT_SVINSSEQ, RIGHT_SVINSSEQ, SOMATICSCORE, PR_ALT, SR_ALT, CONCORDANT)), quote = F, row.names = F, sep = "\t", file = "HighConf_SVs_FFPE_trios.tsv")

# Total unique SVs
length(unique(result_filt$KEY))  # 51
# Table of concordance vs SVTYPE (note that concordant SV number is half of the number listed)
table(result_filt$CONCORDANT, result_filt$SVTYPE)
# Discordant SVs by sample type
table(result_filt[result_filt$CONCORDANT == 0,]$FF, result_filt[result_filt$CONCORDANT ==0,]$SVTYPE)
# Total FF SVs
table(result_filt$FF)
# Total concordant
table(result_filt$CONCORDANT, result_filt$FF)
# Total patients with any high confidence SV calls
length(unique(result_filt$PATIENT_ID))  # 19


### Put all data together (QC, SV)

# Summarize SVs (SV count total, count FF only, count FFPE only, count overlap, recall, precision - per PATIENT_ID)
SV_summary <- data.frame(
  PATIENT_ID=(result_filt %>% group_by(PATIENT_ID) %>% summarise(n()))$PATIENT_ID,
  TOTAL=(result_filt %>% group_by(PATIENT_ID) %>% summarise(n()))$`n()`,
  TOTAL_FF=table(result_filt$PATIENT_ID, result_filt$FF)[,1],
  TOTAL_FFPE=table(result_filt$PATIENT_ID, result_filt$FF)[,2],
  OVERLAP=table(result_filt$PATIENT_ID, result_filt$CONCORDANT)[,2]/2)

SV_summary$FF_ONLY <- SV_summary$TOTAL_FF - SV_summary$OVERLAP
SV_summary$FFPE_ONLY <- SV_summary$TOTAL_FFPE - SV_summary$OVERLAP
SV_summary$RECALL <- SV_summary$OVERLAP / SV_summary$TOTAL_FF
SV_summary$PRECISION <- SV_summary$OVERLAP / SV_summary$TOTAL_FFPE

  
# Add QC details to the SV summary table
QC_table <- QC_portal_trios %>% filter(SAMPLE_TYPE == "FFPE") %>% dplyr::select(PATIENT_ID, CENTER_CODE, TumorType, SAMPLE_WELL_ID, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER)
names(QC_table)[4:11] <- paste0("FFPE_", names(QC_table)[4:11])
QC_table <- left_join(QC_table, (QC_portal_trios %>% filter(SAMPLE_TYPE == "FF") %>% dplyr::select(PATIENT_ID, SAMPLE_WELL_ID, LIBRARY_TYPE, TUMOUR_PURITY, GC_DROP, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER)), by = "PATIENT_ID")
names(QC_table)[12:20] <- paste0("FF_", names(QC_table)[12:20])
SV_summary <- left_join(SV_summary, QC_table, by = "PATIENT_ID")

# Write the full table
write.csv(SV_summary, file = paste0("Full_SV_summary_", today, ".csv"), row.names = F, quote = F)



#########  Create plots and tables #########  


# Barplots of overlapping and unique SVs (normalized) for all 19 trios with SV calls
# First recast the data (each of 3 bp values in separate row, with PATIENT_ID, with indexes 1-2-3), needs package "reshape"
SV_summary_m <- as.data.frame(t(SV_summary %>% dplyr::select(PATIENT_ID, OVERLAP, FF_ONLY, FFPE_ONLY)))
names(SV_summary_m) <- SV_summary_m[1,]
SV_summary_m <- SV_summary_m[2:4,]
SV_summary_m <- melt(cbind(SV_summary_m, ind = rownames(SV_summary_m)), id.vars = c('ind'))
# Plot (needs package "scales")
ggplot(SV_summary_m, aes(x = variable, y = value, fill = ind)) + geom_bar(position = "fill",stat = "identity") + scale_y_continuous(labels = percent_format()) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1), axis.title = element_blank()) + theme(legend.title=element_blank()) + labs(x = "Patient ID", y = element_blank()) + blank



# Recall of FF vs recall of FFPE (recall vs precision)
ggplot(SV_summary[complete.cases(SV_summary),], aes(x = RECALL, y = PRECISION)) + geom_jitter() + geom_smooth(method = "lm", se = T) + labs(x = "Percent Recall of FF", y = "Percent Recall of FFPE") 
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$PRECISION, method = "spearman")  # 0.2496266

# RECALL of FF by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, colored by center
ggplot(SV_summary, aes(x = RECALL, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$FFPE_AT_DROP, method = "spearman")  # 0.1051732

ggplot(SV_summary, aes(x = RECALL, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # 0.102608

ggplot(SV_summary, aes(x = RECALL, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$FFPE_CHIMERIC_PER, method = "spearman")  # -0.2978906

ggplot(SV_summary, aes(x = RECALL, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Percent Recall of FF", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$FFPE_MAPPING_RATE_PER, method = "spearman")  # 0.1513467

# Recall by fragment size
ggplot(SV_summary, aes(x = RECALL, y = FFPE_AV_FRAGMENT_SIZE_BP, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent SV Recall", y = "FFPE Fragment Size") + theme(legend.title=element_blank()) + regr_line
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$FFPE_AV_FRAGMENT_SIZE_BP, method = "spearman")  # -0.4822574



# Recall of FF by AT dropout of FFPE, colored by tumor type (not very informative)
# ggplot(SV_summary, aes(x = RECALL, y = FFPE_AT_DROP, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
# ggplot(SV_summary, aes(x = RECALL, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())

# Recall by tumor purity (not very informative)
ggplot(SV_summary, aes(x = RECALL, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SV_summary[complete.cases(SV_summary),]$RECALL, SV_summary[complete.cases(SV_summary),]$FFPE_TUMOUR_PURITY, method = "spearman")  # -0.01158164

# ggplot(SV_summary, aes(x = RECALL, y = FF_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FF Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
# 
# # Recall by high vs low FF and FFPE purity (both > 50) (not very informative)
# ggplot(SV_summary, aes(x = RECALL, y = FFPE_AT_DROP, col = factor((FFPE_TUMOUR_PURITY > 50) & (FF_TUMOUR_PURITY > 50)))) + geom_jitter() + labs(x = "Percent Recall of FF", y = "FFPE AT Dropout") + theme(legend.title=element_blank())

# # Divide samples into good and bad FF recall (>50%), see which variables make a difference (not very informative)
# ggplot(data=SV_summary, aes(x=(RECALL>=0.5), y=FFPE_AT_DROP, fill=(RECALL>=0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High SV Recall (>=50%)", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
# ggplot(data=SV_summary, aes(x=(RECALL>=0.5), y=FFPE_COVERAGE_HOMOGENEITY, fill=(RECALL>=0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High SV Recall (>=50%)", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
# ggplot(data=SV_summary, aes(x=(RECALL>=0.5), y=FFPE_TUMOUR_PURITY, fill=(RECALL>=0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High SV Recall (>=50%)", y = "FFPE Tumour Purity") + theme(legend.title=element_blank())
# ggplot(data=SV_summary, aes(x=(RECALL>=0.5), y=FF_TUMOUR_PURITY, fill=(RECALL>=0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High SV Recall (>=50%)", y = "FF Tumour Purity") + theme(legend.title=element_blank())

# # Recall by % chimeric fragments
# ggplot(data=SV_summary, aes(x=(RECALL>=0.5), y=FFPE_CHIMERIC_PER, fill=(RECALL>=0.5))) + geom_boxplot() + geom_jitter() + labs(x = "High SV Recall (>=50%)", y = "FFPE Chimeric Frag") + theme(legend.title=element_blank())
# # Recall by tumor type (not enough data...)
# ggplot(data=SV_summary, aes(x=TumorType, y=RECALL)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "Tumor Type", y = "Percent SV Recall") + theme(legend.title=element_blank())
# # Recall by GMC, color by tumor type
# ggplot(data=SV_summary, aes(x=CENTER_CODE, y=RECALL)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "GMC", y = "Percent SV Recall") + theme(legend.title=element_blank())
# # Recall by FF library type
# ggplot(data=SV_summary, aes(x=FF_LIBRARY_TYPE, y=RECALL)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "FF Library Type", y = "Percent SV Recall") + theme(legend.title=element_blank())




# PRECISION by AT dropout/coverage homogeneity, chim reads, mapping rate of FFPE, colored by center
ggplot(SV_summary, aes(x = PRECISION, y = FFPE_AT_DROP, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE AT Dropout") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$PRECISION, SV_summary[complete.cases(SV_summary),]$FFPE_AT_DROP, method = "spearman")  # -0.007730393

ggplot(SV_summary, aes(x = PRECISION, y = FFPE_COVERAGE_HOMOGENEITY, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Unevenness of Coverage") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$PRECISION, SV_summary[complete.cases(SV_summary),]$FFPE_COVERAGE_HOMOGENEITY, method = "spearman")  # -0.07988073

ggplot(SV_summary, aes(x = PRECISION, y = FFPE_CHIMERIC_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE % Chimeric") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$PRECISION, SV_summary[complete.cases(SV_summary),]$FFPE_CHIMERIC_PER, method = "spearman")  # -0.472073

ggplot(SV_summary, aes(x = PRECISION, y = FFPE_MAPPING_RATE_PER, col = factor(CENTER_CODE))) + geom_jitter() + regr_line + labs(x = "Precision", y = "FFPE Mapping Rate") + theme(legend.title=element_blank())
cor(SV_summary[complete.cases(SV_summary),]$PRECISION, SV_summary[complete.cases(SV_summary),]$FFPE_MAPPING_RATE_PER, method = "spearman")  # -0.4097108

# Precision by tumor purity (not informative)
ggplot(SV_summary, aes(x = PRECISION, y = FFPE_TUMOUR_PURITY, col = factor(TumorType))) + geom_jitter() + labs(x = "Precision", y = "FFPE Tumour Purity")  + theme(legend.title=element_blank()) + regr_line
cor(SV_summary[complete.cases(SV_summary),]$PRECISION, SV_summary[complete.cases(SV_summary),]$FFPE_TUMOUR_PURITY, method = "spearman")  # 0.05429203
# PRECISION by GMC, color by tumor type (not informative)
ggplot(data=SV_summary, aes(x=CENTER_CODE, y=PRECISION)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "GMC", y = "Precision") + theme(legend.title=element_blank())
# PRECISION by FF library type (not informative)
ggplot(data=SV_summary, aes(x=FF_LIBRARY_TYPE, y=PRECISION)) + geom_boxplot() + geom_jitter(aes(col = TumorType)) + labs(x = "FF Library Type", y = "Precision") + theme(legend.title=element_blank())
# PRECISION by FFPE fragment size
ggplot(SV_summary, aes(x = PRECISION, y = FFPE_AV_FRAGMENT_SIZE_BP, col = factor(TumorType))) + geom_jitter() + labs(x = "Precision", y = "FFPE Fragment Size") + theme(legend.title=element_blank()) + regr_line
cor(SV_summary[complete.cases(SV_summary),]$PRECISION, SV_summary[complete.cases(SV_summary),]$FFPE_AV_FRAGMENT_SIZE_BP, method = "spearman")  # -0.4822574


### Read support for high confidence SVs

# Median read support by SV type
result_filt %>% group_by(SVTYPE) %>% summarise(median(SR_ALT), median(SR_REF), median(PR_ALT), median(PR_REF))
# Concordant SVs, support in FF and FFPE
ggplot((result_filt %>% filter(CONCORDANT == 1)), aes(x=factor(SVTYPE), y=PR_ALT)) + geom_boxplot(aes(colour = factor(FF))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Paired Reads")
ggplot((result_filt %>% filter(CONCORDANT == 1)), aes(x=factor(SVTYPE), y=SR_ALT)) + geom_boxplot(aes(colour = factor(FF))) + blank + theme(axis.title.x = element_blank()) + ggtitle("Supporting Split Reads")




#########  Overlap with Domain 1, 2 genes #########  

# Read tables containing Domain 1 (actionable) and Domain 2 (cancer census) genes
domain1 <- read.csv("GENOMONCOLOGY_SOLID_TUMOUR.SV.v1.3.csv")
domain2 <- read.csv("CANCER_CENSUS_GENES.v1.2.csv")

# Remove duplicates and functional consequence column to get a unique set
domain1 <- domain1 %>% select(-(alteration))
domain1 <- unique(domain1)

# Check and examine duplicated IDs (note that some genes appear duplicated since they exist on alt haplotypes; they will have different gene_ID)
sum(duplicated(domain1$gene_name)) # 1 (one gene exists on a contig as well, with different gene/transcripts ID, hence double entry)
sum(duplicated(domain1$transcript_ID))  # 0
sum(duplicated(domain1$gene_ID))  # 0
sum(duplicated(domain2$gene_name)) # 46
sum(duplicated(domain2$transcript_ID))  # 0
sum(duplicated(domain2$gene_ID_ID))  # 0

# Read the gene/transcript annotations from the VCF file (HPC)
genesFromVCF <- function(patientID){
  
  # Get FF and FFPE VCF paths for specified patient ID and read into a Large CollapsedVCF object (VariantAnnotation package)
  # Need to read QC_portal_trios separately, clean and add VCF paths
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
  # Remove filtered entries, keep Manta only, keep only ID and CSQT
  Manta_ff <- SVinfo_ff %>% filter(FILTER == "PASS", Application == "Manta") %>% dplyr::select(ID, CSQT)
  
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
  # Remove filtered entries, keep Manta only, keep only ID and CSQT
  Manta_ffpe <- SVinfo_ffpe %>% filter(FILTER == "PASS", Application == "Manta") %>% dplyr::select(ID, CSQT)
  
  # Merge FF and FFPE into one table, keep info on SV source (FF or FFPE?)
  Manta_ff$FF <- "FF"
  Manta_ffpe$FF <- "FFPE"
  ff_ffpe_merged <- rbind(Manta_ff, Manta_ffpe)
  
  # Add PATIENT_ID
  ff_ffpe_merged$PATIENT_ID <- patientID
  
  # Write out the R object
  saveObject(ff_ffpe_merged, file = paste0(patientID, "_SV_genes", ".RData"))
  #save(ff_ffpe_merged, file = paste0(patientID, "_SV_genes", ".RData"))
  
}

# Run the function
patientIDs <- unique(QC_portal_trios$PATIENT_ID)

# Run SV function for each patient ID
sapply(patientIDs, genesFromVCF)

# Get all data together
result2 <- data.frame()
result2 <- lapply(patientIDs, function(x){
  table <- loadObject(paste0(x, "_SV_genes", ".RData"))
  result2 <- rbind(result2, table)
})

# Merge into one data frame (HPC)
result2 <- merge_recurse(result2)

# Write out the result as an R object (HPC)
saveObject(result2, file = paste0(today, "_SV_genes_all", ".RData"))

# Read the table with genes
trio_genes <- loadObject("2017-02-28_SV_genes_all.RData")

# Add gene info to unfiltered SV results
# NOTE that some SVs share IDs across patients; best to join using ID AND PATIENT_ID
sum(duplicated(trio_genes$ID)) # 13
sum(duplicated(result$ID))  # 11
result <- left_join(result, trio_genes) # this will join by all common columns!

# Add the Domain1 and Domain2 transcript counts to SVs 
# (from annotation in the VCF; note that large insertion VCF entries will have ALL genes listed there)
result$CSQT <- as.character(result$CSQT )
result[result$CSQT == "character(0)",]$CSQT <- ""
# Domain 1
result$Domain1_count <- sapply(1:dim(result)[1], function(x){
  sum(sapply(1:dim(domain1)[1], function(y){
    grepl(domain1$transcript_ID[y], result$CSQT[x])
  }))
})
# Domain 2
result$Domain2_count <- sapply(1:dim(result)[1], function(x){
  sum(sapply(1:dim(domain2)[1], function(y){
    grepl(domain2$transcript_ID[y], result$CSQT[x])
  }))
})

# # Add Domain1 and Domain2 gene names to the results table --- unfinished code, leaving it for now
# result$Domain1_genes <- sapply(1:dim(result)[1], function(x){
#   as.character(domain1$gene_name[sapply(1:dim(domain1)[1], function(y){
#     grepl(domain1$transcript_ID[y], result$CSQT[x])
#   })])
# })
# result$Domain2_genes <- sapply(1:dim(result)[1], function(x){
#   as.character(domain2$gene_name[sapply(1:dim(domain2)[1], function(y){
#     grepl(domain2$transcript_ID[y], result$CSQT[x])
#   })])
# })


### Overiew of Domain 1,2 overlap in filtered-out SVs
# Add duplicated flag to make it easier to plot and count unique SV (now FF and FFPE concordant are listed both)
result$DuplicatedSV <- as.numeric(duplicated(result$KEY))
# Overview of Domain 1 overlap
table(result[result$DuplicatedSV == 0,]$Domain1_count != 0, result[result$DuplicatedSV == 0,]$SVTYPE)
# Overview of Domain 2 overlap
table(result[result$DuplicatedSV == 0,]$Domain2_count != 0, result[result$DuplicatedSV == 0,]$SVTYPE)

### Count the unfiltered SVs if only FF are taken into account (count BND as half!!!)
table(result[((result$DuplicatedSV == 0) & (result$FF == "FF")),]$Domain1_count != 0, result[((result$DuplicatedSV == 0) & (result$FF == "FF")),]$SVTYPE)
table(result[((result$DuplicatedSV == 0) & (result$FF == "FF")),]$Domain2_count != 0, result[((result$DuplicatedSV == 0) & (result$FF == "FF")),]$SVTYPE)


### Overlap of filtered SVs with Domain1/2 genes
sum(duplicated(result_filt$ID))  # 0
sum(duplicated(result[result$FILTERED == 0,]$ID))  #0
result_filt <- left_join(result_filt, (result %>% filter(FILTERED == 0) %>% dplyr::select(PATIENT_ID, ID, Domain1_count, Domain2_count)))
# Add duplicated flag to make it easier to plot and count unique SV (now FF and FFPE concordant are listed both)
result_filt$DuplicatedSV <- as.numeric(duplicated(result_filt$KEY))

### Overview of Domain1,2 overlap in filtered SVs
result_filt %>% filter(Domain2_count != 0, DuplicatedSV == 0) %>% dplyr::select(KEY, SVLEN, PATIENT_ID, TumorType, Domain1_count, Domain2_count)
# Overview of Domain 1 overlap
table(result_filt[result_filt$DuplicatedSV == 0,]$Domain1_count != 0, result_filt[result_filt$DuplicatedSV == 0,]$SVTYPE)
# Overview of Domain 2 overlap
table(result_filt[result_filt$DuplicatedSV == 0,]$Domain2_count != 0, result_filt[result_filt$DuplicatedSV == 0,]$SVTYPE)


####### Write out SV results table ####### 

result$MATE_KEY <- as.character(result$MATE_KEY)
write.csv(result, file = paste0("FFPE_trios_unfiltered_SVs_", today, ".csv"), row.names = F, quote = F)
write.table((result_filt[(!duplicated(result_filt$KEY)),] %>% dplyr::select(PATIENT_ID, TumorType, FF, KEY, MATE_KEY, ID, SVLEN, CIGAR, HOMLEN, LEFT_SVINSSEQ, RIGHT_SVINSSEQ, SOMATICSCORE, PR_ALT, SR_ALT, CONCORDANT, Domain1_count, Domain2_count)), quote = F, row.names = F, sep = "\t", file = "HighConf_SVs_FFPE_trios.tsv")

