##********************************
## Hannah Glowacki - Assignment 4
##
## BINF*6210
##
## Topic: Microbiome of the Colon in Colon Cancer Patients
## 
## #Q1: Within the colorectal cancer (CRC) environment, do tumor tissues show different biodiversity compared to adjacent normal tissues? #Q2 (Secondary): Are patterns consistent across both studies?
##********************************

##_Packeges Used ------

##_Original script Used -----

library(dada2)
library(DECIPHER)
library(phangorn)

library(phyloseq)
library(picante)
library(vegan)
library(ggplot2)
library(effectsize)
library(ape)

set.seed(123)

#Used  Sequence Read Archive (SRA) toolkit to downloads the FASTQ DNA sequence files (zip files located in fastq folder) from the Sequence Read Archive (SRA) on NCBI website. Before downloading, files were screened for 16S sequence data (V3-V4 region) from biopsies of colon cancer patients. To ensure sequence quality and validity, sequence data was only considered if this was in an already published studies within 2020-now time frame (see storyboard for details). Found 2 studies meeting the above criteria.

#For this project: Randomly taken 10 sequences; 5 from tumor biopsy and 5 from normal adjecent tissue, from each study (total of 20 sequences). 


#_Step 1: Quality Control using DADA2 (2 hrs to run) -----
# Set working directory to BINF_Assignment_4, then define path to raw zipped FASTQ files 
path <- "FASTQ_Data/"
list.files(path)

fnFs <- sort(list.files(path = path, pattern = "_1.fastq.gz", full.names = TRUE))

fnRs <- sort(list.files(path = path, pattern = "_2.fastq.gz", full.names = TRUE))

# Extract SRA sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspect read quality profiles
plotQualityProfile(fnFs[c(1,6,11,16)])
plotQualityProfile(fnRs[c(1,6,11,16)])

## Filter & Trim (based on quality profiles)
#Forward reads: Look good - will truncate forward reads at 280 bp
#Reverse reads: Worse then forward (common in Illumina sequencing) - will truncate reverse reads at 220 bp.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

#Name the filtered files
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Trim
#V3-V4 is typically a ~460bp amplicon: Overlap = 280 + 220 - 460 = 40bp overlap (minimum acceptable)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280, 220), trimLeft=c(10, 10), maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE, compress=TRUE, verbose = TRUE)
head(out)


#Calculate percentage of reads retained - good retention between 80-90%
out <- as.data.frame(out)
out$percent_retained <- (out$reads.out / out$reads.in) * 100
out

## Dereplication 
#collapsing sequences that are exactly the same.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the dereplicated objects
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Estimate error parameters (Error Rates)
set.seed(123)#Match seed at the top
errFs <- learnErrors(filtFs, multithread = FALSE)
errRs <- learnErrors(filtRs, multithread = FALSE)

#Sanity Check
plotErrors(errFs, nominalQ=TRUE)

## Infer sample composition
dadaFs <- dada(derepFs, err = errFs, multithread = FALSE)
dadaRs <- dada(derepRs, err = errRs, multithread = FALSE)
print(dadaFs)
print(dadaRs)

## Merging paired reads 
merger <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

## Construct sequence table 
seqtab <- makeSequenceTable(merger)
dim(seqtab)

# Inspect distribution of sequence lengths
sort(table(nchar(getSequences(seqtab))), decreasing = TRUE)
#Top sequences by abundance: 440bp: 3,565 sequences
#Distribution is good (tight clustering around 440-446bp)! This is exactly where V3-V4 amplicons should be = will keep my current parameters of  truncLen=c(280, 220)

#Successful merging - instead of anticipated 40 bp overlap, there is 280 + 220 - 440 = 60 bp overlap (even better) 


## Remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) 
#Identified 12373 bimeras out of 14834 input sequences.
#12,373 chimeras out of 14,834 sequences = 83% chimeric sequences

dim(seqtab.nochim) #Out of 20 samples → 2,461 unique ASVs 

sum(seqtab.nochim)/sum(seqtab) #[1] 0.7611049 retained 76% of reads


## Track reads through the pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger, getN), rowSums(seqtab.nochim))

#Make output readable
colnames(track) <- c("input", "filtered", "percent_retained", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Clean up to save space
rm(derepFs, derepRs, dadaFs, dadaRs, merger, seqtab)
gc()

#----- SAVE CHECKPOINT - DADA2 COMPLETE -----

# Create directory for saved objects if it doesn't exist
if(!dir.exists("RDS_objects")) dir.create("RDS_objects")

message("Saving DADA2 pipeline results...")

# Save essential objects for next steps
saveRDS(seqtab.nochim, "RDS_objects/seqtab_nochim.rds")
saveRDS(track, "RDS_objects/track.rds")
saveRDS(sample.names, "RDS_objects/sample_names.rds")

# Save the ASV sequences separately (needed for taxonomy assignment)
asv_seqs <- colnames(seqtab.nochim)
saveRDS(asv_seqs, "RDS_objects/asv_sequences.rds")

# Save filtering stats for reference
saveRDS(out, "RDS_objects/filtering_stats.rds")

#_Step 2: Assign Taxonomy with Decipher -----
#Start up to skip 2 hr DADA2 pipeline in Step 1.
seqtab.nochim <- readRDS("RDS_objects/seqtab_nochim.rds")
sample.names <- readRDS("RDS_objects/sample_names.rds")
asv_seqs <- readRDS("RDS_objects/asv_sequences.rds")
track <- readRDS("RDS_objects/track.rds")


#Downloaded SILVA_SSU_r138_2_2024.RData from: https://www2.decipher.codes/Downloads.html

# Load the training set
silva_path <- "SILVA/SILVA_SSU_r138_2_2024.RData"
load(silva_path)

# Assign taxonomy
dna <- DNAStringSet(getSequences(seqtab.nochim)) #Takes 10 mins
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=TRUE) #Takes 35 mins

#Check - IdTaxa ran correctly 
ids[[1]]
#everything appears to be working - Taxonomy assigned, output structure is correct

# Extract taxonomy
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
taxa <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxa) <- getSequences(seqtab.nochim)

# Check results
head(taxa)
dim(taxa)  # Should be [2461 rows x 7 columns]

# Clean up
rm(trainingSet, ids, dna)
gc()

#----- SAVE CHECKPOINT - TAXONOMY ASSIGNMENT COMPLETE -----

# Save taxonomy table
saveRDS(taxa, "RDS_objects/taxa.rds")

#_Step 3: Create Phylogenetic Tree (Neighbor joining only) using Decipher & phangorn -----
#Load required objects if not already in memory 
seqtab.nochim <- readRDS("RDS_objects/seqtab_nochim.rds")
sample.names <- readRDS("RDS_objects/sample_names.rds")
track <- readRDS("RDS_objects/track.rds")
taxa <- readRDS("RDS_objects/taxa.rds")

#Extract ASV sequences
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs

# Quick alignment - takes 20 mins
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = TRUE)

# Fast tree using phangorn package (NJ only, skip optimization) this takes about 10 mins
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
tree <- NJ(dm) #unrooted phylogenetic tree with 2461 tips and 2459 internal nodes.

# Clean up
rm(alignment, phang.align, dm)
gc

#----- SAVE CHECKPOINT - PHYLOGENETIC TREE COMPLETE -----

saveRDS(tree, "RDS_objects/tree.rds")

#_Step 4: Create Metadata & Phyloseq -----

# Load required objects
seqtab.nochim <- readRDS("RDS_objects/seqtab_nochim.rds")
taxa <- readRDS("RDS_objects/taxa.rds")
tree <- readRDS("RDS_objects/tree.rds")
sample.names <- readRDS("RDS_objects/sample_names.rds")

#Check order sample names (back in Step 1: Using Quality Control using DADA2) & cross reference with tissue/study using "SraAccessionList" excel in the FASTQ_Data folder. 
print(sample.names)

# Create metadata - MANUALLY MATCHED to your print(sample.names) output
metadata <- data.frame(
  sample_id = sample.names,
  tissue_type = c(
    # Study 2 samples (positions 1-10)
    "tumor",   # [1] SRR25305896
    "tumor",   # [2] SRR25305897
    "tumor",   # [3] SRR25305898
    "tumor",   # [4] SRR25305899
    "normal",  # [5] SRR25305917
    "normal",  # [6] SRR25305918
    "normal",  # [7] SRR25305919
    "normal",  # [8] SRR25305920
    "normal",  # [9] SRR25305921
    "tumor",   # [10] SRR25374795
    # Study 1 samples (positions 11-20)
    "normal",  # [11] SRR32497531
    "normal",  # [12] SRR32497532
    "normal",  # [13] SRR32497533
    "normal",  # [14] SRR32497534
    "normal",  # [15] SRR32497535
    "tumor",   # [16] SRR32497540
    "tumor",   # [17] SRR32497541
    "tumor",   # [18] SRR32497542
    "tumor",   # [19] SRR32497543
    "tumor"    # [20] SRR32497544
  ),
  study_id = c(
    rep("study2", 10),  # First 10 samples are Study 2
    rep("study1", 10)   # Last 10 samples are Study 1
  ),
  stringsAsFactors = FALSE
)
rownames(metadata) <- sample.names

print(table(metadata$tissue_type, metadata$study_id))

# Build phyloseq
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  phy_tree(tree),
  sample_data(metadata)
)

print(ps)

#----- SAVE CHECKPOINT - Metadata COMPLETE -----

saveRDS(metadata, "RDS_objects/metadata.rds")
saveRDS(ps, "RDS_objects/ps.rds")

#_Step 5: Calculate Diversity metrics using phyloseq, vegan, picante-----
#Load required fields
metadata <- readRDS("RDS_objects/metadata.rds")
ps <- readRDS("RDS_objects/ps.rds")

# Extract components from phyloseq
otu_table_ps <- otu_table(ps)
tree_obj <- phy_tree(ps)
sample_meta <- sample_data(ps)

# Prepare OTU matrix (samples as rows, ASVs as columns)
otu_mat <- as.matrix(otu_table_ps)

message(paste("OTU matrix dimensions:", nrow(otu_mat), "samples x", ncol(otu_mat), "ASVs")) # OTU matrix dimensions: 20 samples x 2461 ASVs
message(paste("Tree has", length(tree_obj$tip.label), "tips")) # Tree has 2461 tips

# Verify alignment between OTU table and tree
if(!all(colnames(otu_mat) %in% tree_obj$tip.label)) {
  stop("ERROR: Some ASVs in OTU table are not in the phylogenetic tree!")
}

message("✓ OTU table and tree are properly aligned")
#✓ OTU table and tree are properly aligned

#Calculate Alpha Diversity Metrics ----

# 1 Species Richness - using vegan  
richness <- vegan::specnumber(otu_mat)

# 2 Shannon Diversity - using vegan
shannon <- vegan::diversity(otu_mat, index = "shannon")

# 3 Faith's Phylogenetic Diversity - using picante
library(picante)
faith_pd_result <- pd(otu_mat, tree_obj, include.root = FALSE)
faith_pd <- faith_pd_result$PD

## Verify all vectors have same length (20 samples)
message(paste("Richness length:", length(richness))) #20
message(paste("Shannon length:", length(shannon))) #20
message(paste("Faith's PD length:", length(faith_pd))) #20

# Check if all are length 20
if(!all(c(length(richness), length(shannon), length(faith_pd)) == 20)) {
  stop("ERROR: Not all diversity metrics have 20 values!")
}

## Combining ALL Results 
results <- data.frame(
  sample_id = rownames(otu_mat),
  richness = richness,
  shannon = shannon,
  faith_pd = faith_pd,
  stringsAsFactors = FALSE
)

# Convert sample_data to regular data frame and merge
metadata_df <- data.frame(sample_meta)
metadata_df$sample_id <- rownames(metadata_df)
results <- merge(results, metadata_df, by = "sample_id")

# Reorder columns for clarity
results <- results[, c("sample_id", "tissue_type", "study_id", "richness", "shannon", "faith_pd")]

print(head(results))

### Summary of Diversity Metrics ------
#Mean diversity metrics by tissue type
print(aggregate(cbind(richness, shannon, faith_pd) ~ tissue_type, data = results, FUN = mean))

#Mean diversity metrics by study
print(aggregate(cbind(richness, shannon, faith_pd) ~ study_id, data = results, FUN = mean))

#----- SAVE CHECKPOINT - DIVERSITY METRICS COMPLETE -----
saveRDS(results, "RDS_objects/diversity_metrics_results.rds")

#_Step 6: Statistics & Data Visualization-----
#For this section will use: 
#library(ggplot2)
#library(vegan)
#library(phyloseq)
theme_set(theme_bw())

# Load your diversity metrics if not already in memory
results <- readRDS("RDS_objects/diversity_metrics_results.rds")
ps <- readRDS("RDS_objects/ps.rds")  # If needed for ordination plots

print(head(results))

# Part A: Statistical Tests -----
#Q1: Wilcoxon Test: Tumor vs Normal (for Q1) ----
#"Are tumor and normal biopsy tissues different?" 
#Use Wilcoxon because we only have 10 samples per group (small sample size)

#Test Richness
richness_test <- wilcox.test(richness ~ tissue_type, data = results)
cat("  p-value =", richness_test$p.value, "\n")

# p-value (0.7912601) > 0.05 
# Result: Not significantly different 

#Test Shannon
shannon_test <- wilcox.test(shannon ~ tissue_type, data = results)
cat("  p-value =", shannon_test$p.value, "\n")

# p-value (0.9117972) > 0.05
# Result: Not significantly different 

# Test Faith's PD
faithpd_test <- wilcox.test(faith_pd ~ tissue_type, data = results)
cat("  p-value =", faithpd_test$p.value, "\n")

# p-value (0.9705125 )) > 0.05
# Result: Not significantly different 

#Q2: Consistency across 2 studies -----
#"Do both studies show the same pattern?"

# Split by study and compare
study1_data <- subset(results, study_id == "study1")
study2_data <- subset(results, study_id == "study2")

# Study 1: Tumor vs Normal
study1_rich <- wilcox.test(richness ~ tissue_type, data = study1_data)
cat("  p-value =", study1_rich$p.value, "\n") 
# p value = 0.42

#Study 2: Tumor vs Normal
study2_rich <- wilcox.test(richness ~ tissue_type, data = study2_data)
cat("  p-value =", study2_rich$p.value, "\n")
# p value = 0.84

# Calculate mean differences to see direction
#Study 1 mean richness 
cat("Study 1 mean richness: Tumor =", mean(study1_data$richness[study1_data$tissue_type == "tumor"]),
    "vs Normal =", mean(study1_data$richness[study1_data$tissue_type == "normal"]), "\n")
#Tumor = 222 vs Normal = 265.8 

#Study 2 mean richness 
cat("Study 2 mean richness: Tumor =", mean(study2_data$richness[study2_data$tissue_type == "tumor"]),
    "vs Normal =", mean(study2_data$richness[study2_data$tissue_type == "normal"]), "\n")
#Tumor = 173.2 vs Normal = 164.2 

# Part B: Create 5 figures -----
## Figure 1: Richness by Tissue Type ----
fig1 <- ggplot(results, aes(x = tissue_type, y = richness, fill = tissue_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("normal" = "#4ECDC4", "tumor" = "#FF6B6B")) +
  labs(
    title = "Bacterial Richness: Tumor vs Normal Tissue",
    subtitle = paste("p =", round(richness_test$p.value, 4)),
    x = "Tissue Type",
    y = "Species Richness (# of unique taxa)"
  ) +
  theme(legend.position = "none")

print(fig1)
ggsave("Figure1_Richness_TissueType.png", fig1, path = "Figures/", width = 6, height = 5, dpi = 300)


## Figure 2: Shannon Diversity by Tissue Type ----
fig2 <- ggplot(results, aes(x = tissue_type, y = shannon, fill = tissue_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("normal" = "#4ECDC4", "tumor" = "#FF6B6B")) +
  labs(
    title = "Shannon Diversity: Tumor vs Normal Tissue",
    subtitle = paste("p =", round(shannon_test$p.value, 4)),
    x = "Tissue Type",
    y = "Shannon Diversity Index"
  ) +
  theme(legend.position = "none")

print(fig2)
ggsave("Figure2_Shannon_TissueType.png", fig2, path = "Figures/", width = 6, height = 5, dpi = 300)

## Figure 3: Faith's PD by Tissue Type ----
fig3 <- ggplot(results, aes(x = tissue_type, y = faith_pd, fill = tissue_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("normal" = "#4ECDC4", "tumor" = "#FF6B6B")) +
  labs(
    title = "Faith's Phylogenetic Diversity: Tumor vs Normal",
    subtitle = paste("p =", round(faithpd_test$p.value, 4)),
    x = "Tissue Type",
    y = "Faith's PD (phylogenetic diversity)"
  ) +
  theme(legend.position = "none")

print(fig3)
ggsave("Figure3_FaithPD_TissueType.png", fig3, path = "Figures/", width = 6, height = 5, dpi = 300)

## Figure 4: Richness by Study (Grouped Boxplot) ----
# This shows if both studies have the same pattern
fig4 <- ggplot(results, aes(x = study_id, y = richness, fill = tissue_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("normal" = "#4ECDC4", "tumor" = "#FF6B6B")) +
  labs(
    title = "Richness Consistency Across Studies",
    subtitle = "Checking if both studies show the same tumor vs normal pattern",
    x = "Study",
    y = "Species Richness",
    fill = "Tissue Type"
  ) +
  theme(legend.position = "right")

print(fig4)
ggsave("Figure4_Richness_ByStudy.png", fig4, path = "Figures/", width = 7, height = 5, dpi = 300)

## Figure 5: NMDS Ordination (Community Composition) ----
# This shows if tumor and normal samples cluster separately based on their bacteria

# Extract OTU table
otu_mat <- as.matrix(otu_table(ps))

# Run NMDS (this looks at overall bacterial community differences)
set.seed(123)
nmds <- metaMDS(otu_mat, distance = "bray", k = 2, trymax = 100)

# Extract NMDS coordinates
nmds_points <- data.frame(
  NMDS1 = nmds$points[,1],
  NMDS2 = nmds$points[,2],
  sample_id = rownames(nmds$points)
)

# Add metadata
nmds_data <- merge(nmds_points, results, by = "sample_id")

# Create plot
fig5 <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2, 
                              color = tissue_type, shape = study_id)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(aes(group = tissue_type), level = 0.95, linetype = 2) +
  scale_color_manual(values = c("normal" = "#4ECDC4", "tumor" = "#FF6B6B")) +
  labs(
    title = "Bacterial Community Clustering (NMDS)",
    subtitle = paste("Stress =", round(nmds$stress, 3)),
    color = "Tissue Type",
    shape = "Study"
  ) +
  theme(legend.position = "right")

print(fig5)
ggsave("Figure5_NMDS_Communities.png", fig5, path = "Figures/", width = 7, height = 5, dpi = 300)

# Test if communities are significantly different (PERMANOVA)
# This asks: "Are the bacterial communities really different between tumor and normal?"
permanova <- adonis2(otu_mat ~ tissue_type, data = results, method = "bray", permutations = 999)
print(permanova)

cat("\n✓ All figures saved successfully!\n")
cat("✓ Statistical testing complete!\n")
