##********************************
## Hannah Glowacki - Assignment 4
##
## BINF*6210
##
## Topic: Microbiome of the Colon in Colon Cancer Patients
## 
## Topic: Microbiome of the Colon in Colon Cancer Patients
## 
## #Q1: Do tumor tissues show different phylogenetic diversity than adjacent normal tissues? #Q2 (Secondary): Are patterns consistent across both studies?
##********************************

##_Packeges Used ------
library(dada2)
library(DECIPHER)
library(phyloseq)

library(picante)
library(vegan)
library(ggplot2)
library(effectsize)
library(ape)

set.seed(123)

#Used  Sequence Read Archive (SRA) toolkit to downloads the FASTQ DNA sequence files (zip files located in fastq folder) from the Sequence Read Archive (SRA) on NCBI website. Before downloading, files were screened for 16S sequence data (V3-V4 region) from biopsies of colon cancer patients. To ensure sequence quality and validity, sequence data was only considered if this was in an already published studies within 2020-now time frame (see storyboard for details). Found 2 studies meeting the above criteria.

#For this project: Randomly taken 10 sequences; 5 from tumor biopsy and 5 from normal adjecent tissue, from each study (total of 20 sequences). 


#_Step 1: Quality Control using DADA2 -----
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


#Calculate percentage of reads retained - good retention 
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
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Clean up to save space
rm(derepFs, derepRs, dadaFs, dadaRs, merger, seqtab)
gc()

#_Step 2: Assign Taxonomy with Decipher -----

#Downloaded SILVA_SSU_r138_2_2024.RData from: https://www2.decipher.codes/Downloads.html

# Load the training set
silva_path <- "R/SILVA_SSU_r138_2_2024.RData"
load(silva_path)

# Assign taxonomy
dna <- DNAStringSet(getSequences(seqtab.nochim))
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=TRUE)

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

#_Step 3: Create Phylogenetic Tree (Neighbor joining only) using Decipher & phangorn -----
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs

# Quick alignment
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = FALSE)

# Fast tree (NJ only, skip optimization)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
tree <- NJ(dm)

# Clean up
rm(alignment, phang.align, dm)
gc

#_Step 4: Create Metadata & Phyloseq -----

#Check order sample names (back in Step 1: Using Quality Control using DADA2)
print(sample.names)

# Create METADATA - adjusted to match sample order
metadata <- data.frame(
  sample_id = sample.names,
  tissue_type = rep(c(rep("tumor", 5), rep("normal", 5)), 2),  # Adjust!
  study_id = c(rep("study1", 10), rep("study2", 10))          # Adjust!
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

saveRDS(ps, "phyloseq_object.rds")
