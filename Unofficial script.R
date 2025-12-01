##********************************
## Hannah Glowacki - Assignment 4
##
## BINF*6210
##
## Topic: Microbiome of the Colon in Colon Cancer Patients
## 
## Research Q1: In patients with colon cancer, do tumor samples show reduced phylogenetic diversity compared to adjacent normal tissue?
## Research Q2: In patients with colon cancer, does the magnitude of diversity change correlate across different diversity metrics?
##********************************

##_Packeges Used ------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)

library(rentrez)
library(seqinr)

library(dada2)
library(ShortRead)
library(Biostrings)

library(phyloseq)

#Used  Sequence Read Archive (SRA) toolkit to downloads the FASTQ DNA sequence files (zip files located in fastq folder) from the Sequence Read Archive (SRA) on NCBI website. Before downloading, files were screened for 16S sequence data (V3-V4 region) from biopsies of colon cancer patients. To ensure sequence quality and validity, sequence data was only considered if this was in an already published studies within 2020-now time frame (see storyboard for details). Found 2 studies meeting the above criteria.

#For this project: Randomly taken 10 sequences; 5 from tumor biopsy and 5 from normal adjecent tissue, from each study (total of 20 sequences). 


#_Step 1: Quality Control using DADA2 -----
# Set working directory to BINF_Assignment_4, then define path to raw zipped FASTQ files 
path <- "Fastq/"
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
set.seed(108)
derepFs <- derepFastq(filtFs, verbose=TRUE)

derepRs <- derepFastq(filtRs, verbose=TRUE)

class(derepFs)

derepFs

# Name the dereplicated objects
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Estimate error parameters (Error Rates)
errFs <- learnErrors(derepFs, multithread = FALSE)
errRs <- learnErrors(derepRs, multithread = FALSE)

#Sanity Check
plotErrors(errFs, nominalQ=TRUE)

## Infer sample composition
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
print(dadaFs)
print(dadaRs)

## Merging paired reads 
merger <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Construct sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## Remove chimeras 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

## Track reads through the pipeline 

#_

#In order to import data from NCBI, perform quality control and assign taxonomy in dada2 package- We need to install Sequence Read Archive (SRA) toolkit. 
#This can be done several ways (e.g. use SRAdb Bioconductor package in R  - however this process will require massive amount of disk storage, 35 GB, to download SRA database file). 
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("SRAdb")

#library(SRAdb)

#We will install SRA toolkit outside of R, on the computer. See folder called "Software" for zip file (windows computer). Downloaded November 26, 2025 from https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit#ncbi-sra-toolkit

#Add path - so your computer knows where to find it
#"Edit the system environment variables" - "Environment Variables" - "New" - "Copy & paste the location of file" 

#Step 2 on GitHub: Open command prompt (cdm) on computer 
#write "cd" then paste location of file [enter]
#write "cd bin" then [enter]
#test that the toolkit is functional 

#Step 3 on GitHub:Toolkit configuration 


