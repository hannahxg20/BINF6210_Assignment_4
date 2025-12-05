# **BINF6210_Assignment_4**

Conducting a small *exploratory* project that related to my field of interest: the Gut Microbiome! Specifically, my topic will be looking at the microbiome in colorectal cancer (CRC) patients. 

**Brief Background:** 
For this project I randomly took 10 samples - 5 from tumor biopsies and 5 from normal adjecent tissues - from 2 published studies (20 samples total). I used Sequence Read Archive (SRA) Toolkit to downloads the FASTQ DNA sequence files (zip files located in fastq folder) from the Sequence Read Archive (SRA) on NCBI website. Before downloading, I screened the files to ensure they contained 16S rRNA sequencing data (V3–V4 region) from colon cancer biopsy samples. To ensure data quality and reliability, I only included studies published between 2020 and the present (see storyboard for details). 2 studies met the above criteria.

Through this project I hope to explore the following research questions: 
  
#Q1:**Within the colorectal cancer (CRC) environment, do tumor tissues show different biodiversity compared to adjacent normal tissues?**
#Q2 (Secondary): **Are these patterns consistent across both studies?**

**Important Note**
It takes about 4 hrs to run the entire script from beginning to end (Steps 1: Quality Control using DADA2 - Step 6: Statistical Testing & Visualization). To account for this, I have created "RDS_objects" folder, containing all the files needed to move throughout the script to help avoid the long wait time and make it easier to start/stop at any step. In my script, I have created "Save Checkpoints" at the end of each step for RDS files that are required for the following steps. I also wrote code at the beginning of each step to load required RDS files if they are not already in memory. 
