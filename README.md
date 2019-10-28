# ddRAD_Analysis
This repository contains shell scripts used for ddRAD analysis.
The workflow is broken up into 5 sections, with each breakpoint used for validation or to provide modularity.
Required programs: fastqc Trimmomatic bwa samtools picard-tools GATK v4
Goal of this study is to explore chickadee hybridization relative to metrics of human disturbance.

## Bioinformatics  Objectives:
1. Generate Hybrid Indeces for all invididuals 
2. Calculate Heterozygosity for all individuals 

### Step 1: Demultiplex using script
Script for demultiplexing ddRAD data from Illumina using STACKS
Output files named with sample ID and placed in provided output directory

Required components for script:

1. List of index primers and samples
2. List of illumina adapters for index group
 
 
### Step 2: Trim and QC
 
 
#Quality control
#Remove PhiX sequences
#Adapter trimming
#Quality trimming of reads
#Quality assessment
 
### Step 3: Allign to BCCH and MOCH reference genomes
 
### Step 4: Process allignments for variant calling
 
### Step 5: Variant calling
 
### Step 6: Filter variants
 
### Step 7: Generate hybrid indeces for all individauls
 
### Step 8: Calculate heterozygosity for all individauls
 
### Step 9: Graph Hybrid index v. heterozygosity


