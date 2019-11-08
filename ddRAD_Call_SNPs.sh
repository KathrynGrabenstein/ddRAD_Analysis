#!/bin/bash

#call SNPs in samtools and filter variants

ref="/data0/chickadeeddRAD_workdir/BCCH_pseudochr_v2.fasta"
bamdir="/data0/chickadeeddRAD_workdir/sorted_bam_files/"

ID="ddRAD_chickadees" # This will be used as a prefix for the output file


echo "Making a pileup file for "$ID >> call_filter_log.txt

#pipe mpileup to call to avoid creating intermediate bcf file
#Call genotypes from mpileup results
bcftools mpileup -Ou -f $ref --ignore-RG -a AD,ADF,DP,SP,INFO/AD,INFO/ADF \
"$bamdir"*.bam | bcftools call -mv > "$ID"_raw_variants.vcf
