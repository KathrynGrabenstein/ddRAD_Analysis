#!/bin/bash

# shell script to assemble a bam file, then sort, mark duplicates, and index

if [ $# -lt 1 ]
  then
    echo "Aligns fastq reads to a reference genome using bwa mem.
    Bam file is then sorted, duplicates are marked, and file is indexed using
    samtools and picard-tools.
    [-i] Sample list
    [-r] Reference genome
    [-p] Path to trimmed fastqs - the default is a directory called 'fastqs' as
         produced from the initial sorting
    OPTIONAL ARGUMENTS
    [-t] Number of threads to use
    [-b] Output directory for bam files - default is to make a directory
         called 'bam_files'
    [-s] Output directory for sorted bam files - default is to make a
         directory called 'sorted_bam_files'"
         
else
  while getopts i:r:t:p:b:s: option
  do
  case "${option}"
  in
  i) seqs=${OPTARG};;
  r) ref=${OPTARG};;
  t) threads=${OPTARG};;
  p) fastqs_path=${OPTARG};;
  b) bamoutdir=${OPTARG};;
  s) sortbamoutdir=${OPTARG};;
  esac
  done
  
  threads="${threads:-1}"
  fastqs_path="${fastqs_path:-fastqs/}"
  bamoutdir="${bamoutdir:-bam_files/}"
  sortedbamoutdir="${sortedbamoutdir:-sorted_bam_files/}"

    if [ $sortedbamoutdir == sorted_bam_files/ ]
      then
        mkdir sorted_bam_files
    fi

    #create sorted bam files
    while read -r ID; do
    echo "Sorting BAM file "$ID >> sort_markdup_log.txt

    # -@ is number of threads, here 6
    # -T temp, write temporary files with TEMP/xxxxx
    #output file location and type
    #input file
    samtools sort -T temp -@ 6 \
    -o "$sortedbamoutdir"/"$ID"_sorted.bam \
    "$bamoutdir"/"$ID".bam

    echo "Adding replacing reading groups for "$ID >> sort_markdup_log.txt
    #input BAM file
    #output BAM file
    picard-tools AddOrReplaceReadGroups \
    I="$sortedbamoutdir"/"$ID"_sorted.bam \
    O="$sortedbamoutdir"/"$ID"_sorted_RGadded.bam \
    RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$ID

    echo "Marking duplicates for "$ID >> sort_markdup_log.txt
    picard-tools MarkDuplicates \
    I="$sortedbamoutdir"/"$ID"_sorted_RGadded.bam \
    O="$sortedbamoutdir"/"$ID"_sorted_RGadded_dupmarked.bam \
    M="$sortedbamoutdir"/"$ID"_dupmarked_metrics.txt

    #create samfile index
    samtools index \
    "$sortedbamoutdir"/"$ID"_sorted_RGadded_dupmarked.bam

    #remove intermediate files
    #keep only sorted, reading group added, duplicates marked BAM file
    rm "$sortedbamoutdir"/"$ID"_sorted.bam
    rm "$sortedbamoutdir"/"$ID"_sorted_RGadded.bam

    done<"$seqs"
fi
