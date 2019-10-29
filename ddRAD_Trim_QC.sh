#!/bin/bash

#script for pre-trim QC, trim, and then post-trim QC
#pre-trim and post-trim fastqc files go in separate directories


if [ $# -lt 1 ]
  then
    echo "Runs a pre- and post-trim qc using fastqc for samples provided by an input file. Both a pre- and post-trim directory will be created. Then trims sequences using TrimmomaticPE. Changes to trimming parameters should be made if necessary.
    Following trimming, a post-trim qc will also be run using fastqc.
    [-i] Sample list - sample name should match the prefix of the fastq files
    [-p] Path to fastq files
    [-a] File with adapter sequences used by Trimmomatic.
    OPTIONAL ARGUMENTS
    [-t] Number of threads to use"

  else
    while getopts i:p:a:t: option
    do
    case "${option}"
    in
    i) filename=${OPTARG};;
    p) fastqs_path=${OPTARG};;
    a) adapters=${OPTARG};;
    t) threads=${OPTARG};;
    esac
    done

    threads="${threads:-1}"

    echo "Trim and QC for "$filename >> trim_and_QC_log.txt
    mkdir pre_trim_QC_files
    mkdir post_trim_QC_files

    while read -r ID; do
    echo "Beginning pre-trim QC for "$ID >> trim_and_QC_log.txt
    fastqc -t $threads \
    "$fastqs_path"/"$ID".fq \
    --outdir pre_trim_QC_files/
    echo $ID" pre-trim QC done" >> trim_and_QC_log.txt

    echo "Beginning trimming for "$ID >> trim_and_QC_log.txt
    TrimmomaticSE -threads $threads\
    -phred33\
    -trimlog trim_and_QC_log.txt\
    "$fastqs_path"/"$ID".fq \
    "$fastqs_path"/"$ID"_trimmed.fq.gz \
    ILLUMINACLIP:$adapters:1:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90 >> trim_and_QC_log.txt

    echo "Beginning post-trim QC for "$ID >> trim_and_QC_log.txt
    fastqc -t $threads \
    "$fastqs_path"/"$ID"_trimmed.fq.gz \
    --outdir post_trim_QC_files/
    echo $ID" post-trim QC done" >> trim_and_QC_log.txt

    done<"$filename"
fi
