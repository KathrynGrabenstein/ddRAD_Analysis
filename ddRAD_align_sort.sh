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
    if [ $bamoutdir == bam_files/ ]
      then
        mkdir bam_files
    fi

    echo "Beginning alignment for "$seqs >> bwa_alignment_log.txt
    echo "Indexing reference" >> bwa_alignment_log.txt

    #index reference BCCH genome
    bwa index $ref

    while read -r ID; do
    echo "Aligning " $ID >> bwa_alignment_log.txt

    #Align reads to reference genome
    #pipe sam file to samtools (since don't need to save sam file)
    bwa mem -t $threads $ref "$fastqs_path"/"$ID"_trimmed.fq.gz | \
    #convert sam file to bam file
    samtools view -b -o "$bamoutdir"/"$ID".bam -S

    echo $ID" SAM file piped into SAMtools view to convert to .bam" >> bwa_alignment_log.txt

    done < "$seqs"

    if [ $sortedbamoutdir == sorted_bam_files/ ]
      then
        mkdir sorted_bam_files
    fi

    #create sorted bam files
    while read -r SAMPLE; do

    # -@ is number of threads, here 6
    # -T temp, write temporary files with TEMP/xxxxx
    samtools sort -T temp -@ 6 \
    #output file location and type
    -o "$sortedbamoutdir"/"$SAMPLE"_sorted.bam \
    #input file
    "$bamoutdir"/"$SAMPLE".bam

    #
    picard-tools AddOrReplaceReadGroups \
    #input BAM file
    I="$sortedbamoutdir"/"$SAMPLE"_sorted.bam \
    #output BAM file
    O="$sortedbamoutdir"/"$SAMPLE"_sorted_RGadded.bam \
    RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE

    picard-tools MarkDuplicates \
    I="$sortedbamoutdir"/"$SAMPLE"_sorted_RGadded.bam \
    O="$sortedbamoutdir"/"$SAMPLE"_sorted_RGadded_dupmarked.bam \
    M="$sortedbamoutdir"/"$SAMPLE"_dupmarked_metrics.txt

    #create samfile index
    samtools index \
    "$sortedbamoutdir"/"$SAMPLE"_sorted_RGadded_dupmarked.bam

    #remove intermediate files
    #keep only sorted, reading group added, duplicates marked BAM file
    rm "$sortedbamoutdir"/"$SAMPLE"_sorted.bam
    rm "$sortedbamoutdir"/"$SAMPLE"_sorted_RGadded.bam

    done<"$seqs"
fi
