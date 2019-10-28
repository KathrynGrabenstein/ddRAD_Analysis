#!/bin/bash

#Script for demultiplexing ddRAD data from Illumina using STACKS
#Output files named with sample ID and placed in output directory
#required components for looping through index groups to make beforehand

#list of index primers and samples
#List of illumina adapters for index group

# if variables entered into script call line are less than 1, echo error message
if [ $# -lt 1 ]
  then
    echo "Demultiplexes single-end ddRAD sequence data from Illumina using STACKS.
    [-a] List of index groups with adapter sequences
    [-o] Path to output directory
        OPTIONAL ARGUMENTS
    [-t] Number of threads to use"

  else
    while getopts p:i:o:a:t: option
    do
    case "${option}"
    in
    a) adapters=${OPTARG};;
    o) output_dir=${OPTARG};;
    t) threads=${OPTARG};;
    esac
    done

#default options for
    threads="${threads:-1}"

    # if doesn't exist, make directory for where demultiplexed fastq files go

    #count number of rows in index file to get range for i
    Num_Index_Groups=$(cat "$adapters" | wc -l)


    for i in $(seq "$Num_Index_Groups"); do #number of index groups to loop through
        #locate fastq files in input directory
        Input_file=$"/data0/chickadeeddRAD_workdir/untrimmed_fastq/*KCG_index_"$i"_*.fastq"
        echo "Input file is ""$Input_file" >> Demultiplex_log.txt
        
        #assign barcodes file for proper index group
        Barcodes=$"KCG_Index_"$i".txt"
        #print demulitplexing for index group to log file
        echo "Demultiplexing ""$Barcodes" >> Demultiplex_log.txt

        illumina_adapt_index=$(awk -F "," -v index_group=$i 'NR == index_group {print $2}' "$adapters")
        echo "Illumina adapter is ""$illumina_adapt_index" >> Demultiplex_log.txt
 
        #for each index group, run STACKS
        #path to files is given in flag -p
        #barcodes and sample names given in flag -b
        #output directory for demultiplexed files given in flag -o
        #all other flags here are likely project specific, will need to update

        /home/kagr3234/programs/stacks-2.41/process_radtags -f $Input_file -b "$Barcodes" -o "$output_dir" -e sbfI -c -q -E phred33 --inline_null -i fastq --adapter_1 "$illumina_adapt_index" --adapter_mm 1 --filter_illumina

        echo "$ID"" demultiplexing done" >> Demultiplex_log.txt
        #done < "$Barcodes"
    done
fi


#takes ~ 1 hr to run for 22 index groups
