#!/bin/sh

#  PCA_plot.sh
#  
#
#  Created by Kathryn Grabenstein on 11/7/19.
#

#shell script to generate a PCA plot from a VCF file
#requires file automatedPCA.R to be located in home directory

if [ $# -lt 1 ]
  then
    echo "Generates two PCA plots, one defining points by population and another defining points by individual.
    [-v] Path to vcf file
    [-o] Output directory for files
    [-p] Population file
    [-s] Sample file
    [-n] Project name, prefix for output files"

  else
    while getopts v:o:p:n:s: option
    do
    case "${option}"
    in
    v) dataset=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pops=${OPTARG};;
    n) name=${OPTARG};;
    s) sample=${OPTARG};;
    esac
    done
  Rscript /data0/chickadeeddRAD_workdir/PCA/automatedPCA.R $dataset $outDir $pops $name

fi
