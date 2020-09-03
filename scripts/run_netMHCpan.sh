#!/bin/bash
# Script for running netMHCpan in parallel

# Set parameters
TOT_LEN=2300000*2
SPLITS=$1
L=$((TOT_LEN/SPLITS))
echo "Each file has "$L" lines"

# Clear out temporary directory
SEARCH_DIR='./TEMP'
rm $SEARCH_DIR/*

ALLELE=$2 
ALLELE_=$( echo $ALLELE | tr ':' '_' ) 

ITR=0
# Split the input file
split -l $L  ../tim_data/MOD_frameshiftPeptidesComplete.fasta ./TEMP/

# Run netMHCpan on each split file
for SEQ_FILE in "$SEARCH_DIR"/*
do
  echo "$ITR": "$SEQ_FILE""$ALLELE"
  ../netMHCpan-4.1/netMHCpan $SEQ_FILE -l 8,9,10,11,12 -BA -xls -t 25 -a $ALLELE -xlsfile ./results/$ALLELE_"_"$ITR.xls >./results/$ALLELE_"_"$ITR & 
  ITR=$((ITR + 1))
done
echo "Waiting"
wait
