#!/bin/bash

set -euo pipefail  #Set to terminate the script in-case there is any error in the script

INPUT="/path/to/bcl_files/raw_data_folder"
OUTPUT="/path/to/output/bcl2fastq"

mkdir -p "$OUTPUT"  #Creates the output folder if it does not exist

#We use the bcl-convert tool to convert the base calling files into fastq files for usage

bcl-convert \
    --bcl-input-directory "$INPUT" \
    --output-directory "$OUTPUT" \
    --sample-sheet /path/to/sample_sheet_rRNA_bcl2fastq_run.csv \
    --bcl-num-compression-threads 20 \  #threading option for numerous samples
    --force
    
echo "Conversion complete"
