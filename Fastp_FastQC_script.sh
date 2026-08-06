#!/bin/bash

INPUT="/path/to/raw_fastqs"
OUTPUT="/path/to/trimmed_fastqs/output"
OUT_QC_TRIM="/path/to/trimmed_fastqs/report/output"


#Creating the output directories if these do not exist
mkdir -p "$OUTPUT"
mkdir -p "$OUT_QC_TRIM"

#A for loop to read the paired-end fastq files
for f1 in "$INPUT"/*_R1_*.fastq.gz;
do
   [[ "$f1" == *Undetermined* ]] && continue  #avoiding the Undetermined_fastqs

       f2="${f1/_R1_/_R2_}"

       if [[ ! -f "$f2" ]]; then
           echo "Missing pair for $f1"  #this is to flag if any read is missing its pair file
           continue
       fi

       sample=$(basename "${f1%%_L001_R*}")  #this extracts only the base name from the entire file-path

    echo "Processing $sample"

#Using fastp tool to perform the universal adapter trimming
   fastp \
    -i "$f1" \
    -I "$f2" \
    -o "$OUTPUT/${sample}_trimmed_L001_R1.fastq.gz" \
    -O "$OUTPUT/${sample}_trimmed_L001_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --thread 20 \
    -h "$OUTPUT/${sample}_report.html" \
    -j "$OUTPUT/${sample}_report.json"

   echo "Running Trimmed QC..."

#Running fastqc on the trimmed files that were generated (paired-end)
   fastqc \
    -t 10 \
    -o "$OUT_QC_TRIM" \
    "$OUTPUT/${sample}_trimmed_L001_R1.fastq.gz" \
    "$OUTPUT/${sample}_trimmed_L001_R2.fastq.gz"

done

echo "Workflow complete!"
