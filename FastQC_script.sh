#!/bin/bash

INPUT="/path/to/raw_fastq_files"
OUTPUT="/path/to/fastqc_reports/output"

mkdir -p "$OUTPUT"  #Creates the output if the output folder does not exist

# A for loop to identify the fastq files that are supposed to be used for the report generation (individual reports are generated for all paired-end reads)
for file in "$INPUT"/*.fastq.gz;
do
   if [[ "$file" == *"Undetermined"* ]]; then  #An if condition to avoid reading the fastq files (paired-end) for Undetermined fastqs that are generated during bcl to fastq conversion
        continue
   fi

   fastqc -t 10 "$file" -o "$OUTPUT"
done

echo "FastQC reports generated and saved in the "$OUTPUT" directory"
