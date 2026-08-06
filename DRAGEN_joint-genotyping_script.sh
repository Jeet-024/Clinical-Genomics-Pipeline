# This script is written to comprehend the DRAGEN-FPGA based pipeline to generate a jointly-called or jointly-genotyped variant calling file
# The reason we use DRAGEN is that it takes less time to run a larger cohort of samples (as we deal with 100+ samples); we do not need to run alignment and variant calling separately.

#!/bin/bash

HASH_TABLE_DIR="/path/to/reference/human/hg38_hash_table"  #A hash-table is required as a reference, besides the fasta file of the human genome
TARGET_BED="/path/to/target.bed/file" #since we are dealing with Whole Exome sequencing, we need a target-bed file to target alignment against those regions (Twist exome bed can be used, Illumina site)
REF="/path/to/reference/hg38.fa" #FASTA reference file
FASTQ_DIR="/path/to/trimmed/fastq/files"
OUTPUT_DIR="/path/to/Dragen_Output"
TEMP_DIR="/path/to/temp_dir"


mkdir -p "$OUTPUT_DIR"  #Creating the output directory
mkdir -p "$TEMP_DIR"    #Creating the temporary directory (we do this to store the temporary files created during running the DRAGEN pipeline)

GVCF_LIST="$OUTPUT_DIR/gvcf_list.txt"
> "$GVCF_LIST"

echo "Starting single-sample gVCF generation..."

for f1 in "$FASTQ_DIR"/GRIP*_R1_001.fastq.gz;
do
        filename=$(basename "$f1")
        base_name="${filename%_R1_001.fastq.gz}"
        f2="${FASTQ_DIR}/${base_name}_R2_001.fastq.gz"

        echo "Processing sample: $base_name"
        echo "Read 1 : $f1"
        echo "Read 2 : $f2"
        sample_out_dir="$OUTPUT_DIR/$base_name"

        mkdir -p "$sample_out_dir"

        dragen -f \
        -r "$HASH_TABLE_DIR" \
        -1 "$f1" \
        -2 "$f2" \
        --RGID "$base_name" \
        --RGSM "$base_name" \
        --output-directory "$sample_out_dir" \
        --output-file-prefix "$base_name" \
        --intermediate-results-dir "$TEMP_DIR" \
        --enable-map-align-output true \
        --enable-variant-caller true \
        --enable-duplicate-marking true \
        --enable-bam-indexing true \
        --vc-target-bed-padding 100 \
        --vc-emit-ref-confidence gVCF \
        --vc-target-bed "$TARGET_BED"

        echo "${sample_out_dir}/${base_name}.hard-filtered.gvcf.gz" >> "$GVCF_LIST"
done

echo "gVCF generation complete. Starting Joint Genotyping..."

dragen -f \
-r "$HASH_TABLE_DIR" \
--ht-reference "$REF" \
--enable-joint-genotyping true \
--enable-map-align-output true \
--variant-list "$GVCF_LIST" \
--output-directory "$OUTPUT_DIR/cohort_joint_output" \
--output-file-prefix all_samples \
--vc-target-bed-padding 100 \
--vc-target-bed "$TARGET_BED"

echo "Workflow finished successfully!"
