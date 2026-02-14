#!/bin/bash

################################################################################
# GBStrim Wrapper Script
# Processes a single FASTQ file
# Used by parallel to process multiple files
################################################################################

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <input_fastq.gz> <output_directory> <R1/R2>"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_DIR=$2
READ=$3

BASENAME=$(basename "${INPUT_FILE}" .fastq.gz)
OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_trimmed.fastq"

echo "Processing: ${BASENAME}"

/dados01/jorge/TATIANA/bin/gbstrim/gbstrim.pl \
    --enzyme2 psti \
    --enzyme1 mspi \
    --read $READ \
    --croplength 70 \
    --removecutsite \
    --fastqfile "${INPUT_FILE}" \
    --outputfile "${OUTPUT_FILE}" \
    2>&1 > ${OUTPUT_DIR}/${BASENAME}.gbstrim.log

echo "Completed: ${BASENAME}"
