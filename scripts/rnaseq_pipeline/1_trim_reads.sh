#!/bin/bash

# Trimming low quality reads with Trimmomatic
# Input: raw_r1.fastq.gz, raw_r2.fastq.gz
# Output: trimmed_r1.fastq.gz, trimmed_r2.fastq.gz, unpaired_r1.fastq.gz, unpaired_r2.fastq.gz, trimmomatic.log

# Define paths
RAW_R1="path_to_raw_r1.fastq.gz"
RAW_R2="path_to_raw_r2.fastq.gz"
TRIMMED_R1="path_to_trimmed_r1.fastq.gz"
UNPAIRED_R1="path_to_unpaired_r1.fastq.gz"
TRIMMED_R2="path_to_trimmed_r2.fastq.gz"
UNPAIRED_R2="path_to_unpaired_r2.fastq.gz"
TRIM_LOG="path_to_trimmomatic.log"

# Run Trimmomatic
trimmomatic PE $RAW_R1 $RAW_R2 \
    $TRIMMED_R1 $UNPAIRED_R1 \
    $TRIMMED_R2 $UNPAIRED_R2 \
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25 -phred33 2> $TRIM_LOG
