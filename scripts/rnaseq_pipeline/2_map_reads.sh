#!/bin/bash

# Aligning reads to reference genome using HISAT2 and sorting with samtools
# Input: trimmed_r1.fastq.gz, trimmed_r2.fastq.gz, mtbrf.fa
# Output: sorted_bam.bam

# Define paths
TRIMMED_R1="path_to_trimmed_r1.fastq.gz"
TRIMMED_R2="path_to_trimmed_r2.fastq.gz"
REFERENCE_SEQUENCE="mtbrf.fa"
SORTED_BAM="path_to_sorted.bam"

# Indexing the reference genome with HISAT2 and samtools
hisat2-build $REFERENCE_SEQUENCE $REFERENCE_SEQUENCE
samtools faidx $REFERENCE_SEQUENCE

# Aligning reads and sorting with HISAT2 and samtools
hisat2 -x $REFERENCE_SEQUENCE -1 $TRIMMED_R1 -2 $TRIMMED_R2 | samtools sort -o $SORTED_BAM
