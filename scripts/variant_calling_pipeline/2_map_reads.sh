#!/bin/bash

# Mapping with BWA-MEM
# Input: trimmed_r1.fastq.gz, trimmed_r2.fastq.gz, reference.fasta
# Output: sorted.bam

# Define paths
REFERENCE="path_to_reference.fasta"
TRIMMED_R1="path_to_trimmed_r1.fastq.gz"
TRIMMED_R2="path_to_trimmed_r2.fastq.gz"
SORTED_BAM="path_to_sorted.bam"

# Index the reference genome
bwa index $REFERENCE

# Run BWA-MEM and sort with Samtools
bwa mem -aM $REFERENCE $TRIMMED_R1 $TRIMMED_R2 | samtools sort -o $SORTED_BAM -
