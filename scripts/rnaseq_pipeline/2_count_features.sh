#!/bin/bash

# Merge BAM files and calculate feature counts
# Input: sorted_bam1.bam, sorted_bam2.bam, sorted_bam3.bam, mtbrf.gff
# Output: feature_counts.txt

# Define paths for merging BAM files
SORTED_BAM1="path_to_sorted_bam1.bam"
SORTED_BAM2="path_to_sorted_bam2.bam"
SORTED_BAM3="path_to_sorted_bam3.bam"
MERGED_BAM="path_to_merged_bam.bam"

# Run samtools merge
samtools merge $MERGED_BAM $SORTED_BAM1 $SORTED_BAM2 $SORTED_BAM3

# Define paths for feature counts
REFERENCE_ANNOTATION="path_to_mtbrf.gff"
FEATURE_TYPE="exon"
ATTRIBUTE_TYPE="gene_id"
COUNTS_FILE="path_to_feature_counts.txt"

# Run featureCounts
featureCounts -p -t $FEATURE_TYPE -g $ATTRIBUTE_TYPE -a $REFERENCE_ANNOTATION -o $COUNTS_FILE $MERGED_BAM
