#!/bin/bash

# Marking duplicates with Picard
# Input: sorted.bam
# Output: dupmarked.bam, dupmetrics.txt

# Define paths
SORTED_BAM="path_to_sorted.bam"
DUPMARKED_BAM="path_to_dupmarked.bam"
DUPMETRICS="path_to_dupmetrics.txt"

# Run Picard MarkDuplicates
picard MarkDuplicates I=$SORTED_BAM O=$DUPMARKED_BAM M=$DUPMETRICS
