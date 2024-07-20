#!/bin/bash

# Concatenating and normalizing VCF files
# Input: list_of_vcfs.txt, reference.fasta
# Output: concatenated.vcf, normalized.vcf

# Define paths
LIST_OF_VCFS="path_to_list_of_vcfs.txt"
REFERENCE="path_to_reference.fasta"
CONCAT_VCF="path_to_concatenated.vcf"
NORM_VCF="path_to_normalized.vcf"

# Concatenate VCF files
bcftools concat -f $LIST_OF_VCFS -o $CONCAT_VCF

# Normalize VCF
bcftools norm -f $REFERENCE -o $NORM_VCF -d all $CONCAT_VCF
