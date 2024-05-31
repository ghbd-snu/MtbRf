#!/bin/bash

# Calling variants with bcftools
# Input: dupmarked.bam, reference.fasta
# Output: variants.vcf.gz

# Define paths
REFERENCE="path_to_reference.fasta"
DUPMARKED_BAM="path_to_dupmarked.bam"
VCF_GZ="path_to_variants.vcf.gz"

# Run bcftools mpileup and call
bcftools mpileup -f $REFERENCE $DUPMARKED_BAM | bcftools call --ploidy 1 -mv -Oz -o $VCF_GZ
