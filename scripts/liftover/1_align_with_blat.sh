#!/bin/bash

# Align query to reference using BLAT
# Input: query.fa, reference.fa
# Output: psl file

# Define paths
QRY_FILE="path_to_query.fa" # Complete genome
REF_FILE="path_to_reference.fa" # Backbone genome (H37Rv-1)
OUT_PATH="path_to_output_directory"

REF_2BIT="$OUT_PATH/$(basename $REF_FILE .fa).2bit"
QRY_2BIT="$OUT_PATH/$(basename $QRY_FILE .fa).2bit"
REF_CHROM_SIZES="$OUT_PATH/$(basename $REF_FILE .fa).chrom_sizes"
QRY_CHROM_SIZES="$OUT_PATH/$(basename $QRY_FILE .fa).chrom_sizes"
PSL_FILE="${QRY_2BIT%.2bit}.psl"

# Convert FASTA to 2bit and generate chrom sizes
faToTwoBit $REF_FILE $REF_2BIT
faToTwoBit $QRY_FILE $QRY_2BIT
twoBitInfo $REF_2BIT $REF_CHROM_SIZES
twoBitInfo $QRY_2BIT $QRY_CHROM_SIZES

# Align query to reference using BLAT
blat $QRY_2BIT $REF_FILE $PSL_FILE -tileSize=11 -minScore=30 -minIdentity=90
