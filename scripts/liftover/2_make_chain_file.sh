#!/bin/bash

# Generate chain file from BLAT results
# Input: psl file, 2bit files
# Output: final chain file

# Define paths
QRY_FILE="path_to_query.fa" # Complete genome
REF_FILE="path_to_reference.fa" # Backbone genome (H37Rv-1)
OUT_PATH="path_to_output_directory"

REF_2BIT="$OUT_PATH/$(basename $REF_FILE .fa).2bit"
QRY_2BIT="$OUT_PATH/$(basename $QRY_FILE .fa).2bit"
PSL_FILE="${QRY_2BIT%.2bit}.psl"
CHAIN_FILE_1="${QRY_2BIT%.2bit}.chain_1"
NET_FILE="${QRY_2BIT%.2bit}.net"
CHAIN_FILE_2="${QRY_2BIT%.2bit}.chain_2"
FINAL_CHAIN="${QRY_2BIT%.2bit}.chain"

# Create initial chain file
axtChain -linearGap=medium -psl $PSL_FILE $QRY_2BIT $REF_2BIT $CHAIN_FILE_1

# Generate net file
chainNet $CHAIN_FILE_1 $QRY_CHROM_SIZES $REF_CHROM_SIZES $NET_FILE /dev/null

# Subset the net chains
netChainSubset $NET_FILE $CHAIN_FILE_1 $CHAIN_FILE_2

# Swap the chain file orientation
chainSwap $CHAIN_FILE_2 $FINAL_CHAIN
