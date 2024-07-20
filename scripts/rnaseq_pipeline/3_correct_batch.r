#! /usr/bin/env Rscript
# args[1]: count_matrix
# args[2]: batch_list
# args[2]: output_file

args = commandArgs(trailingOnly=TRUE)

# Load library
library("sva")

######################################
# Load data
######################################
# Load uncorrected data
uncorrected_data = read.table(
  args[1], header=TRUE, sep=',', check.names=FALSE
)
# Load batch list
batches <- read.table(
  args[2]
)
batches <- as.vector(batches[['V1']])

######################################
# Run correction
######################################
corrected_data = ComBat_seq(
  counts = as.matrix(uncorrected_data[,2:length(uncorrected_data)]),
  batch = batches
)


######################################
# Save
######################################
# Concatenate
corrected_data = cbind(uncorrected_data[,c(1)], corrected_data)
# Change columns name
first_col_name = names(uncorrected_data)[1]
colnames(corrected_data)[1] <- first_col_name
# Save
write.csv(corrected_data, file = args[3], row.names = FALSE)