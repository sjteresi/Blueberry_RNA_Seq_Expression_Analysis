library(WGCNA)
library(tidyverse)
setwd("/home/scott/Documents/Uni/Research/Projects/Blueberry_RNA_Seq_Expression_Analysis/WGCNA");

# Load the WGCNA package

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Read in the blueberry data set
# Each row is a gene and each column is a library.
BlueberryData = read.csv("Blueberry_TPM.csv", header=TRUE, row.names='Gene_Name')

# Subset data based on whether a single entry is greater than 10, take that row
# 10 for testing purposes, 0.1 for general guidelines
small_BlueberryData = filter_all(as_tibble(BlueberryData), any_vars(. > 0.1))
small_BlueberryData = filter_all(as_tibble(BlueberryData), any_vars(. > 10.0))


# Tutorial 2
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(small_BlueberryData, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")