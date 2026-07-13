#!/usr/bin/env Rscript
##This script reads Beta values of methylation arrays and computes blocks of dna regions with many DMPs,
##Input = beta values of pre-processed methylation data
##Input = metadata information with sample ids and group information (categories: case, control or different groups)
##output = csv files with p values of blocks

library(dplyr)
library(readr)
library(ChAMP)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd, cell_composition_correction.R and rem_conf_probes_adj_age.R
###Load necessary data
bVals <- read_csv("$bVALS_SNPPROBES")
metadata <- read.csv("$extensive_metadata", header = FALSE) %>%
    filter(V1 %in% c(colnames(bVals)))
####Match metadata order to columns in bVals
metadata <- metadata[match(c(colnames(bVals)), metadata\$V1),]

##Choose the samples field
##Choosing only the field with information of the category of each samples (can be multiple)
Samples <- metadata\$V1
Class <- metadata\$V2

##Choose array type: ("EPIC" or "450K")
ARRAY = "EPICv2"

###Choose number of minimum number of regions per block
Min = 5

###Choose maximum gap between 2 clusters in a block
Max = 250000

###Computing blocks
blocks <- champ.Block(
    beta = bVals[,Samples],
    pheno = Class,
    arraytype = ARRAY,
    maxClusterGap = Max,
    minNum = Min
)

write_csv(blocks, "blocks_champ.csv")

# Dump versions
pkgs <- c("dplyr","readr","ChAMP")
pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NA")
rver <- paste(R.version\$major, R.version\$minor, sep=".")
lines <- c(
  sprintf('"%s":', "${task.process}"),
  sprintf('  R: "%s"', rver),
  sprintf('  %s: "%s"', pkgs, vapply(pkgs, pkg_ver, character(1)))
)
writeLines(lines, "versions.yml")
