#!/usr/bin/env Rscript
##This script reads M-values and Beta values of methylation arrays and corrects for different natural confounders (batch effects can be included as well)
##Use after visualizing different confounders and deciding about the one that need to be corrected for
##input = beta values, metadata with samples column under the name "sample_id", and names of the columns from metadata that you want to adjust the data for.
##output = beta values and M values adjusted for the confounders

library(ChAMP)
library(lumi)
library(readr)
library(dplyr)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd)
bVals <- read_csv("$bVALS_SNPPROBES")
metadata <- read_csv("$extensive_metadata")

####Match metadata order to columns in bVals
metadata <- as.data.frame(metadata[match(setdiff(colnames(bVals), "probe"), metadata\$sample_id), ])


##Correct for bmi and age using ChAMP
bVals <- bVals %>% tibble::column_to_rownames("probe")
bVals <- champ.runCombat(
    beta = bVals,
    pd = metadata,
    variablename = "pheno_sex",
    batchname = c("pheno_dummyColumn"), # NOTE: Best if defined here as otherwise will take all if variable is passed for development
    logitTrans = TRUE  # Change to FALSE if you are using M-values
)

###Log transformation to obtain M-values
mVals <- beta2m(bVals)

bVals <- tibble::rownames_to_column(data.frame(bVals), "probe")
mVals <- tibble::rownames_to_column(data.frame(mVals), "probe")

##Save all necessary R objects for later use
write.csv(mVals, "mVals.cell_comp.cor_bmi_age.csv", row.names=FALSE)
write.csv(bVals, "bVals.cell_comp.cor_bmi_age.csv", row.names=FALSE)

# Dump versions
pkgs <- c("ChAMP","lumi","readr","dplyr")
pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NA")
rver <- paste(R.version\$major, R.version\$minor, sep=".")
lines <- c(
  sprintf('"%s":', "${task.process}"),
  sprintf('  R: "%s"', rver),
  sprintf('  %s: "%s"', pkgs, vapply(pkgs, pkg_ver, character(1)))
)
writeLines(lines, "versions.yml")
