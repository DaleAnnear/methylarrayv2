#!/usr/bin/env Rscript
##This script reads M-values and Beta values of methylation arrays and corrects for cell composition (samples from whole blood)
##Use after visualizing cell composition in the different samples and deciding if you want to correct for that effect
##Input = Beta values
##Output = M-values and beta values corrected for cell commpositions

# Debug intersection of probes
suppressPackageStartupMessages({library(readr); library(tibble); library(ChAMP); library(ChAMPdata)})
df <- read_csv("$bVALS_SNPPROBES", show_col_types=FALSE)
b  <- as.matrix(column_to_rownames(df, "probe"))
storage.mode(b) <- "double"
rownames(b) <- trimws(rownames(b))
cat("b dim:", dim(b)[1], "x", dim(b)[2], "\n")
cat("b range:", paste(range(b, na.rm=TRUE), collapse=" .. "), "\n")
data("CellTypeMeans450K", package="ChAMPdata")
cat("ref dim:", dim(CellTypeMeans450K)[1], "x", dim(CellTypeMeans450K)[2], "\n")
dmrs <- intersect(rownames(CellTypeMeans450K), rownames(b))
cat("overlap:", length(dmrs), "\n")

library(readr)
library(ChAMP)
library(lumi)

##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd)
bVals <- as.matrix(read_csv("$bVALS_SNPPROBES") %>% tibble::column_to_rownames("probe"))

##Correct for cell composition using ChAMP
bVals_corrected <- champ.refbase(beta = bVals, arraytype = "450K")
bVals <- bVals_corrected\$CorrectedBeta
mVals <- beta2m(bVals)

##Save all necessary R objects for later use
write_csv(as.data.frame(mVals), "cmVals.cell_comp.csv")
write_csv(as.data.frame(bVals), "cbVals.cell_comp.csv")

# Dump versions
pkgs <- c("readr","ChAMP","lumi")
pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NA")
rver <- paste(R.version\$major, R.version\$minor, sep=".")
lines <- c(
  sprintf('"%s":', "${task.process}"),
  sprintf('  R: "%s"', rver),
  sprintf('  %s: "%s"', pkgs, vapply(pkgs, pkg_ver, character(1)))
)
writeLines(lines, "versions.yml")
