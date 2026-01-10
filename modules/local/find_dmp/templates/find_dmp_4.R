#!/usr/bin/env Rscript
##This script reads Beta values of methylation arrays and computes DMPs, Differentially methylated positions
##Input = beta values of pre-processed methylation data
##Input = metadata information with sample ids and group information (categories: case, control or different groups)
##output = csv files with p-values for the DMPs in each comparison of categories

##Computations of DMPs using ChAMP and minfi

library(dplyr)
library(readr)
library(minfi)
library(tibble)
library(ChAMP)
library(rio)


##Load previously saved data (RData objects, for more details, please look at pre-processing.Rmd, cell_composition_correction.R and rem_conf_probes_adj_age.R
bVals <- read_csv("$bVALS_SNPPROBES")

if ('probe' %in% colnames(bVals)) {
    bVals <- bVals %>% tibble::column_to_rownames(var = 'probe')
}

metadata <- read.csv("$extensive_metadata", header = FALSE) %>%
            filter(V1 %in% c(colnames(bVals)))
####Match metadata order to columns in bVals
metadata <- metadata[match(c(colnames(bVals)), metadata\$V1),]

##Choose the samples field
##Choosing only the field with information of the category of each samples (can be multiple)
Samples <- metadata\$V1
Class <- metadata\$V2

##Choose the adjustment method: can be "BH",
Method = "BH"

##Choose array type: ("EPIC" or "450K")
ARRAY = "EPIC"

###Choose adjusted P value
P = 0.05
P = 1 # NOTE: for development

###Computing DMPs
dmp_data <- champ.DMP(
    beta = bVals[,Samples],
    pheno = Class,
    adjPVal = P,
    adjust.method = Method,
    arraytype = ARRAY
)

export_list(dmp_data, file = "dmp_champ.%s.csv")

###Finding DMPs with another method (required binary categories of "Class", it will work with multiple categories but the do not specify which comparison it is)
print(bVals[,Samples])
dmp_minfi <- dmpFinder(
        as.matrix(bVals[,Samples]),
        Class,
        type = "categorical") %>%
    filter(pval < P) %>%
    arrange(pval)

write_csv(dmp_minfi, file = "dmp_minfi.csv")

# Dump versions
pkgs <- c("dplyr","readr","minfi","tibble","ChAMP","rio")
pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NA")
rver <- paste(R.version\$major, R.version\$minor, sep=".")
lines <- c(
  sprintf('"%s":', "${task.process}"),
  sprintf('  R: "%s"', rver),
  sprintf('  %s: "%s"', pkgs, vapply(pkgs, pkg_ver, character(1)))
)
writeLines(lines, "versions.yml")
