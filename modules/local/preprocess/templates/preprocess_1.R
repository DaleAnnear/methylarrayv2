#!/usr/bin/env Rscript

# Load necessary libraries
library(minfi)
library(IlluminaHumanMethylationEPICv2manifest)

# Unzip the idat files if they are gzipped
gz_files <- list.files(getwd(), pattern="gz\$", recursive=TRUE, full.names=TRUE)
for (file in gz_files) {
    out <- sub(".gz\$", "", file)
    system2("gzip", args = c("-dc", shQuote(file)), stdout = out)
}

# Mapping table
mapping <- read.csv("${sample_mapping}", header = FALSE)

# Set constants
P = 0.1
Norm_method = "preprocessQuantile"  # Choose normalization method: "preprocessFunnorm" or "preprocessQuantile"
Sample_Name <- "Sample_Name"        # Column name that encodes the sample names

# Create rgSet
rgSet <- read.metharray.exp(
    targets = NULL,
    base = getwd(),
    verbose = TRUE
)

# Assign significant names to the columns (samples)
sampleNames(rgSet) <- mapping\$V2[match(sampleNames(rgSet), mapping\$V1)]

# QC: Detection p-values of the signal quality
detP <- detectionP(rgSet)
keep <- colMeans(detP) < P
rgSet <- rgSet[, keep]

# Remove poor quality samples from targets and detection p-value table
detP <- detP[, keep]

# Normalize the data
if (Norm_method == "preprocessQuantile") {
    mSetSq <- preprocessQuantile(rgSet)
} else if (Norm_method == "preprocessFunnorm") {
    mSetSq <- preprocessFunnorm(rgSet)
}

# Filtering: Put probes in the same order in mSetSq and detP
detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]

# Filter probes failing in 1 or more samples (change threshold if needed)
keep <- rowSums(detP < P) == ncol(mSetSq)
mSetSqFlt <- mSetSq[keep, ]

# Calculate M and Beta values
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

# Save results
print(paste0('No of probes: ', nrow(bVals)))
write.csv(as.data.frame(mVals), "mVals.csv")
write.csv(as.data.frame(bVals), "bVals.csv")
save(mSetSqFlt, file = "mSetSqFlt.RData")
save(rgSet, file = "rgSet.RData")

# Dump versions
pkgs <- c("minfi","IlluminaHumanMethylationEPICv2manifest")
pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error=function(e) "NA")
rver <- paste(R.version\$major, R.version\$minor, sep=".")
lines <- c(
  sprintf('"%s":', "${task.process}"),
  sprintf('  R: "%s"', rver),
  sprintf('  %s: "%s"', pkgs, vapply(pkgs, pkg_ver, character(1)))
)
writeLines(lines, "versions.yml")
