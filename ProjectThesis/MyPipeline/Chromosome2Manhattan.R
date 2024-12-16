##############
# File adds chromosomes to SHAP/GWAS and creates manhattan plots for the thesis
##############

library(qqman)
library(tidyr)
library(dplyr)
library(stringr)
library(grid)
library(gridGraphics)
library(ggplot2)
library(dplyr)
library(feather)
library(arrow)

setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Figures")
shap_vals <- read_feather("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/shap_mass_70k.feather")

# Some prepping is needed to match the SHAP values with chromosome locations
# NED A MAP FILE, gives where each SNP is loccated
map_path <- "C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Data/Helgeland_01_2018.map"
# map_path <- "data/raw/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
chro <- read.table(map_path, header = F)

names(chro) <- c("chr", "snpID", "value", "bp")
#head(chro)
sort(unique(chro$chr))

snp_to_discard <- chro[chro$chr %in% c(0, 30, 32), "snpID"]

SNP_cols <- names(shap_vals)
tmpdf <- data.frame(origsnp = SNP_cols, snpID = NA)


# Need to split the SNP name (in the shap df) to get the basepair
get_bp <- function(snp) {
  snp <- str_split(snp, "[_]")[[1]][1]
  return(snp[[1]])
}

tmpdf$snpID <- lapply(tmpdf$origsnp, get_bp)
#head(tmpdf)

# make df instead of list
tmpdf2 <- as.data.frame(lapply(tmpdf, unlist))

# Not sure what next line does, tmpdfchro not used?
# merge the map-file stuff with the SNPname from the SNP columns in shap
tmpdfchro <- merge(tmpdf2, chro[, c("chr", "snpID", "bp")], by = "snpID")
#head(tmpdfchro)

snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]
#head(chro)
#head(tmpdf)
# we now have a mapping between SHAP SNPs and the chromosome location
sum(chro$snpID %in% tmpdf$snpID)


# remove some SNP on chromosomes not desired
mean_shap <- shap_vals[, c(snp_to_keep)]

#dim(mean_shap)
#dim(shap_vals)

rshap <- t(mean_shap)
rshap_2 <- data.frame(origsnp = rownames(rshap), shap = rshap[, 1], row.names = NULL)

sort(unique(tmpdfchro$chr))
length(tmpdfchro$chr)
# merge shap with location
manhattan_df <- merge(rshap_2, tmpdfchro, by = "origsnp")
# NOW WE CAN PLOT
#pdf(paste("reports/figures/SHAP", save_name, ".pdf", sep = ""))
#pdf(paste("Figures/SHAP", "_manhattan", ".pdf", sep = ""))
pdf(file="MassShap70k.pdf", width = 7.20, height = 7.0)
# 5.81,7.20

par(mar = c(5, 6, 4, 1) + .1)
manhattan(manhattan_df, chr = "chr", bp = "bp", p = "shap", snp = "snpID", logp = FALSE, ylim = c(0, 0.013), cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5, ylab = "")
title(ylab = "mean |SHAP|", line = 4, cex.lab = 2.3)
dev.off()
