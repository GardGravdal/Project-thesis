##############
# File adds chromosomes to SHAP/GWAS and creates "correlation" plot for the two models
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

setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline")
# Getting results from SHAP:
shap_vals <- read_feather("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/shap_mass_70k.feather")
resultGemma <- read_table("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/gwas_mass_70k.txt")

# Getting results from GEMMA:
# Some reformatting on the stored data:
# ---------------------------------
# Gemma stores lines in GWASresults across two lines... Need to fix
txt <- readLines("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/gwas_mass_70k.txt")
str(txt)


txt <- txt[-1]

# merge every 2 subsequent lines into one to form a row of final dataframe
idx <- seq(1, length(txt), by=2)
txt[idx] <- paste(txt[idx], txt[idx+1])
txt <- txt[-(idx+1)]

# final data
resultGemma <- read.table(text=txt, col.names=colnames(resultGemma))

# -------------------------------------------

# Some prepping is needed to match the SHAP values with chromosome locations
# NED A MAP FILE, gives where each SNP is loccated
map_path <- "C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Data/Helgeland_01_2018.map"
# map_path <- "data/raw/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
chro <- read.table(map_path, header = F)

names(chro) <- c("chr", "snpID", "value", "bp")
#head(chro)
sort(unique(chro$chr))

snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

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
#head(tmpdf2)

# merge the map-file stuff with the SNPname from the SNP columns in shap
tmpdfchro <- merge(tmpdf2, chro[, c("chr", "snpID", "bp")], by = "snpID")
#head(tmpdfchro)

# Find SNPs in tmpdf$snpID that are not in resultGemma$rs
snps_not_in_result <- setdiff(tmpdf$snpID, resultGemma$rs)

# Make sure that GWAS and SHAP contains the same SNPs for comparison
snp_to_discard <- unique(c(snp_to_discard, snps_not_in_result))

snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]

# we now have a mapping between SHAP SNPs and the chromosome location
sum(chro$snpID %in% tmpdf$snpID)


# remove some SNP on chromosomes not desired for both SHAP and GEMMA
mean_shap <- shap_vals[, c(snp_to_keep)]
resultGemma <- resultGemma[!(resultGemma$chr %in% c(0, 16,30,32)),]

# They are now of same dimension
#dim(mean_shap)
#dim(shap_vals)

rshap <- t(mean_shap)
rshap_2 <- data.frame(origsnp = rownames(rshap), shap = rshap[, 1], row.names = NULL)


sort(unique(tmpdfchro$chr))
# merge shap with location
manhattan_df <- merge(rshap_2, tmpdfchro, by = "origsnp")
# NOW WE CAN PLOT

# Ensure manhattan_df is ordered to match resultGemma$rs
manhattan_df <- manhattan_df[match(resultGemma$rs, manhattan_df$snpID), ]

# Compute -log10(p-values)
resultGemma$log_p <- -log10(resultGemma$p_lrt)

# Fit a linear model
lm_model <- lm(manhattan_df$shap ~ resultGemma$log_p)

# Define the range of x-values for the line (positive x-values only)
x_vals <- seq(0, length(resultGemma$log_p), length.out = length(resultGemma$p_lrt))

# Compute the corresponding y-values using the regression slope
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals

# Make the scatter plot
setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Figures")
pdf(file="MassShapGWAS.pdf", width = 7.20, height = 7.0)
par(mar = c(5, 6, 4, 1) + .1)  
plot(
  resultGemma$log_p, manhattan_df$shap,  
  xlab = expression(-log[10](italic(p))),  # Styled x-axis label
  ylab = "mean |SHAP|",  # Y-axis label
  cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5,  # Adjust text sizes
  pch = 20, col = "black"  # Scatter point style
)
# Add the regression line to the plot
lines(x_vals, y_vals, col = "red", lwd = 2)
dev.off()
# -------------------------------------------
setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline")



# Results for wing trait

# Getting results from SHAP:
shap_vals <- read_feather("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/shap_wing_70k.feather")
resultGemma <- read_table("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/gwas_wing_70k.txt")


# Getting results from GEMMA:
# Some reformatting on the stored data:
# ---------------------------------
# Gemma stores lines in GWASresults across two lines... Need to fix
txt <- readLines("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/gwas_wing_70k.txt")
str(txt)


txt <- txt[-1]

# merge every 2 subsequent lines into one to form a row of final dataframe
idx <- seq(1, length(txt), by=2)
txt[idx] <- paste(txt[idx], txt[idx+1])
txt <- txt[-(idx+1)]

# final data
resultGemma <- read.table(text=txt, col.names=colnames(resultGemma))

# -------------------------------------------

# Some prepping is needed to match the SHAP values with chromosome locations
# NED A MAP FILE, gives where each SNP is loccated
map_path <- "C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Data/Helgeland_01_2018.map"
# map_path <- "data/raw/combined_200k_70k_helgeland_south_corrected_snpfiltered_2024-02-05.map"
chro <- read.table(map_path, header = F)

names(chro) <- c("chr", "snpID", "value", "bp")
#head(chro)
sort(unique(chro$chr))

snp_to_discard <- chro[chro$chr %in% c(0, 16, 30, 32), "snpID"]

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
#head(tmpdf2)

# merge the map-file stuff with the SNPname from the SNP columns in shap
tmpdfchro <- merge(tmpdf2, chro[, c("chr", "snpID", "bp")], by = "snpID")
#head(tmpdfchro)

# Find SNPs in tmpdf$snpID that are not in resultGemma$rs
snps_not_in_result <- setdiff(tmpdf$snpID, resultGemma$rs)

# Make sure that GWAS and SHAP contains the same SNPs for comparison
snp_to_discard <- unique(c(snp_to_discard, snps_not_in_result))

snp_to_keep <- tmpdf[!(tmpdf$snpID %in% snp_to_discard), "origsnp"]

# we now have a mapping between SHAP SNPs and the chromosome location
sum(chro$snpID %in% tmpdf$snpID)


# remove some SNP on chromosomes not desired for both SHAP and GEMMA
mean_shap <- shap_vals[, c(snp_to_keep)]
resultGemma <- resultGemma[!(resultGemma$chr %in% c(0, 16,30,32)),]

# They are now of same dimension
#dim(mean_shap)
#dim(shap_vals)

rshap <- t(mean_shap)
rshap_2 <- data.frame(origsnp = rownames(rshap), shap = rshap[, 1], row.names = NULL)


sort(unique(tmpdfchro$chr))
# merge shap with location
manhattan_df <- merge(rshap_2, tmpdfchro, by = "origsnp")
# NOW WE CAN PLOT

# Ensure manhattan_df is ordered to match resultGemma$rs
manhattan_df <- manhattan_df[match(resultGemma$rs, manhattan_df$snpID), ]

# Store neg log_p values
resultGemma$log_p <- -log10(resultGemma$p_lrt)

# Fit a linear model
lm_model <- lm(manhattan_df$shap ~ resultGemma$log_p)

# Define the range of x-values for the line (positive x-values only)
x_vals <- seq(0, length(resultGemma$log_p), length.out = length(resultGemma$p_lrt))

# Compute the corresponding y-values using the regression slope
y_vals <- coef(lm_model)[1] + coef(lm_model)[2] * x_vals


# Make the scatter plot
setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Figures")

pdf(file="WingShapGWAS.pdf", width = 7.20, height = 7.0)
par(mar = c(5, 6, 4, 1) + .1)  
plot(
  resultGemma$log_p, manhattan_df$shap,  
  xlab = expression(-log[10](italic(p))),  
  ylab = "mean |SHAP|",  
  cex.axis = 1.5, cex.lab = 2.3, cex.main = 2.5,  # Adjust text sizes
  pch = 20, col = "black"  # Scatter point style
)
# Add the regression line to the plot
lines(x_vals, y_vals, col = "red", lwd = 2)
dev.off()

