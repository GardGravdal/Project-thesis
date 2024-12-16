library(dplyr)
library(stringr)
library(ggplot2)
library(qqman)
library(readr)
library(data.table)
library(tidyr)
library(smplot2)
library(latex2exp)
setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results")
resultGemma <- read_table("gwas_wing_70k.txt")


# Some reformatting on the stored data:
# ---------------------------------
# Gemma stores lines in GWASresults across two lines... Need to fix
txt <- readLines("gwas_wing_70k.txt")
str(txt)


txt <- txt[-1]

# merge every 2 subsequent lines into one to form a row of final dataframe
idx <- seq(1, length(txt), by=2)
txt[idx] <- paste(txt[idx], txt[idx+1])
txt <- txt[-(idx+1)]

# final data
resultGemma <- read.table(text=txt, col.names=colnames(resultGemma))

# Exclude certain chromosomes in plot
resultGemma <- resultGemma[!(resultGemma$chr %in% c(0, 16,30,32)),]
# -------------------------

#compute the Bonferroni threshold
bonferroni<- -log10(0.05/ nrow(resultGemma))

# png(paste("GWAS_",PHENOTYPE,".png",sep=""))
setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Figures")
pdf(file="WingGemma70k.pdf", width = 7.20, height = 7.0)
manhattan(resultGemma,chr="chr",bp="ps",p="p_lrt",snp="rs", ylim = c(0, 7), col = c("gray10", "gray50"),
          cex.axis = 1.5, cex.lab=2.3, cex.main=2.5, ylab="", genomewideline = FALSE, suggestiveline = FALSE)
abline(h = bonferroni, col = "red", lty = 2)
title(ylab = TeX(r"($-log_{10}(p)$)"),line=2, cex.lab=2.3)
dev.off()


