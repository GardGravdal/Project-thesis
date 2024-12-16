# Script for generating nice multiple box-plots
# Load ggplot2 library
library(ggplot2)
setwd("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Figures")
massBV <- read.csv("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/corr_BV_mass.csv")
wingBV <- read.csv("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/corr_BV_wing.csv")
xgb_mass <- read.csv("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/corr_xgb_mass.csv")
xgb_wing <- read.csv("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Results/corr_xgb_wing.csv")

data <- data.frame(
  Category = factor(rep(c("GEBV mass", "XGB mass", "GEBV wing", "XGB wing"), each = 10),
                    levels = c("GEBV mass", "XGB mass", "GEBV wing", "XGB wing")),
  Values = c(massBV$corr,
             xgb_mass$corr,
             wingBV$corr,
             xgb_wing$corr)
)


# Create the box plot
#pdf(file="Predictions.pdf", width = 7.20, height = 7.0)
ggplot(data, aes(x = Category, y = Values, fill = Category)) +
  geom_boxplot(outlier.color = "red", outlier.shape = 19, outlier.size = 2) +
  scale_fill_manual(values = c("GEBV mass" = "aquamarine", "XGB mass" = "darkcyan", "GEBV wing" = "aquamarine", "XGB wing" = "darkcyan")) +
  labs(
    title = "Prediction of genetic contribution for different models",
    y = expression(Corr(bar(y), hat(y)^"*")),
    x = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 21, hjust = 0.5),   # Increase title size
    axis.title.x = element_text(size = 21),             # Increase x-axis label size
    axis.title.y = element_text(size = 21),             # Increase y-axis label size
    axis.text.x = element_text(size = 21),              # Increase x-axis tick label size
    axis.text.y = element_text(size = 21)               # Increase y-axis tick label size
  )
#dev.off()

