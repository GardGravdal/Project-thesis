# Load required libraries
library(ggplot2)

# Set seed for reproducibility
set.seed(12)

# Case 1: Large effect size with high variance
n_1 <- 24  # sample size
beta_large <- 1.5  # large effect size
x_large <- rnorm(n_1)
epsilon_large <- rnorm(n_1, sd = 3)  # high variance in error term
y_large <- beta_large * x_large + epsilon_large  # linear model

# Case 2: Small effect size with low variance
n_2 <- 130
beta_small <- 0.2  # small effect size
x_small <- rnorm(n_2)
epsilon_small <- rnorm(n_2, sd = 1)  # low variance in error term
y_small <- beta_small * x_small + epsilon_small  # linear model

# Fit linear models
model_large <- lm(y_large ~ x_large)
model_small <- lm(y_small ~ x_small)

# Extract coefficients and p-values
beta_large_hat <- coef(model_large)[2]
ci_large <- confint(model_large)[2, ]
pval_large <- summary(model_large)$coefficients[2, 4]

beta_small_hat <- coef(model_small)[2]
ci_small <- confint(model_small)[2, ]
pval_small <- summary(model_small)$coefficients[2, 4]

# Print the p-values
cat("P-value for large effect, high variance:", pval_large, "\n")
cat("P-value for small effect, low variance:", pval_small, "\n")

# Prepare data for plotting
plot_data <- data.frame(
  Effect = c("Large Effect, High Variance", "Small Effect, Low Variance"),
  Beta = c(beta_large_hat, beta_small_hat),
  CI_Lower = c(ci_large[1], ci_small[1]),
  CI_Upper = c(ci_large[2], ci_small[2]),
  P_Value = c(pval_large, pval_small)
)

# Plot the results
ggplot(plot_data, aes(x = Effect, y = Beta, color = Effect)) +
  geom_point(size = 4) +  # Point for estimated beta
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, size = 1.5) +  # Confidence interval
  labs(
    y = "Estimated effect size and CI",
    x = NULL
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 19, hjust = 0.5),   # Increase title size
        axis.title.y = element_text(size = 16),             # Increase y-axis label size
        axis.text.x = element_text(size = 15),              # Increase x-axis tick label size
        axis.text.y = element_text(size = 15)
        ) +
  annotate("text", x = 1.4, y = max(plot_data$Beta) + 0.5, 
           label = expression(paste("Two cases of p" %~~% 0.01)), 
           size = 7, hjust = 0.02, vjust = 0.02, color = "black")  # Add text in top-right corner


