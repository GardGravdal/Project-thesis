library(rpart)
library(rpart.plot)
library(ggplot2)

# Simulate some data
set.seed(1)
x <- runif(100, 0, 10) #sim u(0,10), 100 samples
y <- 2 * x + rnorm(100, sd = 1)  # y \sim 2x + epsilon, epsilon \sim z
data <- data.frame(x = x, y = y)

# Fit a regression tree
tree_model <- rpart(y ~ x, data = data)

# Plot the decision tree
# type: How much information is displayed (see doc)
# digits: how many digits of x/y displayed
rpart.plot(tree_model, type = 3, digits = 2, fallen.leaves = TRUE)


# Extract the split points
split_points <- unique(tree_model$splits[, "index"])

# Predict values using the tree
data$predicted <- predict(tree_model)
# Plot the data and the regression tree predictions
ggplot(data, aes(x = x, y = y)) +
  # geom_point creates scatterplots. alpha -> transparecy degree 
  geom_point(alpha = 0.6) +
  geom_step(aes(y = predicted), color = "blue", size = 1.2) +
  labs(
    x = "x",
    y = "y"
  ) +
  theme_minimal() + theme(
    legend.position = "none",
    plot.title = element_text(size = 19, hjust = 0.5),   # Increase title size
    axis.title.x = element_text(size = 16),             # Increase x-axis label size
    axis.title.y = element_text(size = 16),             # Increase y-axis label size
    axis.text.x = element_text(size = 15),              # Increase x-axis tick label size
    axis.text.y = element_text(size = 15)               # Increase y-axis tick label size
  )

