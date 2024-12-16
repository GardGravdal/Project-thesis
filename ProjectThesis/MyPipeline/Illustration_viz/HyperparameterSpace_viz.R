# Load necessary library
library(plotly)

# Define the hyperparameter space
param1 <- seq(0, 1, length.out = 50) # Hyperparameter 1
param2 <- seq(0, 1, length.out = 50) # Hyperparameter 2

# Create a grid of hyperparameter combinations
grid <- expand.grid(param1 = param1, param2 = param2)

# Define a simple objective function
# Objective = -(param1^2 + param2^2) + 100 (a peak at (0,0))
grid$objective <- (grid$param1 - 0.7)^2 + (grid$param2 - 0.3)^2# Convert to matrix form for visualization
objective_matrix <- matrix(grid$objective, nrow = length(param1), ncol = length(param2))

# Create a 3D plot
plot_ly(
  x = ~param1, y = ~param2, z = ~objective_matrix,
  type = "surface",
  colorscale = "Viridis"
) %>%
  layout(
    title = "Hyperparameter Space",
    scene = list(
      xaxis = list(title = "Hyperparameter 1"),
      yaxis = list(title = "Hyperparameter 2"),
      zaxis = list(title = "Objective Value")
    )
  )

