# Base code to create estimates of above and below ground respiration for ENF and BDF
#Author: Dave Moore
# Purpose - support an analysis by Dave Bowling and Kenny Smith

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# Simulate daily temperatures for T_air
set.seed(42) # For reproducibility
days <- 1:365
mean_annual_temp <- 15 # Mean annual temperature in Celsius

# Simulate seasonal temperature variations for T_air
T_air <- mean_annual_temp + 10 * sin(2 * pi * (days - 80) / 365) + rnorm(365, 0, 5* 0.3)

# Simulate T_soil with 30% of the variance of T_air
T_soil <- mean_annual_temp + 7 * sin(2 * pi * (days - 80) / 365) + rnorm(365, 0, 5*0.3 )

# Create a dataframe with the temperature data
temperature_data <- data.frame(
  Day = days,
  T_air = T_air,
  T_soil = T_soil
)

# Define the varying C_leaf for each ecosystem type
C_leaf_DBF <- ifelse(days < 70 | days > 270, 0, 100 * sin(pi * (days - 70) / 200))
C_leaf_ENF <- 100 + (100 * 0.03) * (0.3 * (days / 365) + 0.7 * pmax(0, (days - 200) / 50))

# Define initial values for other carbon pools (in arbitrary units) for each ecosystem type
initial_values <- list(
  DBF = list(C_stem = 200, C_root = 150, C_litter = 50, C_microbial = 30, C_CWD = 70, C_fast_soil = 120, C_slow_soil = 400),
  ENF = list(C_stem = 220, C_root = 170, C_litter = 60, C_microbial = 35, C_CWD = 80, C_fast_soil = 130, C_slow_soil = 420)
)

# Define Q10 function for respiration
Q10_function <- function(C_pool, T, R_base, Q10) {
  R_base * C_pool * Q10^((T - 15) / 10)
}

# Set parameters for respiration calculation (arbitrary units) for each ecosystem type
params <- list(
  DBF = list(
    leaf = list(R_base = 0.01, Q10 = 2),
    stem = list(R_base = 0.02, Q10 = 2.1),
    root = list(R_base = 0.015, Q10 = 2.2),
    litter = list(R_base = 0.03, Q10 = 1.9),
    microbial = list(R_base = 0.05, Q10 = 2.3),
    CWD = list(R_base = 0.025, Q10 = 2),
    fast_soil = list(R_base = 0.02, Q10 = 2.5),
    slow_soil = list(R_base = 0.01, Q10 = 2.6)
  ),
  ENF = list(
    leaf = list(R_base = 0.011, Q10 = 2.2),
    stem = list(R_base = 0.021, Q10 = 2.15),
    root = list(R_base = 0.016, Q10 = 2.25),
    litter = list(R_base = 0.031, Q10 = 1.95),
    microbial = list(R_base = 0.051, Q10 = 2.35),
    CWD = list(R_base = 0.026, Q10 = 2.1),
    fast_soil = list(R_base = 0.021, Q10 = 2.55),
    slow_soil = list(R_base = 0.011, Q10 = 2.65)
  )
)

# Function to calculate respiration for each time step and ecosystem type
calculate_respiration <- function(T_air, T_soil, ecosystem, C_leaf) {
  C_pools <- initial_values[[ecosystem]]
  param <- params[[ecosystem]]
  
  R_leaf <- Q10_function(C_leaf, T_air, param$leaf$R_base, param$leaf$Q10)
  R_stem <- Q10_function(C_pools$C_stem, T_air, param$stem$R_base, param$stem$Q10)
  R_root <- Q10_function(C_pools$C_root, T_soil, param$root$R_base, param$root$Q10)
  R_litter <- Q10_function(C_pools$C_litter, T_air, param$litter$R_base, param$litter$Q10)
  R_microbial <- Q10_function(C_pools$C_microbial, T_soil, param$microbial$R_base, param$microbial$Q10)
  R_CWD <- Q10_function(C_pools$C_CWD, T_air, param$CWD$R_base, param$CWD$Q10)
  R_fast_soil <- Q10_function(C_pools$C_fast_soil, T_soil, param$fast_soil$R_base, param$fast_soil$Q10)
  R_slow_soil <- Q10_function(C_pools$C_slow_soil, T_soil, param$slow_soil$R_base, param$slow_soil$Q10)
  
  # Calculate R_below and R_above
  R_below <- R_root + R_microbial + R_fast_soil + R_slow_soil
  R_above <- R_leaf + R_stem + R_litter + R_CWD
  
  list(R_leaf = R_leaf, R_stem = R_stem, R_root = R_root, R_litter = R_litter,
       R_microbial = R_microbial, R_CWD = R_CWD, R_fast_soil = R_fast_soil,
       R_slow_soil = R_slow_soil, R_below = R_below, R_above = R_above)
}

# Apply the respiration calculation to each row of the temperature dataframe for both ecosystems
respiration_results_DBF <- temperature_data %>%
  rowwise() %>%
  mutate(
    Ecosystem = "DBF",
    C_leaf = C_leaf_DBF[Day],
    Resp = list(calculate_respiration(T_air, T_soil, "DBF", C_leaf))
  ) %>%
  unnest_wider(Resp)

respiration_results_ENF <- temperature_data %>%
  rowwise() %>%
  mutate(
    Ecosystem = "ENF",
    C_leaf = C_leaf_ENF[Day],
    Resp = list(calculate_respiration(T_air, T_soil, "ENF", C_leaf))
  ) %>%
  unnest_wider(Resp)

# Combine results
respiration_results <- bind_rows(respiration_results_DBF, respiration_results_ENF)

# Calculate the sum of R_above and R_below for each ecosystem type
total_respiration <- respiration_results %>%
  group_by(Ecosystem) %>%
  summarise(
    Sum_R_above = sum(R_above),
    Sum_R_below = sum(R_below)
  ) %>%
  mutate(Ratio_R_above_to_R_below = Sum_R_above / Sum_R_below)

# Plot the variation of R_above and R_below by Days for both ecosystems in stacked plots
plot_respiration <- ggplot(respiration_results, aes(x = Day)) +
  geom_line(aes(y = R_above, color = "R_above")) +
  geom_line(aes(y = R_below, color = "R_below")) +
  facet_wrap(~Ecosystem, ncol = 1) +
  labs(title = "Daily Variation of R_above and R_below by Ecosystem Type",
       x = "Day of Year",
       y = "Respiration Rate (arbitrary units)",
       color = "Respiration") +
  theme_minimal()

# Plot the ratio of R_above to R_below for the two ecosystems
plot_ratio <- ggplot(total_respiration, aes(x = Ecosystem, y = Ratio_R_above_to_R_below, fill = Ecosystem)) +
  geom_bar(stat = "identity") +
  labs(title = "Ratio of R_above to R_below",
       x = "Ecosystem",
       y = "Ratio of R_above to R_below") +
  theme_minimal()

# Combine the two plots side by side
grid.arrange(plot_respiration, plot_ratio, ncol = 2)
