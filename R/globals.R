# R/globals.R
utils::globalVariables(c(
  # Existing variables
  "study", "treatment", "time", "survival", "max_time", "n_times",
  "lower", "upper", "probability",

  # Add these new ones from the NOTE
  "x", "y", "x1", "x2", "y1", "y2",
  "weight", "tau", "lower_50", "upper_50",
  "iteration", "value", "chain"
))
