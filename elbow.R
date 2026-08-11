library(readr)
library(dplyr)
library(ggplot2)

# Usage: Rscript elbow.R [max_K]
# Example: Rscript elbow.R 12
# Reads nmf_error_curve.csv and nmf_out_k<max_K>/W.populations_by_K.csv
args <- commandArgs(trailingOnly = TRUE)
max_K <- if (length(args) >= 1) as.integer(args[1]) else stop("Usage: Rscript elbow.R <max_K>")

# --- Elbow plot ---
e <- read_csv("nmf_error_curve.csv", show_col_types = FALSE)

p1 <- ggplot(e, aes(K, error)) +
  geom_line() +
  geom_point(size = 3) +
  theme_bw() +
  ggtitle("NMF reconstruction error vs K")

print(p1)

# --- Check sparsity at max K ---
W_path <- paste0("nmf_out_k", max_K, "/W.populations_by_K.csv")
if (file.exists(W_path)) {
  W <- read_csv(W_path, show_col_types = FALSE)
  cat("Fraction of populations with > 0.4 weight per component:\n")
  print(apply(W[, -1], 2, function(x) mean(x > 0.4)))
}

# --- Marginal error reduction ---
e <- e %>%
  mutate(drop = lag(error) - error)

p2 <- ggplot(e, aes(K, drop)) +
  geom_col() +
  theme_bw() +
  ggtitle("Marginal error reduction per added component")

print(p2)
