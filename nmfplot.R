library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)


# Usage: Rscript nmfplot.R [K]
# Example: Rscript nmfplot.R 6
args <- commandArgs(trailingOnly = TRUE)
K <- if (length(args) >= 1) as.integer(args[1]) else stop("Usage: Rscript nmfplot.R <K>")


# ---- load W (population x K proportions) ----
W <- read_csv(paste0("nmf_out_k",K,"/W.populations_by_K.csv"))


# Ensure population is first column name
# If needed:
# W <- W %>% rename(population = 1)

W_mat <- W %>% select(-population) %>% as.matrix()
rownames(W_mat) <- W$population

# ---- cosine distance + clustering ----
# cosine distance = 1 - cos_sim
cos_sim <- function(A) {
  A <- A / sqrt(rowSums(A^2))
  A %*% t(A)
}
D <- 1 - cos_sim(W_mat)
hc <- hclust(as.dist(D), method = "average")
pop_order <- hc$labels[hc$order]

# ---- long format + plot ----
W_long <- W %>%
  pivot_longer(-population, names_to = "component", values_to = "prop") %>%
  mutate(
    population = factor(population, levels = pop_order),
    component  = factor(component, levels = unique(component))
  )

p <- ggplot(W_long, aes(x = population, y = prop, fill = component)) +
  geom_col(width = 0.9) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Mixture proportion", fill = "Component",
       title = "NMF admixture proportions (populations ordered by similarity)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

print(p)

# ---- save ----
ggsave(paste0("nmf_admixture_stackedbars.", K, ".png"), p, width = 12, height = 6, dpi = 300)
