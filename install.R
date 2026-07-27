# R package dependencies for figure scripts
# Run from crypt/: Rscript install.R

install.packages(
  c("ggplot2", "dplyr", "tidyr", "stringr",
    "cowplot", "ggraph", "igraph", "ape", "scales"),
  repos = "https://cloud.r-project.org"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("ggtree")
