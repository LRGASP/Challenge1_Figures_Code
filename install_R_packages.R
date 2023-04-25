library(devtools)

# functions to install a package if it's missing
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}
install_github_if_missing <- function(package, repo) {
  if (!requireNamespace(package, quietly = TRUE)) {
    devtools::install_github(repo)
  }
}

packages <- c(
        "ComplexUpset",
        "MetBrewer",
        "UpSetR",
        "dplyr",
        "fmsb",
        "forcats",
        "gg.gap",
        "ggcorrplot",
        "ggplot2",
        "ggplot2movies",
        "ggpmisc",
        "ggpubr",
        "ggrepel",
        "grid",
        "gridExtra",
        "gt",
        "hrbrthemes",
        "huxtable",
        "patchwork",
        "reshape2",
        "scales",
        "tidyr",
        "viridis"
)

lapply(packages, install_if_missing)

install_github_if_missing("ggbreak", "YuLab-SMU/ggbreak")
install_github_if_missing("RColorConesa", "ConesaLab/RColorConesa")

# Load the packages
lapply(packages, library, character.only = TRUE)
