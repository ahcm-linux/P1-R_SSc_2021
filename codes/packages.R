message("Installing and loading required packages ...")

suppressMessages({
  load_lib <- c(
    # data wrangling
    "reshape2",
    "plyr",
    "dplyr",
    # FA and SEM
    "psych",
    "GPArotation",
    "lavaan",
    # CCA
    "mixedCCA",
    # LDA
    "HiDimDA",
    # regression
    "gamlss",
    # data visualization
    "ggplot2",
    "gridExtra",
    "scales",
    "pheatmap",
    # missing data
    "mice",
    # data description
    "MVN",
    "DescTools"
  )
  
  install_lib <- load_lib[!load_lib %in% installed.packages()]
  for (lib in install_lib) install.packages(lib, dependencies = TRUE)
  loading <- sapply(load_lib, require, character = TRUE)
})

if (!any(loading)) {
  message("Packages not working: ",  load_lib[!loading])
  stop("The above packages were not properly installed.")
} else {
  message("Packages successfully loaded")
}