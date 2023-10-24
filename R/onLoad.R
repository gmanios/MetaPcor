
.onLoad <- function(libname, pkgname) {

  # List of CRAN packages to install
  cran_packages <- c("readxl", "igraph", "visNetwork",
                     "plotly", "GGally", "ggrepel",
                     "gprofiler2", "matrixcalc",
                     "ff", "bit", "DescTools",
                     "forcats", "stringr",
                     "purrr", "readr",
                     "tidyr", "tibble", "ggplot2",
                     "tidyverse", "dplyr", "plyr",
                     "metafor", "metadat", "Matrix",
                     "meta", "GeneNet", "fdrtool",
                     "longitudinal", "data.table", "corpcor",
                     "devtools", "BiocManager") # also include devtools and BiocManager

  # Identify packages not yet installed
  new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]

  # Install missing CRAN packages
  if(length(new_cran_packages)) {
    install.packages(new_cran_packages, dependencies = TRUE)
  }

  # Load devtools package for installing from GitHub
  if("devtools" %in% installed.packages()[,"Package"]) {
    library(devtools)

    # GitHub packages
    github_packages <- list(
      c("sctyner", "geomnet"),
      c("cran", "space")
    )

    # Check and install missing GitHub packages
    for(pkg_info in github_packages) {
      owner <- pkg_info[1]
      repo <- pkg_info[2]
      if (!is.element(repo, installed.packages()[,"Package"])) {
        install_github(paste0(owner, "/", repo))
      }
    }
  }

  # Load BiocManager for installing from Bioconductor
  if("BiocManager" %in% installed.packages()[,"Package"]) {
    library(BiocManager)

    # Bioconductor packages
    bioc_packages <- c("GEOquery")

    # Check and install missing Bioconductor packages
    for(pkg in bioc_packages) {
      if (!is.element(pkg, installed.packages()[,"Package"])) {
        BiocManager::install(pkg)
      }
    }
  }
  options(repos = c(CRAN = "https://cran.rstudio.com/"))

  # library(GEOquery)
  # library(corpcor)
  # library(space)
  # library(data.table)
  # library(GeneNet)
  # library(meta)
  # library(metafor)
  # library(dplyr)
  # library(tidyverse)
  # library(purrr)
  # library(DescTools)
  # library(ggplot2)
  # library(ff)
  # library(matrixcalc)
  # library(gprofiler2)
  # library(visNetwork)
  # library(geomnet)
  # library(igraph)
  # library(plotly)
  message("MetaPcor has been loaded.")

}
