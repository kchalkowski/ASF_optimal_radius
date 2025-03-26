# _targets.R

# Targets setup --------------------
setwd(this.path::this.dir())

#load libraries for targets script
library(targets)
library(tarchetypes)
#install.packages("geotargets", repos = c("https://njtierney.r-universe.dev", "https://cran.r-project.org"))
library(geotargets)


# here::here() returns an absolute path, which then gets stored in tar_meta and
# becomes computer-specific (i.e. /Users/andrew/Research/blah/thing.Rmd).
# There's no way to get a relative path directly out of here::here(), but
# fs::path_rel() works fine with it (see
# https://github.com/r-lib/here/issues/36#issuecomment-530894167)
#here_rel <- function(...) {fs::path_rel(here::here(...))}

# Set the _targets store so that scripts in subdirectories can access targets
# without using withr::with_dir() (see https://github.com/ropensci/targets/discussions/885)
#
# This hardcodes the absolute path in _targets.yaml, so to make this more
# portable, we rewrite it every time this pipeline is run (and we don't track
# _targets.yaml with git)
tar_config_set(
  store = file.path(this.path::this.dir(),("_targets")),
  script = file.path(this.path::this.dir(),("_targets.R"))
)

#Source functions in pipeline
lapply(list.files(file.path("1_Scripts","Functions"), full.names = TRUE, recursive = TRUE), source)

#Load packages
tar_option_set(packages = c("Rcpp",
                            "R.matlab",
                            "pracma",
                            "rdist",
                            "tidyverse",
                            "RcppArmadillo",
                            "RcppParallel",
                            "stringr",
                            "dplyr"))

# Pipeline ---------------------------------------------------------

list(
  
  ## Input raw data files -----  
  #Example:
  #tar_target(predlands_path,
  #           "/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Not_Pipeline/Setup/Pipeline_SSF_Weekly/4_Output/indiv_plands",
  #           format="file"),
  
  ## Read data -----  
  #Example:
  #tar_terra_sprc(plands_sprc, ReadLands(predlands_path)), 
  #tar_target(pdisp,ReadRDS(preddisp_path)),
  
  ## Format data -----
  #Need to average probs across season for each grid for nnd calcs
  
  
)



