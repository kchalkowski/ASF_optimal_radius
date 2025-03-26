# _targets.R

# Targets setup --------------------
setwd(this.path::this.dir())

#load libraries for targets script
#install.packages("geotargets", repos = c("https://njtierney.r-universe.dev", "https://cran.r-project.org"))
library(targets)
library(tarchetypes)
library(geotargets)

# This hardcodes the absolute path in _targets.yaml, so to make this more
# portable, we rewrite it every time this pipeline is run (and we don't track
# _targets.yaml with git)
tar_config_set(
  store = file.path(this.path::this.dir(),("_targets")),
  script = file.path(this.path::this.dir(),("_targets.R"))
)

#Source functions in pipeline
lapply(list.files(file.path("Scripts","Functions"), full.names = TRUE, recursive = TRUE), source)

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
  
  ## Read input data -----  
  #Examples:
  #tar_terra_sprc(plands_sprc, ReadLands(predlands_path)), 
  #tar_target(pdisp,ReadRDS(preddisp_path)),
  
  ## Input cpp scripts as files to enable tracking -----  
  #Example:
  #tar_target(cpp_func,
  #           "rel/path/to/script.cpp",
  #           format="file"),
  
  ## Initialize model -----
  ### Initialize grid(s): ---------------
    #Method for class 'SpatRaster'
        #Initialize_Grids(object)
    #Method for class 'SpatRasterCollection': 
        #Initialize_Grids(object)
    #Method for class 'numeric'
        #Initialize_Grids(object,grid.opts="")
    #Arguments
        #type- 
          #"homogenous" or "heterogeneous", if object is SpatRaster or SpatRasterCollection, type is automatically 'ras'
            #ras- 
              #use input raster to create grid
            #homogeneous- 
              #creates grid with 7 columns, no land class variables
            #heterogeneous-
              #creates a neutral random landscape model with X lc variables
    #tar_target(land_grids,InitializeGrids(lands_sprc)),
  
  
  
  )



