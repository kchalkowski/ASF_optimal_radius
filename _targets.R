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
  
  ### Input parameters file: -----------
  #tar_target(parameters_txt,
  #           "/rel/path/parms.txt",
  #           format="file"),
  
  ## Read and format input data -----  
  #Examples:
  #tar_terra_sprc(plands_sprc, ReadLands(predlands_path)), 
  #tar_target(pdisp,ReadRDS(preddisp_path)),
  
  ### Read and format parameters file: -----------
  #tar_target(parameters,
              #FormatParameters(parameters_txt)),
  
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
        #Initialize_Grids(object,grid.opts="homogenous")
    #Arguments
      #type- 
        #"homogenous" or "heterogeneous", default for class numeric is "homogenous" and default for type SpatRaster and SpatRasterCollection is "ras"
          #ras- 
            #use input raster to create grid
          #homogeneous- 
            #creates grid with 7 columns, no land class variables
          #heterogeneous-
            #creates a neutral random landscape model with X lc variables
    #Value
      #a nested list of grid parameters
        
  #tar_target(land_grid_list,InitializeGrids(lands_sprc)),
    
  ### Run Model: ---------------
  #Use tar_force format here because otherwise will only run if code has been updated
  #tar_force(x,RunSimulation(lands_sprc), force=TRUE)

  
  )



