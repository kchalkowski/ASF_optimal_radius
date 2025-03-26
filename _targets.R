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
lapply(list.files(file.path("Scripts","R_Functions"), full.names = TRUE, recursive = TRUE), source)

#Load packages
tar_option_set(packages = c("Rcpp",
                            "R.matlab",
                            "pracma",
                            "rdist",
                            "tidyverse",
                            "RcppArmadillo",
                            "RcppParallel",
                            "stringr",
                            "dplyr",
                            "sf",
                            "terra"))

# Pipeline ---------------------------------------------------------

list(
  
  ## Input raw data files -----  

  ### Input parameters file: -----------
  tar_target(parameters_txt,
             file.path("Parameters.txt"),
             format="file"),
  
  ### Input landscapes directory: -----------
  #tar_target(lands_path,
  #           "/rel/path/ras_dir",
  #           format="file"),

  ## Read and format input data -----  
  #Examples:
  #tar_terra_sprc(plands_sprc, ReadLands(predlands_path)), 
  #tar_target(pdisp,ReadRDS(preddisp_path)),
  
  ### Read and format parameters file: -----------
  tar_target(parameters,FormatSetParameters(parameters_txt)),
  
  ### Read and format landscapes: -----------
  #tar_terra_sprc(lands_sprc, ReadLands(predlands_path)), 
  
  ## Input cpp scripts as files to enable tracking -----  
  tar_target(Fast_FOI_Matrix_script,
            file.path("Scripts","cpp_Functions","Fast_FOI_Matrix.cpp"),
            format="file"),
  tar_target(FindCellfromCentroid_script,
             file.path("Scripts","cpp_Functions","FindCellfromCentroid.cpp"),
             format="file"),
  tar_target(Movement_Fast_Generalized_script,
             file.path("Scripts","cpp_Functions","Movement_Fast_Generalized.cpp"),
             format="file"),
  tar_target(Movement_Fast_RSFavail_script,
             file.path("Scripts","cpp_Functions","Movement_Fast_RSFavail.cpp"),
             format="file"),
  tar_target(SpatialZones_fast_script,
             file.path("Scripts","cpp_Functions","SpatialZones_fast.cpp"),
             format="file")#,
  
  ## Initialize model -----
  ### Initialize grid(s): ---------------
    #Method for class 'SpatRaster'
        #Initialize_Grids(object)
    #Method for class 'SpatRasterCollection': 
        #Initialize_Grids(object)
    #Method for class 'numeric'
        #Initialize_Grids(object,grid.opts="homogenous")
    #Arguments
      #grid.opts- 
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
    
  ## Run Model: ---------------
  #Use tar_force format here because otherwise will only run if code has been updated
  #iniitalize output objects
  #tar_force(x,RunSimulation(land_grid_list, parameters, movement), force=TRUE)
      #add nrep
  
  #qs for Madison
    #help dev parms setup file? txt file good idea?--
      #surveillance options laid out
    #grid setup good? sep pipeline for formatting the grids? collab on that?
    #how to handle changing parms file while still tracking? gitignore?
    #thinking of removing grid.type=ML
    #experience making own grid with DIY opt, addl switches needed?
    #other switches/error catching identified?
  
  #updates needed
    #gitignore parameters.txt
    #time lag between first detect and culling
  )



