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
                            "pracma",
                            "rdist",
                            "tidyverse",
                            "RcppArmadillo",
                            "RcppParallel",
                            "stringr",
                            "dplyr",
                            "sf",
                            "raster",
                            "terra",
														"NLMR"))

# Pipeline ---------------------------------------------------------

list(
  
  ## Input raw data files -----  

  ### Input parameters file: -----------
  tar_target(parameters_txt,
             file.path("Parameters.txt"),
             format="file"),
  
  ### Input landscapes directory: -----------
  tar_target(lands_path,
             file.path("Input","lands"),
             format="file"),
  
  ### Input landscape predictions: -----------
  tar_target(landmat_path,
             file.path("Input","ldsel.rds"),
             format="file"),
  
  ### Input sample design: ---------
  # Only need if sample = 1
  #tar_target(sampling,
  #           file.path("Input", "sampling_scheme.csv"),
  #           format="file"),
  

  ## Read and format input data -----  
  tar_terra_sprc(plands_sprc, ReadLands(lands_path)), 
  tar_target(landmat,ReadRDS(landmat_path)),
  
  ### Read and format parameters file: -----------
  tar_target(parameters,FormatSetParameters(parameters_txt)),
  
  ### Read and format landscapes: -----------

  #tar_terra_sprc(lands_sprc, ReadLands(predlands_path)), 
  
  ### Read and format sampling design scheme: ---------
  # Only need if sample = 1
  #tar_target(sample.design,PrepSurveillance(sampling)),
	
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
             format="file"),
  
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
  tar_target(land_grid_list,InitializeGrids(plands_sprc,"heterogeneous")),
  #tar_target(land_grid_list,InitializeGrids(plands_sprc[1],"heterogeneous"))#,
  #tar_target(land_grid_list,InitializeGrids(c(parameters$len,parameters$inc),parameters$grid.opt))#,
  
  ## Run Model: ---------------
  #Use tar_force format here because otherwise will only run if code has been updated
  #initialize output objects
  tar_force(x,RunSimulationReplicates(land_grid_list, parameters), force=TRUE)
      #add nrep

  )



