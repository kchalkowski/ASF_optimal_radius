# _targets.R

#Friday May 16, get parameters file set up to be able to loop
	#create function to make matrices with a row for each parameter setting
		#inputs: parameters, list of inputs, structured so that each unit is its own setting
			#check that list of inputs, names match parms with 'input' label.. otherwise do warning/stop
			#create matrix with all combinations of parameters
			#output variable parameter matrix

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

#set options
options(clustermq.scheduler="multicore")

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
														"NLMR",
														"EnvStats",
														"clustermq"))

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
  
  ## Read and format input data -----  
  tar_terra_sprc(plands_sprc, ReadLands(lands_path)), 
  tar_target(landmat,ReadRDS(landmat_path)),
  
  ### Read and format parameters file: -----------
  tar_target(parameters0,FormatSetParameters(parameters_txt)),
	tar_target(variables,SetVarParms(parameters0,
		inputs=list(
	"shape_rate"=data.frame(
		"shape"=c(0.7515,0.5657),
		"rate"=c(0.3550,1.9082)
		),
	"density_ss_B1_B2"=data.frame(
		"density"=c(0.1,0.3,0.5),
		"ss"=c(2,4,6),
		"B1"=c(0.009,0.009,0.009),
		"B2"=c(0.009*2,0.009*2,0.009*2)
		))
		)),
	tar_target(parameters00,RemoveRedundantParms(parameters0)),
  
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
  
  ## Initialize surface -----
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
	#multiple landscapes:
  tar_target(land_grid_list,InitializeGrids(plands_sprc,"heterogeneous")),
	
	#single landscape:
	#tar_target(land_grid_list,InitializeGrids(plands_sprc[1],"heterogeneous"))#,
  
	#homogenous grid:
	#tar_target(land_grid_list,InitializeGrids(c(parameters$len,parameters$inc),parameters$grid.opt))#,
	
	### Get surface parameters: ---------------
	tar_target(parameters,GetSurfaceParms(parameters00,plands_sprc[1])),
	
  ## Run Model ---------------
  #Use tar_force format here because otherwise will only run if code has been updated
  #tar_force(
	tar_target(
  	out.list,
  	RunSimulationReplicates(
  		land_grid_list=land_grid_list, 
  		parameters=parameters,
			variables=variables,		
  		cpp_functions=
  			list(Fast_FOI_Matrix_script,
  			Movement_Fast_Generalized_script),
  		reps=3
  		)
  	)
  )

  


