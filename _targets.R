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
library(crew)

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
														"clustermq",
														"deSolve"),
														error = 'stop') # for troubleshooting

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
		inputs=list() ## changed to align parameter values with states (they were mixed)
#     "state_basis" = data.frame("state"=c("FL","SC")),
# 	"density_ss"=data.frame(
# 		"density"=c(1.5,3,5), ## were these the wrong values? different from what's in Parameters.txt
# 		"ss"=c(2,4,6)#,
#         "B1"=c(0.9,0.4,0.2,0.009,0.004,0.002) ## B2 is calculated afterwards according to "B2_B1_factor" parameter
# 		),
# 	,"Radius"=data.frame( # comma at the beginning allows easy commenting out of pieces you don't want to use
# 		"Rad"=c(5,10,15,20)
# 		)
# 			)
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
#   tar_target(land_grid_list,InitializeGrids(plands_sprc,"heterogeneous")),

	#single landscape:
	#tar_target(land_grid_list,InitializeGrids(plands_sprc[1],"heterogeneous"))#,

	#homogenous grid:
# 	tar_target(land_grid_list,InitializeGrids(c(parameters00$len,parameters00$inc),parameters00$grid.opt)),

    ## makes the grid.opt parameters functional to choose lands variation... shifted grid.opts = "heterogeneous" to be randomized landscape (was "random" before, but unlisted)
    ## ugly, but functional:
	tar_target(land_grid_list, {if (parameters00$pop_init_grid_opts == 'homogeneous'){
                                  if(parameters00$grid.opts != 'ras'){ ## if grid.opts is homogeneous or heterogeneous
                                    ## make a grid either uniform or random with even initial pig locations
                                    InitializeGrids(c(parameters00$len,parameters00$inc),parameters00$grid.opt)
                                  } else if (parameters00$grid.opts == 'ras'){ ## if there is an input raster
                                    InitializeGrids(plands_sprc, parameters00$grid.opts == 'ras')
                                  }
                                } else if (parameters00$pop_init_grid_opts == 'heterogeneous'){
                                  ## make a grid with uneven pig initial locations...
                                    if (parameters00$grid.opts == 'homogeneous') {
                                      ## can't do neutral plane with random pig distribution
                                      stop('Cannot run homogeneous grid.opts with heterogeneous pop_init_grid_opts')
                                    } else if (parameters00$grid.opts == 'heterogeneous'){
                                      ## random pig distribution with random landscape
                                      InitializeGrids(c(parameters00$len,parameters00$inc),parameters00$grid.opt)
                                    } else if (parameters00$grid.opts == 'ras'){
                                      ## random pig distribution with raster landscape
                                      InitializeGrids(plands_sprc, parameters00$grid.opts == 'ras')
                                    }
                                }
                              }),
	
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
  		reps=2
  		)
  	)#,
	
	## Run MF Model ---------------
	#tar_target(
 # 	out.list,
 # 	RunMFModel(
 # 		Ct=Ct, 
#			Beta=Beta,	
# 		parameters=parameters,
#			variables=variables		
#  		)
# 	)
	
  )

  


