
### Initialize grid(s): ---------------
    #Method for class 'list'
        #Initialize_Grids(object, parameters=parameters, movement=c(parameters$shape,parameters$rate))
    #Arguments
				#movement
					#default is vector with gamma distribution shape and rate fed from parameters file

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
RunSimulationReplicates<-function(land_grid_list, 
																	parameters, 
																	variables,
																	cpp_functions,
																	reps){
	
	#Pull needed parms from parameters for all reps
			list2env(parameters, .GlobalEnv)
	
#Need nested loops:
	#1, loop through all landscapes
	#2, loop through all parameter settings

		#later, loop through this
		centroids=land_grid_list[[1]]$centroids
		grid=land_grid_list[[1]]$grid

		#later, loop through this
		vars=variables[1,]
		names(vars)[3]="dens"
		vars=as.list(vars)
		list2env(vars, .GlobalEnv)

		#calc vals based on variables
		N0=dens*area
		K=N0*1.5

	pop=InitializeSounders(centroids,grid,c(N0,ss),pop_init_type="init_pop",pop_init_grid_opts="ras")
	outputs=Initialize_Outputs(parameters)
	pop=InitializeInfection(pop,centroids,grid,parameters)
	
	for(r in 1:reps){
	out.list=SimulateOneRun(outputs,pop,centroids,grid,parameters,cpp_functions,K)
	
	}

return(out.list)
}










