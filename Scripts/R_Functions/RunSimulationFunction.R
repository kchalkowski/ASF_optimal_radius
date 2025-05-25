
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
for(v in 1:nrow(variables)){
	#v=1
		print(paste0("variable ",v))
		#later, loop through this
		vars=variables[v,]
		names(vars)[3]="dens"
		vars=as.list(vars)
		list2env(vars, .GlobalEnv)

		#calc vals based on variables
		N0=dens*area
		K=N0*1.5
		
		#loop through landscapes
		for(l in 1:length(land_grid_list)){
		print(paste0("landscape: ",l))
		centroids=land_grid_list[[l]]$centroids
		grid=land_grid_list[[l]]$grid

	pop=InitializeSounders(centroids,grid,c(N0,ss),pop_init_type="init_pop",pop_init_grid_opts="ras")
	outputs=Initialize_Outputs(parameters)
	pop=InitializeInfection(pop,centroids,grid,parameters)
	
	for(r in 1:reps){
	print(r)
	#Do simulations
	out.list=SimulateOneRun(outputs,pop,centroids,grid,parameters,cpp_functions,K)
	#print("test")
	#Handle outputs
	
		#Handle effective removal rate
	Ct.r=out.list$Ct
	Ct.r=cbind(1:thyme,Ct.r)
	Ct.r=cbind(rep(l,times=nrow(Ct.r)),Ct.r)
	Ct.r=cbind(rep(v,times=nrow(Ct.r)),Ct.r)
	colnames(Ct.r)=c("var","land","thyme","Ct")

		#Handle sounderlocs
	solocs.r=sounderlocsSummarize(out.list$sounderlocs,1)
	solocs.r=solocs.r$SEIRCZ_total
	solocs.r=cbind(rep(l,times=nrow(solocs.r)),solocs.r)
	solocs.r=cbind(rep(v,times=nrow(solocs.r)),solocs.r)
	colnames(solocs.r)[c(1,2)]=c("var","land")
	if(r==1&l==1&v==1){
		solocs=solocs.r
		Cto=Ct.r
	} else{
		solocs=rbind(solocs,solocs.r)
		Cto=rbind(Cto,Ct.r)
	}
	}
		
		}
}
		return(list("solocs"=solocs,"Cto"=Cto))
}










