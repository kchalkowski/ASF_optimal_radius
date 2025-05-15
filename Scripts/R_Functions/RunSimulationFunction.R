
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
RunSimulationReplicates<-function(land_grid_list, parameters){
	
	#Pull needed parms from parameters for initialize sounders
			len=parameters$len
			inc=parameters$inc
			density=parameters$density
			centroids=land_grid_list[[1]]$centroids
			grid=land_grid_list[[1]]$grid
			ss=parameters$ss
			
	#calc implicit parameters from input parameters
			km_len=len*inc
			area=len^2
			N0=parameters$density*area

pop=InitializeSounders(centroids,grid,c(N0,ss),pop_init_type="init_pop",pop_init_grid_opts="ras")

return(pop)
#out.list=SimulateOneRun(Pcr,Pir,Pbd,death,F1,F2_int,F2_B,F2i_int,F2i_B,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop,out.opts,"homogeneous",rep,DetP)
#return(out.list)
}

