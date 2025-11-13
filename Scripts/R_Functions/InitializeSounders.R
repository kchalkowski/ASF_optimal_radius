# Purpose ------

#Initializes the population matrix for spatial meta-population model
#Also includes option to initialize individual/group in single location on grid
#Option to incorporate heterogeneous landscape preference is in progress

# Function ------

#Inputs: 
#centroids: numeric matrix, x y coordinates of centroids of every cell in grid
#grid: matrix with coordinates for all cells in grid
#pop_init_args:
  #initialization arguments
  #for pop_init_type="init_pop", need a vector with N0 (initial pop size) and ss (sounder size)
  #for pop_init_type="init_single", need a vector with init_loc (cell number to initialize group/individual) and n (number of individuals to initialize)
#pop_init_type: string, "init_pop" or "init_single"
#pop_init_grid_opts: string, "homogeneous" or "ras" or "heterogeneous"
InitializeSounders<-function(centroids,grid,pop_init_args,pop_init_grid_opts="homogeneous",pop_init_type="init_pop",RSF0_lc=NULL){ ## changed argument order to allow easier pop_init_grid_opts without defining pop_init_type
	
	## Parse input args and check input formatting ------------
	
  #Check main input args
  if(pop_init_type=="init_pop"){
    if(length(pop_init_args)<2){
      stop("Missing input args to initialize population")
    }
    N0=pop_init_args[1]
    ss=pop_init_args[2]
    
  } else{
    
    #for init_pop, init_args needs to be vector with init_locs, n
    if(pop_init_type=="init_single"){
      init_locs=pop_init_args[1]
      n=pop_init_args[2]
      if(init_locs!=round(init_locs)|n!=round(n)){
        stop("Initial cell number and number to initialize need to be integers")
      }
    
    } else{
      stop("Unrecognized population initialization type")
    }
    
  }
	
  #Check grid input args
    if(pop_init_grid_opts=="heterogeneous"&ncol(grid)<8){
      stop("Specified heterogeneous land class preference option, but land class values missing in grid.")
    }
    if(pop_init_grid_opts!="heterogeneous" & 
       pop_init_grid_opts!="homogeneous" & 
       pop_init_grid_opts!="ras"){
      stop("Unrecognized grid input option")
    }


	## Initialize population ----------------------

  if(pop_init_type=="init_pop"){ 
  
    #initialize needed objects
    cells=nrow(centroids)
  
    #Get the initializing number of sounders
    sn_i<-N0/ss 
  
    #default option, randomly assign sounders to cells
    if(pop_init_grid_opts=="homogeneous"|
       pop_init_grid_opts=="ras"){
    assigns<-rbinom(cells,1,sn_i/cells) 
    } else{

    if(pop_init_grid_opts=="heterogeneous"){
    	
      #use this to weight preference so that still end up with N0 size population
      pref.wt=sum(grid[,8])/cells 
      #assign to cells with weighted preference according to column 8 values
      assigns=rbinom(cells,1,grid[,8]*((sn_i/cells)/pref.wt))
    }
    
    }

  
  sn<-sum(assigns) #generated number of rows (sounders), size from random assignment of sounders to cells
  init_locs<-which(assigns==1) #get the locations where sounders have been initialized
  
  #Initialize the sounder population matrix
  #each row is a sounder
  pop<-matrix(nrow=sn,ncol=13)
  pop[,1]=rpois(sn,ss) #sounder size with avg as lambda in a poisson
  pop[,3]=init_locs #this will be grid location (row number)
  pop[,4]=0 #this will be assigned movement distance
  pop[,5]=centroids[pop[,3],1] #present location X 
  pop[,6]=centroids[pop[,3],2] #present location Y
  pop[,7]=0 #previous location (grid row number)
  pop[,8]=pop[,1] #number of S status in sounder
  pop[,9]=0 #number of E status in sounder
  pop[,10]=0 #number of I status in sounder
  pop[,11]=0 #number of R status in sounder
  pop[,12]=0 #number of C status in sounder
  pop[,13]=0 #number of Z status in sounder

  #for homogeneous grid, pref col is just uniform 0
  if(pop_init_grid_opts=="homogeneous"){pop[,2]=0} 

  #for heterogeneous or ras grid, pref col indicates preference val of current cell
  if(pop_init_grid_opts=="heterogeneous"|pop_init_grid_opts=="ras"){pop[,2]=grid[pop[,3],8]} 
  
  if(any(pop[,3]>nrow(centroids))){
    stop("agents initialized off the grid")
  }
  
  ## Deal with raster placements, keep realistic ------------
  #RSF0_lc value is lc with 0 probability to move to (e.g., water body)

  if((pop_init_grid_opts == "heterogeneous" | 
      pop_init_grid_opts == "ras") & !missing(RSF0_lc)){
      #stop condition. will need to do some checking before runs to look at center of each lc 
      #for initializing infected individual.
      #stop("individual initialized in unsuitable location (RSF probability = 0)")
    
      #check if any initialized in lc=0
      #while any initialized in lc=0
        #subset those and randomly sample new location not already occupied by another sounder
      if(any(pop[,2]==RSF0_lc)){
        RSF0.r=which(pop[,2]==RSF0_lc)
        
        for(r in 1:length(RSF0.r)){
          cellsq=1:cells
          cellsq=cellsq[-pop[,3]]
          if(all(centroids[cellsq,3]==RSF0lc)){
            stop("No suitable habitat available")
          }
          cellsq=cellsq[-which(centroids[cellsq,3]==RSF0lc)]
          new.cell=sample(cellsq,1)
          pop[RSF0.r[r],3]=new.cell
          pop[RSF0.r[r],5]=centroids[new.cell,1] #x
          pop[RSF0.r[r],6]=centroids[new.cell,2] #y
          pop[RSF0.r[r],2]=centroids[new.cell,3] #lc
          }
        
      }

    }

  }

	
	## Initialize single group/individual ---------------------

  if(pop_init_type=="init_single"){
   
    if((pop_init_grid_opts == "heterogeneous" | 
        pop_init_grid_opts == "ras") & !missing(RSF0_lc)){
      
    	#RSF0_lc value is lc with 0 probability to move to (e.g., water body)
      if(centroids[init_locs,3]==RSF0_lc){
        #stop condition. will need to do some checking before runs to look at center of each lc 
        #for initializing infected individual.
        stop("individual initialized in unsuitable location (RSF probability = 0)")
      }
      
    }

  #for initializing initial infected individual introduction
  pop<-matrix(nrow=1,ncol=13)
  pop[,1]=n #sounder size with avg as lambda in a poisson
  pop[,2]=0 
  pop[,3]=init_locs #this will be grid location (row number)
  pop[,4]=0 #this will be assigned movement distance
  pop[,5]=centroids[pop[,3],1] #present location X 
  pop[,6]=centroids[pop[,3],2] #present location Y
  pop[,7]=0 #previous location (grid row number)	
  pop[,8]=n #number of S status in sounder
  pop[,9]=0 #number of E status in sounder
  pop[,10]=0 #number of I status in sounder
  pop[,11]=0 #number of R status in sounder
  pop[,12]=0 #number of C status in sounder
  pop[,13]=0 #number of Z status in sounder
  
  	}

	
	## Tidying outputs -----------
	
  #add column names
  colnames(pop)=c("Nlive","pref","cell","dist","ctrx","ctry","pcell","S","E","I","R","C","Z")
  
  #remove rows with no pigs
  pop=pop[pop[,1]>0,,drop=FALSE]

  #error catches
  if(any(pop[,3]>nrow(centroids))){
    stop("agents initialized off the grid")
  }

  return(pop)

  }

