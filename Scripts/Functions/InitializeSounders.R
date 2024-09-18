
#########################
######## Purpose ######## 
#########################

#Initializes the population matrix for spatial meta-population model
#Also includes option to initialize individual/group in single location on grid
#Option to incorporate heterogeneous landscape preference is in progress

##########################
######## Function ######## 
##########################

#Inputs: 
#centroids: numeric matrix, x y coordinates of centroids of every cell in grid
#grid: matrix with coordinates for all cells in grid
#pop_init_args:
  #initialization arguments
  #for pop_init_type="init_pop", need a vector with N0 (initial pop size) and ss (sounder size)
  #for pop_init_type="init_single", need a vector with init_loc (cell number to initialize group/individual) and n (number of individuals to initialize)
#pop_init_type: string, "init_pop" or "init_single"
#pop_init_grid_opts: string, "homogeneous" or "heterogeneous"
InitializeSounders<-function(centroids,grid,pop_init_args,pop_init_type,pop_init_grid_opts){

  #########################################################################
  ############## Parse input args and check input formatting ############## 
  #########################################################################  
  
  
  if(missing(pop_init_type)){
    pop_init_type="init_pop"
  }
  
  if(missing(pop_init_grid_opts)){
    pop_init_grid_opts="homogeneous"
  }
  
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
    if(pop_init_grid_opts=="heterogenous"&ncol(grid)<8){
      stop("Specified heterogeneous land class preference option, but land class values missing in grid.")
    }
    if(pop_init_grid_opts!="heterogeneous" & pop_init_grid_opts!="homogeneous"){
      stop("Unrecognized grid input option")
    }
  

  ###################################################
  ############## Initialize Population ############## 
  ###################################################

  if(pop_init_type=="init_pop"){ 
  
    #initialize needed objects
    cells=nrow(centroids)
  
    #Get the initializing number of sounders
    sn_i<-N0/ss 
  
    #default option, randomly assign sounders to cells
    if(pop_init_grid_opts=="homogeneous"){
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

  #for homogenous grid, pref col is just uniform 0
  if(pop_init_grid_opts=="homogeneous"){pop[,2]=0} 

  #for heterogenous grid, pref col indicates preference val of current cell
  if(pop_init_grid_opts=="heterogeneous"){pop[,2]=grid[pop[,3],8]} 

  if(any(pop[,3]>nrow(centroids))){
    stop("agents initialized off the grid")
  }

  }

  ################################################################
  ############## Initialize Single Group/Individual ############## 
  ################################################################
  
  if(pop_init_type=="init_single"){
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

  
  #############################################
  ############## Tidying outputs ############## 
  #############################################
  
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

