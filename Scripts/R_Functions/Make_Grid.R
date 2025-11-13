
#########################
######## Purpose ######## 
#########################

#Initializes grid to run spatially-explicit meta-population disease spread model
#Includes options for either homogeneous, uniform grid, or grid with land class preference probabilities
  #Currently only option is random landscape, alternative clustering options will be available in future updates.

##########################
######## Function ######## 
##########################

#Inputs: 
  #len of each side (in cells);
  #resolution of grid (in km2); 
  #grid.opt="homogeneous","random", or "ras"
  #ras, raster object, optional input

#Outputs: 
#1-grid matrix: 
  #col1- list of cells 1:ncell
  #col2- left-most x coord of cell
  #col3- upper-most y coord of cell
  #col4- right-most x coord of cell
  #col5- lower-most y coord of cell
  #col6- center x coordinate of cell
  #col7- center y coordinate of cell
  #if grid.opt="heterogeneous":
  #col8- preference probability
#2- cells, integer of number of cells in grid
#3- centroids, two col x/y of just centroids of grid of each cell
#Make_Grid(len,inc,grid.opt)


#input options
#required: object
  #raster, or c(len, inc)
#optional: grid.opt, 
#Make_Grid(object)

Make_Grid<-function(object,grid.opt="homogeneous",sample=0,sample.design=NULL){
  require(terra)
  require(NLMR)
  
  if(grid.opt=="homogenous"){
    grid.opt="homogeneous"
  }

  if(class(object)=="numeric"){
    len=object[1]
    inc=object[2]
  }
  
  if(class(object)=="SpatRaster"){
    ras=object
    len=dim(ras)[1]
    inc=res(ras)[1]/1000
    grid.opt="heterogeneous"
    
    if(dim(ras)[1]!=dim(ras)[2]){stop("raster needs to be square")}
      if(dim(ras)[1]!=len){
        stop("dimensions of raster do not match input len")
      }
  
  }

  
  #get number of cells in grid
  cells=len^2
  
  #if grid homogeneous-- will later enter ability to alter LULC
  if("homogeneous"%in%grid.opt & sample != 1){
    #initialize empty grid matrix
    grid=matrix(nrow=round(cells),ncol=7)
  } else {
    #initialize empty grid matrix
    grid=matrix(nrow=round(cells),ncol=8)
  }
  
  #first column is just cell indices
  grid[,1]=1:cells
  
  #Top left X coordinate of each cell
  grid[,2]=rep(seq(0,((inc*len)-inc),inc),times=len)

  #Top left Y coordinate of each cell
  grid[,3]=rep(seq(0,((inc*len)-inc),inc),each=len)
  
  #Top right X coordinate of each cell
  grid[,4]=rep(seq(inc,(inc*len),inc),times=len)

  #Top right Y coordinate of each cell
  grid[,5]=rep(seq(inc,(inc*len),inc),each=len)
  
  #Center X coordinate of each cell
  grid[,6]=rep(seq(((0+inc)/2),(((inc*len)-inc)+(inc*len))/2,inc),times=len)

  #Center Y coordinate of each cell
  grid[,7]=rep(seq(((0+inc)/2),(((inc*len)-inc)+(inc*len))/2,inc),each=len)

  #get centroids-only object
  centroids=grid[,c(6,7)]

  if(!("homogeneous"%in%grid.opt & sample != 1)){
    if(class(object)=="SpatRaster"){
      #need to get values from ras
      grid[,8]=round(values(ras),2)
      #assign to centroids
      centroids=cbind(centroids,grid[,8])
      grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
#     } else if("random"%in%grid.opt){ ## added to test for random grid.opt only if there is not "SpatRaster" object
    } else if("heterogeneous"%in%grid.opt){ ## instead of random? seems like this would be more updated

    #simulates a spatially random neutral landscape model with values drawn from a uniform distribution
    #values rescaled to range from 0-1
      r=NLMR::nlm_random(len,len,inc,rescale=TRUE)
      grid[,8]=round(values(r),2)
      centroids=cbind(centroids,grid[,8])
      grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids,"r"=r)

    }

  } else {
    grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
  }

  return(grid.list)
  
}


