
#########################
######## Purpose ######## 
#########################

#Initializes grid to run spatially-explicit meta-population disease spread model
#Includes options for either homogeneous, uniform grid, or grid with land class preference probabilities
  #Currently only option is random landscape, alternative clustering options will be available in future updates.

##########################
######## Function ######## 
##########################

#Inputs: area of grid (in kilometers2); resolution of grid (in km2); type.opt="homogenous" or "clustered"
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
Make_Grid<-function(len,inc,grid.opt){
  require(raster)
  require(NLMR)
  
  if(missing(grid.opt)){
    grid.opt=="homogeneous"
  } 
  
  #get number of cells in grid
  cells=len^2
  
  #if grid homogenous-- will later enter ability to alter LULC
  if("homogenous"%in%grid.opt){
  
  #initialize empty grid matrix
  grid=matrix(nrow=round(cells),ncol=7)
  
  }
  
  if(!("homogenous"%in%grid.opt)){
    
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

  if(!("homogenous"%in%grid.opt)){
    
    #simulates a spatially random neutral landscape model with values drawn from a uniform distribution
    #values rescaled to range from 0-1
    if("random"%in%grid.opt){
    r=NLMR::nlm_random(len,len,inc,rescale=TRUE)
    grid[,8]=round(values(r),2)
    
    
    centroids=cbind(centroids,grid[,8])
    grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids,"r"=r)
    
    }
    
  } else{
    grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
  }
  
  return(grid.list)
  
}


set.seed(1984)
m = matrix(sample.int(25,25), 5)

col_indexes = c(1, 3)
col_indexes
min_loc_r = which(m[, col_indexes] == min(m[, col_indexes]), arr.ind = TRUE)




