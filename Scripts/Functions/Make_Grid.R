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
#2- cells, integer of number of cells in grid
#3- centroids, two col x/y of just centroids of grid of each cell

#area=80^2
#inc=0.4
#0.4*200

#Note: stuck on running clustered simulated landscape with NLMR due to pkg install issues
# install.packages("remotes")
#remotes::install_github("cran/RandomFieldsUtils")
#remotes::install_github("cran/RandomFields")
#remotes::install_github("ropensci/NLMR")

#require(NLMR)
#Make_Grid(100,0.4,"homogenous")
#len=100
#inc=0.4
#grid.opt="homogenous"
###
#troubleshooting:
#grid<-readMat(paste0(home,"/Input/Grid_80x80_0pt4km.mat"))
#gridml<-grid$grid
#centroidsml=gridml[,6:7]
###

Make_Grid<-function(len,inc,grid.opt){
  
  #get number of cells in grid
  cells=len^2
  
  #if grid homogenous-- will later enter ability to alter LULC
  if(grid.opt=="homogenous"){
  
  #initialize empty grid matrix
  grid=matrix(nrow=round(cells),ncol=7)
  
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
  }

  
  #if(grid.opt=="test_pref_binary"){
  #  grid=matrix(nrow=round(cells),ncol=8)
  #  grid[,1]=1:cells
  #  grid[,2]=rep(seq(0,((inc*nrow.grid)-inc),0.4),nrow.grid)
  #  grid[,3]=rep(seq(0,((inc*nrow.grid)-inc),0.4),nrow.grid)
  #  grid[,4]=rep(seq(inc,(inc*nrow.grid),0.4),nrow.grid)
  #  grid[,5]=rep(seq(inc,(inc*nrow.grid),0.4),nrow.grid)
  #  grid[,6]=rep(seq(((0+inc)/2),(((inc*nrow.grid)-inc)+(inc*nrow.grid))/2,0.4),nrow.grid)
  #  grid[,7]=rep(seq(((0+inc)/2),(((inc*nrow.grid)-inc)+(inc*nrow.grid))/2,0.4),nrow.grid)
  #  grid[,8]=rbinom(round(cells),1,0.2)
  #  centroids=grid[,c(6,7)]
  #}
  
  grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
  return(grid.list)
  
}








