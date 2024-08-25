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
#Sent ticket to IT on Apr 9, follow up if no resolution in a week
# install.packages("remotes")
#remotes::install_github("cran/RandomFieldsUtils")
#remotes::install_github("cran/RandomFields")
#remotes::install_github("ropensci/NLMR")

#require(NLMR)
#Make_Grid(100,0.4,"homogenous")

Make_Grid<-function(area,inc,grid.opt){
  cells=(area/(inc^2))
  nrow.grid=area/inc
  if(grid.opt=="homogenous"){
  grid=matrix(nrow=round(cells),ncol=7)
  grid[,1]=1:cells
  grid[,2]=rep(seq(0,((inc*nrow.grid)-inc),inc),nrow.grid)
  grid[,3]=rep(seq(0,((inc*nrow.grid)-inc),inc),nrow.grid)
  grid[,4]=rep(seq(inc,(inc*nrow.grid),inc),nrow.grid)
  grid[,5]=rep(seq(inc,(inc*nrow.grid),inc),nrow.grid)
  grid[,6]=rep(seq(((0+inc)/2),(((inc*nrow.grid)-inc)+(inc*nrow.grid))/2,inc),nrow.grid)
  grid[,7]=rep(seq(((0+inc)/2),(((inc*nrow.grid)-inc)+(inc*nrow.grid))/2,inc),nrow.grid)
  centroids=grid[,c(6,7)]
  }

  
  if(grid.opt=="test_pref_binary"){
    grid=matrix(nrow=round(cells),ncol=8)
    grid[,1]=1:cells
    grid[,2]=rep(seq(0,((inc*nrow.grid)-inc),0.4),nrow.grid)
    grid[,3]=rep(seq(0,((inc*nrow.grid)-inc),0.4),nrow.grid)
    grid[,4]=rep(seq(inc,(inc*nrow.grid),0.4),nrow.grid)
    grid[,5]=rep(seq(inc,(inc*nrow.grid),0.4),nrow.grid)
    grid[,6]=rep(seq(((0+inc)/2),(((inc*nrow.grid)-inc)+(inc*nrow.grid))/2,0.4),nrow.grid)
    grid[,7]=rep(seq(((0+inc)/2),(((inc*nrow.grid)-inc)+(inc*nrow.grid))/2,0.4),nrow.grid)
    grid[,8]=rbinom(round(cells),1,0.2)
    centroids=grid[,c(6,7)]
  }
  
  grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
  return(grid.list)
  
}








