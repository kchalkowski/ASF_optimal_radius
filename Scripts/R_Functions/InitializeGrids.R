#Updates needed here:
  #way inputs are put into Make Grid is clunky... want similar format of object, grid.opt
  #add sprc and ras input options
  #change checks in MakeGrid to look for SpatRaster format raster than raster
    #for now just going to convert types back/forth. shouldnt be too bad bc just run need to run this once
  #add way to trigger multiple homogenous or heterogeneous grids
  #pass len and inc to grid.list


#object=c(parameters$len,parameters$inc)
#grid.opt=parameters$grid.opts
InitializeGrids<-function(object,grid.opt){
  if(missing(grid.opt)){
    if(class(object)=="numeric"){grid.opt="homogenous"}
    #add missing switch for spatraster type
  }
    
    #vector of len, inc
    if(class(object)=="numeric"){
      len=object[1]
      inc=object[2]
      
      grid.list=Make_Grid(len,inc,grid.opt)

    }
    
    return(grid.list)
  
}



