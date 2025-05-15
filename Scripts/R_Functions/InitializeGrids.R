#Updates needed here:
  #way inputs are put into Make Grid is clunky... want similar format of object, grid.opt
  #add sprc and ras input options
  #change checks in MakeGrid to look for SpatRaster format raster than raster
    #for now just going to convert types back/forth. shouldnt be too bad bc just run need to run this once
  #add way to trigger multiple homogenous or heterogeneous grids
  #pass len and inc to grid.list


#object=c(parameters$len,parameters$inc)
#grid.opt=parameters$grid.opts
InitializeGrids<-function(object,grid.opt="homogeneous"){
	
    if(class(object)=="SpatRaster"){
      tar.grid.list=vector(mode="list",length=1)
      
      grid.list=Make_Grid(object,grid.opt)
      tar.grid.list[[1]]=grid.list
    }
  
    if(class(object)=="SpatRasterCollection"){
      tar.grid.list=vector(mode="list",length=length(object))
      
      for(g in 1:length(object)){
      ras=object[g]
      grid.list=Make_Grid(ras,grid.opt)
      tar.grid.list[[g]]=grid.list
      }
      
    }

    #vector of len, inc
    if(class(object)=="numeric"){
      if(length(object)==3){
      tar.grid.list=vector(mode="list",length=object[3])
      for(g in 1:object[3]){
      grid.list=Make_Grid(c(object[1],object[2]),grid.opt)
      tar.grid.list[[g]]=grid.list
      }
      } else{
        tar.grid.list=vector(mode="list",length=1)
        grid.list=Make_Grid(c(object[1],object[2]),grid.opt)
        tar.grid.list[[1]]=grid.list
      }
    }
    
    return(tar.grid.list)
  
}



