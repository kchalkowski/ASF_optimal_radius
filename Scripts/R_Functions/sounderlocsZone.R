
#############################
########## Purpose ##########
#############################

#takes `sounderlocs` output from simulation run, finds chull around infected individuals for each timestep
#includes option to also get cell number (needed for cell-level summary options)

#locs-data.frame output from out.opts="sounderlocs" of simulation run (9 columns, nrow=every location with pig at every time step of one replicate)
#buffer-double, km, width of buffer around convex hull desired
#centroids-matrix of centroid x and y coordinates of grid used in simulation, rows=ncell x cols=x/y
#zonetype-option to draw zone around live or dead infected individuals
#findCells-boolean, option for outputting cell-level summary, with cell IDs

#Updates needed:
  #currently have bypass for chull function, throws error when xmin/ymin are same e.g.
  #Need to incorporate more elegant bypass into cpp func, currently just create false pt which not ideal

sounderlocsZone<-function(locs,buffer,centroids,zonetype,findCells){

  #Process
    #formatting data:
    #convert locs to matrix
    ml=as.matrix(locs)
    timesteps=unique(ml[,3])
    
    ml=cbind(ml,0)
    
    if(findCells){
    ml=cbind(ml,0)
    colnames(ml)[c(10,11)]<-c("inZone","cell")
    } else{
    colnames(ml)[10]<-c("inZone")
    }
    

    if(zonetype=="I"){
    ml.i=ml[ml[,6]>0,,drop=FALSE]
    } else{
      if(zonetype=="C"){
        ml.i=ml[ml[,8]>0,,drop=FALSE]
      } else{
        if(zonetype=="IC"){
          ml.i=ml[ml[,6]>0|ml[,8]>0,,drop=FALSE]
        }
        
      }
      
    }
    
  
    

    areas=matrix(nrow=length(timesteps),ncol=2)
    colnames(areas)=c("timestep","area")
    for(t in 1:length(timesteps)){
    #print(paste0("timestep ",t))
      ml.t=ml[ml[,3]==timesteps[t],,drop=FALSE]
      ml.t.i=ml.i[ml.i[,3]==timesteps[t],,drop=FALSE]
      
      #ml.t.i=ml.t[ml.t[,6]>0,,drop=FALSE]
      
      if(nrow(ml.t.i)>0){
      mlt.xy=ml.t.i[,c(1,2),drop=FALSE]
      #mlt.xy=cbind(mlt.xy,0,0)
      #print(2)
      if(t==1){
        mlt.xy_all=mlt.xy
      } else{
        mlt.xy_all=unique(rbind(mlt.xy_all,mlt.xy))
      }
      }
      
      #print(3)
      #are any of these same?
      {
      minx=which(mlt.xy_all[,1]==min(mlt.xy_all[,1]))[1]
      maxx=which(mlt.xy_all[,1]==max(mlt.xy_all[,1]))[1]
      miny=which(mlt.xy_all[,2]==min(mlt.xy_all[,2]))[1]
      maxy=which(mlt.xy_all[,2]==max(mlt.xy_all[,2]))[1]
      
      if(minx==miny|minx==maxy|maxx==miny|maxx==maxy){
        if(minx==miny){
          #make a fake point to be miny
          #add 0.4 to x, subtract 0.4 from y
          newpt=c(mlt.xy_all[miny,1]+0.4,mlt.xy_all[miny,2]-0.4)
          mlt.xy_all=rbind(mlt.xy_all,newpt)
        }
        
        if(minx==maxy){
          #make a fake point to be miny
          #add 0.4 to x, subtract 0.4 from y
          newpt=c(mlt.xy_all[maxy,1]+0.4,mlt.xy_all[maxy,2]+0.4)
          mlt.xy_all=rbind(mlt.xy_all,newpt)
        }
        
        if(maxx==miny){
          #make a fake point to be miny
          #add 0.4 to x, subtract 0.4 from y
          newpt=c(mlt.xy_all[miny,1]-0.4,mlt.xy_all[miny,2]-0.4)
          mlt.xy_all=rbind(mlt.xy_all,newpt)
        }
        
        if(maxx==maxy){
          #make a fake point to be miny
          #add 0.4 to x, subtract 0.4 from y
          newpt=c(mlt.xy_all[maxy,1]-0.4,mlt.xy_all[maxy,2]+0.4)
          mlt.xy_all=rbind(mlt.xy_all,newpt)
        }
        
      }
      }
      
      chull=FindChullcpp(mlt.xy_all)
      #print(4)
      
      chullsf=st_as_sf(as.data.frame(chull),coords=c(1,2))
      chullpoly=st_combine(chullsf) %>% st_cast("POLYGON")

      areas[i,1]=i
      areas[i,2]=st_area(chullpoly)
      
      if(nrow(mlt.xy_all)>1){
      zone=BufferChullcpp(mlt.xy_all,chull,buffer)
      } else{
        #else zone is just plus minus the point
        if(nrow(mlt.xy_all)==1){
          zone=matrix(nrow=4,ncol=2)
          zone[1,]=c(mlt.xy_all[1,1]-1,mlt.xy_all[1,2])
          zone[2,]=c(mlt.xy_all[1,1]+1,mlt.xy_all[1,2])
          zone[3,]=c(mlt.xy_all[1,1],mlt.xy_all[1,2]-1)
          zone[4,]=c(mlt.xy_all[1,1],mlt.xy_all[1,2]+1)
          
        }
      }
      
      #print(5)
      inZone=determineInsidePoints(zone,ml.t[,c(1,2)])
      #ml.t[,10]=inZone
      #print(6)
      #ml.t=cbind(ml.t,0)
      #colnames(ml.t)[11]<-"cellnum"
      #print(7)
    
      #cbind(ml[ml[,3]==timesteps[t],],inZone)
      
      ml[ml[,3]==timesteps[t],10] <- inZone

      if(findCells){
        cellnums=FindCellfromCentroid(ml.t[,c(1,2)],centroids)
        ml[ml[,3]==timesteps[t],11] <- cellnums
        
      }
      
    }
    
    out.list=list("ml"=as.data.frame(ml),"areas"=areas)
    return(out.list)
    
}


