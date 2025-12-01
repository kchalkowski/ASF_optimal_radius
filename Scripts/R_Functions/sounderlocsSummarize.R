
#input: sounderlocs
#option 1- typical SEIRCZ summary per timestep (used in main script for diagnostics)
#option 2- zone-only SEIRCZ summary per timestep
#option 3- zone-only SEIRCZ summary per timestep, summarized on cell level
#options 1,2 and 3, with apparent prevalence output    
#specify above at 'out.opts' and output list of results with names, like in getoutputs functions

#Updates needed:
  #currently zone gets calc'd twice if ask for both zone and zonecells outputs, need to compress if statements so that only done once

sounderlocsSummarize<-function(sounderlocs,r,sl.summary.opts=NULL,DetP=NULL,N=NULL){
  if(missing(sl.summary.opts)){
    #default will be just basic SEIRCZ summary per timestep
    SEIRCZ.only=sounderlocs[,c('timestep','S','E','I','R','C','Z')] ## to be compatible with homogeneous landscapes, have to choose by name b/c they skip a column in centroids
#     SEIRCZ.only=sounderlocs[,4:10]
    SEIRCZ.rep=SEIRCZ.only %>% 
      dplyr::group_by(timestep) %>% 
      dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
    SEIRCZ.rep$rep=r
    SEIRCZ_total=SEIRCZ.rep[,c(8,1:7)] #want rep in front
    list.out=list("SEIRCZ_total"=SEIRCZ_total)
    
  } else{
    #else, do whatever outputs specified
    #need get totals in order to get apparent prev outs
    
    if("SEIRCZ_total"%in%sl.summary.opts|
       "SEIRCZ_total_apparent"%in%sl.summary.opts){
      SEIRCZ.only=sounderlocs[,3:9]
      SEIRCZ.rep=SEIRCZ.only %>% 
        dplyr::group_by(timestep) %>% 
        dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
      SEIRCZ.rep$rep=rep
      SEIRCZ_total=SEIRCZ.rep[,c(8,1:7)] #want rep in front
      
    }
    
    if("SEIRCZ_zone"%in%sl.summary.opts|
       "SEIRCZ_zone_apparent"%in%sl.summary.opts|
       "SEIRCZ_zone_cells"%in%sl.summary.opts|
       "SEIRCZ_zone_cells_apparent"%in%sl.summary.opts
       ){
      zoneoutlist=sounderlocsZone(sounderlocs,buffer,centroids,"I",FALSE)
      zone_out_I=zoneoutlist$ml
      areas=zoneoutlist$areas
      
    }
    
    
    if("SEIRCZ_zone"%in%sl.summary.opts|
       "SEIRCZ_zone_apparent"%in%sl.summary.opts){
    
      
      #then, subset to only sounders inZone
      sounderlocs_zone=zone_out_I[zone_out_I$inZone==1,]
      #then get prevalence summaries from that
      SEIRCZ.rep_zone=sounderlocs_zone %>% 
        dplyr::group_by(timestep) %>% 
        dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
      SEIRCZ.rep_zone$rep=rep
      SEIRCZ_zone=SEIRCZ.rep_zone[,c(8,1:7)] #want rep in front
      
    }
    #column 'cell' is not found
    if("SEIRCZ_zone_cells"%in%sl.summary.opts|
       "SEIRCZ_zone_cells_apparent"%in%sl.summary.opts){
      #speed test using sample sounderlocs:
      #slow version, uses R functions
      #{
      #tic()
      #zone_out_I=findZone_sim(sounderlocs,"I")
      #slowfunc=toc()
      #took 360 seconds, or 6 minutes
      #this would be per-replicate in simulation run
      #} 
      #fast version:
      #uses optimized functions with Rcpp
      zoneoutlist=sounderlocsZone(sounderlocs,buffer,centroids,"I",TRUE)
      zone_out_I_cells=zoneoutlist$ml
      #zone_out_I_cells=zone_out_I
      #tic()
    
      #fastfunc=toc()
      #took 11.57 seconds
      #229 times faster!!
      
      
      #subset to only ones in zone
      sounderlocs_zone_cells=zone_out_I_cells[zone_out_I_cells$inZone==1,]
      SEIRCZ_zone_cells=sounderlocs_zone_cells %>% 
        dplyr::group_by(timestep) %>% 
        dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
      
      
      if("SEIRCZ_zone_cells"%in%sl.summary.opts){
      #then summarize across cells, by timestep
      SEIRCZ.rep_zone_cells_summarized=
        SEIRCZ_zone_cells %>%
        dplyr::group_by(timestep) %>% 
        dplyr::summarise(Sme=median(S),
                         Eme=median(E),
                         Ime=median(I),
                         Rme=median(R),
                         Cme=median(C),
                         Zme=median(Z),
                         S25=quantile(S,0.25),
                         E25=quantile(E,0.25),
                         I25=quantile(I,0.25),
                         R25=quantile(R,0.25),
                         C25=quantile(C,0.25),
                         Z25=quantile(Z,0.25),
                         S75=quantile(S,0.75),
                         E75=quantile(E,0.75),
                         I75=quantile(I,0.75),
                         R75=quantile(R,0.75),
                         C75=quantile(C,0.75),
                         Z75=quantile(Z,0.75)) %>% as.data.frame()
      
      SEIRCZ.rep_zone_cells_summarized$rep=rep
      SEIRCZ_zone_cells=SEIRCZ.rep_zone_cells_summarized[,c(20,1:19)] #want rep in front
      }
    }
    
    if("SEIRCZ_total_apparent"%in%sl.summary.opts &
       !missing(DetP) &
       !missing(N)){
      #SEIRCZ.rep
      SEIRCZ_total_apparent=ExtApparentPrev(SEIRCZ.rep,DetP,N,rep)
      
    }
    
    if("SEIRCZ_zone_apparent"%in%sl.summary.opts &
       !missing(DetP) &
       !missing(N)){
      #SEIRCZ.rep_zone
      SEIRCZ_zone_apparent=ExtApparentPrev(SEIRCZ.rep_zone,DetP,N,rep)
      
    }
    
    #invalid times argument
    if("SEIRCZ_zone_cells_apparent"%in%sl.summary.opts &
       !missing(DetP) &
       !missing(N)){
      
      sounderlocs_zone_cells$rep=rep
      ###problem here
      SEIRCZ_zone_cells=sounderlocs_zone_cells[,c(12,1:11)]
      
      SEIRCZ_zone_cells_apparent=ExtApparentPrev(SEIRCZ_zone_cells,DetP,N,rep)
      
    }
    
    
    
    #return all objects in list
    list.out=vector(mode="list",length=length(sl.summary.opts))
    names(list.out)=sl.summary.opts
    
    for(i in 1:length(sl.summary.opts)){
    list.out[[i]]=eval(parse(text=sl.summary.opts[i]))
    }
      
  }
  
  return(list.out)
  
}


