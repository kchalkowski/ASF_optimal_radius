FastMovement=function(pop,centroids,shift,inc,mv_pref){
  
  #get distances from gamma distribution
  pop[,4]=rgamma(nrow(pop),shape=shift[1],scale=shift[2])
  
  #set those less than inc to 0
  pop[pop[,4]<inc,][,4]=0 
  
  #set present locations to previous locations
  pop[,7]=pop[,3]
  
  #convert abundance/locs vector into long format
  abund.mat=matrix(0,nrow=40000,ncol=1)
  abund.df=data.frame("abund"=pop[,1],"cell"=pop[,3])
  abund.df=abund.df %>% dplyr::group_by(cell) %>% dplyr::summarize("abund"=sum(abund)) %>% as.data.frame()
  cells=data.frame("cell"=1:40000)
  abund.df=left_join(cells,abund.df,by="cell")
  abund.df$abund[is.na(abund.df$abund)]<-0
  abund.mat[,1]=abund.df$abund
  #pop[,4]=0
  #pop[1,4]=0.6347192
  m1=parallelMovementRcpp_portion(pop,abund.mat[,1,drop=FALSE],pop[,3,drop=FALSE],centroids,mv_pref)
  #abund.mat[1,1]
  
  
  pop[,3]=m1
  
  #if old locs same as new locs, print
  #if(all(pop[,3]==pop[,7])){
  #warning("all old locs same as new locs")
  #}
  
  #if spop
  #if stop function here.. if no cells to move to
  #any(pop[,3]==nrow(centroids)+1000)
  if(any(pop[,3]==nrow(centroids)+1000)) {
    stop("No cells to move to! This shouldn't happen")
  }
  
  #if stop function here.. 
  #if all sounders with dist equals zero NOT contained in rows for which prev locs=present locs
  if(all(!(which(pop[,4]==0)%in%which(pop[,3]==pop[,7])))){
    stop("All sounders with distance=0 should have same prev. and present locations")
  }
  
  return(pop)
}
