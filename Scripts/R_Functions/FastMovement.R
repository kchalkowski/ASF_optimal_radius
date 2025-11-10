FastMovement=function(pop,centroids,shape,rate,inc,mv_pref,RSF_mat=NULL,RSF_mat0=NULL){

  #run checks to make sure objects input correctly
  if(mv_pref!=3){
    if(!missing(RSF_mat0)|!missing(RSF_mat)){
      message("mv_pref not set to RSF-availability (3), but RSF_mat supplied. Ignoring RSF_mat/RSF_mat0 input.")
    }
  } else{
    if(mv_pref==3){
      if(missing(RSF_mat0)|missing(RSF_mat)){
      stop("mv_pref set to RSF-availability movement (3), but RSF matrices not supplied")
      }
    }
  }
  
  
  #get distances from gamma distribution
  pop[,4]=rgamma(nrow(pop),shape=shape,rate=rate)
  ## to fix the zombies (pig corpses with a movement distance; see note in last error catches below)
  pop[,4][pop[,1]==0] <- 0
  
  #set those less than inc to 0
#   pop[pop[,4]<inc,][,4]=0
  pop[,4][pop[,4] < inc] = 0 # need to do it this way to deal with the single row matrix case (drop=FALSE doesn't work for assigning values)

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

  #different functions for each movement pref version, temporary workaround.
  #mv_pref=3 in separate function because not sure how to use Nullify class input to rcpp parallel for optional args. 
  if(mv_pref!=3){
    m1=parallelMovementRcpp_portion(pop,abund.mat[,1,drop=FALSE],pop[,3,drop=FALSE],centroids,mv_pref)
  } else{
    if(mv_pref==3){
      m1=parallelMovement_RSFavail(pop,abund.mat[,1,drop=FALSE],pop[,3,drop=FALSE],centroids,RSF_mat0,RSF_mat,mv_pref)
    }
  }

  pop[,3]=m1
  
  
  if(mv_pref==2|mv_pref==3){
    #update lc vals of current new cell to pop
    pop[,2]=centroids[pop[,3],3]
  }
  
  #if old locs same as new locs, print
  #if(all(pop[,3]==pop[,7])){
  #warning("all old locs same as new locs")
  #}
  
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
