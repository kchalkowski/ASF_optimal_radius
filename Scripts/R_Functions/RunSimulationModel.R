#Runs each replicate of simulation
###%%%%%%%%%
#source(paste0(home,"/Scripts/SetParameters.R")) #uncomment this to set parms if not running in RunASFSimReplicates.R

#########################
####Choose State Rules
#########################
#define movement characteristics of the population (slow=FL, fast=SC)
RunSimulationModel<-function(rep){

######################
####Initialize Population
#####################
pop=InitializeSounders(centroids,grid,c(N0,ss),"init_pop",pop_init_grid_opts)
#FOI_cutoff=round(max(rgamma(2*nrow(pop)*thyme,shape=shift[1],scale=shift[2])),2)
out.list=SimulateOneRun(Pcr,Pir,Pbd,death,F1,F2_int,F2_B,F2i_int,F2i_B,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop,out.opts,"homogeneous",rep,DetP)
return(out.list)
}

