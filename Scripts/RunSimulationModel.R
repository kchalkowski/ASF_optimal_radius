#Runs each replicate of simulation
###%%%%%%%%%
#source(paste0(home,"/Scripts/SetParameters.R")) #uncomment this to set parms if not running in RunASFSimReplicates.R

#########################
####Choose State Rules
#########################
#define movement characteristics of the population (slow=FL, fast=SC)

RunSimulationModel<-function(){

######################
####Initialize Population
#####################
pop<-InitializeSounders(N0,ss,cells,centroids,0,0,0)

######################
####RunModel
#####################

out.list=SimulateOneRun(Pcr,Pir,Pbd,death,F1,F2_int,F2_B,F2i_int,F2i_B,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop)
return(out.list)
}

