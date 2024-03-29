##The purpose of this script is to run a single rep of the ASF control optimization model

SimulateOneRun<-function(Pcr,Pir,Pbd,death,F1,F2_int,F2_B,F2i_int,F2i_B,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop){

###########################################
######## Initialize Output Objects ######## 
###########################################

Nall=matrix(nrow=thyme) #track total abundance
BB=matrix(nrow=thyme) #track births

POSlive=as.list(rep(0,thyme)) #Positive cases observed and removed from landscape
POSdead=as.list(rep(0,thyme))#Positive carcasses observed and removed from landscape
NEGlive=as.list(rep(0,thyme)) #Negative tests of detected carcasses that are removed from landscape
NEGdead=as.list(rep(0,thyme)) #Negative tests of carcasses that are removed from landscape

POSlive_locs<-as.list(rep(0,thyme))
POSdead_locs<-as.list(rep(0,thyme))

idZONE=matrix(nrow=1,ncol=3) #optimzing #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
#idZONE<-as.list(rep(NA,thyme)) #original
Tculled=matrix(0,nrow=thyme) #total number culled at each time step
ZONEkm2=matrix(0,nrow=thyme) 
Carea=matrix(0,nrow=thyme) #area of culling zone at each time step
Spread=matrix(0,nrow=thyme, ncol=3) #number of infectious individuals, area of infection, max distance between any two cases
Incidence=matrix(0,nrow=thyme) #store new cases for each time step
I_locs=vector("list",thyme)
C_locs=vector("list",thyme)
removalcells=vector("list",thyme)
I_locs[1:thyme]<-0
C_locs[1:thyme]<-0
Isums<-matrix(0,nrow=thyme)
Csums<-matrix(0,nrow=thyme)
out=matrix(c(0,0,0),nrow=thyme,ncol=3)
ICtrue=matrix(0,nrow=thyme,ncol=1)

######################################
######## Initialize Infection ######## 
######################################

#num_inf_0=1 #how many pigs to infect starting off

#find the midpoint of the grid
id=which(centroids[,1]>=midpoint[1]&centroids[,2]>=midpoint[2])[1] #location on grid closest to midpoint

#infected<-InitializeSounders(N0,ss,cells,centroids,1,id,1)
infected<-InitializeSounders(N0,ss,cells,centroids,num_inf_0,id,1)
infected[,8]<-0
infected[,10]<-1

#combine infected pig with pop matrix
pop<-rbind(pop,infected)

#track first infection in Incidence matrix
Incidence[1]<-num_inf_0

##################################
######## Start simulation ######## 
##################################

#start the timestep loop
##i=1
##detectday=i
for(i in 1:thyme){
if (any(pop[,9,drop=FALSE]!=0|pop[,10,drop=FALSE]!=0|pop[,12,drop=FALSE]!=0)){

#for(i in 1:21){ #for manual troubleshooting of loop, in place of 1:thyme
#for(i in (detectday):thyme){ #for manual troubleshooting of loop, in place of 1:thyme
#for(i in 22:thyme){ #for manual troubleshooting of loop, in place of 1:thyme
#for(i in 22:24){

print(i)
#####################################
######## Track I/C locations ######## 
#####################################

if(nrow(pop[pop[,10]>0,,drop=FALSE])>0){
Isums[i]<-nrow(pop[pop[,10]>0,,drop=FALSE])
} else{Isums[i]=0}

if(nrow(pop[pop[,12]>0,,drop=FALSE])>0){
Csums[i]<-nrow(pop[pop[,12]>0,,drop=FALSE])
} else{Csums[i]=0}


if(any(pop[,10]>0)){		
I_locs[[i]]<-rep(pop[pop[,10]>0,3],pop[pop[,10]>0,10])
} else{
I_locs[[i]]<-pop[pop[,10]>0,3]
}

if(any(pop[,12]>0)){		
C_locs[[i]]<-rep(pop[pop[,12]>0,3],pop[pop[,12]>0,12])
} else{
C_locs[[i]]<-pop[pop[,12]>0,3]
}

##########################
######## Movement ######## 
##########################
	
pop<-FastMovement(pop,centroids,shift,inc)

###############################
######## State Changes ######## 
###############################
#births, natural deaths, disease state changes (exposure, infection, recovery, death), carcass decay
st.list<-StateChanges(pop,centroids,cells,Pbd,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B,K,death,Pcr,Pir,Incidence,BB,i)
pop<-st.list[[1]]
Incidence<-st.list[[2]]
BB<-st.list[[3]]



###**start on day of first detection
###################################
######## Initiate Response ######## 
###################################

#if it's detect day, and there are infected pigs to detect, and Rad>0
if(i==detectday&sum(pop[,c(9,10,12)])>0&Rad>0){
fd.list<-FirstDetect(pop,i,POSlive,POSdead,POSlive_locs,POSdead_locs)
pop=fd.list[[1]]
POSlive=fd.list[[2]]
POSdead=fd.list[[3]]
POSlive_locs=fd.list[[4]]
POSdead_locs=fd.list[[5]]
}

###**start day after day of first detection
#######################################
######## Initiate Culling Zone ######## 
#######################################

#if(i>1){
#print(paste0("POSlive ",i,": ",POSlive[[i-1]]))
#print(paste0("POSlive_locs ",i,": ",POSlive_locs[[i-1]]))
#print(paste0("POSdead ",i,": ",POSdead[[i-1]]))
#print(paste0("POSdead_locs ",i,": ",POSdead_locs[[i-1]]))
#}
#print(paste0("actual num. EIC ",i,": ",sum(colSums(pop)[c(9,10,12)])))

#if it is at least day after detect day, and Rad>0
if(i > detectday & Rad > 0){
	print("Entering Culling if")
	#new detections from last step, bc day lag 
	#(either from initial detection or last culling period)
	#get locations in grid for detections
	idNEW=c(unique(POSlive_locs[[i-1]]),unique(POSdead_locs[[i-1]])) #checked
	
	#remove NA/0 (may get NAs/zeroes if no live/dead detected)
	idNEW<-idNEW[idNEW>0&!is.na(idNEW)] 
	print(paste0("idNEW: ", idNEW))

	#get all unique grid cells of idZONE
	#idZONE_amal=unique(do.call(rbind,idZONE)) #original
	
	#idZONE_amal <- idZONE_amal[complete.cases(idZONE_amal),,drop=FALSE] #original
	
	
	#if there were detections in previous time steps, only get newly detected infected grid cells
	#"infected grid cell"=grid cell where there was an infected pig or carcass
	#if(nrow(idZONE_amal)>0){ #original
		if(!all(is.na(idZONE))){
	#remove any zero values
	#idZONE_amal=idZONE_amal[idZONE_amal[,1]>0,] #original
	idZONE=idZONE[idZONE[,1]>0,]

	#determine which cell ids are new, not already in zone from prev. timesteps
	uniqueidNEW<-which(!(idNEW %in% unique(idZONE[,1])))
	idNEW<-idNEW[uniqueidNEW]
	} else{idNEW=idNEW}
#idZONE_amal

	#Culling process
	#idZONE in S1R needs to contain all paired cells in zone from previous time zones
	#input to cullingonerun can be rbinded/unique version of this
	output.list<-CullingOneRun(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detectday,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead)

	POSlive[[i]]<-output.list[[1]]
	POSdead[[i]]<-output.list[[2]]
	POSlive_locs[[i]]<-output.list[[3]]
	POSdead_locs[[i]]<-output.list[[4]]
	NEGlive[[i]]<-output.list[[5]]
	NEGdead[[i]]<-output.list[[6]]
	#idZONE[[i]]<-output.list[[7]] #original
	idZONE<-output.list[[7]]
	removalcells[[i]]<-output.list[[8]]
	culled<-output.list[[9]]
	ZONEkm2[i,]<-output.list[[10]]

	pop<-output.list[[11]]
	#Total number culled at each timestep
	Tculled[i]=culled
} #if greater than detectday closing bracket



#############################
####Track true spatial spread
#############################
#if any infected individuals
if(nrow(pop[pop[,9,drop=FALSE]>0|pop[,10,drop=FALSE]>0|pop[,12,drop=FALSE]>0,,drop=FALSE])>0){
out[i,]<-areaOfinfection(pop,centroids,inc)
} else{out[i,]=c(0,0,0)}

#############################
####Summarize infections
#############################

#sum all infectious cases (I,C,E) at each timestep
#ICtrue = sum(I + C,2); sum of all infectious cases over time
if(i==1){ICtrue[i]=num_inf_0}
#print(ICtrue)
#print(length(which(ICtrue!=0)))
if(i==detectday){
ICtrue[i]<-(sum(colSums(pop)[c(9,10,12)])+1) #account for having removed that first detected
} else{
	ICtrue[i]<-sum(colSums(pop)[c(9,10,12)])
}

#} #for manual testing of loop

#comment brackets below for manual testing
} else{print("Exiting loop, no infections")} #if any infected closing bracket/else
	} #for timestep closing bracket

#############################
#############################

	
	
out.list<-GetOutputs(pop,Incidence,Tculled,ICtrue,out,detectday,I_locs,Clocs,POSlive_locs,POSdead_locs,Isums,Csums)

return(out.list)

	} #function closing bracket



