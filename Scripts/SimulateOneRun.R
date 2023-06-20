##The purpose of this script is to run a single rep of the ASF control optimization model

SimulateOneRun<-function(Pcr,Pir,Pbd,death,F1,F2,F2i,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop){

###########################################
######## Initialize Output Objects ######## 
###########################################

Nall=matrix(nrow=thyme) #track total abundance
BB=matrix(nrow=thyme) #track births

POSlive=as.list(rep(0,thyme)) #Positive cases observed and removed from landscape
POSdead=as.list(rep(0,thyme))#Positive carcasses observed and removed from landscape
NEGlive=as.list(rep(0,thyme)) #Negative tests of hunted carcasses that are removed from landscape
NEGdead=as.list(rep(0,thyme))#Negative tests of carcasses that are removed from landscape

POSlive_locs<-as.list(rep(0,thyme))
POSdead_locs<-as.list(rep(0,thyme))

idZONE=matrix(nrow=1,ncol=3) #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
#idZONE<-as.list(rep(as.integer(0),thyme))
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
Etrue=matrix(0,nrow=thyme,ncol=1)
IConDD=0
ICatEnd=0
#out[i,]<-areaOfinfection(pop,centroids,inc)
#print(out)

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

#initialize indices for weekly movement
indWM=seq(1,72,7)
indP=1:length(indWM)

#start the timestep loop
##i=1
##detectday=i
for(i in 1:thyme){
if (any(pop[,9,drop=FALSE]!=0|pop[,10,drop=FALSE]!=0|pop[,12,drop=FALSE]!=0)){
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

I_locs[[i]]<-pop[pop[,10]>0,3]
C_locs[[i]]<-pop[pop[,12]>0,3]
	
##########################
######## Movement ######## 
##########################
	
pop<-FastMovement(pop,centroids,shift,inc)

###############################
######## State Changes ######## 
###############################
#births, natural deaths, disease state changes (exposure, infection, recovery, death), carcass decay
st.list<-StateChanges(pop,centroids,cells,Pbd,B1,B2,F1,F2,Fi,K,death,Pcr,Pir,Incidence,BB,i)
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

#if it is at least day after detect day, and Rad>0
if(i > detectday & Rad > 0){
	
	#new detections from last step, bc day lag 
	#(either from initial detection or last culling period)
	#get locations in grid for detections
	idNEW=c(POSlive_locs[[i-1]],POSdead_locs[[i-1]])
	
	#remove NA/0 (may get NAs/zeroes if no live/dead detected)
	idNEW<-idNEW[idNEW>0&!is.na(idNEW)]

	#keep only new grid cells that weren't already identified in previous time steps
	if(any(idZONE[[i]])!=0){
	uniqueidNEW<-which(!(idNEW %in% idZONE[[i-1]][,1]))
	idNEW<-idNEW[uniqueidNEW]
	} else{idNEW=idNEW}

	#Culling process
	output.list<-CullingOneRun(pop,idNEW,idZONE[[i-1]],Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead)

	POSlive[[i]]<-output.list[[1]]
	POSdead[[i]]<-output.list[[2]]
	POSlive_locs[[i]]<-output.list[[3]]
	POSdead_locs[[i]]<-output.list[[4]]
	NEGlive[[i]]<-output.list[[5]]
	NEGdead[[i]]<-output.list[[6]]
	idZONE[[i]]<-output.list[[7]]
	removalcells[[i]]<-output.list[[8]]
	culled<-output.list[[9]]
	ZONEkm2[i,]<-output.list[[10]]

	pop<-output.list[[11]]
	#Total number culled at each timestep
	Tculled[i]=culled
} #if greater than detectday closing bracket

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################


#############################
####Track true spatial spread
#############################
#print("before aoi")
#print(pop[rowSums(is.na(pop)) > 0,])
#if any infected individuals
#areaOfinfection()
#print("before area of infection")
if(nrow(pop[pop[,9,drop=FALSE]>0|pop[,10,drop=FALSE]>0|pop[,12,drop=FALSE]>0,,drop=FALSE])>0){
#print("inside if statement")
out[i,]<-areaOfinfection(pop,centroids,inc)
#print("after area of infection output")
} else{out[i,]=c(0,0,0)}
#print("after area of infection")
#Plotting function in matlab goes here

#############################
####Generate Outputs
#############################

#Incidence is matrix with one column, nrow timesteps
#needs to be total exposures at each timestep
#first timestep one infected, incidence
#rest of the timesteps will be informed by Eep
#Incidence
#I_locs
#C_locs

#total exposures at each timestep
#this will just be inherent in new way of storing incidence
#IncOverTime = sum(Incidence,2);

#% find the last infectious case (true)
#sum all infectious cases (I,C) at each timestep
#ICtrue = sum(I + C,2); % sum of all infectious cases over time
ICtrue[i]<-sum(colSums(pop)[c(10,12)])

#same as above but with one removed at detection added back into total
#if Intensity > 0; ICtrue(detectday) = ICtrue(detectday) + 1; end % add back the one we removed at deteection
if(Intensity > 0){
ICtrue[detectday]=ICtrue[detectday]+1
	}

#% number of I, C, and E on detection day (add one to account for the one we subtract)
#IConDD = sum(I(detectday,:),2)+sum(C(detectday,:),2)+sum(E(detectday,:),2)+1; 
if(i==detectday){
IconDD=sum(colSums(pop)[c(9,10,12)]+1)
	}


#this just sums all of the exposures up until current time step
#Tinc = sum(IncOverTime); % true total new cases over all time
Tinc=sum(Incidence)

#this is sum of all exposures starting from day after detection day
#TincFromDD = sum(IncOverTime(detectday+1:length(IncOverTime))); % true total new cases from the time of detection
TincFromDD<-sum(Incidence[(detectday+1):length(Incidence)])
	
	
#this is sum of all exposures through detection day
#TincToDD = sum(IncOverTime(1:detectday));
TincToDD<-sum(Incidence[1:detectday])
	
#print("after TinctoDD")

#find last day there was an infectious individual
#idT = find(ICtrue > 0,1,'last'); % find the last day there was an infectious individual
idT=0
pl=1
#if last day there is infectious individual, just set that day 
#otherwise gets stuck in while loop
if(ICtrue[i]>0){idT=i} else {
while(pl>0){
idT=idT+1
pl=ICtrue[idT]	
}
idT=idT-1
}

#% how many detections?
#total number of detections made
#DET = sum(sum(POSlive+POSdead));
#print(POSlive[[i]])
#print(POSdead[[i]])
DET=sum(POSlive[[i]],POSdead[[i]])
#print("after DET")

#Spread = zeros(time,3); % number of infectious individuals, area of infection, max distance between any two cases
#%Mspread = [max(Spread(:,1)) max(Spread(:,2)) max(Spread(:,3))];
#want second column, area of infection
Mspread<-max(out[,2])

iCatEnd=sum(colSums(pop)[c(9,10,12)])
#print("after iCatEnd")

#% Reduced version as this is all we are tracking for now
#All = [Tinc, sum(Tculled), idT, max(Spread(:,2)), IConDD, ICatEnd, TincFromDD, TincToDD, DET];
list.all<-vector(mode="list",length=13)
list.all[[1]]=Tinc #sum of all exposures over simulation 
list.all[[2]]=Tculled #total number culled at each timestep 
list.all[[3]]=idT #last day there is an infectious individual
list.all[[4]]=Mspread #max spread of infection TROUBLESHOOT
list.all[[5]]=IConDD #number of I, C, and E on detection day TROUBLESHOOT
list.all[[6]]=ICatEnd #number of I,C,E on last day 
list.all[[7]]=TincFromDD #sum of all exposures starting day after detection day
list.all[[8]]=TincToDD #sum of all exposures up until detection day
list.all[[9]]=DET #total number of detections
list.all[[10]]=I_locs
list.all[[11]]=C_locs
list.all[[12]]=Isums
list.all[[13]]=Csums
#print("after create outputs")
} else{print("Exiting function, no infections")} #if any infected closing bracket/else
	} #for timestep closing bracket

######################################################################
#%find how many infections are present at the end
#% number of I, C, and E on the last day
#ICatEnd = sum(I(time,:),2)+sum(C(time,:),2)+sum(E(time,:),2);


return(list.all)
#return(pop)

	} #function closing bracket



