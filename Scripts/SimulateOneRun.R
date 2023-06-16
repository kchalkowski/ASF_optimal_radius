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

#idZONE=matrix(nrow=1,ncol=3) #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
idZONE<-as.list(rep(as.integer(0),thyme))
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

#######################################
######## Births/Natural Deaths ######## 
#######################################

#determine disease state change probabilities
Pse<-FOI(pop,centroids,cells,B1,B2,F1,F2,Fi) #force of infection
Pei=1-exp(-1/rpois(nrow(pop),4)/7) #transitions exposure to infected
Pic=1-exp(-1/rpois(nrow(pop),5)/7) #transitions infected to either dead or recovered

#conduct the state changes
#pop<-StateChanges(pop,Pse,Pei,)

######
idN=pop[pop[,8,drop=FALSE]>0|pop[,9,drop=FALSE]>0|pop[,10,drop=FALSE]>0|pop[,11,drop=FALSE]>0,] #get all sounder sets with live individuals; subset
liveind<-sum(colSums(pop)[8:11]) #N live individuals
liverows<-which(pop[,8,drop=FALSE]>0|pop[,9,drop=FALSE]>0|pop[,10,drop=FALSE]>0|pop[,11,drop=FALSE]>0) #rownums with live indiv

#        if id <= length(idN)
########
#    liveind = sum(N);
#    Brate = Pbd*liveind*(1-liveind/K); %density-dependent birth rate; seasonally varying
Brate=Pbd*liveind*(1-liveind/K)
#    Tbirths = poissrnd(Brate); % Total births
Tbirths=rpois(1,Brate)
#    BB(i) = Tbirths; % record total births this time step
BB[i]=Tbirths
#    a = round(randi([0,min(Tbirths,10)], 1, cells)); b = cumsum(a); % Assign births to cells with a max of 10 per cell
#a=round(runif(cells,min=0, max=min(Tbirths,10)))
#a=round(runif(nrow(pop),min=0, max=min(Tbirths,10)))
#    id = find(b >= Tbirths,1,'first'); % Pick out enough numbers that sum to Tbirths - this will determine how many cells get births
id<-0
n=1
while (sum(id) < Tbirths) {
    birthset_i <- round(runif(1,min=0, max=min(Tbirths,10)))
    id[n]<-birthset_i
    n=n+1
}

#    Sdpb = zeros(1,cells);
#Sdpb=matrix(nrow=cells,ncol=1)

#    if isempty(id) < 1
if(length(id)>1){
#        idN = find(N > 0);

if(length(id)<=nrow(idN)){
#            id2 = randsample(idN,id); % pick which cells with pigs will get the births
id2=sample(1:nrow(idN),length(id))
#        else
} else {
	#            id2 = randsample(idN,id,true); % if there are more births than cells only add births cells where the pigs are (so fewer births will be happening)
id2=sample(1:length(liverows),length(liverows))
#        end
}
}


Sdpb=matrix(nrow=nrow(pop),ncol=1)
Sdpb[,1]=0

for(j in 1:length(id)){
Sdpb[id2[j],1]<-id[j]
	}

###################################
######## State Transitions ######## 
###################################

#natural deaths
Sdpd<-matrix(nrow=nrow(pop),ncol=1)
Edpd<-matrix(nrow=nrow(pop),ncol=1)
Idpd<-matrix(nrow=nrow(pop),ncol=1)
Rdpd<-matrix(nrow=nrow(pop),ncol=1)
Sdpd[,1]=0
Edpd[,1]=0
Idpd[,1]=0
Rdpd[,1]=0

#disease state change recording
Eep=matrix(nrow=nrow(pop),ncol=1)
Eep[,1]=0
Iep=matrix(nrow=nrow(pop),ncol=1)
Iep[,1]=0
Rep=matrix(nrow=nrow(pop),ncol=1)
Rep[,1]=0
Cep=matrix(nrow=nrow(pop),ncol=1)
Cep[,1]=0

#Carcass decay recording
Ccd=matrix(nrow=nrow(pop),ncol=1)
Ccd[,1]=0
Zcd=matrix(nrow=nrow(pop),ncol=1)
Zcd[,1]=0

#Conduct the state changes
#############################

for(k in 1:nrow(pop)){

#operations on Susceptible individuals
if(pop[k,8]>0){
#print(pop[k,3])
#print(pop[k,8])
Sdpd[k]<-sum(rbinom(pop[k,8],1,death))
Eep[k]<-sum(rbinom(pop[k,8],1,Pse[pop[k,3]])) #Exposure (S -> E) infection based on probability using their location
#print("popk")
#print(pop[k,8])
#print(Pse[pop[k,3]])
#print(Eep[k])
}	

#operations on Exposed individuals
if(pop[k,9]>0){
Edpd[k]<-sum(rbinom(pop[k,9],1,death))
Iep[k]<-sum(rbinom(pop[k,9],1,Pei))
}

#operations on Infected individuals	
if(pop[k,10]>0){
Idpd[k]<-sum(rbinom(pop[k,10],1,death))
Rep[k]<-sum(rbinom(pop[k,10],1,Pir*Pic))
Cep[k]<-sum(rbinom(pop[k,10],1,(1-Pir)*(Pic))) 
}	

#operations on Recovered individuals
if(pop[k,11]>0){
Rdpd[k]<-sum(rbinom(pop[k,11],1,death))
}	

#operations on Carcasses (infected)
if(pop[k,12]>0){
Ccd<-sum(rbinom(pop[k,12],1,Pcr))	
}	

#operations on Carcasses (uninfected)
if(pop[k,13]>0){
#operations on Uninf carcass individuals	
Zcd<-sum(rbinom(pop[k,13],1,Pcr))	
	}	
}

Incidence[i]<-Incidence[i]+sum(Eep)

#update states in pop matrix
###################################
pop[,8]=pop[,8]-Eep+Sdpb-Sdpd #S
pop[,9]=pop[,9]-Iep+Eep-Edpd #E
pop[,10]=pop[,10]-Rep-Cep+Iep-Idpd#I
pop[,11]=pop[,11]+Rep-Rdpd #R
pop[,12]=pop[,12]+Idpd+Cep-Ccd #C
pop[,13]=pop[,13]+Sdpd+Rdpd+Edpd-Zcd #Z

#sometimes end up with negative numbers 
#(i.e. all pigs in sounders chosen for natural mort and disease mort)
#just set anything below zero to zero
pop[which(pop[,8]<0),8]<-0
pop[which(pop[,9]<0),9]<-0
pop[which(pop[,10]<0),10]<-0
pop[which(pop[,11]<0),11]<-0
pop[which(pop[,12]<0),12]<-0
pop[which(pop[,13]<0),13]<-0

#move dead individuals (C or Z) into their own rows
#pop[,12] and pop[,13] > 0
deadguys<-pop[pop[,12]>0|pop[,13]>0,,drop=FALSE]

#if there are deadguys....
if(length(deadguys)!=0){
#remove abundance and all live guy counts from deadguy set
deadguys[,1]=0
deadguys[,8]=0
deadguys[,9]=0
deadguys[,10]=0
deadguys[,11]=0

#set all deadguys in other pop to zero
pop[which(pop[,12]>0),12]<-0
pop[which(pop[,13]>0),13]<-0
pop<-rbind(pop,deadguys)

}

#Update abundance numbers (live individuals only count in abundance)
pop[,1]=rowSums(pop[,8:11])

########This block begins at detection day
############################################################################
############################################################################
############################################################################

#############################
####Initiate Response based on day of first detection
#############################
#detection is the row of the pop of the infected pig that was detected
#POSlive is a matrix with a row for each timestep
#column one of poslive is the number of infected pigs detected at that timestep
#remov this column? column two of poslive is the grid cell ID of the infected pigs
######################################################################################################
if((i==detectday)&(sum(pop[,9]+pop[,10]+pop[,12])>0)){
  detection<-as.integer(sample(as.character(which(pop[,9]>0|pop[,10]>0|pop[,12]>0)),1))
  POSlive[[i]]<-min(pop[detection,9]+pop[detection,10],1)
  POSdead[[i]]<-min(pop[detection,12],0)
  
  
  #update the surveillance data
  if(POSlive[[i]]>0){POSlive_locs[[i]]<-pop[detection,3]}
  if(POSdead[[i]]>0){POSdead_locs[[i]]<-pop[detection,3]}
  
  #Store the pop rows of the detected pigs
  detected<-pop[detection,,drop=FALSE]
  
  #Remove the detected case
  #if an infected live or infected carcass is discovered and removed,
  #remove from disease status column and remove one from general number of pigs in cell column
  if(pop[detection,9]>0|pop[detection,10]>0|pop[detection,12]>0){
    pop[detection,1]<-pop[detection,1]-1
    pop[detection,9]<-max(pop[detection,9]-1,0)
    pop[detection,10]<-max(pop[detection,10]-1,0)
    pop[detection,11]<-max(pop[detection,11]-1,0)
  }
  #if an infected live AND infected carcass is discovered and removed,
  #remove from disease status column and remove two from general number of pigs in cell column
  if(pop[detection,9]>0|pop[detection,10]>0&pop[detection,12]>0){
    pop[detection,1]<-pop[detection,1]-2
    pop[detection,9]<-max(pop[detection,9]-1,0)
    pop[detection,10]<-max(pop[detection,10]-1,0)
    pop[detection,11]<-max(pop[detection,11]-1,0)
  }
  
  #if every pig in cell with infected pigs removed, remove the row from the population
  if(pop[detection,1]==0){pop<-pop[-detection,]}
  #print("after initiate response")
} 

#This block after detection step
############################################################################
############################################################################
############################################################################

#############################
####Response: Culling Zone
#############################

#% RESPONSE: CULLING ZONE (ALL CULLED PIGS ARE TESTED)
if(i > detectday & Rad > 0){
idNEW=c(POSlive_locs[[i-1]],POSdead_locs[[i-1]])
idNEW<-idNEW[idNEW>0&!is.na(idNEW)]

if(any(idZONE[[i]])!=0){uniqueidNEW<-which(!(idNEW %in% idZONE[[i-1]][,1]))
idNEW<-idNEW[uniqueidNEW]
} else{idNEW=idNEW}

output.list<-CullingOneRun(pop,idNEW,idZONE[[i-1]],Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead)

#############################
####Update surveillance from culling zone
#############################

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



