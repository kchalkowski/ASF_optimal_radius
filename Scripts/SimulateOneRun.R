##The purpose of this script is to run a single rep of the ASF control optimization model

SimulateOneRun<-function(Pcr,Pir,Pbd,death,F1,F2,F2i,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop){

#outputs:
#All,Nall,BB,IncOverTime	
	
#############################
####Initialize Output Objects
#############################
Nall=matrix(nrow=thyme) #track total abundance
BB=matrix(nrow=thyme) #track births

#############################
####Initialize Disease Surveillance/Culling Variables
#############################

POSlive=matrix(0, nrow=thyme,ncol=1) #Positive cases observed and removed from landscape
POSdead=matrix(0, nrow=thyme,ncol=1) #Positive carcasses observed and removed from landscape
NEGlive=matrix(0, nrow=thyme,ncol=1) #Negative tests of hunted carcasses that are removed from landscape
NEGdead=matrix(0, nrow=thyme,ncol=1) #Negative tests of carcasses that are removed from landscape

#empty_list <- vector(mode = "list", length = desired_length)
POSlive_locs<-vector(mode="list",length=thyme) #put this in beginning at initialization
POSdead_locs<-vector(mode="list",length=thyme) #put this in beginning at initialization


#############################
####Initialize Data Recording
#############################

idZONE=matrix(nrow=1,ncol=3) #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
Tculled=matrix(0,nrow=thyme) #total number culled at each time step
Carea=matrix(0,nrow=thyme) #area of culling zone at each time step
Spread=matrix(0,nrow=thyme, ncol=3) #number of infectious individuals, area of infection, max distance between any two cases
Incidence=matrix(0,nrow=thyme) #store new cases for each time step
I_locs=vector("list",thyme)
C_locs=vector("list",thyme)
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

#############################
####Initialize Infection
#############################

#find the midpoint of the grid
id=which(centroids[,1]>=midpoint[1]&centroids[,2]>=midpoint[2])[1] #location on grid closest to midpoint

infected<-InitializeSounders(N0,ss,cells,centroids,1,id,1)
infected[,8]<-0
infected[,10]<-1
#infected[,9]<-1
#combine infected pig with pop matrix
pop<-rbind(pop,infected)
Incidence[1]<-1

#############################
####Simulate pop over time
#############################

#initialize indices for weekly movement
indWM=seq(1,72,7)
indP=1:length(indWM)

#start the timestep loop
##i=1
##detectday=i
for(i in 1:thyme){
if (any(pop[,9]!=0|pop[,10]!=0|pop[,12]!=0)){
#print("top of loop")
#print(i)
#print(pop[rowSums(is.na(pop)) > 0,])
#print(pop[pop[,3]==0,])
#print(head(pop))	
#print(paste0("any infections?:",any(pop[,10]!=0|pop[,12]!=0)))


#############################
####Movement Function
#############################
#print(pop[,3])
#print("pre movement prints")
#print(unique(sapply(pop,class)))
#pop<-Movement(pop,centroids,shift,inc) #This function will be translated to RCPP
pop<-FastMovement(pop,centroids,shift,inc)
#print("after movement")
#print(pop[rowSums(is.na(pop)) > 0,])
#print(pop[pop[,3]==0,])
#print(unique(sapply(pop,class)))
#print(unique(sapply(pop, anyNA)))
#print(pop[anyNA(pop),])
#print("post movement prints end")
#print(pop[,7])
# watch daily movement changes, plotting
#check = find(temp>0); %ids of cells with pigs at previous timestep	
#scatter(centroids(check,1),centroids(check,2),'MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor',[0.9 0.9 0.9],'sizedata',4); hold on;
#set(gca,'fontsize',10,'fontname','times','xtick',[],'ytick',[]); %xlabel('km east-west'); ylabel('km north-south');
#text(2,gridlen-5,sprintf('Week = %d',i-1),'fontweight','bold','fontsize',12,'fontname','times');
#Icheck = find(I(i-1,:) > 0); scatter(centroids(Icheck,1),centroids(Icheck,2),'r','filled','sizedata',10);
#Ccheck = find(C(i-1,:) > 0); scatter(centroids(Ccheck,1),centroids(Ccheck,2),'b','filled','sizedata',15);
#axis([0 gridlen 0 gridlen]);

#may choose to plot at each timestep later (may save time/memory)
#but for now, just store movement records in a matrix
#only need to record the locations of the infected ones

#print(pop[pop[,10]>0,,drop=FALSE])
#print(pop[pop[,10]>0,,drop=FALSE][7])
#print(pop[pop[,10]>0,7,drop=FALSE])

#print(pop[pop[,12]>0,,drop=FALSE])
#print(pop[pop[,12]>0,][7])
#print(pop[pop[,12]>0,7,drop=FALSE])
#print(nrow(pop[pop[,10]>0,,drop=FALSE]))
if(nrow(pop[pop[,10]>0,,drop=FALSE])>0){
Isums[i]<-nrow(pop[pop[,10]>0,,drop=FALSE])
} else{Isums[i]=0}

if(nrow(pop[pop[,12]>0,,drop=FALSE])>0){
Csums[i]<-nrow(pop[pop[,12]>0,,drop=FALSE])
} else{Csums[i]=0}

I_locs[[i]]<-pop[pop[,10]>0,7]
C_locs[[i]]<-pop[pop[,12]>0,7]
#############################
####Births and Natural Deaths
#############################

#% Births and natural deaths
#%tic;
#    N = sum([S(i-1,:); E(i-1,:); R(i-1,:)],1);
#N=nrow(pop) 
########
idN=pop[pop[,8]>0|pop[,9]>0|pop[,10]>0|pop[,11]>0,] #get all sounder sets with live individuals; subset
liveind<-sum(colSums(pop)[8:11]) #N live individuals
liverows<-which(pop[,8]>0|pop[,9]>0|pop[,10]>0|pop[,11]>0) #rownums with live indiv

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
#Sdpb[id2[i],1]<-id[i]	
#which(pop[,8]>0|pop[,9]>0|pop[,10]>0|pop[,11]>0)
Sdpb[id2[j],1]<-id[j]
	}

#        Sdpb(id2) = a(1:id); % assign the births to cells where they will happen 
#for(i in 1:length(id2)){
#pop[id2[i],1]=pop[id2[i],1]+id[i] #add to total number of pigs in sounder
#pop[id2[i],8]=pop[id2[i],1]+id[i] #add to total number of susceptible pigs in sounder
#	}
#    end



#############################
####State Transitions
#############################
##########Notes: currently doing state transitions
#need Sdpd, Edp etc to be lists of 1/0/2/whatever but created such that they will only go in columns
#where there are relevant individuals (i.e. cant have recoveries where no infected)

#Natural Deaths
#    Sdpd = binornd(S(i-1,:),death); % natural death of S (a constant)
#    Edp = binornd(E(i-1,:),death); % natural death of E (a constant)
#    Rdp = binornd(R(i-1,:),death); % natural death of R (a constant)
#Sdpd<-rbinom(nrow(pop),1,death)
#Edpd<-rbinom(nrow(pop),1,death)
#Idpd<-rbinom(nrow(pop),1,death)
#Rdpd<-rbinom(nrow(pop),1,death)
#pop[liverows,]
#
Sdpd<-matrix(nrow=nrow(pop),ncol=1)
Edpd<-matrix(nrow=nrow(pop),ncol=1)
Idpd<-matrix(nrow=nrow(pop),ncol=1)
Rdpd<-matrix(nrow=nrow(pop),ncol=1)
Sdpd[,1]=0
Edpd[,1]=0
Idpd[,1]=0
Rdpd[,1]=0
#Force of Infection
#Pse = FOI(I(i-1,:),C(i-1,:),B1,B2,F1,F2,F2i,centroids,cells,S(i-1,:));
#Pei = 1-exp(-1./(poissrnd(4,1,cells)./7)); %weekly scale %1-exp(-1./max(1,poissrnd(4,1,cells))); % incubation period (fix as mean of truncated Poisson)
#Pic = 1-exp(-1./(poissrnd(5,1,cells)./7)); %weekly scale % infectious period until death (fix as mean of truncated Poisson)
#QUESTION: say in notes truncated poisson. is this old? poissrnd is regular poisson right?
#print("which is na pop2")
#print(pop[which(is.na(pop)),])

#int("any infectious individuals?:")
#print(length(pop[pop[,10]>0|pop[,12]>0,1])>0)
Pse<-FOI(pop,centroids,cells,B1,B2,F1,F2,Fi)
#print("after FOI")
#print(pop[rowSums(is.na(pop)) > 0,])
#print(which(Pse==1))
#print(pop[which(Pse==1),])
Pei=1-exp(-1/rpois(nrow(pop),4)/7)
Pic=1-exp(-1/rpois(nrow(pop),5)/7)
#print(head(Pse))
#Disease transitions
#    Eep = binornd(S(i-1,:),Pse); % Exposure
#    Iep = binornd(E(i-1,:),Pei); % transitions from E to I
#    Rep = binornd(I(i-1,:),Pir.*Pic); % transitions from I to R, natural deaths of R
#    Cep = binornd(I(i-1,:),(1-Pir).*Pic); % transitions from I to C
Eep=matrix(nrow=nrow(pop),ncol=1)
Eep[,1]=0
#Eep<-rbinom(nrow(pop),1,Pse[pop[,3],]) #Exposure (S -> E)
Iep=matrix(nrow=nrow(pop),ncol=1)
Iep[,1]=0
#Iep<-rbinom(nrow(pop),1,Pei) #Disease progression (E -> I)
Rep=matrix(nrow=nrow(pop),ncol=1)
Rep[,1]=0
#Rep<-rbinom(nrow(pop),1,Pir*Pic) #Recovery (I -> R)
Cep=matrix(nrow=nrow(pop),ncol=1)
Cep[,1]=0
#Cep<-rbinom(nrow(pop),1,(1-Pir)*(Pic)) #Death (I -> C)

#Carcass decay
Ccd=matrix(nrow=nrow(pop),ncol=1)
Ccd[,1]=0
#Ccd<-rbinom(nrow(pop),1,Pcr) #prob of removal from landscape
Zcd=matrix(nrow=nrow(pop),ncol=1)
Zcd[,1]=0
#Zcd<-rbinom(nrow(pop),1,Pcr) #prob of removal from landscape
#print(head(pop))
for(k in 1:nrow(pop)){
#print(i)
if(pop[k,8]>0){
#operations on Susceptible individuals
#print(pop[k,3])
#print(pop[k,8])
Sdpd[k]<-sum(rbinom(pop[k,8],1,death))
Eep[k]<-sum(rbinom(pop[k,8],1,Pse[pop[k,3]])) #Exposure (S -> E)
#print("popk")
#print(pop[k,8])
#print(Pse[pop[k,3]])
#print(Eep[k])
}	
	#print(i)
if(pop[k,9]>0){
#operations on Exposed individuals
Edpd[k]<-sum(rbinom(pop[k,9],1,death))
Iep[k]<-sum(rbinom(pop[k,9],1,Pei))
}
	#print(i)
if(pop[k,10]>0){
#operations on Infected individuals	
Idpd[k]<-sum(rbinom(pop[k,10],1,death))
Rep[k]<-sum(rbinom(pop[k,10],1,Pir*Pic))
Cep[k]<-sum(rbinom(pop[k,10],1,(1-Pir)*(Pic))) 
}	
	#print(i)
if(pop[k,11]>0){
#operations on Recovered individuals
Rdpd[k]<-sum(rbinom(pop[k,11],1,death))
}	
	#print(i)
if(pop[k,12]>0){
#operations on Inf carcass individuals	
Ccd<-sum(rbinom(pop[k,12],1,Pcr))	
}	
	#print(i)
if(pop[k,13]>0){
#operations on Uninf carcass individuals	
Zcd<-sum(rbinom(pop[k,13],1,Pcr))	
	}	
}

Incidence[i]<-Incidence[i]+sum(Eep)
#print("after assignments")
#print(pop[rowSums(is.na(pop)) > 0,])
#############################
####Update States based on Demographic and Epidemiological Processes
#############################
#print(paste0("new exposures:",any(Eep>0)))

pop[,8]=pop[,8]-Eep+Sdpb-Sdpd #S
pop[,9]=pop[,9]-Iep+Eep-Edpd #E
pop[,10]=pop[,10]-Rep-Cep+Iep#I
pop[,11]=pop[,11]+Rep-Rdpd #R
pop[,12]=pop[,12]+Cep-Ccd #C
pop[,13]=pop[,13]+Rdpd+Sdpd+Edpd-Zcd #Z

pop[which(pop[,8]<0),8]<-0
pop[which(pop[,9]<0),9]<-0
pop[which(pop[,10]<0),10]<-0
pop[which(pop[,11]<0),11]<-0
pop[which(pop[,12]<0),12]<-0
pop[which(pop[,13]<0),13]<-0

#Make Nall output here

#############################
####Response: Culling Zone
#############################

#% RESPONSE: CULLING ZONE (ALL CULLED PIGS ARE TESTED)
#    if i > detectday && Rad > 0 % now that the initial detection is made, start testing all culled individuals; only do if there are pigs in the zone to cull
if(i > detectday & Rad > 0){
##    idNEW = find(POSlive(i-1,:) + POSdead(i-1,:) > 0); % find cells that had positives in the previous time step from detectday(assume 1-day lag between sampling and test results)
#idNEW=c(POSlive[i-1,2],POSdead[i-1,2])
idNEW=c(POSlive_locs[[i-1]],POSdead_locs[[i-1]])
#POSlivei
idNEW<-idNEW[idNEW>0&!is.na(idNEW)]
#% determine which are newly identified grid cells
#    [~,temp] = intersect(idNEW, unique(idZONE(:,1))); % gives IDs in idNEW that overlap with pre-existing positive IDs
#    idNEW(temp) = []; % remove the id's that overlap so we are just looking at the unique ids
#idZONE=matrix(nrow=1,ncol=3) #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
uniqueidNEW<-which(!(idNEW %in% idZONE[,1]))
idNEW<-idNEW[uniqueidNEW]
#    X = [S(i,:);E(i,:);I(i,:);R(i,:);C(i,:);Z(i,:)]; % get current status 
#    [idZONE,ids,Plive,Pdead,culled,areaC] = CullingOneRun(X,idNEW,idZONE,Intensity,alphaC,centroids,Rad,cullstyle,inc,i,stplot,enplot); % update numbers that were culled and tested, 

#CullingOneRun(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detected)
output.list<-CullingOneRun(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead)
#    if isempty(culled) == 0; Tculled(i) = culled; else; Tculled(i) = 0; end
#    if isempty(areaC) == 0; Carea(i) = areaC; else; Carea(i) = 0; end
#print("exited CullingOneRun")

#############################
####Update surveillance from culling zone
#############################

POSlive<-output.list[[1]]
POSdead<-output.list[[2]]
POSlive_locs<-output.list[[3]]
POSdead_locs<-output.list[[4]]
NEGlive<-output.list[[5]]
NEGdead<-output.list[[6]]
idZONE<-output.list[[7]]
removalcells<-output.list[[8]]
culled<-output.list[[9]]
ZONEkm2<-output.list[[10]]

#print("after update surveillance")
#print(pop[rowSums(is.na(pop)) > 0,])
#############################
####Update States based on Management Processes
#############################
#print(paste("nrow pop before removals:",nrow(pop)))
pop<-output.list[[11]]
#print(paste("nrow pop after removals:",nrow(pop)))


#Total number culled at each timestep
Tculled[i]=culled

} #if greater than detectday closing bracket

    

#############################
####Initiate Response based on day of first detection
#############################

#% INITIATE RESPONSE ON THE DAY OF FIRST DETECTION (PRE-DETERMINED)
#    if i == detectday && sum(I(i,:)+C(i,:)+E(i,:)) > 0 && Rad > 0 % day the the first detection is made

#detection is the row of the pop of the infected pig that was detected
#POSlive is a matrix with a row for each timestep
#column one of poslive is the number of infected pigs detected at that timestep
#remov this column? column two of poslive is the grid cell ID of the infected pigs
if((i==detectday)&(sum(pop[,9]+pop[,10]+pop[,12])>0)){
detection<-as.integer(sample(as.character(which(pop[,9]>0|pop[,10]>0|pop[,12]>0)),1))
POSlive[i,1]<-min(pop[detection,9]+pop[detection,10],1)
POSdead[i,1]<-min(pop[detection,12],0)
#update the surveillance data
#instead of second column of POSlive/POSdead, create a list with vector for each timestep

if(POSlive[i,1]>0){POSlive_locs[[i]]<-pop[detection,3]}
if(POSdead[i,1]>0){POSdead_locs[[i]]<-pop[detection,3]}

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

} 
#############################
####Track true spatial spread
#############################
#print("before aoi")
#print(pop[rowSums(is.na(pop)) > 0,])
#if any infected individuals
#areaOfinfection()
if(nrow(pop[pop[,9]>0|pop[,10]>0|pop[,12]>0,,drop=FALSE])>0){
out[i,]<-areaOfinfection(pop,centroids,inc)
} else{out[i,]=c(0,0,0)}

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
DET=sum(POSlive,POSdead)

#Spread = zeros(time,3); % number of infectious individuals, area of infection, max distance between any two cases
#%Mspread = [max(Spread(:,1)) max(Spread(:,2)) max(Spread(:,3))];
#want second column, area of infection
Mspread<-max(out[,2])

iCatEnd=sum(colSums(pop)[c(9,10,12)])

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

} else{print("Exiting function, no infections")} #if any infected closing bracket/else
	} #for timestep closing bracket

######################################################################
#%find how many infections are present at the end
#% number of I, C, and E on the last day
#ICatEnd = sum(I(time,:),2)+sum(C(time,:),2)+sum(E(time,:),2);


#return(list.all)
return(pop)

	} #function closing bracket



