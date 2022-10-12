#The purpose of this script is to initialize the sounder pop matrix
#This will replace the various state variable matrices (ie S, I)
#instead, values will be assigned to indicate statuses and location

InitializeSounders<-function(N0,ss,cells,centroids,type,init_locs,n){
if(type==0){
sn_i<-N0/ss #Get the initializing number of sounders
assigns<-rbinom(cells,1,sn_i/cells) #randomly assign sounders to cells
#sn<-sum(assigns) #generated num sounders from random assignment of sounders to cells
initsounds<-assigns*rpois(cells,ss)
init_locs<-which(assigns==1) #get the locations where sounders have been initialized
N0=sum(initsounds)

#Initialize the sounder population matrix
pop<-matrix(nrow=N0,ncol=7)
pop[,1]=1:N0 #used to be sounder size, could change this back if want to play with this
pop[,2]=0 #this will be disease status: S=0; E=1; I=2; R=3; C=4; Z=5
pop[,3]=0 #this will be present grid location (row number)
pop[,4]=0 #this will be assigned movement distance
pop[,5]=0#centroids[pop[,3],1] #present location X 
pop[,6]=0#centroids[pop[,3],2] #present location Y
pop[,7]=0 #previous location (grid row number)

wheresounders<-which(initsounds>0)
soundnums<-initsounds[wheresounders]

for(i in 1:length(wheresounders)){
#pop[,3]=init_locs #this will be present grid location (row number)	
#pop[1:4]<-wheresounders[1]	
#print(i)
if(i==1){
pop[(1):(soundnums[i]),3]<-wheresounders[i]
}
if(i==2){
pop[(1+soundnums[i-1]):(soundnums[i]+soundnums[i-1]),3]<-wheresounders[i]	
#pop[(villages$Dogs[i-1]+1):(villages$Dogs[i-1]+villages$Dogs[i]),18]<-villages$object_[i]

}
if(i>2){
pop[(sum(soundnums[1:(i-1)])+1):((sum(soundnums[1:(i-1)]))+soundnums[i]),3]<-wheresounders[i]
#pop[(sum(villages$Dogs[1:(i-1)])+1):(sum(villages$Dogs[1:(i-1)])+villages$Dogs[i]),18]<-villages$object_[i]
}
}

pop[,5]=centroids[pop[,3],1] #present location X 
pop[,6]=centroids[pop[,3],2] #present location Y

} else{ 
	
pop<-matrix(nrow=n,ncol=7)
pop[,1]=1:n #
pop[,2]=0 #this will be disease status: S=0; E=1; I=2; R=3; C=4; Z=5
pop[,3]=init_locs #this will be present grid location (row number)
pop[,4]=0 #this will be assigned movement distance
pop[,5]=centroids[pop[,3],1] #present location X 
pop[,6]=centroids[pop[,3],2] #present location Y
pop[,7]=0 #previous location (grid row number)	
	
	}
	
return(pop)

}
