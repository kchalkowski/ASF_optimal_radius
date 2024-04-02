#The purpose of this script is to initialize the sounder pop matrix
#This will replace the various state variable matrices (ie S, I)
#instead, values will be assigned to indicate statuses and location

#N0-population size, determined in InitializeASFModel using given density and area of grid
#ss-average sounder size, setting manually
#centroids-center coordinates of each cell
#type- if 0, initialize population for start of population; else, use init_locs and total number to initialize new births in population
#init_locs- used for births
#n-number of births
InitializeSounders<-function(N0,ss,cells,centroids,type,init_locs,n){
if(type==0){
sn_i<-N0/ss #Get the initializing number of sounders
assigns<-rbinom(cells,1,sn_i/cells) #randomly assign sounders to cells
sn<-sum(assigns) #generated number of rows (sounders), size from random assignment of sounders to cells
init_locs<-which(assigns==1) #get the locations where sounders have been initialized

#Initialize the sounder population matrix
#each row is a sounder
pop<-matrix(nrow=sn,ncol=13)
pop[,1]=rpois(sn,ss) #sounder size with avg as lambda in a poisson
pop[,2]=0 #
pop[,3]=init_locs #this will be grid location (row number)
pop[,4]=0 #this will be assigned movement distance
pop[,5]=centroids[pop[,3],1] #present location X 
pop[,6]=centroids[pop[,3],2] #present location Y
pop[,7]=0 #previous location (grid row number)
pop[,8]=pop[,1] #number of S status in sounder
pop[,9]=0 #number of E status in sounder
pop[,10]=0 #number of I status in sounder
pop[,11]=0 #number of R status in sounder
pop[,12]=0 #number of C status in sounder
pop[,13]=0 #number of Z status in sounder

} else{ 
	
pop<-matrix(nrow=type,ncol=13)
pop[,1]=n #sounder size with avg as lambda in a poisson
pop[,2]=0 #this will be disease status: S=0; E=1; I=2; R=3; C=4; Z=5
pop[,3]=init_locs #this will be grid location (row number)
pop[,4]=0 #this will be assigned movement distance
pop[,5]=centroids[pop[,3],1] #present location X 
pop[,6]=centroids[pop[,3],2] #present location Y
pop[,7]=0 #previous location (grid row number)	
pop[,8]=n #number of S status in sounder
pop[,9]=0 #number of E status in sounder
pop[,10]=0 #number of I status in sounder
pop[,11]=0 #number of R status in sounder
pop[,12]=0 #number of C status in sounder
pop[,13]=0 #number of Z status in sounder
	}

pop=pop[pop[,1]>0,,drop=FALSE]

return(pop)

}

