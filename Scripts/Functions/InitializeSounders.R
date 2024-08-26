#The purpose of this script is to initialize the sounder pop matrix
#This will replace the various state variable matrices (ie S, I)
#instead, values will be assigned to indicate statuses and location

#N0-population size, determined in InitializeASFModel using given density and area of grid
#ss-average sounder size, setting manually
#centroids-center coordinates of each cell
#type- if 0, initialize population for start of population; else, use init_locs and total number to initialize new births in population
#init_locs- used for births
#n-number of births

#optional inputs
#RSF mats...
#first col-cell class, double
#second col-RSF pref

#testing:
#RSF_mat=matrix(nrow=2,ncol=2)
#RSF_mat[,1]=c(0.0,1.0)
#RSF_mat[,2]=c(0.1,0.9)

InitializeSounders<-function(N0,ss,cells,centroids,type,init_locs,n,grid.opts){

#for initializing whole population
if(type==0){ 
sn_i<-N0/ss #Get the initializing number of sounders
  
  if(grid.opts=="homogenous"){
  assigns<-rbinom(cells,1,sn_i/cells) #randomly assign sounders to cells
  }
  
#testing
#if just do prob RSF * sn_i/cells, could end up with lower density than chose according to pop size
#want to choose overall pop size, then just have aggs different dep on landscape
#sum(prefs/cells) #averaged proportion of 1/0
#need like, X*0.2=(sn_i/cells
#X=(sn_i/cells)/0.2 as main probability
#then that probability is applied to each RSF prob per cell

#RSF_mat is matrix with col1=land class (ie, zero or one for a binary)
#col2=RSF preference (based on proportion of 0-1, with 1=highly preferred)
  if(grid.opts!="homogenous"){
    #if not homogenous, is preference of some kind
    #For now, rule will be higher number in cell=more preferred. 
    #For now is 1/0 for "test_pref_binary"
    #later want something flexible to enter in classes as input
    #prefs=grid[,8] #vector of lscape classes in grid
    #prefs[prefs==1]=RSF_mat[RSF_mat[,1]==1,2] #convert to RSF_prefs
    #prefs[prefs==0]=RSF_mat[RSF_mat[,1]==0,2] #convert to RSF_prefs
    pref.wt=sum(grid[,8])/cells #use this to weight preference so that still end up with N0 size population
    assigns=rbinom(cells,1,grid[,8]*((sn_i/cells)/pref.wt))
  }
  
sn<-sum(assigns) #generated number of rows (sounders), size from random assignment of sounders to cells
init_locs<-which(assigns==1) #get the locations where sounders have been initialized
  
#Initialize the sounder population matrix
#each row is a sounder
pop<-matrix(nrow=sn,ncol=13)
pop[,1]=rpois(sn,ss) #sounder size with avg as lambda in a poisson
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

if(grid.opts=="homogenous"){pop[,2]=0} #use this for pref.. just zero when homogenous
if(grid.opts!="homogenous"){pop[,2]=grid[pop[,3],8]} 

if(any(pop[,3]>nrow(centroids))){
  stop("agents initialized off the grid")
}

} else{ 
#for initializing initial infected individual introduction
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

if(any(pop[,3]>nrow(centroids))){
  stop("agents initialized off the grid")
}

return(pop)

}

