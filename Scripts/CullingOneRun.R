CullingOneRun<-function(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead){

  #profvis({
######################
###### Get Zone ######
######################
	
#if there were infected pigs detected in the last time step
if(length(idNEW)>0){

	#initiate vector to store pairs of each infected grid cell ID with every other cell ID
	#nrow for each combo, 3 cols
	#col 1-each grid cell ID with detected infection
	#col 2-each paired grid cell
	#col 3-distance between infected grid cell ID with detection and each grid cell
	pairedIDs<-matrix(nrow=cells*length(idNEW),ncol=3)

	#get matrix of each infected grid cell ID paired with every other cell ID
	if(length(idNEW)==1){
		pairedIDs[1:cells,1]<-idNEW
		pairedIDs[1:cells,2]<-1:cells
	} else{ if(length(idNEW)>1){
	
			for(j in 1:length(idNEW)){
				if(j==1){
					pairedIDs[1:cells,1]<-idNEW[j]	
					pairedIDs[1:cells,2]<-1:cells
				}
	
				if(j>1){
					cells=nrow(centroids)
					pairedIDs[((1)+((j-1)*cells)):((j*cells)),1]<-idNEW[j]
					pairedIDs[((1)+((j-1)*cells)):((j*cells)),2]<-1:cells
				}
			}
		}
	}
	
	#get distance between each infected grid cell ID paired with every other cell ID
	#below results in matrix with one column of each grid cell id and the second column the distance to each infected
	pairedIDs[,3]<-sqrt((centroids[pairedIDs[,1],1]-centroids[pairedIDs[,2],1])^2 + (centroids[pairedIDs[,1],2]-centroids[pairedIDs[,2],2])^2)

	#get all grid cells where the distance between that grid cell and an infected detected pig is less than the pre-determined radius
	idout<-pairedIDs[pairedIDs[,3]<=Rad,]

	} else{ #else, if no infections detected in last time step...no new grid cells added
		idout<-matrix(nrow=0,ncol=3)
	} 

#remove NAs
#idZONE_i <- idZONE[!is.na(idZONE)]

#fullZONE contains all grid cells with detected infections, from last time step and all before	
#if(nrow(idZONE)>1){
#idZONE_i=do.call(rbind,idZONE_i)	
fullZONE=rbind(idZONE,idout)
#} else{
#		if(length(idZONE_i)==1){
#			fullZONE=rbind(idZONE_i[[1]],idout)	
#			} else{
#				fullZONE=idout #this happens for first detection
#	}
#	}

#get all unique grid cells in the zone		
allINzone=unique(fullZONE[,2])

#get total number of grid cells in the zone	
Uall=length(allINzone)

#get total area of the zone	
ZONEkm2=Uall*inc^2
	
#get which rows in the zone	
#total number of sounders=soundINzone
soundINzone<-which(pop[,3]%in%allINzone) 

#get total number of pigs in zone
#(col 1 abundance of all live, plus num carcasses in 12, 13)
pigsinzone<-sum(pop[soundINzone,1],pop[soundINzone,12],pop[soundINzone,13])

#total number of infected pigs in zone
EICinzone<-sum(pop[soundINzone,9],pop[soundINzone,10],pop[soundINzone,12])

#get total number of pigs outside the zone
pigsoutzone<-sum(pop[-soundINzone,1],pop[-soundINzone,10],pop[-soundINzone,12])

#get total number of infected pigs outside the zone
EICoutzone<-sum(pop[-soundINzone,9],pop[-soundINzone,10],pop[-soundINzone,12]) 

#get total number of  individuals (inside and outside zone)	
totalpigs=sum(pop[,1],pop[,10],pop[,12])

#get total number of infected individuals (inside and outside zone)	
totalEIC=EICinzone+EICoutzone

#####################################
###### Begin Culling Algorithm ######
#####################################

#if there are pigs to cull... 
if(pigsinzone>0){
	
#get number of pigs for each grid cell in zone
#and get their status, SEIRCZ
#initiate empty matrix, nrow for each grid cell, 7 for each of SEIRCZ
#tic()
SEIRCZpigs<-matrix(0,nrow=nrow(fullZONE),ncol=7)
fullZONEpigs<-cbind(fullZONE,SEIRCZpigs)
popINzone<-pop[soundINzone,]
for(u in 1:nrow(popINzone)){
u_row<-which(fullZONEpigs[,2]==popINzone[u,3])
fullZONEpigs[u_row,4]<-popINzone[u,1] #total number of pigs
fullZONEpigs[u_row,5]<-popINzone[u,8] #total number susceptible pigs
fullZONEpigs[u_row,6]<-popINzone[u,9] #total number exposted pigs
fullZONEpigs[u_row,7]<-popINzone[u,10] #total number infected pigs
fullZONEpigs[u_row,8]<-popINzone[u,11] #total number recovered pigs
fullZONEpigs[u_row,9]<-popINzone[u,12] #total number infected carcasses
fullZONEpigs[u_row,10]<-popINzone[u,13] #total number uninfected carcasses
}
#toc()

#####speed test, alt
#tic()
#SEIRCZpigs<-matrix(0,nrow=nrow(fullZONE),ncol=7)
#fullZONEpigs<-cbind(fullZONE,SEIRCZpigs)
#popINzone<-pop[soundINzone,]
#  r_row<-which(fullZONEpigs[,2]%in%popINzone[u,3])
#  fullZONEpigs2=fullZONEpigs[r_row,]
#  for(u in 1:nrow(popINzone)){
#  u_row<-which(fullZONEpigs2[,2]==popINzone[u,3])
#  fullZONEpigs2[u_row,4]<-popINzone[u,1] #total number of pigs
#  fullZONEpigs2[u_row,5]<-popINzone[u,8] #total number susceptible pigs
#  fullZONEpigs2[u_row,6]<-popINzone[u,9] #total number exposted pigs
#  fullZONEpigs2[u_row,7]<-popINzone[u,10] #total number infected pigs
#  fullZONEpigs2[u_row,8]<-popINzone[u,11] #total number recovered pigs
#  fullZONEpigs2[u_row,9]<-popINzone[u,12] #total number infected carcasses
#  fullZONEpigs2[u_row,10]<-popINzone[u,13] #total number uninfected carcasses
#  }
#  toc()

#####

#remove rows from fullZONEpigs without pigs
fullZONEpigs<-fullZONEpigs[fullZONEpigs[,1]>0,,drop=FALSE]

#Cullstyle start in, start with closest pigs from detections
fullZONEpigs<-as.matrix(arrange(as.data.frame(fullZONEpigs),fullZONEpigs[,3]))

#%density of all live and dead pigs in the zone
Dr=pigsinzone/ZONEkm2

#determine density-dependent capture probability in this radius  	
cprob=1-(1/(1+alphaC)^Dr)
	
#get total number culled/removed/sampled in the zone
numb=rbinom(pigsinzone,1,cprob*Intensity)

#get cumulative sum of targeted pigs
#cpigs = cumsum(tpigs);
cpigs=sum(numb)

#determine how far down the list to remove pigs from cells
removals=0 #total number of removals, go through loop until first time it is equal to or greater than cpigs
incr=0 #row number where culling stops
while(removals<cpigs){
incr=incr+1
#if(i==62){
 # print(fullZONEpigs)
  #print(removals)
  #print(cpigs)
  #print(numb)
  #print(pigsinzone)
  #print(soundINzone)
  #print(pop)
  #}
removals<-removals+fullZONEpigs[incr,4]
	}

#determine which pigs culled
culled=removals[[1]]

#% list of cells that pigs will be eliminated from (column index was 1 in old version)		
removalcells<-fullZONEpigs[1:incr,2]

#get which pigs culled
removalpigs<-fullZONEpigs[1:incr,,drop=FALSE]

######################################
###### Update surveillance data ######
######################################


#POSlive_i is a matrix with a row for each timestep
#column one of poslive is the number of exposed/infected pigs detected at that timestep
#sum removalpigs column 6,7
POSlive_i<-sum(removalpigs[,7],removalpigs[,6])

#POSdead
#POSdead is a matrix with a row for each timestep
#column one of poslive is the number of infected carcasses detected at that timestep
#sum removalpigs column 9
POSdead_i<-sum(removalpigs[,9])

#POSlive_locs
#list of length thyme, each timestep is vector of grid cell locations where liive infected pigs detected at that ts
#removalpigs col 2 where column 6 or 7 >0 (need check that should be E Ior just I)
if(POSlive_i>0){
lll<-length(removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2])
#POSlive_locs_i=vector(mode="integer",length=lll)
#POSlive_locs[[i]]<-removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2]
POSlive_locs_i<-removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2]

#POSlive_locs_detect

#print(paste("number of removalpigs rows", length(removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2])))
#print(paste("POSlive_locs[i]",POSlive_locs[i]))
#print(paste("POSlive_locs[i]",POSlive_locs_i))
#print(paste("removalpigs locations",removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2]))
} else {POSlive_locs_i<-0}

#POSdead_locs
#list of length thyme, each timestep is vector of grid cell locations where dead infected pigs detected at that ts
#removalpigs col 2 where column 9>0
if(POSdead_i>0){
#POSdead_locs[[i]]<-removalpigs[removalpigs[,9]>0,2]
POSdead_locs_i<-removalpigs[removalpigs[,9]>0,2]
#lld<-length(removalpigs[removalpigs[,9]>0,2])
#POSdead_locs_i=vector(mode="integer",length=lld)
#POSdead_locs_i[[i]]<-removalpigs[removalpigs[,9]>0,2]
} else {POSdead_locs_i<-0}

#% Nlive(ids) = sum(X([1 2 4],ids),1); % count of S,E,R removed
#vector of nrow timestop, count of total SR removed
#sum removalpigs column 5,8
NEGlive_i<-sum(removalpigs[,5],removalpigs[,8])

#% Ndead(ids) = X(6,ids); % count of Z removed
#vector of nrow timestep, count of total Z removed
#sum removalpigs column 10
NEGdead_i<-sum(removalpigs[,10])

#idZONE:
#grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
idZONE=fullZONE 

#remove removed sounders from pop
removalrows<-which(pop[,3] %in% removalcells)
removedpop<-pop[-removalrows,,drop=FALSE]

#send updated objects to output list
output.list<-vector(mode="list",length=11)
output.list[[1]]<-POSlive_i
output.list[[2]]<-POSdead_i
output.list[[3]]<-POSlive_locs_i
output.list[[4]]<-POSdead_locs_i
output.list[[5]]<-NEGlive_i
output.list[[6]]<-NEGdead_i
output.list[[7]]<-idZONE
output.list[[8]]<-removalcells
output.list[[9]]<-culled
output.list[[10]]<-ZONEkm2
output.list[[11]]<-removedpop
} else{
output.list<-vector(mode="list",length=11)
output.list[[1]]<-0
output.list[[2]]<-0
output.list[[3]]<-0
output.list[[4]]<-0
output.list[[5]]<-0
output.list[[6]]<-0
output.list[[7]]<-idZONE
output.list[[8]]<-0
output.list[[9]]<-0
output.list[[10]]<-0
output.list[[11]]<-pop	
	
	}
#}) #profvis closer
  
return(output.list)

} #function closing bracket

