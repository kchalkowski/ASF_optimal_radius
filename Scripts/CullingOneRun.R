#CullingOneRun<-function(X,idNEW,idZONE,Intensity,alphaC,CT,Rad,inc,i){
CullingOneRun<-function(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead){
#print("inside cullingonerun")
#print(length(idNEW))
if(length(idNEW)>0){
#idout = []; ids = []; Plive = []; Pdead = []; culled = []; #initiating objects
#print("starting CullingOneRun")
#print(POSlive[i])
#print(NEGdead[i])
#print(idNEW)
#get list of each infected grid cell ID paired with every other cell ID
#pairedIDs = [reshape(repelem(idNEW,size(CT,1),1),length(idNEW)*size(CT,1),1) repmat((1:size(CT,1))',length(idNEW),1)];
pairedIDs<-matrix(nrow=cells*length(idNEW),ncol=3)
if(length(idNEW)==1){
pairedIDs[1:cells,1]<-idNEW
pairedIDs[1:cells,2]<-1:cells
}

if(length(idNEW)>1){
for(i in 1:length(idNEW)){
if(i==1){
pairedIDs[1:cells,1]<-idNEW[i]	
pairedIDs[1:cells,2]<-1:cells
	}
if(i>1){
	cells=nrow(centroids)
#pairedIDs[(cells+1):nrow(pairedIDs),1]<-idNEW[2]
#pairedIDs[(cells+1):nrow(pairedIDs),2]<-1:cells
	#print(i)
	#print(length(((1)+((i-1)*cells)):((i*cells))))
	#print(dim(pairedIDs))
	#print(paste("length paired id subset:",length(pairedIDs[((1)+((i-1)*cells)):((i*cells)),1])))
	#print(paste("length 1:cells:",length(1:cells)))
pairedIDs[((1)+((i-1)*cells)):((i*cells)),1]<-idNEW[i]
pairedIDs[((1)+((i-1)*cells)):((i*cells)),2]<-1:cells
}
	}
	}

#get distance between each infected grid cell ID paired with every other cell ID
#below results in matrix with one column of each grid cell id and the second column the distance to each infected
#dist = [pairedIDs(:,2) sqrt((CT(pairedIDs(:,1),1)-CT(pairedIDs(:,2),1)).^2 + (CT(pairedIDs(:,1),2)-CT(pairedIDs(:,2),2)).^2)]; % list of grid cell ids with distance to each infected
pairedIDs[,3]<-sqrt((centroids[pairedIDs[,1],1]-centroids[pairedIDs[,2],1])^2 + (centroids[pairedIDs[,1],2]-centroids[pairedIDs[,2],2])^2)

#get all grid cells where the distance between that grid cell and an infected detected pig is less than the pre-determined radius
#temp = find(dist(:,2) <= Rad);
#3 columns: grid cell id for infection, each grid cell id in the zone, distance of each grid cell id in the zone from the infection	
#idout =[pairedIDs(temp,1) pairedIDs(temp,2) dist(temp,2)]; % Store: grid cell id for infection, each grid cell id in the zone, distance of each grid cell id in the zone from the infection
idout<-pairedIDs[pairedIDs[,3]<=10,]

} else{idout<-matrix(nrow=0,ncol=3)}

#print("starting Cull")
#% full set of ids in culling zone with distances to detections (newly infected ID, ID of all others in zone, distance)
#essentially binding idout from above with idZONE
#fullZONE = [idZONE; idout]; 
fullZONE=rbind(idZONE,idout)
#remove NAs
fullZONE<-fullZONE[!is.na(fullZONE[,1]),,drop=FALSE]

#get vector of total pigs in each cell
#pigs = sum(X,1); % vector of total pigs in each cell
#already have this information stored in pop

#get all unique grid cells in the zone		
#allINzone = unique(fullZONE(:,2));
allINzone=unique(fullZONE[,2])
#nrow(fullZONE)
#length(allINzone)

#total number of grid cells in the zone	
#Uall = length(allINzone); % total number of grid cells in the zone
Uall=length(allINzone)

#get total area of the zone	
#% Multiply number of grid cells in zone by area of one grid cell (0.2 * 0.2 = 0.04 km2; 0.4 * 0.4 = 0.16 km2)
#ZONEkm2 = Uall*inc^2; 
ZONEkm2=Uall*inc^2
	
#get number of total number infected pigs in the zone	
#EICinzone = sum(pigs(allINzone)); % current total number of pigs in the zone
soundINzone<-which(pop[,3]%in%allINzone)

#get total number of pigs and infected pigs inside zone
pigsinzone<-sum(pop[soundINzone,1]) #total number of pigs in zone
EICinzone<-sum(pop[soundINzone,9],pop[soundINzone,10],pop[soundINzone,12]) #total number of infected pigs in zone

#get total number of pigs outside the zone
#get total number of infected pigs outside the zone
#EICoutzone
pigsoutzone<-sum(pop[-soundINzone,1]) #total number of pigs in zone
EICoutzone<-sum(pop[-soundINzone,9],pop[-soundINzone,10],pop[-soundINzone,12]) #total number of infected pigs in zone

#get total number of  individuals (inside and outside zone)	
#get total number of infected individuals (inside and outside zone)	
#EIC = sum(X([2 3 5],:),1); % sum all the individuals that have virus
totalpigs=sum(pop[,1])
totalEIC=EICinzone+EICoutzone

###########Culling algorithm	
#if EICinzone > 0 #if there are pigs to cull... 
#print("pigs in zone")
print(pigsinzone)
if(pigsinzone>0){
	
#get number of pigs for each grid cell in zone
#temp2 = [fullZONE pigs(fullZONE(:,2))']; %add number of pigs per grid cell on fullZONE
#% id infected, id of cell in zone, distance between 1 and 2, number of pigs in grid cell (remove grid cells without pigs)	
#temp2 = temp2(temp2(:,4) > 0,:); 
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

#remove rows from fullZONEpigs without pigs
fullZONEpigs<-fullZONEpigs[fullZONEpigs[,4]>0,,drop=FALSE]

#Cullstyle start in, start with closest pigs from detection
#% sort by grid id and then distance (ascending order of total distance) 
#temp2 = sortrows(temp2,[3 -4]);
fullZONEpigs<-as.matrix(arrange(as.data.frame(fullZONEpigs),fullZONEpigs[,3]))

#get total pigs sorted in order of grid cells that are being targeted	
#tpigs = temp2(:,4); % get total pigs sorted in order of grid cells that are being targeted
	
#%density of all live and dead pigs in the zone
#Dr = TotalPIGS/ZONEkm2; 
Dr=pigsinzone/ZONEkm2

#determine density-dependent capture probability in this radius  	
#prob = 1-(1./(1+alphaC).^Dr);
cprob=1-(1/(1+alphaC)^Dr)
	
#get total number culled/removed/sampled in the zone
#numb = binornd(TotalPIGS,prob.*Intensity); %get the total number that will be culled/removed/sampled in the zone:
numb=rbinom(pigsinzone,1,cprob*Intensity)

#get cumulative sum of targeted pigs
#cpigs = cumsum(tpigs);
cpigs=sum(numb)

#determine how far down the list to remove pigs from cells
#id = find(cpigs >= numb,1,'first')
removals=0 #total number of removals, go through loop until first time it is equal to or greater than cpigs
incr=0 #row number where culling stops
while(removals<cpigs){
incr=incr+1
removals<-removals+fullZONEpigs[incr,4]
	}

culled=removals[[1]]

#% list of cells that pigs will be eliminated from (column index was 1 in old version)		
#ids = temp2(1:id,2);	
removalcells<-fullZONEpigs[1:incr,2]

#get total number of pigs culled
#culled = sum(pigs(ids)); %total pigs culled
#is just removals, already done in while loop above
removalpigs<-fullZONEpigs[1:incr,,drop=FALSE]

#}

#%Update surveillance data

#outputs needed
#idZONE,ids,POSlive,POSdead,culled,areaC]
#print(paste("POSlive_locsc174:",POSlive_locs[i-1]))
#print(paste("POSdead_locsc175:",POSdead_locs[i-1]))
#POSlive
#POSlive is a matrix with a row for each timestep
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
POSlive_locs_i=vector(mode="numeric",length=lll)
POSlive_locs[i]<-removalpigs[removalpigs[,6]>0|removalpigs[,7]>0,2]
} else {POSlive_locs_i<-0}

#POSdead_locs
#list of length thyme, each timestep is vector of grid cell locations where dead infected pigs detected at that ts
#removalpigs col 2 where column 9>0
if(POSdead_i>0){
POSdead_locs[[i]]<-removalpigs[removalpigs[,9]>0,2]

lld<-length(removalpigs[removalpigs[,9]>0,2])
POSdead_locs_i=vector(mode="numeric",length=lld)
POSdead_locs_i[i]<-removalpigs[removalpigs[,9]>0,2]
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
#assign fullzone to IDzone, is just without NA
idZONE=fullZONE
#print(idZONE)
#ids:
#removalcells, list of cells where pigs removed from

#culled: 
#total number of pigs culled

#areaC: 
#area of culling zone, ZONEkm2

#Last step: remove selected pigs from population
#removing all pigs in grid cells in removal cells
#assuming removing all pigs in those cells
#find rows in pop equal to those grid cells
#remove those rows from the population
#print(paste("number culled:",culled))
#print(paste("nrow pop before removals:",nrow(pop)))

removalrows<-which(pop[,3] %in% removalcells)
removedpop<-pop[-removalrows,,drop=FALSE]
#print(paste("nrow pop after removals:",nrow(removedpop)))
#print(paste("POSlive_locsc236:",POSlive_locs[[i]]))
#print(paste("POSdead_locsc237:",POSdead_locs[[i]]))

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
output.list[[7]]<-0
output.list[[8]]<-0
output.list[[9]]<-0
output.list[[10]]<-0
output.list[[11]]<-pop	
	
	}

return(output.list)

} #function closing bracket

