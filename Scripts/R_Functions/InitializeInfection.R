InitializeInfection<-function(pop,centroids,grid,parameters){
	
######################################
######## Initialize Infection ######## 
######################################
#num_inf_0=1 #how many pigs to infect starting off
num_inf_0=parameters$num_inf_0
	
	
#find the midpoint of the grid
midpoint=c(median(centroids[,1]),median(centroids[,2]))
id=which(centroids[,1]>=midpoint[1]&centroids[,2]>=midpoint[2])[1] #location on grid closest to midpoint
infected=InitializeSounders(centroids,grid,c(id,num_inf_0),pop_init_type="init_single",pop_init_grid_opts="homogeneous")
infected[,8]<-0
infected[,10]<-1

#combine infected pig with pop matrix
pop<-rbind(pop,infected)

return(pop)

}
