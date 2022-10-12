FastMovement=function(pop,centroids,shift,inc){

#get distances from gamma distribution
pop[,4]=rgamma(nrow(pop),shape=shift[1],scale=shift[2])

#set those less than inc to 0
pop[pop[,4]<inc,][,4]=0 

#set present locations to previous locations
pop[,7]=pop[,3]

#get new locations with movement func
pop[,3]<-parallelMovementRcpp_portion(pop,pop[,1,drop=FALSE],pop[,3,drop=FALSE],centroids)	
#parallelMovementRcpp_portion(pop[,5,drop=FALSE],pop[,6,drop=FALSE],pop[,4,drop=FALSE],pop[,1,drop=FALSE],pop[,3,drop=FALSE],centroids)

return(pop)
}
#pop[,4]
#parallelMovementRcpp_portion(pop[,5,drop=FALSE],pop[,6,drop=FALSE],pop[,4,drop=FALSE],pop[,1,drop=FALSE],pop[,3,drop=FALSE],centroids)

#for(i in 1:100){
#FastMovement(pop,centroids,shift,inc)
#print(pop[pop[,9]>0&pop[,10]>0,])
#}


#View(pop)
#microbenchmark(FastMovement(pop,centroids,shift,inc))
#77 millseconds