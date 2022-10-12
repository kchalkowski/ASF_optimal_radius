areaOfinfection<-function(pop,centroids,inc){
#output:
#vector out: number of infectious indiv., area of infection, max distance between any two cases
#pop[1:4,9]<-1 #for testing function
#print(head(pop))
infected<-pop[pop[,9]>0|pop[,10]>0|pop[,12]>0,,drop=FALSE]	

infectcells<-infected[,3]
	
if(length(infectcells)==1){ #get the area of one grid cell

	#out = [1 inc^2 0];	
out=c(1,inc^2,0)

} else if(length(infectcells)==2){ #area of two grid cells, get distance between their centroids

	#dist = sqrt((CT(id(1),1)-CT(id(2),1)).^2 + (CT(id(1),2)-CT(id(2),2)).^2);
dist=sqrt((centroids[infectcells,][1,1]-centroids[infectcells,][2,1])^2 + (centroids[infectcells,][1,2]-centroids[infectcells,][2,2])^2)
out=c(2,2*(inc^2),dist)

} else {

#[~,A] = boundary(CT(id,:));
#out = [length(id) A distmat];
#isolate xy points of centroids of cells with infected individuals
#print("start aoi print statement")
#print(pop[is.na(pop[,3]),])
#print(infectcells)
xy<-centroids[infectcells,]	
#print(xy)
#print(anyNA(infected))
#print(head(infected))
#print(infected)
#print(pop[rowSums(is.na(pop)) > 0,])
#print("end aoi print statement")
#print(xy)
#use chull on object of xy points
#ch<-chull(xy)
ch<-chull(xy)

#subset xy point set to boundary points from chull
#boundary<-xy[ch,]
boundary<-xy[ch,]

#get area within boundary
#abs(polyarea(boundary[,1],boundary[,2]))	
A<-abs(polyarea(boundary[,1],boundary[,2]))

#distmat = max(max(squareform(pdist(CT(id,:)))));
#pdist- gets pairwise distance between pairs of observations
#squareform- puts pairwise distance matrices into a square format showing the pairs
#get the max of these
#install.packages("rdist")
#library(rdist)
maxdist<-max(pdist(xy))

out=c(length(infectcells),A,maxdist)

	} #greater than 2 infected closing bracket
return(out)
	} #function closing bracket



###################
