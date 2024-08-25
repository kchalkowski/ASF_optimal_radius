areaOfinfection<-function(pop,centroids,inc){

#subset out infected rows
infected<-pop[pop[,9]>0|pop[,10]>0|pop[,12]>0,,drop=FALSE]	

#determine which cells infected
infectcells<-unique(infected[,3])
	
if(length(infectcells)==1){ #get the area of one grid cell

out=c(1,inc^2,0)

} else if(length(infectcells)==2){ #area of two grid cells, get distance between their centroids

dist=sqrt((centroids[infectcells,][1,1]-centroids[infectcells,][2,1])^2 + (centroids[infectcells,][1,2]-centroids[infectcells,][2,2])^2)
out=c(2,2*(inc^2),dist)

} else {

#get all centroids xy coords
xy<-centroids[infectcells,,drop=FALSE]	

#compute convex hull of points
ch<-chull(xy)

#get boundary coordinates of convex hull
boundary<-xy[ch,,drop=FALSE]

#get area within boundary
A<-abs(polyarea(boundary[,1,drop=FALSE],boundary[,2,drop=FALSE]))

#distmat = max(max(squareform(pdist(CT(id,:)))));
#pdist- gets pairwise distance between pairs of observations
#squareform- puts pairwise distance matrices into a square format showing the pairs
#get the max of these
#install.packages("rdist")
#library(rdist)
maxdist<-max(pdist(xy))

#number of unique infected cells, area of infection, max distance between infected cells
out=c(length(infectcells),A,maxdist)

	} #greater than 2 infected closing bracket
return(out)
	} #function closing bracket
