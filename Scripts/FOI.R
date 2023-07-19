FOI<-function(pop,centroids,cells,B1,B2,F1,F2,Fi){
#Pse = zeros(1,cells);
Pse=matrix(nrow=cells,ncol=1)	
#F1 contact probability at 0 distance
#F2 contact probability live pigs at other distances (define based on landscape dist matrix)
#F2i contact probability infected carcasses at other distances
#B1 = transmission probability given contact for direct contact (a constant)
#B2 = transmission probability given contact for indirect contact (a constant)
if(length(pop[pop[,10]>0|pop[,12]>0,1])>0){#if there are any infectious agents (2=I; 4=C)
################
	#within cell transmission
#W = [F1.*I; F1.*C]; #contact probability at 0 distance, get this for each infectious individual
#B = zeros(length(id),cells); % creates a matrix with row=num infectious individuals, cols=cells


#between cell transmission
#get the centroid of each infected individual
#X1 = CT(id(i),1); Y1 = CT(id(i),2); %centroid of infected
#for one infected individual, one number each for X1 Y1
#have this already, would be X1=pop[pop[,2]==2,][5]; Y1=X1=pop[pop[,2]==2,][6]
	
#distance between each infected individual, and all other cells
#dist = sqrt((CT(:,1)-X1).^2 + (CT(:,2)-Y1).^2); % distance to all other cells
#is vector with 1 col, 40000 rows (bc 40000 cells)
num_I<-sum(pop[pop[,10]>0,10])
num_C<-sum(pop[pop[,12]>0,12])
Imat<-pop[pop[,10]>0,,drop=FALSE]
Cmat<-pop[pop[,12]>0,,drop=FALSE]
dist_I=matrix(nrow=num_I, ncol=cells)
dist_C=matrix(nrow=num_C, ncol=cells)
prob_I=matrix(nrow=num_I, ncol=cells)
prob_C=matrix(nrow=num_C, ncol=cells)
#B=matrix(nrow=num_C,ncol=cells)

#initialize 3d array for probs/dists for I and C
#dim=c(no rows, no cols, no arrays)
#paste0("I_",1:num_I)
#paste0("cell_",1:cells)

#dimnames_I=list(paste0("I_",1:num_I),paste0("cell_",1:cells),c("dist","prob"))
#dimnames_C=list(paste0("C_",1:num_C),paste0("cell_",1:cells),c("dist","prob"))
#pdI<-array(0,dim=c(num_I,cells,2), dimnames=list(paste0("I_",1:num_I),paste0("cell_",1:cells),c("dist","prob")))
#pdC<-array(0,dim=c(num_C,cells,2), dimnames=list(paste0("I_",1:num_I),paste0("cell_",1:cells),c("dist","prob")))

#print(paste0("nrow Imat:",nrow(Imat)))
#dist = sqrt((CT(:,1)-X1).^2 + (CT(:,2)-Y1).^2)
if(nrow(Imat)>0){
dimnames_I=list(paste0("I_",1:num_I),paste0("cell_",1:cells),c("dist","prob"))
pdI<-array(0,dim=c(num_I,cells,2), dimnames=list(paste0("I_",1:num_I),paste0("cell_",1:cells),c("dist","prob")))

for(i in 1:nrow(Imat)){
pdI[i,,1]=sqrt((centroids[,1]-Imat[,5])^2 + (centroids[,2]-Imat[,6])^2)
pdI[i,,2]=exp(predict(F2,newdata=data.frame(xx=pdI[i,,1])))
#pdI[,pdI[i,,1:2][,1]==0,][[2]]<--F1
}
}

#C matrix###############3
#print(paste0("nrow Cmat:",nrow(Cmat)))
if(nrow(Cmat)>0){
dimnames_C=list(paste0("C_",1:num_C),paste0("cell_",1:cells),c("dist","prob"))
#print(dimnames_C)
#print(length(dimnames_C))
#print(c(num_C,cells,2))
#print(length(c(num_C,cells,2)))
pdC<-array(0,dim=c(num_C,cells,2), dimnames=dimnames_C)

for(i in 1:nrow(Cmat)){
pdC[i,,1]=sqrt((centroids[,1]-Cmat[,5])^2 + (centroids[,2]-Cmat[,6])^2)
#}

#convert vector of distances to contact probabilities 	
#prob = predict(F2,dist); probi = predict(F2i,dist);


#if(nrow(Cmat)>0){
#for(i in 1:nrow(Cmat)){
pdC[i,,2]=exp(predict(F2,newdata=data.frame(xx=pdC[i,,1])))
#}
#}

#predict(F2,newdata=data.frame(xx=dist_I[1,]))

#set the prob at 0 distance = 0 because that is covered in within-group FOI term	
#prob(id(i)) = 0; % set the prob at 0 distance = 0 because that is covered in within-group FOI term	
#W = [F1.*I; F1.*C]; #contact probability at 0 distance, get this for each infectious individual

	
#if(nrow(Cmat)>0){
#for(i in 1:nrow(Cmat)){
#pdC[,pdC[i,,1:2][,1]==0,][[2]]<--F1
#}
}
}
#B- FOI matrix from infected to all other cells
#I- matrix of locations for infected individuals (Imat[,7])
#C- matrix of locations for infected carcasses (Cmat[,7])
#id- row number of grid/centroids where there are infected individuals or carcasses (Imat[,7])(Cmat[,7])
#B1- transmission probability given contact for indirect contact with infected live pigs 
#B2- transmission probability given contact for indirect contact with infected carcasses
#B- matrix of FOI from all infected cells

#B(i,:) = -B1.*I(id(i)).*prob' - B2.*C(id(i)).*probi'; %FOI from infected cell id(i) to all other cells
#B=matrix(nrow=num_C/I,ncol=cells)

#force of infection from each infected pig
#-B1.*I(id(i)).*prob'
#force of infection from infected carcasses
#B2.*C(id(i)).*probi'

if(nrow(Imat)>0){
B_I=matrix(nrow=cells,ncol=dim(pdI)[1])
B_I[,1]=0
for(i in 1:dim(pdI)[1]){
B_I[,i]=(pdI[i,,2]*B1)
}
B_I<-rowSums(B_I)
for(i in 1:nrow(Imat)){
I_direct<-which(pdI[i,,1:2][,1]==0)
B_I[I_direct]=(F1)
}
}

if(nrow(Cmat)>0){
B_C=matrix(nrow=cells,ncol=dim(pdC)[1])
B_C[,1]=0
for(i in 1:dim(pdC)[1]){
B_C[,i]=(pdC[i,,2]*B2)
}
B_C<-rowSums(B_C)
for(i in 1:nrow(Cmat)){
C_direct<-which(pdC[i,,1:2][,1]==0)
B_C[C_direct]=(F1)
}
}

#for(i in 1:dim(pdI)[1]){
#B[,i]=(pdI[i,,2]*B1)+(pdC[i,,2]*B2)
#}


if(nrow(Cmat)>0&nrow(Imat)>0){
	Bsum=B_I+B_C
} else if(nrow(Imat)>0){
	Bsum=B_I
	} else if(nrow(Cmat)>0){Bsum=B_C}


#Bsum<-as.matrix(rowSums(B))

#Pse = 1-exp(-W(1,:) - W(2,:) - sum(B,1));  %indirect transmission from other cells (kernel values > 0 distance)
Pse=1-exp(-Bsum)

#Pse[Pse[,1]>0.001,]
#Pse[Pse<0]<-0
#now, want the Pse for each member of the population

} else {Pse<-matrix(0,nrow=cells,ncol=1)} 	#if any infectious closing bracket
#which(Pse>0.001)
#which(Pse>0.01)
#print("any FOI 1?")
#print(any(is.na(Pse)))
#print(any(is.na(Pse)))

return(Pse)

	} #function closing bracket
