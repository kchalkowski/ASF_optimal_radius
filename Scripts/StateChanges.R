#outputs: pop,Incidence,BB
StateChanges<-function(pop,centroids,cells,Pbd,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B,K,death,Pcr,Pir,Incidence,BB,i){
####################################################################
########### Initialize state change probability matrices ########### 
####################################################################

#births
Sdpb=matrix(nrow=nrow(pop),ncol=1)
Sdpb[,1]=0
	
#natural deaths
Sdpd<-matrix(nrow=nrow(pop),ncol=1)
Edpd<-matrix(nrow=nrow(pop),ncol=1)
Idpd<-matrix(nrow=nrow(pop),ncol=1)
Rdpd<-matrix(nrow=nrow(pop),ncol=1)
Sdpd[,1]=0
Edpd[,1]=0
Idpd[,1]=0
Rdpd[,1]=0

#disease state change recording
Eep=matrix(nrow=nrow(pop),ncol=1)
Eep[,1]=0
Iep=matrix(nrow=nrow(pop),ncol=1)
Iep[,1]=0
Rep=matrix(nrow=nrow(pop),ncol=1)
Rep[,1]=0
Cep=matrix(nrow=nrow(pop),ncol=1)
Cep[,1]=0

#Carcass decay recording
Ccd=matrix(nrow=nrow(pop),ncol=1)
Ccd[,1]=0
Zcd=matrix(nrow=nrow(pop),ncol=1)
Zcd[,1]=0
	
########################################
########### Determine Births ########### 
########################################

#subset sounder sets with live, uninfected individuals
idN=pop[pop[,8,drop=FALSE]>0|pop[,9,drop=FALSE]>0|pop[,11,drop=FALSE]>0,]

#Number of live, uninfected individuals
liveind<-sum(colSums(pop)[c(8,9,11)])

#get row indices of live individuals
liverows<-which(pop[,8,drop=FALSE]>0|pop[,9,drop=FALSE]>0|pop[,11,drop=FALSE]>0) #rownums with live indiv

#density-dependent birth rate
Brate=Pbd*liveind*(1-liveind/K)

#get total births, using Brate as mean in a poisson
Tbirths=rpois(1,Brate)

#record total births this time step
BB[i]=Tbirths

#Pick out enough numbers that sum to Tbirths - this will determine how many cells get births
id<-0
n=1
while (sum(id) < Tbirths) {
    birthset_i <- round(runif(1,min=0, max=min(Tbirths,10)))
    id[n]<-birthset_i
    n=n+1
}

#if there are some births
if(length(id)>1){

#if there are more live sounders than births needed	
if(length(id)<=nrow(idN)){

#%pick which cells with pigs will get the births
id2=sample(1:nrow(idN),length(id))

} else {
#%if there are more births than cells only add births cells where the pigs are (so fewer births will be happening)
id2=sample(1:length(liverows),length(liverows))
}
}

for(j in 1:length(id)){
Sdpb[id2[j],1]<-id[j]
	}


##############################################################
######## Determine disease state change probabilities ######## 
##############################################################

#Pse<-FOI(pop,centroids,cells,B1,B2,F1,F2) #force of infection #R version
Pse<-FOIParallelFull(pop,centroids,cells,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B) #cpp parallel version, 22x faster than R version

Pei=1-exp(-1/(rpois(cells,4)/7)) #transitions exposure to infected
Pic=1-exp(-1/(rpois(cells,5)/7)) #transitions infected to either dead or recovered


###############################################
######## Conduct the State Transitions ######## 
###############################################

for(k in 1:nrow(pop)){

#operations on Susceptible individuals
if(pop[k,8]>0){
#print(pop[k,3])
#print(pop[k,8])
Sdpd[k]<-sum(rbinom(pop[k,8],1,death))
Eep[k]<-sum(rbinom(pop[k,8],1,Pse[pop[k,3]])) #Exposure (S -> E) infection based on probability using their location
#print("popk")
#print(pop[k,8])
#print(Pse[pop[k,3]])
#print(Eep[k])
}	

#operations on Exposed individuals
if(pop[k,9]>0){
Edpd[k]<-sum(rbinom(pop[k,9],1,death))
Iep[k]<-sum(rbinom(pop[k,9],1,Pei[k]))
}

#operations on Infected individuals	
if(pop[k,10]>0){
#Idpd[k]<-sum(rbinom(pop[k,10],1,death)) #idpd not in matlab model! maybe reason for lower fadeout rate
Rep[k]<-sum(rbinom(pop[k,10],1,Pir*Pic[k]))
Cep[k]<-sum(rbinom(pop[k,10],1,(1-Pir)*(Pic[k]))) 
}	

#operations on Recovered individuals
if(pop[k,11]>0){
Rdpd[k]<-sum(rbinom(pop[k,11],1,death))
}	

#operations on Carcasses (infected)
if(pop[k,12]>0){
Ccd<-sum(rbinom(pop[k,12],1,Pcr))	
}	

#operations on Carcasses (uninfected)
if(pop[k,13]>0){
#operations on Uninf carcass individuals	
Zcd<-sum(rbinom(pop[k,13],1,Pcr))	
	}	
}

Incidence[i]<-Incidence[i]+sum(Eep)


###################################
######## Update pop matrix ######## 
###################################

#update states in pop matrix
###################################
pop[,8]=pop[,8]-Eep+Sdpb-Sdpd #S
pop[,9]=pop[,9]-Iep+Eep-Edpd #E
pop[,10]=pop[,10]-Rep-Cep+Iep#-Idpd#I
pop[,11]=pop[,11]+Rep-Rdpd #R
pop[,12]=pop[,12]+Cep-Ccd#+Idpd #C
pop[,13]=pop[,13]+Sdpd+Rdpd+Edpd-Zcd #Z

#sometimes end up with negative numbers 
#(i.e. all pigs in sounders chosen for natural mort and disease mort)
#just set anything below zero to zero
pop[which(pop[,8]<0),8]<-0
pop[which(pop[,9]<0),9]<-0
pop[which(pop[,10]<0),10]<-0
pop[which(pop[,11]<0),11]<-0
pop[which(pop[,12]<0),12]<-0
pop[which(pop[,13]<0),13]<-0

#move dead individuals (C or Z) into their own rows
#pop[,12] and pop[,13] > 0
deadguys<-pop[pop[,12]>0|pop[,13]>0,,drop=FALSE]

#if there are deadguys....
if(nrow(deadguys)!=0){
#remove abundance and all live guy counts from deadguy set
deadguys[,1]=0
deadguys[,8]=0
deadguys[,9]=0
deadguys[,10]=0
deadguys[,11]=0

#set all deadguys in pop rows to zero
pop[which(pop[,12]>0),12]<-0
pop[which(pop[,13]>0),13]<-0

#add deadguys to pop matrix
pop<-rbind(pop,deadguys)

}

#Update abundance numbers (live individuals only count in abundance)
pop[,1]=rowSums(pop[,8:11])

return(list(pop,Incidence,BB))

}
