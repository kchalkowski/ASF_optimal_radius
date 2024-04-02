GetOutputs<-function(pop,Incidence,Tculled,ICtrue,out,detectday){
#List of outputs created here:	
	#Tinc #sum of all exposures over simulation 
	#sum(Tculled)
	#idT #last day there is an infectious individual
	#Mspread #max spread of infection
	#IConDD #number of I, C, and E on detection day
	#ICatEnd #number of I,C,E on last day 
	#TincFromDD #sum of all exposures starting day after detection day
	#TincToDD #sum of all exposures up until detection day
	#DET #total number of detections

#Tinc, this just sums all of the exposures
Tinc=sum(Incidence)

#sum all culled
sumTculled=sum(Tculled)

#Find last day there was an infectious individual
if(any(ICtrue!=0)){
idT=which(ICtrue!=0)[length(which(ICtrue!=0))]
} else{idT=1}

#Find max spread of infection
Mspread<-max(out[,2])

#IConDD #number of I, C, and E on detection day TROUBLESHOOT
IConDD=ICtrue[detectday]

#ICatEnd #number of I,C,E on last day
ICatEnd=ICtrue[idT]

#TincToDD #sum of all exposures up until detection day
TincToDD<-sum(Incidence[idT:detectday])

#TincFromDD #sum of all exposures starting day after detection day
TincFromDD<-sum(Incidence[detectday:idT])

#DET #total number of detections
DET=sum(unlist(POSlive),unlist(POSdead))

#send to list
list.all=list("Tinc"=Tinc,
							"sumTculled"=sumTculled,
							"idT"=idT,
							"Mspread"=Mspread,
							"IConDD"=IConDD,
							"ICatEnd"=ICatEnd,
							"TincToDD"=TincToDD,
							"TincFromDD"=TincFromDD,
							"DET"=DET)

return(list.all)

	}