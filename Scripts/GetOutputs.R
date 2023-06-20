GenerateOutputs<-function(pop,Incidence,ICtrue,out,detectday){
#List of outputs created here:	
	#Tinc #sum of all exposures over simulation 
	#idT #last day there is an infectious individual
	#Mspread #max spread of infection
	#IConDD #number of I, C, and E on detection day TROUBLESHOOT
	#ICatEnd #number of I,C,E on last day 
	#TincFromDD #sum of all exposures starting day after detection day
	#TincToDD #sum of all exposures up until detection day
	#DET #total number of detections

#Tinc, this just sums all of the exposures
Tinc=sum(Incidence)

#Find last day there was an infectious individual
idT=which(ICtrue!=0)[length(which(ICtrue!=0))]

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

list.all=list(Tinc,idT,Mspread,IConDD,ICatEnd,TincToDD,TincFromDD,DET)

return(list.all)

	}