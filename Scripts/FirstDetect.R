FirstDetect<-function(pop,i,POSlive,POSdead,POSlive_locs,POSdead_locs){
#randomly detects a single infected pig (carcass or live)
	
	#Randomly detect an E,I, or C pig
	detection<-as.integer(sample(as.character(which(pop[,9]>0|pop[,10]>0|pop[,12]>0)),1))
	
	#If chose row with combination of E, I, C, pick one
	if(sum(pop[detection,c(9,10,12)])>1){
		found=sample(c(9:12),1)
	} else{ #else, found is just the one
			found=c(9,10,12)[which(pop[detection,c(9,10,12)]>0)]
		}
	
	##Set POSlive, POSdead
	if(found==9|found==10){ #if live one (E,I) detected, record 1 live
		POSlive[[i]]<-1
		POSdead[[i]]<-0
	} else{
		if(found==12){ #if dead one detected (C), record 1 dead
			POSlive[[i]]<-0
			POSdead[[i]]<-1
			}	
		}

  #update the surveillance data
  if(POSlive[[i]]>0){POSlive_locs[[i]]<-pop[detection,3]}
  if(POSdead[[i]]>0){POSdead_locs[[i]]<-pop[detection,3]}
  
  #Store the pop rows of the detected pigs
  detected<-pop[detection,,drop=FALSE]
  
  #Remove the detected case
  pop[detection,found]<-pop[detection,found]-1
  
  #if found was live, update number of pigs in sounder with infected pig/carcass (carcasses already removed from abundance count during state change function)
  if(found==9|found==10){
  	pop[detection,1]<-pop[detection,1]-1
  }
  
  #if every pig in cell with infected pigs removed, including carcasses, remove the row from the population
  if(pop[detection,1]==0&pop[detection,12]==0&pop[detection,13]==0){
  	pop<-pop[-detection,]
  	}

  return(list(pop,POSlive,POSdead,POSlive_locs,POSdead_locs))
  
}