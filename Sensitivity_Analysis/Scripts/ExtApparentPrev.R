###function to mimic apparent prevalence with user-input detection probability
#Will later have to incorporate apparent prev into model structure with control
#this function intended to be used as passive surveilance method, on 'sounderlocs' output

ExtApparentPrev<-function(SEIRCZ.rep,DetP,N,rep){
  require(plyr)
  #convert each timestep (week) into rep lists
  for(r in 1:nrow(SEIRCZ.rep)){
  S<-rep("S",SEIRCZ.rep[r,2]) 
  E<-rep("S",SEIRCZ.rep[r,3]) #Would E have antibodies?
  I<-rep("I",SEIRCZ.rep[r,4]) 
  R<-rep("S",SEIRCZ.rep[r,5]) #Would R have antibodies?
  allpig=c(S,E,I,R)
  samples=sample(allpig,min(N,length(allpig)))
  counts.r=plyr::count(samples)
  counts.r$appar<-0
  if(length(counts.r[counts.r$x=="I",]$freq)>0){
  num.appar.inf<-sum(rbinom(c(1,0),counts.r[counts.r$x=="I",]$freq,prob=c(DetP,(1-DetP))))
  num.add.sus=counts.r[counts.r$x=="I",]$freq-num.appar.inf
  counts.r[counts.r$x=="S",]$appar=counts.r[counts.r$x=="S",]$freq+num.add.sus
  counts.r[counts.r$x=="I",]$appar=num.appar.inf
  } else{
    counts.r$appar=counts.r$freq
  }
  counts.r$timestep=r
  if(r==1){
    counts=counts.r
  } else{
    counts=rbind(counts,counts.r)
  }
  }

  counts$rep=rep
  
  counts=pivot_wider(counts,names_from=x,values_from=c(freq,appar)) %>% as.data.frame
  counts[is.na(counts)]<-0
  counts$freq_prev=counts$freq_I/(counts$freq_I+counts$freq_S)
  counts$appar_prev=counts$appar_I/(counts$appar_I+counts$appar_S)
  
  #move rep to front
  counts=counts[,c(2,1,3:8)]
  
  return(counts)
}
