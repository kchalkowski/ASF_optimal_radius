###function to mimic apparent prevalence with user-input detection probability
#Will later have to incorporate apparent prev into model structure with control
#this function intended to be used as passive surveilance method, on 'sounderlocs' output
dat=Izone
SEIRCZ.rep2=SEIRCZ.rep
type="zonecells"
ExtApparentPrev<-function(dat,DetP,N,rep,type){
  require(plyr)
  
  if(type=="fullzone"){
    SEIRCZ.rep=dat
  }
  
  if(type=="onlyzone"){
    dat=dat[dat$inZone==1,]
    SEIRCZ.rep=dat %>% 
      dplyr::group_by(timestep) %>% 
      dplyr::summarise(S=sum(S),
                       E=sum(E),
                       I=sum(I),
                       R=sum(R),
                       C=sum(C),
                       Z=sum(Z)) %>%
      as.data.frame()
  }
  
  if(type=="zonecells"){
    dat=dat[dat$inZone==1,]
    SEIRCZ.rep=dat %>% 
      dplyr::group_by(timestep,cell) %>% 
      dplyr::summarise(S=sum(S),
                       E=sum(E),
                       I=sum(I),
                       R=sum(R),
                       C=sum(C),
                       Z=sum(Z)) %>%
      as.data.frame()
  }
  
  
  #convert each timestep (week) into rep lists
  for(r in 1:nrow(SEIRCZ.rep)){
    #print(1)
  S<-rep("S",SEIRCZ.rep[r,which(colnames(SEIRCZ.rep)=="S")]) 
  E<-rep("S",SEIRCZ.rep[r,which(colnames(SEIRCZ.rep)=="E")]) #Would E have antibodies? assign as S for now
  I<-rep("I",SEIRCZ.rep[r,which(colnames(SEIRCZ.rep)=="I")]) 
  R<-rep("S",SEIRCZ.rep[r,which(colnames(SEIRCZ.rep)=="R")]) #Would R have antibodies? assign as S for now
  C<-rep("C",SEIRCZ.rep[r,which(colnames(SEIRCZ.rep)=="C")])
  Z<-rep("Z",SEIRCZ.rep[r,which(colnames(SEIRCZ.rep)=="Z")])
  #print(1)
  allpig=c(S,E,I,R)
  samples=sample(allpig,min(N,length(allpig)))
  counts.r=plyr::count(samples)
  counts.r$appar<-0
  #print(1)
  allcarcass=c(C,Z)
  csamples=sample(allcarcass,min(N,length(allcarcass)))
  counts.c=plyr::count(csamples)
  #counts.c$appar<-0
  #print(1)
  if(length(counts.r[counts.r$x=="I",]$freq)>0){
  num.appar.inf<-sum(rbinom(c(1,0),counts.r[counts.r$x=="I",]$freq,prob=c(DetP,(1-DetP))))
  num.add.sus=counts.r[counts.r$x=="I",]$freq-num.appar.inf
  counts.r[counts.r$x=="S",]$appar=counts.r[counts.r$x=="S",]$freq+num.add.sus
  counts.r[counts.r$x=="I",]$appar=num.appar.inf
  } else{
    counts.r$appar=counts.r$freq
  }
  #print(1)
  if(nrow(counts.c)>0){
  if(length(counts.c[counts.c$x=="C",]$freq)>0){
    counts.c$appar=NA
    num.appar.inf<-sum(rbinom(c(1,0),counts.c[counts.c$x=="C",]$freq,prob=c(DetP,(1-DetP))))
    num.add.sus=counts.c[counts.c$x=="C",]$freq-num.appar.inf
    counts.c[counts.c$x=="Z",]$appar=counts.c[counts.c$x=="Z",]$freq+num.add.sus
    counts.c[counts.c$x=="C",]$appar=num.appar.inf
  } else{
    counts.c$appar=counts.c$freq
  }
  
  } else{
    counts.c=data.frame(x="Z",freq=0,appar=0)
  }
  
  counts.r$src="alive"
  counts.c$src="carcass"
  
  counts.rep=rbind(counts.r,counts.c)
  counts.rep$week=r
  
  if(type=="zonecells"){
    counts.rep$cell=SEIRCZ.rep[r,]$cell
  }
  
  #print(1)
  if(r==1){
    counts=counts.rep
  } else{
    counts=rbind(counts,counts.rep)
  }
  }

  counts$rep=rep
  counts=pivot_wider(counts,names_from=x,values_from=c(freq,appar)) %>% as.data.frame
  counts[is.na(counts)]<-0
  
  counts$freq_prev<-0
  counts[counts$src=="alive",]$freq_prev=counts[counts$src=="alive",]$freq_I/(counts[counts$src=="alive",]$freq_I+counts[counts$src=="alive",]$freq_S)
  counts[counts$src=="carcass",]$freq_prev=counts[counts$src=="carcass",]$freq_C/(counts[counts$src=="carcass",]$freq_C+counts[counts$src=="carcass",]$freq_Z)
  counts[is.na(counts)]<-0
  
  
  counts$appar_prev<-0
  counts[counts$src=="alive",]$appar_prev=counts[counts$src=="alive",]$appar_I/(counts[counts$src=="alive",]$appar_I+counts[counts$src=="alive",]$appar_S)
  counts[counts$src=="carcass",]$appar_prev=counts[counts$src=="carcass",]$appar_C/(counts[counts$src=="carcass",]$appar_C+counts[counts$src=="carcass",]$appar_Z)
  
  #move rep to front

  if(type=="zonecells"){
    counts=counts[,c(1,4,2,3,5:14)]
  } else{
    counts=counts[,c(1,3,2,4:13)]
  }
  
  #reorganize to fit outbreak data better
  counts.alive=counts[counts$src=="alive",]
  counts.carcass=counts[counts$src=="carcass",]
  
  #need to make this more general.. indices with logical not great
  if(type=="zonecells"){
    counts.alive=counts.alive[,c(1:6,9,10,13,14)]
    counts.carcass=counts.carcass[,c(1:4,7,8,11,12,13,14)]
  } else{
    counts.alive=counts.alive[,c(1:4,7,8,11,12,13)]
    counts.carcass=counts.carcass[,c(1:3,5,6,9,10,12,13)]
  }
  
  if(type=="zonecells"){
    colnames(counts.alive)[c(5:8)]<-c("freq_Nneg","freq_Npos","appar_Nneg","appar_Npos")
    colnames(counts.carcass)[c(5:8)]<-c("freq_Nneg","freq_Npos","appar_Nneg","appar_Npos")
  } else{
    colnames(counts.alive)[c(4:7)]<-c("freq_Nneg","freq_Npos","appar_Nneg","appar_Npos")
    colnames(counts.carcass)[c(4:7)]<-c("freq_Nneg","freq_Npos","appar_Nneg","appar_Npos")
  }
  
  allcounts=rbind(counts.alive,counts.carcass)
  allcounts[is.na(allcounts)]<-0
  
  if(type=="zonecells"){
    #summarize across cells within zone for each week
    allcounts=
      allcounts %>% 
      dplyr::group_by(src,rep,week) %>% 
      dplyr::summarise(freq_Nneg_md=median(freq_Nneg),
                       freq_Nneg_25=quantile(freq_Nneg,0.25),
                       freq_Nneg_75=quantile(freq_Nneg,0.75),
                       freq_Npos_md=median(freq_Npos),
                       freq_Npos_25=quantile(freq_Npos,0.25),
                       freq_Npos_75=quantile(freq_Npos,0.75),
                       appar_Nneg_md=median(appar_Nneg),
                       appar_Nneg_25=quantile(appar_Nneg,0.25),
                       appar_Nneg_75=quantile(appar_Nneg,0.75),
                       appar_Npos_md=median(appar_Npos),
                       appar_Npos_25=quantile(appar_Npos,0.25),
                       appar_Npos_75=quantile(appar_Npos,0.75),
                       freq_prev_md=median(freq_prev),
                       freq_prev_25=quantile(freq_prev,0.25),
                       freq_prev_75=quantile(freq_prev,0.75),
                       appar_prev_md=median(appar_prev),
                       appar_prev_25=quantile(appar_prev,0.25),
                       appar_prev_75=quantile(appar_prev,0.75)) %>%
      as.data.frame()
  
  }
  
  return(allcounts)
  
}
