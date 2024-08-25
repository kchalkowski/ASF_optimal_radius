#FOI_NEW
FOI_R<-function(pop,centroids,cells,B1,B2,F1,F2,F2i){
  #Create I/C matrices. Just want vector with locations of I/C, with nums of I/C in right place
  I=matrix(0,nrow=cells,ncol=1)
  C=matrix(0,nrow=cells,ncol=1)
  I[pop[which(pop[,10]!=0),3],1]=pop[which(pop[,10]!=0),10]
  C[pop[which(pop[,12]!=0),3],1]=pop[which(pop[,12]!=0),12]
  
  Pse=matrix(0,nrow=cells,ncol=1)
  temp=I+C
  id=which(temp>0)
  W=matrix(0,nrow=cells,ncol=2)
  #test
  #F1=0.7381
  W[,1]=F1*I
  W[,2]=F1*C
  
  B=zeros(length(id),cells)
  
  if(length(id)==1){
  X1=matrix(nrow=length(id),ncol=1)
  Y1=matrix(nrow=length(id),ncol=1)
  
  X1=centroids[id,1]
  Y1=centroids[id,2]
  
  #in real one, remove brackets after X1-- just testing it out
  dist=sqrt((centroids[,1]-X1)^2+(centroids[,2]-Y1)^2)
  dist=as.data.frame(dist)
  colnames(dist)="X"
  prob=predict(F2,dist,type="response")
  probi=predict(F2i,dist,type="response")
  prob[dist==0]=0
  probi[dist==0]=0
  
  B=B1*I[id[1]]*prob + B2*C[id[1]]*probi
  }
  
  if(length(id)>1){
    X1=t(repmat(centroids[id,1],cells,1))
    Y1=t(repmat(centroids[id,2],cells,1))
    X2=repmat(centroids[,1],length(id),1)
    Y2=repmat(centroids[,2],length(id),1)
    
    dist=sqrt((X2-X1)^2+(Y2-Y1)^2)
    #checked dist vals with ML version, look good
    
    #obj=reshape(,nrow(dist),ncol(dist),direction="long")
    obj=c(dist)
    obj=as.data.frame(obj)
    colnames(obj)="X"
    prob.temp=predict(F2,obj,type="response")
    probi.temp=predict(F2i,obj,type="response")
    
    prob=matrix(prob.temp,nrow = length(id),ncol = cells)
    probi=matrix(probi.temp,nrow = length(id),ncol = cells)
    
    prob[dist==0]<-0
    probi[dist==0]<-0
    
    B=B1*t(repmat(I[id],cells,1))*prob+B2*t(repmat(C[id],cells,1))*probi
    #sum(B) in R=351.743
    #sum(B,"all") in ML=351.743
    
  }
  #Pse=1-exp(-W(1,:)-W(2,:)-sum(Bm1));  
  #obj=colSums(B)
    #sum of obj in R and ML- 351.742
  if(length(id)==1){
  Pse=1-exp(-W[,1]-W[,2]-t(B))
  } else{
    if(length(id)>1){
      Pse=1-exp(-W[,1]-W[,2]-t(colSums(B)))
    }
  }
  #sum Pse in ML and R- 364.5351
  return(Pse)
  
}
