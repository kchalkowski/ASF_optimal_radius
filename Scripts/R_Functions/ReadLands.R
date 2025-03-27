ReadLands<-function(path){

  fs=list.files(path,full.names=TRUE)
  plands_list=vector(mode="list",length=length(fs))
  
  for(f in 1:length(fs)){
    #print(paste0("reading land ",f," of ",length(fs)))
    print(fs[f])
    plands_list[[f]]=terra::rast(fs[f])
  }
  
  print("formatting list as sprc")
  plands_sprc<-terra::sprc(plands_list)
  
  return(plands_sprc)
}