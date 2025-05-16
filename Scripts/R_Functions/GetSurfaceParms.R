#GetSurfaceParms(parameters00,plands_sprc[1])
GetSurfaceParms<-function(parameters,land){
	inc=res(land)[1]/1000 #get resolution in km
	km_len=dim(land)[1]*inc
	area=km_len^2
	
	#add to parameters list
	parameters=append(parameters,inc)
	names(parameters)[length(parameters)]<-"inc"
	
	parameters=append(parameters,km_len)
	names(parameters)[length(parameters)]<-"km_len"
	
	parameters=append(parameters,area)
	names(parameters)[length(parameters)]<-"area"
	
	return(parameters)
	
}
