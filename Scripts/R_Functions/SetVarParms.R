#Parameters is list of parameters
#inputs is list of inputs that will vary, must match names in parameters
#xFL=c(0.7515,0.3550)
#xSC=c(0.5657,1.9082)
#inputs=list(
#	"shape_rate"=data.frame(
#		"shape"=c(0.7515,0.5657),
#		"rate"=c(0.3550,1.9082)
#		),
#	"density_ss_B1_B2"=data.frame(
#		"density"=c(1.5,3,5),
#		"ss"=c(2,4,6),
#		"B1"=c(0.009,0.009,0.009),
#		"B2"=c(0.009*2,0.009*2,0.009*2)
#		))

SetVarParms<-function(parameters, inputs){

## this doesn't work -- no "input" parameter? commenting out lets the model run
# 	variable=parameters[which(parameters=="input")]
# 	for(i in 1:length(names(variable))){
# 	if(length(grep(names(variable)[i],names(inputs)))<1){
# 		stop(paste0("Missing parameter: ",names(variable)[i]))
# 		}
# 	}
	
	for(i in 1:(length(inputs)-1)){
		
		if(i==1){
		mat1=as.matrix(inputs[[i]])
		} else{
		mat1=result
		}
		
		mat2=as.matrix(inputs[[i+1]])
		indc=expand.grid(1:nrow(mat1),1:nrow(mat2))
		result <- cbind(mat1[indc[, 1], ], mat2[indc[, 2], ])
		
	}
	
	return(result)
	
}
