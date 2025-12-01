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
## I switched to a character vector in Parameters.txt that defines which of the available variables (currently coded in _targets.R) should be used
# 	variable=parameters[which(parameters=="input")]
	## pull in parameters
	variable = parameters$input
	## get the parameters with more than one value
	variable_messy = parameters[which(lapply(parameters, length)>1)]
	## filter out parameters that are supposed to have more than one value (or that have values connected to values of other parameters)
	variable_messy = variable_messy[names(variable_messy) %in% c('out.opts', 'input',names(variable_messy)[grep('^B1',names(variable_messy))],'ss') == FALSE]
	## get all combinations
	temptab <- expand.grid(variable_messy)
	## build out table of parameters that are defined in sync, e.g. B1, density, ss with specific state
	canonical.params <- data.frame(state = rep(parameters$state, each=length(parameters$density)),
								density=rep(parameters$density, length(parameters$state)),
								ss=rep(parameters$ss, length(parameters$state)),
								B1=unlist(lapply(parameters$state, function(x){
								parameters[paste0('B1_',x)]}))
# 								,B1=rep(c(parameters$B1_FL, parameters$B1_SC), each=length(parameters$state)/length(parameters$density))
)
	## join user defined variables with canonical parameters
	common.columns = names(canonical.params)[names(canonical.params) %in% names(temptab)]
	result <- inner_join(temptab, canonical.params, by=common.columns) ## makes 'noise' without by= argument, but changes depending on user input
	## calculate B2 (relative to B1)
	result$B2 <- result$B1*parameters$B2_B1_factor

	shapes <- parameters[grep('^shape_',names(parameters))]
	rates <- parameters[grep('^rate_',names(parameters))]
	F1s <- parameters[grep('^F1_', names(parameters))]
	shaperate_table <- data.frame(state = unlist(lapply(strsplit(names(rates), '_'), function(x) unlist(x)[2])),
								  shape = unlist(shapes),
								  rate = unlist(rates),
								  F1 = unlist(F1s))

	result <- inner_join(result, shaperate_table, by='state')
# # 	for(i in 1:length(names(variable))){
# 	for(i in seq(variable)){
# # 		if(length(grep(names(variable)[i],names(inputs)))<1){
# 		if(length(grep(variable[i],names(inputs)))<1){
# 			stop(paste0("Missing parameter: ",names(variable)[i]))
# 		}
# 	}
#
# 	## filter only inputs defined in parameters file
# 	inputs <- inputs[variable]
#
# 	for(i in 1:(length(inputs)-1)){
# 		## I feel like there is an easier way to do this, but it works
#
# 		if(i==1){
# 			mat1=as.matrix(inputs[[i]])
# 		} else{
# 			mat1=result
# 		}
#
# 		mat2=as.matrix(inputs[[i+1]])
# 		indc=expand.grid(1:nrow(mat1),1:nrow(mat2))
# 		result <- cbind(mat1[indc[, 1], ], mat2[indc[, 2], ])
#
# 	}
# 	## check names for good measure (i.e. to include singletons like "Rad")
# 	colnames(result) <- unlist(lapply(strsplit(names(unlist(inputs, recursive=FALSE)), '.', fixed=TRUE), function(x) x[2]))
	return(result)
	
}
