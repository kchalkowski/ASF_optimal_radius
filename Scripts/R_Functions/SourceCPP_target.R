cpp_path=tar_read(Fast_FOI_Matrix_script)
sourceCPP_target<-function(cpp_path){
	#readLines(cpp_path)
	#readLines(cpp_path) %>% paste0(collapse="\n") %>% cat
	Rcpp::sourceCpp(cpp_path)
	}