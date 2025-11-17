#input path to parameters_txt file
#parameters_txt="Parameters.txt"

#Future updates:
  #1-haven't tested some formats not currently used, i.e. nv_vars, may break
  #2-need input switch for manual inputs for B1, ss
  #3-stop check for nonsensical combinations of parameters
  #4-stop check if needed parameters missing
  #5-stop check if needed parameters are incorrectly formatted

FormatSetParameters<-function(parameters_txt){
  
  #read text file
  raw=readLines("Parameters.txt")
  
  #remove all but parameters set
  parms=raw[grep("^~",raw)]
  
  #trim ~
  parms=str_sub(parms,3L)
  
  #remove all whitespace
  parms=gsub("\\s","",parms,perl=TRUE)
  
  #Get values
    #get vars using lookahead perl regex
    vars=gsub("(.*)(?=\\=)","",parms,perl=TRUE)
    #trim =
    vars=str_sub(vars,2L)
  
  #Get names
    #get names using lookahead perl regex
    pnames=gsub("(?<=\\=)(.*)","",parms,perl=TRUE)
    #trim =
    pnames=str_sub(pnames,1L,-2L)
    
  #sep character and numeric
    #character
    cvars=vars[grep("\\\"",vars)]
    cvars=cvars[-grep("c\\(\\\"",cvars)]
    cvarnames=pnames[grep("\\\"",vars)]
    cvarnames=cvarnames[grep("\\\"",cvars)]
    
    #character vector
    cv_vars=vars[grep("c\\(\\\"",vars)]
    cv_varnames=pnames[grep("c\\(\\\"",vars)]
    
    #numeric
    nvars0=vars[-grep("\\\"",vars)]
    nvars = nvars0[-grep("c\\(\\d", nvars0)] ## this was grabbing the numeric vectors as well
    nvarnames=pnames[-grep("\\\"",vars)]
    nvarnames = nvarnames[-grep("c\\(\\d", nvars0)] ## need a temporary object to transfer index

    #numeric vector
    nv_vars=vars[grep("c\\(\\d",vars)]
    nv_varnames=pnames[grep("c\\(\\d",vars)]
    
  #handle characters, remove extra \"
    cvars=gsub("\\\"","",cvars)
    
  #handle numeric
    nvars=as.numeric(nvars)
    
  #handle vectors
    if(length(cv_vars) > 1){
      cv_vars=lapply(cv_vars[seq(cv_vars)], function(x) eval(parse(text=x))) # handles more than one character vector (i.e. if there are multiple out.opts AND multiple "input" names for variables)
    } else {
      cv_vars <- eval(parse(text=cv_vars)) # also handles single and missing vectors (as before; preserves the functionality of "is.null" statements for allvars below)
    }
    if(length(nv_vars) > 1){
      nv_vars=lapply(nv_vars[seq(nv_vars)], function(x) eval(parse(text=x)))
    } else {
      nv_vars <- eval(parse(text=nv_vars))
    }

  allvars=c("cvars","cv_vars","nvars","nv_vars")
  allvars=allvars[which(!(c(is.null(cvars),
                    is.null(cv_vars),
                    is.null(nvars),
                    is.null(nv_vars))))]

  #initiate output list
  list.all=vector(mode="list",length=0)
  
  #go through existing parms, append to list
    if("cvars"%in%allvars){
      templist=as.list(cvars)
      names(templist)=cvarnames
      list.all=append(list.all,templist)
    }
  
    if("cv_vars"%in%allvars){
      templist=as.list(cv_vars)
      names(templist)=cv_varnames
      list.all=append(list.all,templist)
    }
  
    if("nvars"%in%allvars){
      templist=as.list(nvars)
      names(templist)=nvarnames
      list.all=append(list.all,templist)
    }
  
    if("nv_vars"%in%allvars){
      templist=as.list(nv_vars)
      names(templist)=nv_varnames
      list.all=append(list.all,templist)
    }
  
  #Now, use vals in list.all to auto-set certain parameters
  #dens=list.all$density
    #ss
    #B2?? automatically calculated as B1*0.5
  
  #if(dens==1.5){
  #  ss=2
  #  B1=0.9
  #}
  #if(dens==3){
  #  ss=4
  #  B1=0.4
  #}
  #if(dens==5){
  #  ss=2
  #  B1=0.2
  #}
  
  #B2=B1*0.5
  
  #add addl calc'd parameters to list
  #templist=list(ss,B1,B2)
  #names(templist)=c("ss","B1","B2")
  #list.all=append(list.all,templist)
  
  return(list.all)
    
}
