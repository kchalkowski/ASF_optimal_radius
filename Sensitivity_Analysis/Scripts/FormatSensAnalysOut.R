
#############################
########## Purpose ##########
#############################

#The purpose of this script is to format output from 
#doSensitivityAnalysis.R and plot against observed Poland 
#ASF prevalence data

#input: sl.summaries, edat/wdat (formatted, de-identified, aggregated poland data)
#output: plot comparisons

##################################
########## Script Setup ##########
##################################

#set home dir where ASF_optimal_radius pipeline lives
home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius"
setwd(home)

#set run_folder where sensitivity analysis rds file saved
run_folder="TestRun_2024_09_13_08_44_35"

#read in sl.summary opts-- needs be same as set in doSensitivityAnalysis.R
sl.summary.opts=c("SEIRCZ_total",
                  "SEIRCZ_zone",
                  "SEIRCZ_zone_cells",
                  "SEIRCZ_total_apparent",
                  "SEIRCZ_zone_apparent",
                  "SEIRCZ_zone_cells_apparent")


#load libraries
library(hrbrthemes)
library(dplyr)
library(ggplot2)

#read in sl.summaries
sl.summaries=readRDS(paste0(home,"/Sensitivity_Analysis/Output/",run_folder,"/sl.summaries.rds"))
names(sl.summaries)=sl.summary.opts

#output datasets
SEIRCZ_total_apparent=sl.summaries$SEIRCZ_total_apparent
SEIRCZ_zone_apparent=sl.summaries$SEIRCZ_zone_apparent

#read ASF outbreak data
edat=readRDS(paste0("Sensitivity_Analysis/Input/datE80km_weekly.rds"))
wdat=readRDS(paste0("Sensitivity_Analysis/Input/datW80km_weekly.rds"))

edatz=readRDS(paste0("Sensitivity_Analysis/Input/E_wk.zone.summary.rds"))
wdatz=readRDS(paste0("Sensitivity_Analysis/Input/W_wk.zone.summary.rds"))

#######################################
########## Summary functions ##########
#######################################
{
#remove unneeded cols from obs
rm_cols_obs=function(dat){
  dat=dat[,-c(which(colnames(dat)=="outbreak"),
              which(colnames(dat)=="startdt"),
              which(colnames(dat)=="N"),
              which(colnames(dat)=="Npos"),
              which(colnames(dat)=="Ncells_Ipa"),
              which(colnames(dat)=="Ncells_Cpa"))]
  return(dat)
}


#use add_vars for both to add missing vars
add_vars <- function(df,vars,val){
  df[vars] <- val
  df
}

#function assumes edat/wdat don't have cols not in sim_outbreak
#sim_outbreak=SEIRCZ_zone_apparent
#edat=edatz
#wdat=wdatz
CombineSimObsOutbreak<-function(sim_outbreak,edat,wdat){
  #assumes that edat2 and wdat2 have same col names
  new_vars=colnames(sim_outbreak)[!(colnames(sim_outbreak)%in%colnames(edat))]
  
  edat=edat |>
    add_vars(df = _,vars = new_vars,val = NA)
  wdat=wdat |>
    add_vars(df = _,vars = new_vars,val = NA)
  
  
  l<-list(observed_east=edat,
          observed_west=wdat,
          simulated=sim_outbreak)
  
  #all(colnames(wdat)%in%colnames(sim_outbreak))
  #all(colnames(sim_outbreak)%in%colnames(wdat))
  
  out.breaks=do.call(rbind, lapply(l, function(x) x[match(names(l[[1]]), names(x))]))
  rownames(out.breaks)=NULL
  out.breaks<-as.data.frame(out.breaks)
  return(out.breaks)
}

#simdata=SEIRCZ_zone_cells_apparent
format_simdata<-function(simdata){
  #match colnames
  if(any(colnames(simdata)=="Ipa_med")){
  colnames(simdata)[
    c(which(colnames(simdata)=="Ipa_med"),
      which(colnames(simdata)=="timestep"),
      which(colnames(simdata)=="Cpa_med"))]<-
    c("Ipa","week","Cpa")
  } else{
    colnames(simdata)[
        which(colnames(simdata)=="timestep")]<-
      c("week")
  }
  
  
  
  simdata=as.data.frame(simdata)
  simdata$type="simulated"
  return(simdata)
}


PlotOutbreaks<-function(outbreaks,plot.opt){
  if(plot.opt=="live"|plot.opt=="carcass"){
    #need finish adding geom_ribbon to all funcs
    
    outsums=
      outbreaks %>% dplyr::group_by(type,week) %>%
      dplyr::summarize(Ipa_med=median(Ipa),
                       Ipa_min=min(Ipa),
                       Ipa_max=max(Ipa),
                       Cpa_med=median(Cpa),
                       Cpa_min=min(Cpa),
                       Cpa_max=max(Cpa),
                       Ipt_med=median(Ipt),
                       Ipt_min=min(Ipt),
                       Ipt_max=max(Ipt),
                       Cpt_med=median(Cpt),
                       Cpt_min=min(Cpt),
                       Cpt_max=max(Cpt)
      ) %>%
      as.data.frame()
    
    colnames(outsums)[c(which(colnames(outsums)=="Ipa_med"),
                        which(colnames(outsums)=="Cpa_med"))]<-c("Ipa","Cpa")
    outsums=outsums[outsums$week<=72,]
    if(plot.opt=="live"){
      p=ggplot(outsums,aes(x=week,y=Ipa,color=type))+
        geom_line()+
        #geom_line(outsums[outsums$type=="simulated",],mapping=aes(x=week,y=Ipt_med),color="black",inherit.aes=FALSE)+
        geom_ribbon(outsums,mapping=aes(x=week,ymin=Ipa_min,ymax=Ipa_max,fill=type),alpha=0.5,inherit.aes=FALSE)+
        theme_ipsum()
    }
    
    if(plot.opt=="carcass"){
      p=ggplot(outsums,aes(x=week,y=Cpa,color=type))+
        geom_line()+
        #geom_line(outsums[outsums$type=="simulated",],mapping=aes(x=week,y=Cpt_med),color="black",inherit.aes=FALSE)+
        geom_ribbon(outsums,mapping=aes(x=week,ymin=Cpa_min,ymax=Cpa_max,fill=type),alpha=0.5,inherit.aes=FALSE)+
        theme_ipsum()
    }
    
  }
  
  if(plot.opt=="area"){
    
  }
  
  return(p)
  
}


}

#############################################################
########## Match formatting and rbind sim/obs data ##########
#############################################################

{
edat<-rm_cols_obs(edat)
wdat<-rm_cols_obs(wdat)
edatz<-rm_cols_obs(edatz)
wdatz<-rm_cols_obs(wdatz)

SEIRCZ_total_apparent=format_simdata(SEIRCZ_total_apparent)
SEIRCZ_zone_apparent=format_simdata(SEIRCZ_zone_apparent)

outbreaks_full=CombineSimObsOutbreak(SEIRCZ_total_apparent,edat,wdat)
outbreaks_zone=CombineSimObsOutbreak(SEIRCZ_zone_apparent,edatz,wdatz)
}

################################
########## Make plots ##########
################################

PlotOutbreaks(outbreaks_full,"live")
PlotOutbreaks(outbreaks_full,"carcass")

PlotOutbreaks(outbreaks_zone,"live")
PlotOutbreaks(outbreaks_zone,"carcass")
