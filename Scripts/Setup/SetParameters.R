
#######################
#### Set Variables ####
#######################
#parameters that may need to vary

state=1 #state switch, FL or SC params. 1=FL, 2=SC
density=1.5 #density per X
detectday=60
ss=2 #sounder size

##############################
#### Grid option settings ####
##############################

#DIY or ML:
  #ML- use original grid made in ML by Kim-- is 80kmx80km, 400m resolution
  #DIY- make grid in R with preset options using Make_Grid.R
  #county - make grid that is proportional to county size and extends just beyond
  grid.type="County"
  len=200 #nrow/ncol of grid
  inc=0.4 #resolution of grid (km)
  #grid.type="ML"

#grid.opt: "homogeneous" or "random", optional input  
  #default is "homogeneous"
  #homogeneous- grid with 7 cols is created, no land class designations
  #random- random neutral landscape model created, with land class values draw from uniform distribution and rescaled from 0 to 1.

  grid.opt="homogeneous"
  #grid.opt="ras"

  #grid.opt="homogeneous"
  
#How to initialize the population
  #pop_init_grid_opts="homogeneous", random distribution
  #pop_init_grid_opts="heterogeneous", distribute pigs according to land class preference
    #if "heterogeneous", initialized grid MUST be something other than "homogeneous" or InitializeSounders will throw error.
  #pop_init_grid_opts="homogeneous"
  pop_init_grid_opts="homogeneous"
  
  # pop_init_type="init_pop"

##################################
#### Movement option settings ####
##################################
#movement preference switch
  #0=random movement preference (distance-only)
  #1=abundance-avoidance movement preference
  #2=rsf preference (only works when pop_init_grid_opts="heterogeneous" and grid_opts!="homogeneous")
mv_pref=0
  
######################################
#### Surveillance option settings ####
######################################
  
DetP=1 
  #probability of infected pigs sampled as positive
  #for now, is just used as an input in some of the sounderlocs output processing steps.
  #later, may be integrated into culling/surveillance process within simulation.

# sampling switch
#0=none
#1=sample (will only work if detectday=0 [need flag])
sample=1
  
# path to sampling design file
sample.design <- read.csv(paste0(home,"/Input/Bladen_county_updated.csv")) # maybe make this more generic so user doesn't have to put in path?
  
# path to county shape file
county_shapefile <- st_read(paste0(home,"/Input/pvs_batch_form_37/partnership_shapefiles_24v2_37017/PVS_24_v2_county_37017.shp"))

#########################
####Load fixed parameters
#########################
#parameters that stay the same

Pir = 0.05 #; %proportion of individuals that recover instead of die
death=0.00639 # probability of natural death
alphaC = 1.1 #% scaling parameter on relationship of effort to capture success
loc=c("FL","SC") #state strings

Rad=10 #culling radius
thyme=52 #total number of time steps. is spelled 'thyme' bc 'time' is a function in stats namespace
intensity=0.05 #
cullstyle="startOUT" #option in original ASF opt rad 2022 model

Intensity = 0.05 #proportion of population targeted for removal per day based on capabilities
num_inf_0=1 #how many pigs to infect starting off
mc_time=0.0027
Pcr=0.23 # something to do with the carcasses
Pbd=0.0192 #From opt rad ms

#set output options
#options: "sounderlocs","idzone","alldetections","incidence"
#if none selected (i.e., out.opts=c()), just get standard outputs
out.opts=c("alldetections")
#out.opts=c()

#sounderlocs:
#subset of pop matrix taken at each timestep: 
#rows where there is a pig in any state, number in each status, and the cell number (location)
#first col is timestep number
#this can get quite large-- not good practice to just merge across 500-1000 reps
#summary options for replicates: sum of num SEIRCZ across each timestep, each replicate (tracking population/state changes)

#idzone
#matrix with timestep in first column, second column is each cell within control zone at each timestep
#This can also get quite large. Mostly for plotting/visualization for checking individual replicates.
#No cross-replicate summary options at the moment.

#alldetections:
#matrix with each row timestep,
#col1 timestep, col2=type (1 is live, 0 is dead), col3 is number of live/dead detected at each timestep, col4 is locations of detections
#cross replicate summary-- just add rep column, rbind for each replicat
#for now only include number.. later can add opt to output locs if/when needed

#incidence:
#matrix with timestep, num exposed, num infected, and assoc. locs.
#can get large, cross-replicate summaries for now ignore locs, just summarise num exp/inf across timestep
#in future, will add num susc to this to calculate different measures of Rt.



