
#######################
#### Set Variables ####
#######################
#parameters that may need to vary

state=1 #state switch, FL or SC params. 1=FL, 2=SC
density=1.5 #density per X
detectday=73
ss=2 #sounder size


#########################
#### Toggle switches ####
#########################
#ML- use original grid made in ML by Kim-- is 80kmx80km, 400m resolution
#DIY- make grid in R with preset options using Make_Grid.R

#Grid settings
grid.type="DIY"
#len=200 #nrow/ncol of grid
len=200
inc=0.4 #resolution of grid (m)

#grid.type="ML"

#Surveillance settings
#DetP=NULL #DetP=NULL if no detection probability, means DetP=1
DetP=0.5 #otherwise, set detection probability-- of all tested pigs that are positive
#this probability will determine how many of those are actually detected
  #Later-- get this probability from sensitivity/specificity (whichever one relates to false negatives)
  #Later-- consider whether should incorporate false positives, or if this matters
  #Also-- for now, not incorporatig DetP into FirstDetect. Point of that is to initiate a detection, and response. 
    #probably better ways to look at effect of late response due to missed detections
  #also assuming DetP is same for live and dead pigs for now
#########################
####Load fixed parameters
#########################
#parameters that stay the same
grid.opts="homogenous"
Pir = 0.05 #; %proportion of individuals that recover instead of die
death=0.00639 # probability of natural death
alphaC = 1.1 #% scaling parameter on relationship of effort to capture success
loc=c("FL","SC") #state strings
Rad=0 #culling radius
thyme=72 #total number of time steps. is spelled 'thyme' bc 'time' is a function in stats namespace
intensity=0.05 #
cullstyle="startOUT" #option in original ASF opt rad 2022 model
inc=0.4 #width of cells in grid 
Intensity = 0.05 #proportion of population targeted for removal per day based on capabilities
num_inf_0=1 #how many pigs to infect starting off
mc_time=0.0027
Pcr=0.23
Pbd=0.0192 #From opt rad ms

#set output options
#options: "sounderlocs","idzone","alldetections","incidence"
#if none selected (i.e., out.opts=c()), just get standard outputs
out.opts=c("sounderlocs")

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



