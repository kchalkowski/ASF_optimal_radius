###Notes
#As model is currently, dead pigs move!
#Need to change this
#before doing move stuff, detect any dead pigs in sounders
#then, move this pig into its own row 
#put if statement for only rows with live pigs

#Need way to remove all pigs infected, remove the row don't just subtract one from infected column
#but, still need record of that pig

#make sure removing rows and not just individuals when doing the births and deaths!!

###Outline of ASF Model in MatLab

#OneRunWeekly.m sets the variables and then runs the main function, SimulateOneRun

#SimulateOneRun activities
	#initializes state variables/agents (try alternative method.. data table with IDs for location in grid
	#for 2-time loop simulates the population over time
		#if function to stop early if no more I,C,E
			#Processes: 
			#Movement 
				#movement function (RCPP)
				#assigns new locs to matrices (RCPP)
				#then plots (create output to plot later?)
			#Births and Deaths (RCPP)
			#Epidemiological State Transitions
				#Calculate force of infection (RCPP)
				#transitions: e to i, i to r, etc
				#Carcass decay (natural removal from landscape)
				#update states based on demographic/epi processes
			#Response: Culling Zone
				#when i=detectday, make initial detection, start testing all culled individuals
				#find new positives (uniques)
				#CullingOneRun function
				#update states based on mgmt process
				#initiate response on day of first detection
					#get all positive samples
					#update surveillance data
					#remove detected case
			#Track true spatial spread (RCPP)
				#area of infection function (RCPP)
				#plot infected cases on top
			#Get outputs
				#incidence sums, etc etc..

#R ASF Model outline....
#InitializeModel.R sets all the parameters, initializes agents, etc.

#To do list
#Generate grid function in R


