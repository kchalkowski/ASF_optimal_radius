##Current troubleshooting needed
#issue with chull function
#initialize sounders creates rows with 0 sounders
#need to add if statement to stop dead pigs from moving


#issue with fastmovement
#think problem is that it is assigning 0 to pop[,3]
#confirmed! did test and that does show up. and all NA pop[,9] 10 etc. are where pop[,3] is 0. 
#Just need to figure out why zero is sometimes chosen and make it not do that.