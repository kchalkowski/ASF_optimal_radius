#######################
####### Purpose #######
#######################

# Processes sampling file that is user inputted indicating
# collection date, longitude, latitude, quantity and acreage of sample effort
  # Currently, only option is to have landscape prediction for county of interest
  # i.e. grid.opt="ras"

#######################
######## Function #####
#######################

# Inputs:
  # sampling- file with sampling scheme
  # inc - grid cell resolution (not sure if this will be true with new MakeGrid implementation)

#Outputs:
#1-sample.design- dataframe containing sampling info to be used in Surveilance.R
  #col1-date
  #col2-lat
  #col3-long
  #col4-quantity
  #col5-acres
  #col6-FY start date
  #col7-week #


PrepSurveillance<-function(sampling){

  # Read sampling file
  sample.design <- read.csv("sampling_scheme.csv")

  # Aggregate by week and store weeks to be sampled in variable
  # Convert "collection_date" to Date format
  names(sample.design)[1] <- "dates"  # Rename the first column to 'dates'
  sample.design$dates <- as.Date(sample.design$dates, format = "%m/%d/%Y") # change to standard date format
  
  # Apply mutate on the entire dataframe, not just on a vector
  sample.design <- sample.design %>%
    mutate(# if collection date is in october or later, fiscal year starts in that year
      # can handle out of order dates from spreadsheet
      fiscal_year_start = as.Date(paste0(ifelse(month(dates) >= 10, year(dates), year(dates) - 1), "-10-01")), 
      # Which fiscal week are we in
      fiscal_week = as.integer(difftime(dates, fiscal_year_start, units = "weeks")) + 1,
      # Make week number stays between 1-52 (handle any overflow)
      fiscal_week = ifelse(fiscal_week < 1, fiscal_week + 52, fiscal_week),  # Handle weeks that are before the fiscal year start
      fiscal_week = ifelse(fiscal_week > 52, fiscal_week %% 52, fiscal_week)  # Ensure the week number stays between 1 and 52
    )
  
  # Decide which cells are sampling locations and add "1" to 8th column of grid
  names(sample.design)[2] <- "latitude"
  names(sample.design)[3] <- "longitude"
  
  sample_coords <- st_as_sf(sample.design, coords = c("longitude", "latitude"), crs=4269) # pull out x and y coords from sample.design, convert to sf object
  sample_sf_transformed <- st_transform(sample_coords, crs = 26917) # transform coordinates to meter based CRS
  sample_coords_transformed <- st_coordinates(sample_sf_transformed) # extract the coordinates only
  sample_coords_transformed <- sample_coords_transformed / scale_factor
      
  return(sample.design)
}