# Creating metrics of fire regimes to delineate Australia's pyroregions

# This script creates the metrics of fire and fuel used to delineate Australia's pyroregions.
# Because file sizes of the satellite datasets are large, we have not provided them on Github, but they will be archived on a public 
# repository upon publication. 

# The output of this script is a stack of rasters summarising spatiotemporal metrics of Australia's fire regimes, which can
# be used in the clustering script to delineate Australia's pyroregions.

# This script is broken into the following sections:
# 1. Organise hotspot data
# 2. Calculate hotspot density
# 3. Timing of peak hotspots
# 4. Fire radiative power
# 5. Patch size
# 6. Length of fire season
# 7. Time since fire metrics (pyrodiversity and contagion)
# 8. Fire emissions
# 9. AVHRR Burned area
# 10. Prop of years with fire
# 11. Stack fire metrics and export
# 12. Fuel proxies: NPP and NDVI

rm(list = ls())

# packages
pacman::p_load(tidyr, tidyverse, sf, stringr, parallel, terra, tidyterra, 
               fasterize, landscapemetrics, exactextractr, leaflet, vegan, ncdf4)

# set working directory
setwd("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Analysis")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 1. ORGANISE HOTSPOT DATA ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# load csv of fire hot spots & format
# dates and times are in UTC
hot <- read.csv("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/2022 - MODIS HOTSPOTS/modis_hotspots.csv") %>%
  # reformat times/dates and add related columns
  mutate(HHMM = as.character(HHMM), YYYYMMDD = as.character(YYYYMMDD)) %>%
  mutate(HHMM = ifelse(nchar(HHMM) == 2, paste("00", HHMM, sep = ""), HHMM)) %>% # e.g. for cases where 00:12 was just shown as "12"
  mutate(HHMM = ifelse(nchar(HHMM) == 3, paste("0", HHMM, sep = ""), HHMM)) %>% # e.g. for cases where 03:12 was just shown as "312"
  mutate(date_time_UTC = paste(str_sub(.$YYYYMMDD, 1, 4), "-", str_sub(.$YYYYMMDD, 5, 6), "-", str_sub(.$YYYYMMDD, 7, 8), " ",
                                                str_sub(.$YYYYMMDD, 1, 2), ":", str_sub(.$YYYYMMDD, 3, 4), sep = "")) %>%
  mutate(date_time_UTC = as.POSIXct(date_time_UTC, tz = "GMT"),
         year = year(date_time_UTC), doy = as.numeric(strftime(date_time_UTC, format = "%j")), month = month(date_time_UTC), week = as.numeric(strftime(date_time_UTC,format="%W")))
         
# turn spatial and omit Nov and Dec of 2000, so that each month has been surveyed the same number of times
hot1 <- hot %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) %>% 
  filter(!(month %in% c(11, 12) & year == 2000))
rm(hot)

# load raster to use as template for calculating fire metrics; ERA = 0.25 deg
aus_land_area <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Grid templates/aus_land_area.tif")
aus_land_binary <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Grid templates/aus_land_binary.tif")
aus_land_prop <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Grid templates/aus_land_prop.tif")

# Extract the cell IDs and area for each point location
area <- extract(aus_land_area, hot1, ID = F, cells = T, xy = T) 
hot1 <- hot1 %>% mutate(
  cellID = area$cell,
  cell_x = area$x,
  cell_y = area$y,
  cellArea = area$area)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 2. HOTSPOT DENSITY ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# temporally static day and night hotspot density
hsDensity_day <- (setNames(rasterize(hot1 %>% filter(dn == "D"), aus_land_binary, fun=sum), "hsDensity_day")) / aus_land_area # divide by land area of cell
hsDensity_day[is.na(hsDensity_day)] <- 0 # add in zeroes for cells where no hotspots were observed
hsDensity_day <- hsDensity_day + aus_land_binary # add NAs back in for ocean (aus_land_binary has land as 0 and ocean as NA)

hsDensity_night <- (setNames(rasterize(hot1 %>% filter(dn == "N"), aus_land_binary, fun=sum), "hsDensity_night")) / aus_land_area # divide by land area of cell
hsDensity_night[is.na(hsDensity_night)] <- 0 # add in zeroes for cells where no hotspots were observed
hsDensity_night <- hsDensity_night + aus_land_binary # add NAs back in for ocean (aus_land_binary has land as 0 and ocean as NA)

# plot
plot(c(hsDensity_day,hsDensity_night))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 3. TIMING OF PEAK HOTSPOT COUNTS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Calculate week of peak hotspot. 
# One aim is to distinguish the northern and southern fire seasons.
# However, because week of year is a circular variables, and because cluster algorithms cannot perform circular analysis, 
# we instead linearise week of year into "weeks from Sept", or summer-ness and winter-ness. Under this approach, 
# spring and autumn fires would not be distunguished. However, we will attempt to capture those differences with a different variable.

# first, loop through the year and visually identify the period when the difference between north and south is greatest 
weeks <- seq(3,52, by = 3)
weekList <- list()
for(i in weeks) {
  print(i)
  referenceWeek = i
  peakHotspotWeek <- hot1 %>% data.frame() %>% 
    group_by(cellID, week, cell_x, cell_y) %>% 
    summarise(hotspotCount = n()) %>%
    group_by(cellID) %>%
    filter(hotspotCount == max(hotspotCount)) %>%
    # clunky way of calculating number of weeks (plus or minus) to the reference week
    mutate(weeksDiff1 = abs(week - referenceWeek),
           weeksDiff2 = abs(week - referenceWeek - 52),
           weeksDiff3 = abs(week + 52 - referenceWeek)) 
  peakHotspotWeek$peakFire_weeksFrom <- unlist(pmap(peakHotspotWeek[,c("weeksDiff1", "weeksDiff2", "weeksDiff3")], min))
  
  peakHotspot_weeksFromReferenceWeek <- peakHotspotWeek %>% 
    filter(!is.na(cell_x)) %>%
    st_as_sf(coords = c("cell_x", "cell_y"), crs = st_crs(4326)) %>%
    rasterize(., aus_land_binary,  field = "peakFire_weeksFrom", fun = mean) 
  names(peakHotspot_weeksFromReferenceWeek) <- paste("peakFire_weeksFrom", i, sep = "")
  weekList[[i]] <- peakHotspot_weeksFromReferenceWeek
}
weekRasts <- rast(weekList)
plot(rast(weekList))

# north-south pattern is most distinct in about week 36, so moving forward with week 36 (Apprpox Sept 3, varies by up to a few days each year)
plot(weekRasts$peakFire_weeksFrom36)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 4. FIRE RADIATIVE POWER ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

FRP_95 <- rasterize(hot1, aus_land_binary, field = "FRP", fun = function(i) quantile(i, probs = 0.95, na.rm = T)); names(FRP_95) <- "FRP_95"
plot(FRP_95)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 5. PATCH SIZE ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Metrics of fire size/shape from MODIS burnt area POLYGONS. I separtely rasterized the burn area polygon files (see "MODIS burn area.R"), which 
# produced a raster where pixels within the daily polygons were given the value corresponding to the area of the polygon.

# list and load all rasters (of MODIS burnt area polygons)
rastList <- list.files("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/2022 - MODIS Burnt Area/Rasters/burntPatchArea/", pattern = "tif")
burnAreaStack_0.01 <- rast(paste("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/2022 - MODIS Burnt Area/Rasters/burntPatchArea/", rastList, sep = ""))
names(burnAreaStack_0.01) <- paste("x", substr(rastList, 10, 16), sep = "")

# find maximum daily burn area in km2 for reporting in text
max(minmax(burnAreaStack_0.01))/1000000 # km2
max(minmax(burnAreaStack_0.01))/10000 # ha
min(minmax(burnAreaStack_0.01))/1000000

# coarsen to same scale (0.25 deg) as other products
burnAreaStack_0.25 <- aggregate(burnAreaStack_0.01, fact = 25, na.rm = T)

# values are in m2; convert to km2
burnAreaStack_0.25_km2 <- burnAreaStack_0.25/1000000

# calculate mean burn size for each pixel
# because NA's are removed, the inference is on the fires that actually burn; i.e. if a fire burns, it is likely to be of these sizes/variability
medianPatchSize <- median(burnAreaStack_0.25_km2, na.rm = T); names(medianPatchSize) <- "medianPatchSize"
maxPatchSize <- max(burnAreaStack_0.25_km2, na.rm = T); names(maxPatchSize) <- "maxPatchSize"

# plot
plot(log(c(medianPatchSize, maxPatchSize)))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 6. LENGTH OF FIRE SEASON ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# for season length, we need a slightly coarser scale to reduce stochasticity; thus, aggregating raster and extracting new cell ID
seasonLengthTemplate <- aus_land_area %>% aggregate(fact = 2)
cellID_sLength <- extract(seasonLengthTemplate, hot1, ID = F, cells = T, xy = T)
hot1$cellID_sLength <- cellID_sLength$cell
hot1$cell_x_sLength <- cellID_sLength$x
hot1$cell_y_sLength <- cellID_sLength$y

# create row for every day for each cell, which will be used for days without any hotspots recorded
zeroes <- hot1 %>% 
  data.frame() %>% 
  dplyr::select(cellID_sLength, cell_x_sLength, cell_y_sLength) %>% 
  unique() %>% 
  slice(rep(1:n(), each = 366)) %>%
  mutate(doy = rep(1:366, times = length(unique(cellID_sLength))), hs_count = 0)

# calculate daily hotspots per cell (including zeroes)
dailyHot <- hot1 %>% data.frame() %>%
  group_by(cellID_sLength, cell_x_sLength, cell_y_sLength, doy) %>%
  summarise(hs_count = n()) %>%
  # add in zeroes
  bind_rows(zeroes) %>%
  # count total hotspots including zeroes
  group_by(cellID_sLength, cell_x_sLength, cell_y_sLength, doy) %>%
  summarise(hotspotCount = sum(hs_count)) %>%
  filter(!is.na(cell_y_sLength)) %>% 
  # turn into sf object
  st_as_sf(coords = c("cell_x_sLength", "cell_y_sLength"), crs = st_crs(4326), remove = F) %>%
  ungroup()

# Function to loop through different start dates of "fire year" to find the shortest continuous period that contains 
# "fireQuantity" proportion of all historical fire activity for each cell. Option to set different "dayIncrement" of days. 
# The function takes an sf data frame of daily hot spot counts for each cell, and returns an sf object for the minimum
# continuous number of days required for each cell to reach "fireQuantity" proportion of the year's fire activity.
findSeasonLength <- function(dailyHotspotCount, dayIncrement, fireQuantity, minHotspots){
  
  # initialise list to store results in 
  seasonLengthList <- list()
  
  # loop through different start days for fire year and calculate fire day of year (fire_doy) relative to the start of the fire year (doy_season_start)
  for(i in seq(1, 365, by = dayIncrement)) {
    print(i)
    
    # set start date of fire year
    dailyHot1 <- dailyHotspotCount %>%
      mutate(doy_season_start = i,
             days_to_end_year = 365 - doy_season_start, 
             fire_doy = doy - doy_season_start,
             fire_doy = ifelse(fire_doy < 0, doy + days_to_end_year, fire_doy)) 
    
    # calculate cumulative fire sum/proportion as the fire year progresses
    cumulativeFire <- dailyHot1 %>% 
      arrange(cellID_sLength, fire_doy) %>%
      group_by(cellID_sLength, cell_x_sLength, cell_y_sLength) %>%
      mutate(cumulativeFire = cumsum(hotspotCount),
             cumulativePropFire = cumulativeFire/max(cumulativeFire),
             totFire = sum(hotspotCount))
    
    # Find the day when cumulative proportion first exceeds thresholds determined by "fireQuantity"
    fireSeason <- cumulativeFire %>% 
      # make sure each cell is arranged from start of year to end of year, and then group by cellID
      arrange(cellID_sLength, fire_doy) %>%
      group_by(cellID_sLength) %>%
      # exclude cells without at least minHotspots; While this threshold is arbitrary, the reason it is needed is because having few hotpots means it's 
      # difficult/impossible to calculate cumulative percentiles reliably -> essentially too little data
      filter(totFire > minHotspots) %>%
      summarise(day_exceed_05 = min(which(cumulativePropFire > (1 - fireQuantity)/2 )),
                day_exceed_95 = min(which(cumulativePropFire > (1 - fireQuantity)/2 + fireQuantity ))) %>%
      mutate(fireSeasonLength = day_exceed_95 - day_exceed_05,
             doy_season_start = i) 
    
    # store results in list
    seasonLengthList[[paste("yrStart_", i, sep = "")]] <- fireSeason
  } 
  
  # bind all runs together, and find the minimum season length to reach the threshold of cumulative fire activity
  seasonLength_ <- seasonLengthList %>% 
    bind_rows() %>% 
    group_by(cellID_sLength) %>%
    summarise(minSeasonLength = min(fireSeasonLength))
  
  # returning a spatial data frame with one row for each pixel, containing the minimum season length
  return(seasonLength_)
}

# loop through every day of year as candidate starting point
seasonLength <- findSeasonLength(dailyHotspotCount = dailyHot, dayIncrement = 1, fireQuantity = 0.9, minHotspots = 50)

# rasterize and resample to final resolution to match other layers
fireSeasonLength <- rasterize(seasonLength, seasonLengthTemplate,  field = "minSeasonLength") %>% 
  resample(aus_land_area) 
fireSeasonLength <- fireSeasonLength + aus_land_binary # zero out ocean
names(fireSeasonLength) <- "fireSeasonLength"

plot(fireSeasonLength)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 7. TIME SINCE FIRE METRICS (PYRODIVERISTY AND CONTAGION) ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Calculate richness of ages
count_ages <- function(x) {
  diversity_count <- length(unique(x))
  return(diversity_count)
}

tsf <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/2022 - MODIS Burnt Area/Derived summaries/timesSinceFire_0_005.tif")
tsf[tsf == 0] <- NA  
ages_richness <- tsf %>% 
  aggregate(fact = 50, fun = count_ages) %>%
  project(aus_land_binary) 
ages_richness <- ages_richness + aus_land_binary # NA out ocean
names(ages_richness) <- "timeSinceFire_richness"
plot(ages_richness, colNA = "blue")

# Calculate contagion
calculate_contagion <- function(x) {
  contag <- lsm_l_contag(x)
  return(contag$value)
}

contag <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/2022 - MODIS Burnt Area/Derived summaries/timesSinceFire_0_005.tif") %>%
  project("epsg:4326") %>%
  focal(w = 9, FUN = calculate_contagion, na.rm = T) %>% 
  aggregate(fact = 50, na.rm = T) %>%
  resample(aus_land_binary)
names(contag) <- "timeSinceFire_contagion"
plot(contag)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 8. FIRE EMISSIONS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# emissions data from here: https://gmd.copernicus.org/articles/15/8411/2022/
emmissions_stack <- list()
netCDFList <- list.files("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Fire emissions/", pattern = ".nc")
for(i in netCDFList) {
  print(i)
  # load netcdf 
  fileName <- paste("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Fire emissions/", i, sep = "")
  our_nc_data <- nc_open(fileName)
  
  # extract lat/long/time and CO2 emissions for above and below ground
  lat <- ncvar_get(our_nc_data, "MOD_CMG025/lat")
  lon <- ncvar_get(our_nc_data, "MOD_CMG025/lon")
  time <- ncvar_get(our_nc_data, "MOD_CMG025/time")
  tunits <- ncatt_get(our_nc_data, "MOD_CMG025/time", "units") 
  aboveGround <- ncvar_get(our_nc_data, "MOD_CMG025/emissions/C_AG_TOT") 
  belowGround <- ncvar_get(our_nc_data, "MOD_CMG025/emissions/C_BG_TOT") 
  
  # turn into a matrix
  lonlattime <- as.matrix(expand.grid(lon,lat,time)) 
  total_emissions <- as.vector(aboveGround) + as.vector(belowGround)

  # organise and rasterize
  emiss <- data.frame(cbind(lonlattime, total_emissions)) %>%
    rename(Long = Var1, Lat = Var2, Date = Var3) %>%
    filter(!is.na(total_emissions),
           Lat > -44.125 & Lat < -9.875,
           Long > 109.875 & Long < 155.125) %>%
    st_as_sf(coords = c("Long", "Lat"), crs = st_crs(4326)) %>%
    rasterize(aus_land_area, field = "total_emissions", fun = mean, na.rm = T)
  
  emmissions_stack[[paste("yr_", substr(i, 28, 31), sep = "")]] <- emiss
}
emmissions_stack <- rast(emmissions_stack) 
plot(emmissions_stack[[1:19]], colNA = "blue")
meanMonthlyEmissions <- mean(emmissions_stack) / aus_land_prop
names(meanMonthlyEmissions) <- "meanMonthlyEmissions"

plot(log1p(meanMonthlyEmissions), colNA = "blue")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 9. AVHRR BURNED AREA ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

load_avhrr_ba <- function(parentDir = "C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/AVHRR Burned Area Global 1982-2018/") {
  folderList <- list.files(parentDir)
  folderList <- folderList[grep(pattern = "txt", folderList, invert = T)] # do not include text files
  yearlyRastList <- list()
  # loop through each yearly folder
  for(i in folderList) {
    print(i)
    
    # find file names in each yearly folder
    fileList <- list.files(paste(parentDir, i, "/", sep = ""))
    
    yearlyRastStack <- rast()
    # load each file
    for(j in 1:length(fileList)) {
      # load each file, and select bured area layer
      r <- rast(paste(parentDir, i, fileList[[j]], sep = "/")) %>%
        dplyr::select(burned_area)    
      yearlyRastStack <- c(yearlyRastStack, r)
    }
    
    yearlyRastList[[i]] <- yearlyRastStack
  }
  return(yearlyRastList)
}

# load rasters. Each year is an element of this list. Within each year, there are 12 monthly rasters for area burned.
avhrr_ba <- load_avhrr_ba()

# calculate annual total for each year
annual_ba <- rast()
for(i in names(avhrr_ba)) {
  print(i)
  annual_ba[[i]] <- sum(avhrr_ba[[i]])
}

# global
plot(mean(annual_ba))

meanGlobal <- mean(annual_ba)
cSizes <- cellSize(meanGlobal, unit="m")
plot(cSizes)

# crop for Aus and resample to match grid
annual_ba_aus <- annual_ba %>% crop(ext(aus_land_area)) %>% resample(aus_land_area)
annual_ba_aus <- annual_ba_aus + aus_land_binary
annual_ba_aus <- annual_ba_aus / 1000000 # get from m2 to km2
plot(mean(annual_ba_aus), colNA = "blue") 

# summarise rasters
avhrr_ba_median <- median(annual_ba_aus); names(avhrr_ba_median) <- "avhrr_ba_median"
avhrr_ba_sd <- stdev(annual_ba_aus)
avhrr_ba_cv <- avhrr_ba_sd / avhrr_ba_mean; names(avhrr_ba_cv) <- "avhrr_ba_cv"
plot(avhrr_ba_median)
plot(avhrr_ba_cv)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### 10. PROP OF YEARS WITH FIRE ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# disaggreate grid to a finer cell size, and extract the cell IDs and area for each point location
fine_disag_fact <- 10
fine_grid <- aus_land_area %>% disagg(fact = fine_disag_fact) %>% mutate(area = area/fine_disag_fact)
fine_grid$zeroes <- fine_grid; fine_grid$zeroes[fine_grid$zeroes > 0] <- 0
area_fine <- extract(fine_grid, hot1, ID = F, cells = T, xy = T) 

hot1$cellID_fine <-  area_fine$cell
hot1$cell_x_fine <- area_fine$x  
hot1$cell_y_fine <- area_fine$y
hot1$cellArea_fine <- area_fine$area

yrs_with_fire <- hot1 %>% 
  data.frame() %>% 
  group_by(cellID_fine, year, cell_x_fine, cell_y_fine, dn) %>% 
  summarise(hotspotBinary = ifelse(n() > 0, 1, 0)) %>%
  group_by(cellID_fine, cell_x_fine, cell_y_fine, dn) %>% 
  summarise(yrs = sum(hotspotBinary)) %>% 
  filter(!is.na(cell_x_fine)) 

# create separate layers for day and night
yrs_with_fire_list <- list()
for(i in unique(yrs_with_fire$dn)) {
  yrs_with_fire_list[[i]] <- yrs_with_fire %>%
    filter(dn == i) %>%
    st_as_sf(coords = c("cell_x_fine", "cell_y_fine"), crs = st_crs(4326)) %>%
    rasterize(., fine_grid$area,  field = "yrs") / fine_grid$area * max(fine_grid$area) # divide by area of cell to account for different opportunity to detect a fire, then multiply by max grid size so that inference on "yrswithfire" is calculated per grid cell
}

yrs_with_fire_rast <- rast(yrs_with_fire_list) 
yrs_with_fire_rast[is.na(yrs_with_fire_rast)] <- 0 # adding zeroes in place of NAs, because NAs over Aus land represent cells where no fire has been observed
yrs_with_fire_rast <- yrs_with_fire_rast + fine_grid$zeroes # inserting NAs for ocean
yrs_with_fire_rast <- aggregate(yrs_with_fire_rast, fact = fine_disag_fact, na.rm = T)
names(yrs_with_fire_rast) <- paste("yrsWithFire", names(yrs_with_fire_rast), sep = "_")
plot(yrs_with_fire_rast)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 11. STACK FIRE RASTERS AND EXPORT ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

fireMetricStack <- rast(list(hsDensity_day, yrs_with_fire_rast$yrsWithFire_D, yrs_with_fire_rast$yrsWithFire_N, FRP_95, medianPatchSize, maxPatchSize, 
                             ages_richness$timeSinceFire_richness,contag$timeSinceFire_contagion,weekRasts$peakFire_weeksFrom36,fireSeasonLength,
                             meanMonthlyEmissions, avhrr_ba_median, avhrr_ba_cv))
fireMetricStack <- (fireMetricStack + aus_land_binary) # ensure ocean is zeroed out for layers calculated at coarser cell size
plot(fireMetricStack)

# some key variables have a few holes in otherwise contiguous areas with information. Thus, borrowing information by smoothing from surrounding cells; 
# if a cell with NA is adjacent to one or more cells that are not NAs, borrow that data by taking the mean within the focal window.
# na.policy = "only" will not change non-NA values
fireMetricStack$medianPatchSize <- focal(fireMetricStack$medianPatchSize, w = 3, fun = mean, na.rm = T, na.policy = "only")
fireMetricStack$maxPatchSize <- focal(fireMetricStack$maxPatchSize, w = 3, fun = mean, na.rm = T, na.policy = "only")
fireMetricStack$timeSinceFire_contagion <- focal(fireMetricStack$timeSinceFire_contagion, w = 3, fun = mean, na.rm = T, na.policy = "only")
fireMetricStack$peakFire_weeksFrom36 <- focal(fireMetricStack$peakFire_weeksFrom36, w = 3, fun = mean, na.rm = T, na.policy = "only")
fireMetricStack$fireSeasonLength <- focal(fireMetricStack$fireSeasonLength, w = 3, fun = mean, na.rm = T, na.policy = "only")
fireMetricStack <- fireMetricStack + aus_land_binary # ensure ocean is zeroed out for layers calculated at coarser cell size
plot(fireMetricStack, colNA = "blue")

# writeRaster(fireMetricStack, "fireMetricStack_20230821.tiff", overwrite = T)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 12. FUEL PROXIES: NPP & NDVI ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Calculate mean NPP from annual rasters of MODIS NPP 
nppFiles <- list.files('C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/MODIS_NPP/Net_PP_Yearly_500m_v6/Npp', pattern = "tif", full.names = T)
npp_yearly <- rast(nppFiles) # load one rast per year
npp_yearly[npp_yearly > 30000] <- NA # water was stored as large values, replacing with NA
npp_mean <- mean(npp_yearly, na.rm = T); plot(npp_mean) # take mean of all years
# save raster
#writeRaster(npp_mean, "npp_mean.tif")

# Calculate NDVI difference bwteen summer and winter months
fileList <- list.files("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/MODIS_Veg Indices/VI_Monthly_005dg_v6/NDVI")
r <- rast(paste("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/MODIS_Veg Indices/VI_Monthly_005dg_v6/NDVI",
                fileList, sep = "/"))
plot(r)
NDVI_diff <- (r$MOD13C2_NDVI_2001_032 - r$MOD13C2_NDVI_2001_213) %>% 
  resample(fireMetricStack1) %>%
  focal(w = 3, na.policy = "only", fun = "mean") # interpolating a handful of interior NAs
NDVI_diff <- NDVI_diff + aus_land_binary
names(NDVI_diff) <- "NDVI_diff"
plot(NDVI_diff)
# save raster
#writeRaster(NDVI_diff, "NDVI_diff.tif")
