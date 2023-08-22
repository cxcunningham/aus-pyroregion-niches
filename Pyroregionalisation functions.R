# CUSTOM FUNCTIONS FOR PYROREGIONALISATION ANALYSIS

# function for plotting rasters
raster_plot_function <- function(rasterForPlotting, varName, title, legendName, trans = "none", breaks_ = c(1), categorical = F) {
  p <- ggplot() +
    theme_map() +
    geom_spatraster(data = rasterForPlotting, aes(fill = .data[[varName]]))  +
    scale_fill_viridis(name = legendName, na.value="grey100") +
    theme(axis.text = element_blank()) +
    labs(title = title) + 
    theme(legend.position = c(0.2,0.05), legend.direction = "horizontal", plot.title = element_text(hjust = 0.5, size = 10))
  
  if(trans != "none") 
    p <- p + scale_fill_viridis(trans = trans, name = legendName, breaks = breaks_, na.value="grey100")
  
  if(categorical == T)
    p <- p + scale_fill_viridis_d(name = legendName, na.value="grey100")
  
  plot(p)
  return(p)
}

# plot layers for multipanel plot
raster_multilayer_plot_function <- function(rasterForPlotting = fireMetricStack_trans, varName, 
                                            trans = NA, begin = 0,
                                            rounding = 0.01, n.breaks_ = 4, breaks_ = NA, legendName = element_blank()) {
  p <- ggplot() +
    theme_map() +
    geom_spatraster(data = rasterForPlotting, aes(fill = .data[[varName]]))  +
    scale_fill_continuous_sequential(palette = "Inferno", begin = begin, rev = FALSE, na.value="grey100", name = element_blank(), n.breaks = 3) +
    theme(axis.text = element_blank()) +
    labs(title = varName) + 
    theme(legend.position = c(0.1,0.05), legend.direction = "horizontal", plot.title = element_text(hjust = 0.5, size = 14),
          legend.key.size = unit(0.4, 'cm'))  +
    geom_sf(data = aus_states, fill = NA, colour = "white", linewidth = 0.4) +
    lims(x = c(113.375, 153.625), y = c(-43.625, -10.875))
  
  if(!is.na(trans)) 
    p <- p + 
      scale_fill_continuous_sequential(trans = trans, palette = "Inferno", begin = begin, rev = FALSE, na.value="grey100", name = element_blank(), n.breaks = 3) 
  
  if (!is.na(breaks_[1])) {
    if (is.expression(legendName)) {
      p <- p + 
        scale_fill_continuous_sequential(trans = trans, palette = "Inferno", begin = begin, rev = FALSE, na.value = "grey100",
                                         name = legendName, breaks = breaks_)
    } else {
      p <- p + 
        scale_fill_continuous_sequential(trans = trans, palette = "Inferno", begin = begin, rev = FALSE, na.value = "grey100",
                                         name = legendName, breaks = breaks_, labels = function(x) format(breaks_, scientific = F))
    }
  }
  plot(p)
  return(p)
}

# function that automatically chooses the best method of normalising data (attempting to make Gaussian) and then scales from 0-1
normalise_rasters <- function(rasters){
  
  rasterList <- list()
  transformList <- list()
  
  for(i in names(rasters)){
    # get data
    print(i)
    rast1 <- rasters[[i]] # select relevant raster
    fireVect <- values(rast1) # get values stored in raster
    
    # attempt to make data Gaussian-distributed
    bestNormaliser <- bestNormalize(values(rast1)) # identify best normaliser
    normalisedValues <- bestNormaliser$chosen_transform$x.t # apply the best transformation
    
    # create new raster with normalised data
    normalisedRaster <- rast1 
    values(normalisedRaster) <- normalisedValues # replace values of initial raster with transformed values,
    normalisedRaster <- climateStability::rescale0to1(normalisedRaster) # scale normalised values between zero and one
    
    rasterList[[i]] <- normalisedRaster
    transformList[[i]] <- bestNormaliser$chosen_transform
    
  }
  
  returnList <- list(rasters = rasterList, normalizers = transformList)
  return(returnList)
}


# function to load, crop, and create mean of outputs of several different climate models
bioclimEnsemble_fun <- function(listOfFileNames, folder = "C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/worldClim/future projections/5 minutes/",
                                varsOfInterest = c("PrecipDriestQuarter", "MeanTempDriestQuarter", "TempSeasonality", "AnnualMeanTemp", "AnnualPrecip")){
  
  # initialise list to store each stack of GCMs in
  rastList <- list()
  
  # load each dataset and select variables of interest
  for(i in names(listOfFileNames)){
    fileLocation <- paste(folder, listOfFileNames[[i]], sep = "")
    print(paste("loading dataset:", i, fileLocation))
    r <- rast(fileLocation) 
    names(r) <- worldClimBioNames
    r <- r[[varsOfInterest]] 
    rastList[[i]] <- r
  }
  
  means_ <- list()
  
  # loop through each variable of interest for each GCM
  for(var in varsOfInterest) {
    
    varList <- list()
    
    print(paste("working on", var))
    
    for(gcm in names(rastList)) {
      
      # store the same variables from diff GCMs
      varList[[paste(gcm, var, sep = "_")]] <- rastList[[gcm]][[var]]
    }
    
    # calculate mean of a given variable
    mn <- mean(rast(varList))
    names(mn) <- var
    print(mn)
    means_[[var]] <- mn
  }
  return(rast(means_))
}

# function to extract specific climate variables
extract_climate_variables <- function(fireRegimes, rasters, scenario, selectedClusters = NULL) {
  extracted <- fireRegimes %>% dplyr::select(clusterLetter) %>%
    mutate(
      scenario = scenario,
      PrecipDriestQuarter = extract(rasters$PrecipDriestQuarter, ., ID = F)[[1]],
      #PrecipSeasonality = extract(rasters$PrecipSeasonality, ., ID = F)[[1]],
      MeanTempDriestQuarter = extract(rasters$MeanTempDriestQuarter, ., ID = F)[[1]],
      TempSeasonality = extract(rasters$TempSeasonality, ., ID = F)[[1]],
      AnnualMeanTemp = extract(rasters$AnnualMeanTemp, ., ID = F)[[1]],
      AnnualPrecip = extract(rasters$AnnualPrecip, ., ID = F)[[1]]) %>% 
    drop_na()
  
  if(!is.null(selectedClusters)){
    extracted <- extracted %>% filter(clusterLetter %in% selectedClusters)
  }
  
  return(extracted)
}

# function to report object sizes, to ensure not keeping large objects that are no longer needed
sort(sapply(ls(),function(x){format(object.size(get(x)), units = "MB")})) 
