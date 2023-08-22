# DELINEATING REGIONS OF AUSTRALIA WITH INTERNALLY CONSISTENT FIRE REGIMES IN A CHANGING CLIMATE

# This script is broken into the following sections:
# 1. Organise preprocessed fire and fuel metrics - see "Metrics of fire regimes_v2.R"
# 2. Fit Gaussian mixture model with Mclust
# 3. Evaluate and map cluster assignment uncertainty
# 4. Create maps of pyroregions
# 5. Random forest to evaluate important features of pyroregions
# 6. Associations with historical climate
# 7. Organise climate projections
# 8. Comparison of climate hypervolumes
# 9. Quantify land use of each pyroregion
# 10. Plot heatmaps of seasonal patterns of hotspots

# load packages
pacman::p_load(ggplot2, tidyr, leaflet, tidyverse, sf, plotly, stringr, parallel, terra, tidyterra, 
               viridis, ggthemes, cowplot, fasterize, landscapemetrics, ggcorrplot, mclust, GGally, missForest, climateStability,
               NbClust, FactoMineR, factoextra, clustree, segmented, MASS, colorspace, bestNormalize, MODISTools, rgee, stars,
               MODIStsp, hypervolume, ggh4x, caret, randomForest)

# set working directory
setwd("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Analysis")

# source custom functions to help with mapping/plotting, data normalisation, and creating ensemble averages of multiple climate models,  
source("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Analysis/Pyroregionalisation functions.R")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 1. ORGANISE PREPROCESSED FIRE AND FUEL METRICS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# load rasters to use as templates for land; using ERA grid = 0.25 deg
aus_land_area <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Grid templates/aus_land_area.tif")
aus_land_binary <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Grid templates/aus_land_binary.tif")
aus_land_prop <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Grid templates/aus_land_prop.tif")

# load shapefile of Aus states
aus_states <- st_read("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Aus states shapefile/STE_2021_AUST_SHP_GDA2020/STE_2021_AUST_GDA2020.shp") %>%
  st_transform(st_crs(aus_land_area)) 
aus_states1 <- aus_states %>% st_simplify(dTolerance = 2000) 

# load raster with fire-only metrics; see "Metrics of fire regimes_v2.R"
fireMetricStack1 <- rast("fireMetricStack_20230821.tiff")

# visualise cells with at least one variable missing data (usually areas where MODIS hotspots have not been observed)
plot(sum(fireMetricStack1), colNA = "blue") 

# normalise the fire rasters and stack
fireMetricStack_normalised <- normalise_rasters(rasters = fireMetricStack1) # see "Pyroregionalisation functions.R"
fireMetricStack_trans <- rast(fireMetricStack_normalised$rasters)

# organise proxies of fuel, starting with NPP
npp_mean <- rast("npp_mean.tif")
npp_mean1 <- npp_mean %>% 
  resample(fireMetricStack1) %>% # resample to be on same grid as other variables
  sqrt(.) %>% # square root transform
  focal(w = 3, fun = mean, na.rm = T, na.policy = "only") # there's a handful of NAs; interpolate these values only, and then reinsert zeroes for ocean
npp_mean1 <- npp_mean1 + aus_land_binary # ensure ocean is zeroed
names(npp_mean1) <- "NPP"; plot(npp_mean1); rm(npp_mean)
fireMetricStack_trans$NPP <- climateStability::rescale0to1(npp_mean1) # add to stack

# NDVI
NDVI_diff <- rast("NDVI_diff.tif") # load preprocessed raster
fireMetricStack_trans$NDVI_diff <- climateStability::rescale0to1(NDVI_diff) # rescale and add to stack

# add new variable for relative difference between day and night fires
fireMetricStack_trans$nightFireDominance <- climateStability::rescale0to1(fireMetricStack_trans$yrsWithFire_N - fireMetricStack_trans$yrsWithFire_D)

# select 16 variables of interest, and reorder variables for plotting in thematic order
fireMetricStack_trans <- fireMetricStack_trans[[c("hsDensity_day", "nightFireDominance", "medianPatchSize", "maxPatchSize", "FRP_95", "meanMonthlyEmissions", 
  "timeSinceFire_richness", "timeSinceFire_contagion", "avhrr_ba_median","avhrr_ba_cv","peakFire_weeksFrom36", "fireSeasonLength", "NPP", "NDVI_diff")]]

# for variables that have intuitive units, plot them on their original scales
names(fireMetricStack_trans)
p1 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "hsDensity_day", begin = 0.1, trans = "log1p", breaks_ = c(0, 2, 6), legendName = expression(hs/km^{2})) + labs(title = "Day hotspot density"); p1
p3 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "medianPatchSize", trans = "log", breaks_ = c(1, 25, 1000), legendName = expression(km^{2})) + labs(title = "Median daily patch area"); p3
p4 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "maxPatchSize", trans = "log", breaks_ = c(0, 1, 50, 3000), legendName = expression(km^{2})) + labs(title = "Max daily patch area"); p4
p5 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "FRP_95", trans = "log", breaks_ = c(0, 5, 100, 3000), legendName = "MW") + labs(title = "Fire radiative power q95"); p5
p6 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "meanMonthlyEmissions", trans = "log1p", breaks_ = c(0, 5, 30), legendName = expression("g/C/m"^{2}~"/month"), begin = 0.1) + labs(title = "Mean monthly emissions") + theme(legend.position = c(0,0.05)); p6
p7 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "timeSinceFire_richness", trans = "identity", breaks_ = c(0,10, 20), legendName = "richness\nof ages") + labs(title = "Pyrodiversity"); p7
p9 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "avhrr_ba_median", trans = "log1p", breaks_ = c(0, 5, 50, 500), legendName = expression(km^{2})) + labs(title = "Median burned area"); p9
p10 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "avhrr_ba_cv", trans = "log1p", breaks_ = c(0.2, 2, 5), legendName = "CV ") + labs(title = "Burned area CV"); p10
p11 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "peakFire_weeksFrom36", trans = "identity", breaks_ = c(0,10, 20), legendName = "weeks") + labs(title = "Peak fire (weeks to Sept 1)"); p11
p12 <- raster_multilayer_plot_function(rasterForPlotting = fireMetricStack1, varName = "fireSeasonLength", trans = "identity", breaks_ = c(10,150,300), legendName = "days") + labs(title = "Fire season length"); p12
p13 <- raster_multilayer_plot_function(rasterForPlotting = (npp_mean1^2)*0.1, varName = "NPP") + labs(title = "Net primary prod.") + scale_fill_viridis_c(option = "D", na.value="grey100", trans = "sqrt", breaks = c(30,500,2000), name = expression("g/C/m"^{2}~year)) + theme(legend.position = c(0,0.05)); p13

# otherwise, plot unit-less variables on their transformed scales
p2 <- raster_multilayer_plot_function(varName = "nightFireDominance", begin = 0.1, breaks = c(0,0.5,1), trans = "identity", legendName = "scaled ") + labs(title = "Night-day hotspot difference"); p2
p8 <- raster_multilayer_plot_function(varName = "timeSinceFire_contagion", breaks = c(0,0.5,1), trans = "identity", legendName = "scaled ")  + labs(title = "Time since fire contagion"); p8
p14 <- raster_multilayer_plot_function(varName = "NDVI_diff") + labs(title = "NDVI summer-winter diffrence") + scale_fill_viridis_c(option = "D", na.value="grey100", n.breaks = 3, name = "scaled "); p14
allPlots <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14, ncol = 5, labels = "AUTO", label_x = -0.01)

# save figures
#pdf(file = "allFireMetrics_20230807.pdf", width = 16, height = 10)
allPlots
dev.off()

#jpeg(file = "allFireMetrics_20230807.jpg", width = 16, height = 10, res = 600, units = "in")
allPlots
dev.off()

# visualise correlation matrix of all 16 metrics
corr_matrix  <- layerCor(fireMetricStack_trans, fun = "pearson", na.rm = T)
corr_matrix  <- corr_matrix$pearson %>% data.frame()
colnames(corr_matrix) <- names(fireMetricStack_trans); row.names(corr_matrix) <- names(fireMetricStack_trans)
ggcorrplot(corr_matrix, hc.order = TRUE, colors = c("#6D9EC1", "white", "#E46726"), lab = T)

# create data frame with one column for each variable for use in cluster analysis
fireVars <- as.data.frame(fireMetricStack_trans[[1]], xy = T) %>% dplyr::select(x, y)
for(i in names(fireMetricStack_trans)){
  print(i)
  dat <- as.data.frame(fireMetricStack_trans[[i]], xy = T)
  fireVars <- left_join(fireVars, dat)
}
clustDat <- fireVars %>% drop_na() # keep for later for easy access of xy coords 
clustDat1 <- clustDat %>% dplyr::select(-x, -y)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 2. FIT GAUSSIAN MIXTURE MODEL WITH MCLUST ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# fit all possible model types with between 1 and 30 clusters
#set.seed(123); mod_all <- Mclust(clustDat1, 1:30, prior=priorControl()) 
# plot(mod_all, "BIC")

# VVV model was clearly best, so proceeding with that

# fit VVV model with between 1 and 30 clusters
set.seed(123); mod1 <- Mclust(clustDat1, 1:30, prior=priorControl(), modelNames = "VVV")
summary(mod1)
plot(mod1, "BIC")

# plot BIC to identify when addition of clusters no longer appreciably improves fit
IC_df <- t(t(mod1$BIC[,1])) %>% data.frame() %>% 
  rename(IC = ".") %>% 
  mutate(ICdiff = IC - lag(IC), components = 1:nrow(.)) 
evalPlot1 <- IC_df %>%
  ggplot(aes(components, IC)) +
  theme_minimal() +
  geom_point() +
  geom_line(colour = "grey50") +
  labs(x = "Number of clusters", y = "BIC") +
  geom_vline(xintercept = 17, linetype = "dashed", colour = "darkred", linewidth = 0.75)
evalPlot2 <- IC_df %>%
  ggplot(aes(components, ICdiff)) +
  geom_point() +
  theme_minimal() +
  geom_line(colour = "grey50") +
  coord_cartesian(ylim = c(-500,15000)) +
  labs(y = "Improvement in BIC\n with addition of one cluster", x = "Number of clusters") +
  geom_hline(yintercept = 0, linetype = "dashed")  +
  geom_vline(xintercept = 17, linetype = "dashed", colour = "darkred", linewidth = 0.75)

# save figure
#jpeg(file = "BIC plot.jpg", width = 6, height = 3, res = 600, units = "in")
plot_grid(evalPlot1, evalPlot2, labels = "AUTO", label_fontface = "plain")
dev.off()

# refit single model based on BIC-informed optimal number of clusters
set.seed(123); mod_selected <- Mclust(clustDat1, 17, prior=priorControl(), modelNames = "VVV")

# add model estimates to data frame
clustFitted <- clustDat1 %>% 
  data.frame() %>%
  bind_cols(clustDat[,c("x","y")]) %>% # add back xy coords
  mutate(cluster = mod_selected$classification, 
         uncertainty = mod_selected$uncertainty,
         probabilities = mod_selected$z) 

# order clusters based from lowest uncertainty to highest uncertainty
clustersOrdered <- clustFitted %>% 
  group_by(cluster) %>% 
  summarise(ordering = median(uncertainty)) %>%
  arrange(ordering) %>%
  mutate(clusterOrdered = 1:nrow(.), clusterLetter = c(letters, "z")[1:nrow(.)]) %>%
  dplyr::select(-ordering)

# add categorical (letter) cluster label
clustFitted <- clustFitted %>% 
  left_join(clustersOrdered) %>% dplyr::select(-cluster) %>% rename(cluster = clusterOrdered)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 3. EVALUATE AND MAP CLUSTER ASSIGNENT UNCERTAINTY ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# too many categories for established colour palettes, so creating a new colour palette by manually adding some colours to "Set3"
#display.brewer.pal(n = 12, name = 'Set3')
#brewer.pal(n = 12, name = 'Set3')
cols_ <- c("#8DD3C7", "gray30", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", 
           "goldenrod1","wheat4","goldenrod2","slategray1", "goldenrod3","#FFFFB3", "purple", "black") 
#pie(rep(1, length(cols_)), col = cols_)

# map uncertainty
uncertainty_ras <- clustFitted %>% 
  st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>%
  mutate(uncertainty = plyr::round_any(uncertainty, 0.05)) %>%
  rasterize(fireMetricStack_trans[[1]],  field = "uncertainty")

p_uncertainty <- ggplot() +
  theme_map() +
  geom_spatraster(data = uncertainty_ras, aes(fill = last))  +
  theme(axis.text = element_blank()) +
  theme(legend.position = c(0.07,0.13), legend.direction = "horizontal",
        legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10)) +
  labs(title = element_blank()) +
  scale_fill_binned_diverging(na.value = "white", 
                              breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6), 
                              labels = c(0,"","","0.3","","",0.6), 
                              name = "Classification\nuncertainty") +
  lims(x = c(113,154), y = c(-44.125, -10.8))+
  geom_sf(data = aus_states1, fill = NA, colour = "grey50") 
p_uncertainty

# for each cluster, visualise uncertainty above a threshold to identify areas where a given cluster is most uncertain 
thresh <- 0.2
uncertaintyMap <- ggplot() +
  theme_map() +
  geom_tile(clustFitted %>% filter(uncertainty <= thresh), mapping = aes(x, y), fill = "grey90") +
  geom_tile(clustFitted %>% filter(uncertainty > thresh), mapping = aes(x, y, fill = factor(clusterLetter))) +
  scale_fill_manual(values = cols_, name = "Fire regime") +
  labs(title = paste("Model uncertainty >", thresh) ) 
uncertaintyMap
# some areas of high uncertainty and interchangability, especially arid west

# report median uncertainty
uncertaintyLabel <- clustFitted %>% group_by(clusterLetter) %>% summarise(medianUncertainty = median(uncertainty)) %>% 
  arrange(medianUncertainty) %>%
  mutate(label = paste("m = ", signif (medianUncertainty, 1), sep = ""))

# histogram of uncertainty
uncertaintyHist <- ggplot(clustFitted, aes(uncertainty, after_stat(ndensity))) + 
  theme_minimal() +
  facet_wrap(~clusterLetter, nrow = 3) +
  geom_histogram(aes(fill = factor(clusterLetter)), alpha = 1) +
  scale_fill_manual(values = cols_, name = "Fire regime") +
  labs(x = "Uncertainty (probability)", y = element_blank()) +
  geom_text(data = uncertaintyLabel, aes(x = 0.35, y = 0.9, label = label), size = 3) 
uncertaintyHist

# export Supplementary Figure
#jpeg(file = "uncertaintyHistograms.jpg", width = 7, height = 5, res = 600, units = "in")
uncertaintyHist
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 4. CREATE MAPS OF PYROREGIONS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# rasterize clusters
cluster_ras <- clustFitted %>% 
  st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>% 
  rasterize(fireMetricStack_trans[[1]],  field = "cluster") 
plot(cluster_ras)

# ad posthoc category for missing data regime, constituting a regime in which fire is very rare
cluster_ras[is.na(cluster_ras)] <-  nrow(unique(cluster_ras$last)) + 1
cluster_ras <- cluster_ras + aus_land_binary # zero out ocean
plot(cluster_ras)

# majority focal window
cluster_ras_smoothed <- cluster_ras %>% focal(w = 5, fun = modal, na.rm = T, na.policy = "omit")
cluster_ras_smoothed <- cluster_ras_smoothed + aus_land_binary
plot(cluster_ras_smoothed)

# set categories
clusts <- data.frame(id=unique(values(cluster_ras_smoothed)[,1])) %>% arrange(id) %>% drop_na() 
clusts$clusters <- letters[clusts$id]
levels(cluster_ras_smoothed) <- clusts
plot(cluster_ras_smoothed)

# turn into polygons, smooth and plot
cluster_poly <- cluster_ras_smoothed %>% 
  as.polygons(dissolve = T) %>% st_as_sf() %>%
  smoothr::smooth(method = "ksmooth", smoothness = 1) 
plot(cluster_poly)

# create label at centroid of each distinct polygon by first breaking up non-contiguous polys
clusterLabels <- cluster_poly %>%
  st_cast(to = 'MULTILINESTRING') %>% 
  st_cast(to = "LINESTRING") %>%
  st_cast(to = "POLYGON") %>%
  rename(id = clusters) %>%  
  mutate(area = as.numeric(st_area(.))/1000000) %>%
  st_point_on_surface() %>%
  #st_centroid() %>%
  mutate(x = st_coordinates(.)[,1], y = st_coordinates(.)[,2]) %>%
  data.frame()

# plot polygons and outline of Aus states
p_map_smoothed <- ggplot(cluster_poly) +
  theme_map() +
  geom_sf(aes(fill = clusters), colour = NA, alpha = .6) +
  scale_color_manual(values = cols_) +
  scale_fill_manual(values = cols_, name = "Fire regime") +
  geom_text(data = clusterLabels %>% filter(area > 15000), aes(x = x, y = y, label = id, colour = id), size = 4, colour = "black") +
  geom_sf(data = aus_states1, fill = NA, colour = "grey50") +
  theme(legend.position = "right", plot.margin = unit(c(5.5,5.5,5.5,100), "pt")) +
  lims(x = c(113.375, 153.625), y = c(-43.625, -10.875))
p_map_smoothed

# export shapefile of pyroregions
#st_write(cluster_poly, "cluster_poly.shp", append = F)

# make multi-panel map of pyroregions and uncertainty
multiMap <- ggdraw(p_map_smoothed) +
  draw_plot(p_uncertainty , x = -.32, y = .27, scale = 0.33, height = 1.1) 

# export map (Fig 2)
#jpeg(file = "clusterMaps_and_uncertainty_aug15.jpg", width = 10, height = 8, res = 600, units = "in")
multiMap
dev.off()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 5. RANDOM FOREST TO EVALUATE IMPORTANT FEATURES OF PYROREGIONS ---- 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# format raster data to fit random forest model
clusters <- clustDat %>%  # orig predictors with lat longs
  st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>% # turn spatial
  mutate(clusterLetter = extract(cluster_ras_smoothed,., ID = F)[[1]]) # extract smoothed pyroregion 

# for reporting in paper, calc median patch sizes of pyroregions
medPatchSizes <- clusters %>% mutate(medianPatchSize = extract(fireMetricStack1$medianPatchSize, .)[,2]) %>%
  data.frame() %>%  group_by(clusterLetter) %>% summarise(median(medianPatchSize))

# summarise each variable using violin plots
selectedVarNames <- names(clustFitted)[1:14]
clustSummaries <- clustFitted %>% 
  group_by(clusterLetter) %>% 
  dplyr::select(clusterLetter, all_of(selectedVarNames)) %>%
  pivot_longer(cols =  selectedVarNames,
               values_to = "metricValue", names_to = "metricType") 

# select data for random forest
rf_dat <- data.frame(clusters) %>% dplyr::select(-geometry) %>% filter(clusterLetter != "r")
rf_dat$clusterLetter <- factor(rf_dat$clusterLetter, sort(unique(rf_dat$clusterLetter)))

# option to split data into training and test sets, or to use k-fold cross-validation instead. 
# In this case, we went with k-fold CV, so retaining all data (p = 1) at this step
set.seed(2021)
inTrain <- createDataPartition(rf_dat$clusterLetter, p = 1, list = FALSE)[,1] # p = 1 because we ended up using k-fold CV 
trainDat <- rf_dat[inTrain,]
#testDat  <- rf_dat[-inTrain,]

# fit global model classifying the pyroregions
trControl_global <- trainControl(method = "cv", number = 10, allowParallel = T)
tuneGrid_global <-  expand.grid(.mtry = seq(2, (ncol(trainDat) - 1),2))
set.seed(825)
rf1 <- train(clusterLetter ~ ., data = trainDat, method = "rf", metric = "Kappa", trControl = trControl_global, 
             tuneGrid = tuneGrid_global, importance = T, ntree = 1000)
#saveRDS(rf1, "rf1_20230712.rds")
#rf1 <- readRDS("rf1_20230712.rds")

# retrieve variable importance values
# see here for blog on class-specific variable importance: https://ellisvalentiner.com/post/variable-importance-bad-behavior/
varImp <- importance(rf1$finalModel, scale = F) %>% # from https://explained.ai/rf-importance/ "Make sure that you don't use the MeanDecreaseGini column in the importance data frame. You want MeanDecreaseAccuracy, which only appears in the importance data frame if you turn on importance=T when constructing the Random Forest. rf <- randomForest(hyper-parameters..., importance=T) imp <- importance(rf, type=1, scale = F) # permutation importances"
  data.frame() %>% 
  dplyr::select(-MeanDecreaseGini) %>% # essential to choose type 1 importance: https://explained.ai/rf-importance/ but see here which shows that the first-class specific importance values are only type 1: https://ellisvalentiner.com/post/variable-importance-bad-behavior/
  mutate(variable = rownames(.)) %>% 
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "cluster", values_to = "importance") %>%
  group_by(cluster) %>%
  mutate(relativeImportance = importance/sum(importance),
         cluster = ifelse(cluster == "MeanDecreaseAccuracy", "Mean", cluster)) 

# order regimes 
varImp$cluster <- factor(varImp$cluster, unique(varImp$cluster))

# add median values for each regime/variable
clusterMedians <- clustSummaries %>% group_by(metricType, clusterLetter) %>% summarise(medianValue = median(metricValue)) %>% rename(cluster = clusterLetter, variable = metricType)
varImp <- varImp %>% left_join(clusterMedians)

# order predictors according to average importance across all pyroregions
varOrder <- varImp %>% filter(cluster == "Mean") %>% arrange(desc(relativeImportance))
varImp$variable <- factor(varImp$variable, varOrder$variable)
varImp$cluster <- factor(varImp$cluster, c(letters[1:length(unique(varImp$cluster))-1], "Mean")) # order so that "mean" is plotted last

# identify each regime's most important features
mostImportant <- varImp %>% group_by(cluster) %>% arrange(cluster, desc(relativeImportance)) %>% slice_head(n = 3) %>% mutate(rank = 1:3)

# change variable names
facet_var_names <- c("NPP" = "NPP", "NDVI_diff" = "NDVI difference","avhrr_ba_median" = "Burned area median",
                     "fireSeasonLength" = "Fire season length", "avhrr_ba_cv" = "Burned area CV",
                     "timeSinceFire_richness" = "Pyrodiversity", "meanMonthlyEmissions" = "Emissions",
                     "hsDensity_day" = "Day hotspot density","timeSinceFire_contagion" = "Contagion",
                     "medianPatchSize" = "Median patch size","maxPatchSize" = "Max patch size",
                     "FRP_95" = "FRP 95","nightFireDominance" = "Night fire","peakFire_weeksFrom36" = "Peak fire week")

# heat map of variable importance, regimes on top (Fig 2A)
featureImpPlot <- ggplot(varImp, aes(cluster, variable)) +
  theme_minimal() +
  geom_point(aes(colour = medianValue, size = relativeImportance)) +
  scale_size("Relative\nimportance", range = c(1,15)) +
  scale_x_discrete(position = "top") +
  theme(legend.position = "bottom", plot.margin = ggplot2::margin(0, 1, 0, 0, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(label = plyr::round_any(relativeImportance, 0.01)), colour = "black", size = 3.5) +
  scale_colour_gradient2(midpoint = 0.5, 
                       high = scales::muted("darkred"), 
                       low = scales::muted("steelblue"),
                       name = "Median\npredictor\nvalue", breaks = c(0.2,0.5,0.8)) +
  scale_y_discrete(limits=rev, "Features of fire regimes",
                   labels = facet_var_names)+
  labs(x = "Fire regime")  +
  theme(legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = 9),
        axis.text = element_text(size = 13), axis.title = element_text(size = 13))
featureImpPlot

# violin plots of the distributions of each variable for each pyroregions
violinSummaries <- clustSummaries %>% 
  mutate(metricType = factor(metricType, levels(varImp$variable))) %>% # order variables according to random forest variable importance
  ggplot(aes(metricValue, y = as.factor(clusterLetter), group = clusterLetter, 
             fill = as.factor(clusterLetter))) +
  theme_minimal()+
  facet_wrap(~metricType, ncol = 3, labeller = labeller(metricType = facet_var_names)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), linewidth = 0.4, alpha = 0.7) +
  theme(legend.position = "none") +
  labs(y = "Fire regime", x = "Fire metric value") +
  scale_x_continuous(breaks = c(0,0.5,1))+ 
  scale_y_discrete(limits=rev) +
  scale_fill_manual(values = cols_)
violinSummaries

#jpeg(file = "violins.jpg", width = 10, height = 13, res = 600, units = "in")
violinSummaries
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 6. ASSOCIATIONS WITH HISTORICAL CLIMATE ---- 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# names of each bioclim variable from worldclim
worldClimBioNames <- c("AnnualMeanTemp","MeanDiurnalRange","Isothermality","TempSeasonality","MaxTempWarmestMonth","MinTempColdestMonth","TempAnnualRange","MeanTempWettestQuarter",
                       "MeanTempDriestQuarter","MeanTempWarmestQuarter","MeanTempColdestQuarter","AnnualPrecip","PrecipWettestMonth","PrecipDriestMonth","PrecipSeasonality","PrecipWettestQuarter",
                       "PrecipDriestQuarter",'PrecipWarmestQuarter',"PrecipColdestQuarter")

# load rasters of historical climate
worldClim_bio <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/worldClim/wc2.1_country/AUS_wc2.1_30s_bio.tif") %>% crop(ext(109.875, 155.125, -44.125, -9.875)) %>% aggregate(fact = 10)
names(worldClim_bio) <- worldClimBioNames; plot(worldClim_bio)

# function to scale between zero and one
rescale_zeroToOne <- function(x) ((x - min(x))/(max(x) - min(x)))

# normalise the extracted data for use in random forest
normalise_vec <- function(dat){
  bestNormaliser <- bestNormalize(dat) # identify best normaliser
  normalisedValues <- bestNormaliser$chosen_transform$x.t # apply the best transformation
  print(bestNormaliser$chosen_transform)
  return(normalisedValues)
}

# extract environmental variables for current and future climates
clusterPoints <- cluster_ras_smoothed %>% as.data.frame(xy = T) %>% st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>% rename(clusterLetter = clusters)
cluster_enviro <- extract_climate_variables(fireRegimes = clusterPoints, rasters = worldClim_bio, scenario = "historic") 
cluster_enviro_scaled <- cluster_enviro %>% 
  data.frame() %>%
  mutate(PrecipDriestQuarter = rescale_zeroToOne(normalise_vec(PrecipDriestQuarter)),
         MeanTempDriestQuarter = rescale_zeroToOne(normalise_vec(MeanTempDriestQuarter)),
         TempSeasonality = rescale_zeroToOne(normalise_vec(TempSeasonality)),
         AnnualMeanTemp = rescale_zeroToOne(normalise_vec(AnnualMeanTemp)),
         AnnualPrecip = rescale_zeroToOne(normalise_vec(AnnualPrecip)))

# summarise enviro variable using violin plots
selectedClimVarNames <- c("PrecipDriestQuarter", "MeanTempDriestQuarter", "TempSeasonality", "AnnualMeanTemp", "AnnualPrecip")

# summarise cluster attributes
climateSummaries <- cluster_enviro_scaled %>% 
  dplyr::select(clusterLetter, all_of(selectedClimVarNames)) %>%
  pivot_longer(cols =  selectedClimVarNames,
               values_to = "metricValue", names_to = "metricType") 

# facet labels
climate_facet_names <- c('AnnualMeanTemp' = "Mean Temperature", 'AnnualPrecip'="Mean Precipitation",'MeanTempDriestQuarter'="Mean Temperature\nDriest Quarter",
                         'PrecipDriestQuarter'="Precipitation\nDriest Quarter",'TempSeasonality'="Temperature\nSeasonality")

# violin plots of historical climate variables
violinSummaries_climate <- climateSummaries %>% 
  ggplot(aes(metricValue, y = as.factor(clusterLetter), group = clusterLetter, 
             fill = as.factor(clusterLetter))) +
  theme_minimal()+
  facet_wrap(~metricType, ncol = 5, labeller = labeller(metricType = climate_facet_names)) +
  geom_violin(draw_quantiles = c(0.25,0.5,0.75), linewidth = 0.4, alpha = 0.7) +
  theme(legend.position = "none") +
  labs(y = "Fire regime", x = "Fire metric value") +
  scale_x_continuous(breaks = c(0,0.5,1))+
  scale_y_discrete(limits=rev)  +
  scale_fill_manual(values = cols_, name = "Fire regime")
violinSummaries_climate

#jpeg(file = "violins_climate.jpg", width = 10, height = 4, res = 600, units = "in")
violinSummaries_climate
dev.off()

# Fit random forest to evaluate relative importance of climate in associating with each pyroregion
# select data for random forest
rf_dat_env <- cluster_enviro_scaled %>% dplyr::select(clusterLetter, PrecipDriestQuarter, MeanTempDriestQuarter, TempSeasonality, AnnualMeanTemp, AnnualPrecip)

# option to split data into training and test sets; here we've retained all data because using k-fold CV instead
set.seed(2021)
inTrain_env <- createDataPartition(rf_dat_env$clusterLetter, p = 1, list = FALSE)[,1] # p = 1 because I ended up using k-fold CV instead of splitting data into train and test
trainDat_env <- rf_dat_env[inTrain_env,]
#testDat_env  <- rf_dat_env[-inTrain_env,] 

# fit model classifying the pyroregions in response to environmental variables
trControl_env <- trainControl(method = "cv", number = 10, allowParallel = T)
tuneGrid_env <-  expand.grid(.mtry = seq(1, (ncol(trainDat_env) - 1),1))
set.seed(825)
rf_env1 <- train(clusterLetter ~ ., data = trainDat_env, method = "rf", metric = "Kappa", 
                 trControl = trControl_env, tuneGrid = tuneGrid_env, importance = T, ntree = 1000)
#saveRDS(rf_env1, "rf_env1_20230712.rds")
#rf_env1 <- readRDS("rf_env1_20230712.rds")
summary(rf_env1)

# class-specific variable importance: https://ellisvalentiner.com/post/variable-importance-bad-behavior/
varImp_env <- importance(rf_env1$finalModel, scale = F) %>% # from https://explained.ai/rf-importance/ "Make sure that you don't use the MeanDecreaseGini column in the importance data frame. You want MeanDecreaseAccuracy, which only appears in the importance data frame if you turn on importance=T when constructing the Random Forest. rf <- randomForest(hyper-parameters..., importance=T) imp <- importance(rf, type=1, scale = F) # permutation importances"
  data.frame() %>% 
  dplyr::select(-MeanDecreaseGini) %>% # essential to choose type 1 importance: https://explained.ai/rf-importance/ but see here which shows that the first-class specific importance values are only type 1: https://ellisvalentiner.com/post/variable-importance-bad-behavior/
  mutate(variable = rownames(.)) %>% 
  pivot_longer(cols = 1:(ncol(.)-1), names_to = "cluster", values_to = "importance") %>%
  group_by(cluster) %>%
  mutate(relativeImportance = importance/sum(importance),
         cluster = ifelse(cluster == "MeanDecreaseAccuracy", "Mean", cluster)) 

# calculate median climate values
climateMedians <- rf_dat_env %>% 
  pivot_longer(cols =  2:6, values_to = "metricValue", names_to = "variable") %>%
  group_by(clusterLetter, variable) %>%
  summarise(medianVal = median(metricValue)) %>% 
  rename(cluster = clusterLetter)

# join variable importance with median values
varImp_env <- varImp_env %>%
  left_join(climateMedians)

# order predictors according to mean influence
varOrder_env <- varImp_env %>% filter(cluster == "Mean") %>% arrange(desc(relativeImportance))
varImp_env$variable <- factor(varImp_env$variable, varOrder_env$variable)

# order regimes
varImp_env$cluster <- factor(varImp_env$cluster, unique(varImp_env$cluster))

# identify each regime's most important features
mostImportant_env <- varImp_env %>% group_by(cluster) %>% arrange(cluster, desc(relativeImportance)) %>% slice_head(n = 3) %>% mutate(rank = 1:3)

# heat map of variable importance
varImpPlot_Clim <- ggplot(varImp_env, aes(cluster, variable)) +
  theme_minimal() +
  geom_point(aes(size = relativeImportance, colour = medianVal)) +
  scale_size("Relative\nimportance", range = c(1,10)) +
  labs(x = element_blank(), y = "Climate associations") +
  scale_x_discrete(position = "top") +
  theme(legend.position = "bottom", plot.margin = ggplot2::margin(0, 1, 0, 0, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(aes(label = plyr::round_any(relativeImportance, 0.01)), colour = "black", size = 3) +
  scale_y_discrete(limits=rev)  +
  scale_colour_gradient2(midpoint = 0.5, 
                         high = scales::muted("darkred"), 
                         low = scales::muted("steelblue"),
                         name = "Median\npredictor\nvalue", breaks = c(0.3,0.5,0.7)) +
  theme(legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = 9),
        axis.text = element_text(size = 13), axis.title = element_text(size = 13))
varImpPlot_Clim

# save multipanel plot of map + variable importance
#jpeg(file = "varImpPlot_circles.jpg", width = 9, height = 9, res = 600, units = "in")
plot_grid(featureImpPlot,
          NULL,
          varImpPlot_Clim,
          ncol = 1, rel_heights = c(1,0.05,0.4), axis = "lr", align = "v")
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 7. ORGANISE CLIMATE PROJECTIONS ----
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Calc the ensemble-average of the outputs of 5 different climate models under 3 SSPs
# SSP3-7.0
ensemble_SSP370 <- bioclimEnsemble_fun(
  listOfFileNames = list(MIROC = "SSP370/AustraliaCropped/Australia_wc2.1_5m_bioc_MIROC6_ssp370_2081-2100.tif", 
                         ACCESS = "SSP370/AustraliaCropped/Australia_wc2.1_5m_bioc_ACCESS-CM2_ssp370_2081-2100.tif",
                         BCC = "SSP370/AustraliaCropped/Australia_wc2.1_5m_bioc_BCC-CSM2-MR_ssp370_2081-2100.tif",
                         GISS = "SSP370/AustraliaCropped/Australia_wc2.1_5m_bioc_GISS-E2-1-G_ssp370_2081-2100.tif",
                         MPI = "SSP370/AustraliaCropped/Australia_wc2.1_5m_bioc_MPI-ESM1-2-HR_ssp370_2081-2100.tif"))
plot(ensemble_SSP370)

# SSP2-4.5
ensemble_SSP245 <- bioclimEnsemble_fun(
  listOfFileNames = list(MIROC = "SSP245/AustraliaCropped/Australia_wc2.1_5m_bioc_MIROC6_ssp245_2081-2100.tif", 
                         ACCESS = "SSP245/AustraliaCropped/Australia_wc2.1_5m_bioc_ACCESS-CM2_ssp245_2081-2100.tif",
                         BCC = "SSP245/AustraliaCropped/Australia_wc2.1_5m_bioc_BCC-CSM2-MR_ssp245_2081-2100.tif",
                         GISS = "SSP245/AustraliaCropped/Australia_wc2.1_5m_bioc_GISS-E2-1-G_ssp245_2081-2100.tif",
                         MPI = "SSP245/AustraliaCropped/Australia_wc2.1_5m_bioc_MPI-ESM1-2-HR_ssp245_2081-2100.tif"))
plot(ensemble_SSP245)

# SSP1-2.6
ensemble_SSP126 <- bioclimEnsemble_fun(
  listOfFileNames = list(MIROC = "SSP126/AustraliaCropped/Australia_wc2.1_5m_bioc_MIROC6_ssp126_2081-2100.tif", 
                         ACCESS = "SSP126/AustraliaCropped/Australia_wc2.1_5m_bioc_ACCESS-CM2_ssp126_2081-2100.tif",
                         BCC = "SSP126/AustraliaCropped/Australia_wc2.1_5m_bioc_BCC-CSM2-MR_ssp126_2081-2100.tif",
                         GISS = "SSP126/AustraliaCropped/Australia_wc2.1_5m_bioc_GISS-E2-1-G_ssp126_2081-2100.tif",
                         MPI = "SSP126/AustraliaCropped/Australia_wc2.1_5m_bioc_MPI-ESM1-2-HR_ssp126_2081-2100.tif"))
plot(ensemble_SSP126)

# calculate climate anomalies and export supplementary figure
meanTempAnom_126 <- resample(ensemble_SSP126$AnnualMeanTemp, worldClim_bio$AnnualMeanTemp) - worldClim_bio$AnnualMeanTemp
meanTempAnom_245 <- resample(ensemble_SSP245$AnnualMeanTemp, worldClim_bio$AnnualMeanTemp) - worldClim_bio$AnnualMeanTemp
meanTempAnom_370 <- resample(ensemble_SSP370$AnnualMeanTemp, worldClim_bio$AnnualMeanTemp) - worldClim_bio$AnnualMeanTemp
anoms <- c(meanTempAnom_126, meanTempAnom_245, meanTempAnom_370)
names(anoms) <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0")
#jpeg(file = "tempAnomalies.jpg", width = 7, height = 3, res = 600, units = "in")
ggplot() +
  geom_spatraster(data = anoms) +
  facet_wrap(~lyr) +
  theme_map() +
  scale_fill_viridis_b(breaks = seq(0, 4.5, by = 0.5), labels = c(0, "", 1, "", 2, "", 3, "", 4, ""), na.value="grey100", name = "Mean\ntemp\nanomaly") +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))
dev.off()  

# extract projected climates for each cell of each pyroregion
cluster_ensemble_ssp126 <- extract_climate_variables(fireRegimes = clusterPoints, rasters = ensemble_SSP126, scenario = "SSP1-2.6_ensemble")
cluster_ensemble_ssp245 <- extract_climate_variables(fireRegimes = clusterPoints, rasters = ensemble_SSP245, scenario = "SSP2-4.5_ensemble")
cluster_ensemble_ssp370 <- extract_climate_variables(fireRegimes = clusterPoints, rasters = ensemble_SSP370, scenario = "SSP3-7.0_ensemble")

# bind datasets together so that rescaling is done jointly
joinedClimates <- cluster_enviro %>% 
  bind_rows(cluster_ensemble_ssp126, cluster_ensemble_ssp245, cluster_ensemble_ssp370) %>%
  ungroup() %>% data.frame() %>%
  mutate(PrecipDriestQuarter = rescale_zeroToOne(PrecipDriestQuarter),
         MeanTempDriestQuarter = rescale_zeroToOne(MeanTempDriestQuarter),
         TempSeasonality = rescale_zeroToOne(TempSeasonality),
         AnnualMeanTemp = rescale_zeroToOne(AnnualMeanTemp),
         AnnualPrecip = rescale_zeroToOne(AnnualPrecip))

# density plots of how each variable is predicted to change in future; comparing two pyroregions
densP <- joinedClimates %>% filter(clusterLetter %in% c("b", "d"), scenario %in% c("historic", "SSP2-4.5_ensemble")) %>%
  mutate(scenario = ifelse(scenario == "SSP2-4.5_ensemble", "SSP2-4.5", scenario),
         clusterLetter = ifelse(clusterLetter == "b", "b: Monsoon Savanna", "d: Wet Temperate Forest")) %>%
  pivot_longer(cols = c(PrecipDriestQuarter, MeanTempDriestQuarter, TempSeasonality, AnnualMeanTemp, AnnualPrecip), 
               names_to = "climateVariable", values_to = "climateValues") %>%
  ggplot(aes(climateValues*100, colour = clusterLetter, fill = clusterLetter, linetype = scenario)) +
  theme_minimal() +
  facet_wrap(~climateVariable, scales = "free_y", labeller = labeller(climateVariable = climate_facet_names), nrow = 3) + #rows = vars(clusterLetter),
  geom_density(linewidth = 0.7, alpha = 0.25)+
  theme( panel.grid.minor = element_blank(), legend.position = c(0.75,0.15)) +
  scale_color_manual(values = c("gray30",  "#FB8072"), name = "Pyroregion") +
  scale_fill_manual(values = c("gray30",  "#FB8072"), name = "Pyroregion") +
  scale_linetype_manual(values = c("solid","11"), name = "Climate") +
  labs(x = "Climate values (scaled)") +
  theme(legend.spacing.y = unit(0.1, "cm"), legend.key.size = unit(0.35, "cm"), legend.title = element_text(size = 9), legend.text = element_text(size = 8)) +
  scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1"))
  
# manually change each facet's axis range; necessary to discern patterns because of differing concentrations of the densities
densP <- densP + 
  facetted_pos_scales(
    y = list(climateVariable == "AnnualMeanTemp" ~ scale_y_continuous(limits = c(0, 0.19)),
             climateVariable == "AnnualPrecip" ~ scale_y_continuous(limits = c(0, 0.06)),
             climateVariable == "MeanTempDriestQuarter" ~ scale_y_continuous(limits = c(0, 0.16)),
             climateVariable == "PrecipDriestQuarter" ~ scale_y_continuous(limits = c(0, 0.65)),
             climateVariable == "TempSeasonality" ~ scale_y_continuous(limits = c(0, 0.045)))
  )

#jpeg(file = "nicheComparison.jpg", width = 3.5, height = 4.5, res = 600, units = "in")
densP
dev.off()

# split scaled datasets
currentDat <- joinedClimates %>% filter(scenario == "historic")
ensemble_ssp126_scaled <- joinedClimates %>% filter(scenario == "SSP1-2.6_ensemble")
ensemble_ssp245_scaled <- joinedClimates %>% filter(scenario == "SSP2-4.5_ensemble")
ensemble_ssp370_scaled <- joinedClimates %>% filter(scenario == "SSP3-7.0_ensemble")
rm(joinedClimates)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 8. COMPARISON OF CLIMATE HYPERVOLUMES ---- 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# function to create hypervolume for each pyroregion
create_hypervol <- function(dataset1, hv_method = "svm", quantile.requested_ = 0.95, svm.gamma_ = 0.5, verbose_ = F, 
                            vars = c("PrecipDriestQuarter", "MeanTempDriestQuarter", "TempSeasonality", "AnnualMeanTemp", "AnnualPrecip")) {
  
  outputList <- list()
  
  # loop through each fire regime
  for(i in unique(dataset1$clusterLetter)) {
    
    print(i)
    
    # filter pyroregion and variables of interest 
    dat1 <- dataset1 %>% filter(clusterLetter == i) %>% dplyr::select(vars)

    # compute hypervolumes for both datasets 
    if(hv_method == "gaussian") {
      set.seed(3); hv1 <- hypervolume_gaussian(dat1, name = unique(dataset1$scenario), verbose = verbose_, quantile.requested = quantile.requested_, sd = 4) # option to modify quantile; hard-coding increase in sd from 3 to 4 to sample a larger potential region and therefore increase accuracy of hypervolume estimate at the expense of computation time
    }
    if(hv_method == "svm") {
      set.seed(3); hv1 <- hypervolume_svm(dat1, name = unique(dataset1$scenario), verbose = verbose_, svm.gamma = svm.gamma_) # option to modify the wrap tightness
    }
    
    outputList[[i]]$hv1 <- hv1
    }
  return(outputList)
}

# function to test whether future climate is within the historical hypervolume; either on a pyroregion-specific basis ("yes") or continental basis ("no")
inclusion_test_fun <- function(hypervol, newdat, vars = c("PrecipDriestQuarter", "MeanTempDriestQuarter", "TempSeasonality", "AnnualMeanTemp", "AnnualPrecip"),
                               pyrome_specific_newclimate = "yes") {
  
  outDat <- list()
  
  if(pyrome_specific_newclimate == "yes"){
    
    # loop through each pyroregion
    for(i in names(hypervol)){
      print(i)
      # select relevant cluster of new data
      nd <- newdat %>% filter(clusterLetter == i) 
      
      # add column for whether new PYROME-SPECIFIC CLIMATE DATA is inside or outside of the same pyroregion's historical hypervolume
      nd$incl <- hypervolume_inclusion_test(hypervol[[i]]$hv1,  nd %>% dplyr::select(all_of(vars)), reduction.factor = 1, fast.or.accurate = "accurate", verbose = T)
      
      # add column for whether the hypervolume was Guassian or SVM
      nd$hypervolume_type <- hypervol[[i]]$hv1@Method
      nd$hypervolume_type <- ifelse(grepl("Gaus", nd$hypervolume_type), "Gaussian", "SVM")
      outDat[[i]] <- nd
    }
    outDat1 <- bind_rows(outDat) %>% mutate(noveltyType = "Pyrome-specific")
  }
  
  if(pyrome_specific_newclimate == "no"){
    
    # loop through each pyroregion
    for(i in names(hypervol)){
      print(i)
      # select relevant cluster of new data
      nd <- newdat 
      
      # add column for whether new CONTINENTAL CLIMATE DATA is inside or outside of a given pyroregion's historical hypervolume
      nd$incl <- hypervolume_inclusion_test(hypervol[[i]]$hv1,  nd %>% dplyr::select(all_of(vars)), reduction.factor = 1, fast.or.accurate = "accurate", verbose = T)
      
      # add column for whether the hypervolume was Guassian or SVM
      nd$hypervolume_type <- hypervol[[i]]$hv1@Method
      nd$hypervolume_type <- ifelse(grepl("Gaus", nd$hypervolume_type), "Gaussian", "SVM")
      outDat[[i]] <- nd
    }
    # summarise how many pyroregion's historical hyperclimates does the future climate fall within; a value of zero indicates it is entirely novel 
    outDat1 <- bind_rows(outDat) %>% group_by(clusterLetter, scenario, geometry, hypervolume_type) %>% summarise(sumOfIncluded = sum(incl)) %>% mutate(noveltyType = "Continental") %>% data.frame()
  }
  return(outDat1)
}


# create hypervolumes using both gaussian and svm methods
climatehypervol_gauss <- create_hypervol(dataset1 = currentDat, hv_method = "gaussian")
climatehypervol_svm <- create_hypervol(dataset1 = currentDat, hv_method = "svm")
#saveRDS(list(climatehypervol_gauss, climatehypervol_svm), "hypervolume_list.rds")
# hypervolList <- readRDS("hypervolume_list.rds")
# climatehypervol_gauss <-  hypervolList[[1]]
# climatehypervol_svm <- hypervolList[[2]]

climatehypervol_svm$d$hv1@Volume/climatehypervol_svm$b$hv1@Volume

# calculate pyrome-specific novelty; i.e., whether the future climate in a pyorome's location is within its historical hypervolume
inclusion_currentTest_gaus <- inclusion_test_fun(climatehypervol_gauss, currentDat)
inclusion_currentTest_svm <- inclusion_test_fun(climatehypervol_svm, currentDat)
inclusion_ssp126_gaus <- inclusion_test_fun(climatehypervol_gauss, ensemble_ssp126_scaled)
inclusion_ssp126_svm <- inclusion_test_fun(climatehypervol_svm, ensemble_ssp126_scaled)
inclusion_ssp245_gaus <- inclusion_test_fun(climatehypervol_gauss, ensemble_ssp245_scaled)
inclusion_ssp245_svm <- inclusion_test_fun(climatehypervol_svm, ensemble_ssp245_scaled)
inclusion_ssp370_gaus <- inclusion_test_fun(climatehypervol_gauss, ensemble_ssp370_scaled)
inclusion_ssp370_svm <- inclusion_test_fun(climatehypervol_svm, ensemble_ssp370_scaled)

# calculate whether climate is novel to the all pyrome's on the continent; i.e., whether the future climate fall within the climate of one or more other pyromes
inclusion_currentTest_gaus_cont <- inclusion_test_fun(climatehypervol_gauss, currentDat, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_currentTest_svm_cont <- inclusion_test_fun(climatehypervol_svm, currentDat, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_ssp126_gaus_cont <- inclusion_test_fun(climatehypervol_gauss, ensemble_ssp126_scaled, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_ssp126_svm_cont <- inclusion_test_fun(climatehypervol_svm, ensemble_ssp126_scaled, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_ssp245_gaus_cont <- inclusion_test_fun(climatehypervol_gauss, ensemble_ssp245_scaled, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_ssp245_svm_cont <- inclusion_test_fun(climatehypervol_svm, ensemble_ssp245_scaled, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_ssp370_gaus_cont <- inclusion_test_fun(climatehypervol_gauss, ensemble_ssp370_scaled, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))
inclusion_ssp370_svm_cont <- inclusion_test_fun(climatehypervol_svm, ensemble_ssp370_scaled, pyrome_specific_newclimate = "no") %>% mutate(incl = ifelse(sumOfIncluded > 0, TRUE, FALSE))

# join all projected hypervolume inclusions together
inclusion_all <- bind_rows(inclusion_currentTest_gaus, inclusion_currentTest_svm, inclusion_ssp126_gaus, inclusion_ssp126_svm, inclusion_ssp245_gaus, inclusion_ssp245_svm, inclusion_ssp370_gaus, inclusion_ssp370_svm,
                           inclusion_currentTest_gaus_cont, inclusion_currentTest_svm_cont, inclusion_ssp126_gaus_cont, inclusion_ssp126_svm_cont, data.frame(inclusion_ssp245_gaus_cont), data.frame(inclusion_ssp245_svm_cont), data.frame(inclusion_ssp370_gaus_cont), data.frame(inclusion_ssp370_svm_cont)) %>%
  # the historic climate is a test to make sure the hypervolume approach correctly classifies newly fed current climate data as being within the niche. Turn on or off as needed.
  filter(scenario != "historic") %>%
  mutate(notIncluded = ifelse(incl == T, 0, 1),
         scenario = substr(scenario, 1, 8)) 
#saveRDS(inclusion_all, "hypervolumeInclusionAll.rds")
inclusion_all <- readRDS("hypervolumeInclusionAll.rds")

# barplot of mean proportions
barPlot_dat <- inclusion_all %>% 
  group_by(scenario, noveltyType, clusterLetter) %>%
  summarise(meanNotIncluded = mean(notIncluded)) %>%
  mutate(noveltyType = ifelse(noveltyType == "Pyrome-specific", "Climate novel for pyroregion", "Climate novel for continent")) %>%
  mutate(noveltyType = factor(noveltyType, levels = c("Climate novel for pyroregion", "Climate novel for continent")))
barMean <- barPlot_dat %>% group_by(scenario, noveltyType) %>% summarise(meanOfMeans = mean(meanNotIncluded))
barPlot_proportions <- ggplot(barPlot_dat, aes(x = clusterLetter, y = meanNotIncluded)) +
  facet_grid(rows = vars(noveltyType), cols = vars(scenario)) +
  theme_minimal() +
  geom_col(width = 0.5) +
  labs(x = "Pyroregion", y = "Novel climate (proportion of pyroregion)") +
  geom_hline(data = barMean, aes(yintercept = meanOfMeans), linetype = "dashed") 
barPlot_proportions
# summarise mean of mean proportions
inclusion_all %>% group_by(noveltyType, scenario) %>% summarise(meanNotIncluded = mean(notIncluded))

# tally inclusions spatially
mapList <- inclusion_all %>% 
  group_by(noveltyType, scenario, geometry) %>% 
  summarise(notIncludedSum = sum(notIncluded)) %>%
  st_as_sf() %>%
  mutate(noveltyType = ifelse(noveltyType == "Pyrome-specific", "Climate novel to pyroregion", "Climate novel to continent")) %>%
  mutate(noveltyType = factor(noveltyType, levels = c("Climate novel to pyroregion", "Climate novel to continent")))

projMaps <- bind_rows(mapList) %>%
  mutate(notIncludedSum = as.factor(notIncludedSum),
         x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>%
  data.frame() %>%
  ggplot(aes(x, y, fill = notIncludedSum)) +
  theme_map() +
  facet_grid(rows = vars(noveltyType), cols = vars(scenario), switch = "y") +
  geom_tile() +
  scale_fill_manual(values = c("steelblue","darkorange", "darkred"), name = "novel\nclimate")  +
  theme(#strip.text.x = element_blank(), 
        legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
        legend.position = c(0.05,0.5)) 
projMaps
  
#jpeg(file = "hypervol_overlaps_pyromeAndContinent.jpg", width = 12, height = 4.6, res = 600, units = "in")
plot_grid(projMaps, barPlot_proportions, ncol = 2, rel_widths = c(1,0.6), labels = c("A", "C"), label_size = 10, label_fontface = "plain") +
  draw_label("B", x = 0.005, y = 0.48, size = 10)
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 9. QUANTIFY LAND USE OF EACH PYROREGION ---- 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# load land use raster and table of categories
landCover <- rast("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Landscape variables/DLCD MODIS land cover/DLCD_v2-1_MODIS_EVI_9_20100101-20111231.tif")
landCov <- read.csv("C:/Users/cxc/OneDrive - University of Tasmania/UTAS/Fire/Australias fire regimes/Data/Landscape variables/DLCD MODIS land cover/LandCoverCategories.csv") %>%
  dplyr::select(LandUseClass, Common.name, ISO.class.descriptor)

clusterPoints_landuse <- cluster_ras_smoothed %>% 
  # disaggregate pyroregions so we can characterise land use at finer scale 
  disagg(fact = 20) %>% as.data.frame(xy = T) %>% st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>% 
  rename(clusterLetter = clusters) %>%
  # extract land use
  mutate(LandUseClass = extract(landCover, ., ID = F)[,1]) %>%
  left_join(landCov) %>%
  # aggregate some groups
  mutate(Common.name = ifelse(Common.name %in% c("Irrigated cropping", "Irrigated pasture", "Irrigated sugar"), "Irrigated cropping/pasture", Common.name),
         Common.name = ifelse(Common.name %in% c("Rain fed cropping", "Rain fed pasture", "Rain fed sugar"), "Rain fed cropping/pasture", Common.name),
         Common.name = gsub(" and ", "/", Common.name))

# less ~1% of rows have NAs, so drop them in tabulate step
sum(is.na(clusterPoints_landuse$Common.name))/nrow(clusterPoints_landuse)

# tabulate
landUseTabulated <- clusterPoints_landuse %>%
  # drop tiny number of NAs, and drop Lakes/Dams and mines/quarries
  filter(!is.na(Common.name), Common.name != "Lakes/dams",  Common.name != "Mines/Quarries", Common.name != "Salt lakes") %>%
  # it turns out that irrigated cropping only had 3% cover "c", 1% cover of "j", and nothing for all other pyroregions, so lumping it with cropping/pasture
  mutate(Common.name = ifelse(Common.name %in% c("Irrigated cropping/pasture", "Rain fed cropping/pasture"), "Cropping/pasture", Common.name)) %>%
  # reorder categories into thematic groups
  mutate(Common.name = factor(Common.name, levels = c("Closed Forest", "Open Forest", "Woodland", "Open Woodland", "Dense Shrubland", "Open Shrubland", "Scattered shrubs/grasses",
                                                         "Closed Tussock Grassland", "Open Tussock Grassland", "Open Hummock Grassland","Alpine meadows","Wetlands", "Cropping/pasture", "Urban areas"))) %>%
  # add nested higher level
  mutate(LandForm = ifelse(Common.name %in% c("Closed Forest", "Open Forest", "Woodland", "Open Woodland"), "Tree", NA),
         LandForm = ifelse(Common.name %in% c("Dense Shrubland", "Open Shrubland", "Scattered shrubs/grasses"), "Shrub", NA),
         LandForm = ifelse(Common.name %in% c("Closed Tussock Grassland", "Open Tussock Grassland", "Open Hummock Grassland","Alpine meadows"), "Grass", NA),
         LandForm = ifelse(Common.name %in% c("Cropping/pasture", "Urban areas"), "Human", NA)) %>%
  data.frame() %>%
  group_by(clusterLetter) %>%
  mutate(ncells = n()) %>%
  group_by(clusterLetter, Common.name, ncells) %>%
  summarise(n_ = n()) %>%
  mutate(prop = n_/ncells)%>%
  ungroup() %>% 
  mutate(textColour = ifelse(prop < 0.35, "white", "black"))

# check everything sums to one
landUseTabulated %>% group_by(clusterLetter) %>% summarise(sum(prop))

# plot; for visualisation purposes, dropping instances where land use comprises <1%
#jpeg(file = "pyromeLandCover.jpg", width = 5, height = 3, res = 600, units = "in")
ggplot(landUseTabulated %>% filter(prop >= 0.01), aes(clusterLetter, Common.name, fill = prop)) +
  theme_minimal() +
  geom_tile()+
  scale_x_discrete(position = "top") +
  labs(x = "Pyroregion", y = "Land cover") +
  geom_text(aes(label = plyr::round_any(prop, 0.01), colour = textColour), size = 2) +
  scale_colour_manual(values = c("black", "white"), guide = "none") +
  scale_fill_viridis(begin = 0, name = "Proportional land cover") +
  theme(legend.title = element_text(size = 7), legend.text = element_text(size = 7), axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"), legend.position = "bottom") +
  scale_y_discrete(limits=rev)
dev.off()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 10. PLOT HEATMAPS OF SEASONAL PATTERNS OF HOTSPOTS ---- 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

hot1$clusterLetter <- extract(cluster_ras_smoothed, hot1)[["clusters"]]
head(hot1)

# because hotspots are presence-only, we also need a file containing zeroes for each cell
zeroes <- expand.grid(clusterLetter = unique(hot1$clusterLetter), cellID = unique(hot1$cellID), week = 1:52, weeklyHotspotCount = 0)
cellAreas <- hot1 %>% data.frame() %>% group_by(cellID) %>% summarise(cellArea = mean(cellArea)) %>% drop_na()

hot1_df <- hot1 %>% data.frame()

# loop through each year and calculate hotspot densities for each week
weeklyHSList <- list()
for(i in 2001:2021) {
  print(i)
  weeklyHSList[[paste("y_", i, sep = "_")]] <- hot1_df %>% 
    filter(year == i) %>%
    group_by(clusterLetter, cellID, week) %>% 
    # calc weekly hotspot count
    summarise(weeklyHotspotCount = n()) %>% 
    # add in zeroes
    bind_rows(zeroes) %>%
    group_by(clusterLetter, cellID, week) %>% 
    summarise(weeklyHotspotCount = sum(weeklyHotspotCount)) %>%
    left_join(cellAreas) %>%
    filter(!is.na(cellArea)) %>%
    mutate(weeklyHotspotDensity = weeklyHotspotCount/cellArea) %>%
    # summarise from cell level to pyroregion level
    group_by(clusterLetter, week) %>%
    summarise(meanHotspotDensity = mean(weeklyHotspotDensity)) %>%
    filter(!is.na(clusterLetter), week > 0 & week < 53) 
}
weeklyHS <- bind_rows(weeklyHSList, .id = "year") %>%
  mutate(year = as.numeric(gsub("y__", "", year)))
rm(zeroes)

# saveRDS(weeklyHS, "weeklyHS_byPyrome.rds")
weeklyHS <- readRDS("weeklyHS_byPyrome.rds")

weeklyHS_mean <- weeklyHS %>% group_by(clusterLetter, week) %>% summarise(meanHS = mean(meanHotspotDensity))

p_hot_barcode <- ggplot(weeklyHS, aes(week, year, fill = meanHotspotDensity*100)) +
  facet_wrap(~clusterLetter, nrow = 3) +
  theme_minimal() +
  geom_tile() +
  scale_fill_gradient2(
                         high = scales::muted("darkred"), 
                         low = scales::muted("steelblue"),
                         name = "Hs/ha/week ", trans = "sqrt", breaks = c(0,0.02,0.06,0.12)) +
  labs(x = "Week", y = "Year", title = "Temporal pattern of hotspot density")+
  scale_x_continuous(breaks = c(0,15,30,45)) +
  theme(panel.grid.minor = element_blank(), legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey30"), plot.title = element_text(hjust = 0.5, size = 10))
p_hot_barcode

#jpeg(file = "hotspotAnnualPatterns_barcodes.jpg", width = 8, height = 3.7, res = 600, units = "in")
p_hot_barcode
dev.off()