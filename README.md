# Code to delineate Australia's pyroregions and their risk of shifting under climate change.

This code accompanies the paper, "Pyrogeography in flux: reorganisation of Australian fire regimes in a hotter world" by Calum X. Cunningham, Grant J. Williamson, Rachael H. Nolan, Lina Teckentrup, Matthias M. Boer, and David M.J.S. Bowman

There are three key scripts involved. 
1. "Pyroregionalisation functions.R" provides custom functions that are sourced by the other scripts. 
2. "Metrics of fire regimes_v2.R" turns satellite datasets into spatiotemporal metrics of fire regimes. Because the satellite datasets are large, we do not store them here. This script produces a stack of geotiffs ("fireMetricStack_20230821.tiff").
3. "Cluster analysis delineating fire regimes v4.R" uses the stack of geotiffs (from 2) to (i) delineate Australia's pyroregions, (ii) characterise their climatic niches, and (iii) evaluate whether projected climates occur within the climatic niches.
