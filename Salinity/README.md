# Soil Pore Salinity

This folder contains files pertaining to modeling soil pore salinity using empirical data.

*Data collection*:  
To quantify soil pore salinity, we used data collected during five years of September post-restoration monitoring
surveys, and two years of spring and summer vegetation survey data collected for the Estuary and Salmon Restoration Program (ESRP). 
Salinity values were obtained by drawing out moisture from the top layer of marsh sediment and analyzing the extracted water using a 
handheld salinity refractometer in the lab.   
  
Associated elevation values for each soil pore salinity sample were derived using Real Time Kinematic Global Positioning System (RTK GPS)
measurements that were collected for vegetation quadrats positioned 0, 20, and 40 m from the edge of the transect slough. For missing 
elevation values, we used transect averages for all other years where RTK data were available. We converted elevation values to inundation
frequency for analysis using the previously-derived inundation model. The shortest path distances to the Nisqually River and the Nisqually
Delta were calculated in ArcGIS 10.2 (ESRI, Redlands, California, USA) using the path distance tool. The associated “cost” raster was a 
slope raster derived from the 2011 DEM. Path distances were calculated in ArcGIS for individual pixels across the delta with the Nisqually 
River or Nisqually Delta (digitized and converted to a raster file) as the target cells.
  
*Data analysis*:  
We tested linear and non-linear models to evaluate the relationship between soil pore salinity and inundation, 
path distance to delta, path distance to river, ratio of river distance to delta distance, and Euclidean distance to tidal channel edge. 
We evaluated several model iterations using a maximum likelihood parameterization process with the “lm” and “nls” functions in R 3.4.1 
(R Core Development Team 2017), and selected a best-fit model using Akaike’s Information Criterion (AIC; Burnham and Anderson 2002). 
In cases where multiple best-fit models occurred (ΔAIC < 3), we used model averaging to derive parameter values.  

For non-linear relationships, we used an exponential decay function to produce an upper bound on the model output with the function:  
  
 ![equation_sps](https://user-images.githubusercontent.com/25207964/33147660-63cf6b62-cf7d-11e7-9de4-b3312335e224.GIF)
  
This was most relevant with respect to distance from the river (or ratio of river distance to delta distance), which exhibited 
a distinct positive correlative relationship with an upper asymptote.
  
*Files*:  
SPS_Params.R - Code for parameterizing soil pore salinity model  
SPS_Rmarkdown.md - R markdown tutorial showing how to run code and model output  
SPS_Rmarkdown_files - folder containing output figures for R markdown document  
