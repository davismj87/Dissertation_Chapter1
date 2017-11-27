# Vegetation Stem Density

This folder contains files pertaining to modeling stem density using empirical data.

*Data collection*:  
We conducted vegetation surveys at 30 monitoring sites in September 2009-2015 and at 12 ESRP sites in spring and summer 
2014-2015. At each site, we collected mean percent cover, maximum canopy height, and stem density data for all species found 
within a 0.25 m^2 sampling quadrat. We sampled three quadrats per-site at locations 0, 20, and 40 m from the edge of an adjacent 
tidal slough. We also noted local habitat characteristics such as water depth and overhead canopy cover. We used RTK GPS to obtain 
high-precision estimates of survey site elevations each year. 
  
*Data analysis*:  
We modeled vegetative stem density as a linear or non-linear function of inundation duration and soil pore salinity. Non-linear 
relationships were derived as a Gaussian curve with respect to inundation duration and an exponential decay function with respect 
to soil pore salinity, with equations:  
  
 ![veg_eq1](https://user-images.githubusercontent.com/25207964/33291986-858368a0-d37c-11e7-95bd-21348cbcc5c9.GIF)
  
Predicted values were bounded such that stem density could never be negative, but were also restricted by an upper density 
threshold. We ran a model parameterization and best-fit selection process for stem density using a maximum likelihood parameterization 
process with the “lm” and “nls” functions in R 3.4.1 (R Core Development Team 2017), and selected a best-fit model using Akaike’s 
Information Criterion (AIC; Burnham and Anderson 2002). In cases where multiple best-fit models occurred (ΔAIC < 3), we used model 
averaging to derive parameter values.
  
*Files*:  
Vegetation_Params.R - Code for parameterizing soil pore salinity model  
Veg_Rmarkdown.md - R markdown tutorial showing how to run code and model output  
Veg_Rmarkdown_files - folder containing output figures for R markdown document  
