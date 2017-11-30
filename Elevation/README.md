# Elevation Change

This folder contains files pertaining to modeling elevation change using empirical data.

*Data collection*:  
We measured elevation change through time using rod Surface Elevation Tables (SETs; Cahoon et al. 2002) located at 14 restoring 
and historic marsh sites. Rod SETs are leveling devices used to measure relative elevation change at a precision of up to 1 mm. 
The SET benchmark is driven into the marsh substrate to resistance, deep enough to withstand processes such as deposition, erosion, 
compaction and freezing and thawing events, making this a consistent and minimally disruptive method for measuring surface elevation. 
We took an initial SET measurement immediately after the restoration in September 2009, and collected yearly measurements every spring 
thereafter through April 2017. Feldspar horizon markers or sediment pins at each site were used to measure corresponding sediment 
accretion rates. 
  
*Data analysis*:  
We modeled yearly elevation change as a function of inundation duration, stem density, distance to the Nisqually River, distance to 
tidal channel edge, and seasonal sediment input. We opted to limit our analysis to linear configurations due to the large number of 
tested parameters, recognizing that an unbounded model might demonstrate greater uncertainty in non-target areas; however, we 
log-transformed some predictor variables to improve model fit. We used a maximum likelihood model parameterization and best-fit model 
selection process with the “lm” and “nls” functions in R 3.4.1 (Burnham and Anderson 2002; R Core Development Team 2017) to determine 
which variables had a significant effect on marsh accretion, and how these variables were related to elevation change. 
  
*Files*:  
Elev_Params.R - Code for parameterizing soil pore salinity model  
Elev_Rmarkdown.md - R markdown tutorial showing how to run code and model output  
Elev_Rmarkdown_files - folder containing output figures for R markdown document  
