# Inundation Duration

This folder contains raw data and code for modelling inundation duration by elevation (m NAVD88).

*Data collection*: To quantify inundation duration based on surface elevation, we used barometrically-compensated 
water level data collected from a Solinst® LTC level logger installed in historic marsh. Data were logged continuously 
at 15 minute intervals between February 2010 and January 2015. We converted water level (m) to elevation (m NAVD88) 
using Real-Time Kinematic Global Positioning System (RTK GPS) validated sensor height. We also used the dataset to derive 
a local tidal datum for the Nisqually River Delta (MHHW = 3.1 m, MHW = 2.8 m, MTL = 1.4 m, MLW = -0.1 m, MLLW = -1 m). 
Because the sensor was positioned above MLW (0.35 m), low tide values were derived using the mean proportional distance 
between MTL and MLW, and MTL and MLLW for ten other nearby sites with active NOAA tide gauges.

*Data analysis*: We calculated the relationship between inundation duration and surface elevation as the total proportion 
of time a cell with a specific elevation x was tidally inundated. We used the “nls” function in R 3.4.1 (R Core 
Development Team 2017) to parameterize a generalized logistic decay function using a maximum likelihood estimation procedure 
on normally-distributed elevation data:  

![equation1](https://user-images.githubusercontent.com/25207964/30882735-ce1dc414-a2be-11e7-8825-f1aa680b0221.PNG) 
  
The function had a lower asymptote (A) of 0 and an upper asymptote (K) of 1 to limit inundation to proportional values. 
To model change in I_x with MTL, we quantified changes in parameters B, v, and Q as MTL increased from 1.4 to 2.5, 
reflecting projected sea-level rise (SLR) estimates of up to a meter or more over the coming decades. This was best reflected in an 
exponential relationship with parameter Q:  

![equation2](https://user-images.githubusercontent.com/25207964/30882738-d04bd5aa-a2be-11e7-907d-822753f3bed0.PNG)

where A_mtl and B_mtl define change in Q with increasing MTL.  

*Files*:
