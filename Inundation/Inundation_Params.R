###Inundation_Params.R###
###Code used to parameterize inundation duration by elevation###

###Clearing the working memory and loading necessary packages###

	rm(list = ls())
	library(bbmle)
	library(nlstools)

###Loading the dataset###

	#Column 1 (DateTime) lists the date and time each data point was collected
	#Column 2 (Year) lists the year each data point was collected
	#Column 3 (Levels) lists the water level elevation (NAVD88)
	#We are interested in columns 2 and 3

		WLL=read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_Inundation.csv")
		head(WLL)

###Preparing the data###

	#The next step is to calculate Inundation Duration (ID%) for each elevation
	#Inundation duration is the proportion of time a cell of elevation x will be covered by the tide
	#Make sure ID% is not calculated below sensor elevation or above high water
	
		Sensor_Elevation = 0.4
		High_Water = max(WLL$Level)
		Ints = as.integer((High_Water-Sensor_Elevation)/0.05)

	#Splitting the dataset by year

		WLL$Year = as.factor(WLL$Year)
			nlevels(WLL$Year)
		Level = split(WLL$Level, WLL$Year)
			Level[6]

	#Creating an array for the data
	
		ID = array(0,dim=c(Ints,1+nlevels(WLL$Year)))

	#Running the for-loop

		for (i in 1:Ints){
			for (j in 1:nlevels(WLL$Year)){
				Threshold = (Sensor_Elevation - 0.05) + i*0.05
				ID[i,1] = Threshold
				ID[i,j+1] = sum(as.data.frame(Level[j]) > ID[i,1])/sum(as.data.frame(Level[j]) > 0)
				}
			}

	#Change headers and restructure

		ID = as.data.frame(ID)
		colnames(ID) = c("Elevation","y2010","y2011","y2012","y2013","y2014","y2015")
		head(ID)

		#Stacking the dataframe

			ID2 = data.frame(ID[1], stack(ID[2:ncol(ID)]))	
				head(ID2)
				plot(ID2$Elevation, ID2$values)	

	#Now we have a new dataframe that lists ID% by elevation and year

		colnames(ID2) = c("El_m","Inund_perc","Yr")
			head(ID2)
			
###Fitting the data to a generalized logistic function###

	#The function

		Inund = function(B, v, Q) {
			I = 1/((1 + Q*exp(-B*El_m))^(1/v))
			}
		
		#Trying starting values for parameters
		
			TestSV=Inund(-0.9, 0.1, 0.017)
			plot(ID2$El_m, TestSV, pch=19)
			points(ID2$El_m, ID2$Inund_perc)

	#Using non-linear least squares to conduct maximum likelihood estimation

		NLS_I=nls(Inund_perc ~ 1/((1 + Q*exp(-B*El_m))^(1/v)), 
			data=ID2, 
			start = list(B=-0.9, v=0.1, Q =0.017),
			lower=c(-10, 0.1, 0.0001),
			upper=c(-0.0001, 1, 1),
			algorithm="port")
		summary(NLS_I)
		confint2(NLS_I)

	#Bootstrapping for parameter confidence intervals

		NLS_I_boot = nlsBoot(NLS_I, niter=500)
			plot(NLS_I_boot, type="boxplot")
			summary(NLS_I_boot)

	#Plotting model output

		plot(ID2$El_m, Inund(coef(NLS_I)[1], coef(NLS_I)[2], coef(NLS_I)[3]), 
			type='p', pch=19, cex=1.5, col="steelblue4", 
			xlab="Elevation (m NAVD88)", ylab="Inundation Duration", ylim=c(0,1), cex.lab=1.35, cex.axis=1.25)
		points(ID2$El_m, ID2$Inund_perc, pch=20, col="black", cex = 1)
		legend(3.25,0.95, c("Predicted","Raw"), pch=c(19,20), col = c("steelblue4","black"))

###Running the model for different MTL###
	#This version allows for adjustments to MTL (sea-level)

	#Create multiple new datasets with incrementally higher water levels

		MTLpresent = WLL$Level
		MTL.2 = WLL$Level + 0.2
		MTL.4 = WLL$Level + 0.4
		MTL.6 = WLL$Level + 0.6
		MTL.8 = WLL$Level + 0.8
		Level_MTL = cbind(MTLpresent,MTL.2,MTL.4,MTL.6,MTL.8)
			head(Level_MTL)

	#Set elevation intervals
	
		Sensor_Elevation = 0.4
		High_Water_MTL = max(Level_MTL)
		Ints_MTL = as.integer((High_Water-Sensor_Elevation)/0.05)

	#Creating an array for the data
	
		ID_MTL = array(0,dim=c(Ints_MTL,6))

	#Running the for-loop

		for (i in 1:Ints_MTL){
			for (j in 1:5){
				Threshold = (Sensor_Elevation - 0.05) + i*0.05
				ID_MTL[i,1] = Threshold
				ID_MTL[i,j+1] = sum(Level_MTL[,j] > ID_MTL[i,1])/sum(Level_MTL[,j] > 0)
				}
			}	

	#Change headers and restructure

		ID_MTL = as.data.frame(ID_MTL)
		colnames(ID_MTL) = c("Elevation","1.34","1.54","1.74","1.94","2.14")
		head(ID_MTL)

		#Removing erroneous values

		for(i in 3:6){
			ID_MTL[,i][ID_MTL[,i] > 0.8] = NA
		}

		#Stacking the dataframe

			ID2_MTL = data.frame(ID_MTL[1], stack(ID_MTL[2:ncol(ID_MTL)]))	
				head(ID2_MTL)
				plot(ID2_MTL$Elevation, ID2_MTL$values)	
				ID2_MTL$ind = as.numeric(as.character(ID2_MTL$ind))

	#Now we have a new dataframe that lists ID% by elevation and year

		colnames(ID2_MTL) = c("El_m","Inund_perc","MTL")
			head(ID2_MTL)

###Fitting the data to a generalized logistic function###

	#The function

		Inund2 = function(B, v, alpha, beta) {
			I = 1/((1 + (alpha*exp(beta*ID2_MTL$MTL))*exp(-B*ID2_MTL$El_m))^(1/v))
			}

		#Trying starting values for parameters
		
			TestSV_MTL=Inund2(-0.9, 0.1, 0.058, -0.9)
			plot(ID2_MTL$El_m, TestSV_MTL, pch=19)
			points(ID2_MTL$El_m, ID2_MTL$Inund_perc)

	#Using non-linear least squares to conduct maximum likelihood estimation

		NLS_MTL=nls(Inund_perc ~ 1/((1 + (alpha*exp(beta*MTL))*exp(-B*El_m))^(1/v)), 
			data=ID2_MTL, 
			start = list(B=-0.9, v=0.1, alpha=0.058, beta=-0.9),
			lower=c(-10, 0.1, 0.0001, -5),
			upper=c(-0.0001, 1, 1, 5),
			algorithm="port")
		summary(NLS_MTL)
		confint2(NLS_MTL)

	#Bootstrapping for parameter confidence intervals

		NLS_MTL_boot = nlsBoot(NLS_MTL, niter=500)
			plot(NLS_MTL_boot, type="boxplot")
			summary(NLS_MTL_boot)

	#Plotting model output

		plot(ID2_MTL$El_m, Inund2(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4]), 
			type='p', pch=19, cex=1.5, col="steelblue4", 
			xlab="Elevation (m NAVD88)", ylab="Inundation Duration", ylim=c(0,1), cex.lab=1.35, cex.axis=1.25)
		points(ID2_MTL$El_m, ID2_MTL$Inund_perc, pch=20, col="black", cex = 1)
		legend(3.25,0.95, c("Predicted","Raw"), pch=c(19,20), col = c("steelblue4","black"))

###PLOT###

	Inund3 = function(B, v, alpha, beta, MTL) {
		I = 1/((1 + (alpha*exp(beta*MTL))*exp(-B*Elev_pred))^(1/v))
		}

	Elev_pred = seq(-1,5,0.01)
	plot(Elev_pred, Inund3(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], MTL=1.34), 
		type='l', lwd=5, col="steelblue4", 
		xlab="Elevation (m NAVD88)", ylab="Inundation Duration", xlim=c(0,5), ylim=c(0,1), cex.lab=1.35, cex.axis=1.25)
		lines(Elev_pred, Inund3(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], MTL=1.54), lwd=5, col="steelblue3")
		lines(Elev_pred, Inund3(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], MTL=1.74), lwd=5, col="steelblue2")
		lines(Elev_pred, Inund3(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], MTL=1.94), lwd=5, col="steelblue1")
		lines(Elev_pred, Inund3(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], MTL=2.14), lwd=5, col="lightsteelblue2")
		points(ID2_MTL$El_m, ID2_MTL$Inund_perc, pch=20, col="black", cex = 1.1)
	legend(3.3,0.99, c("Raw","MTL = 1.34 m", "MTL = 1.54 m", "MTL = 1.74 m", "MTL = 1.94 m", "MTL = 2.14 m"),
		lty = c(3, 1, 1, 1, 1, 1), lwd = 4, col = c("black", "steelblue4", "steelblue3", "steelblue2", "steelblue1", "lightsteelblue2"))


#######################
###CONSOLIDATED CODE###
#######################

rm(list = ls())
library(ggplot2)
library(bbmle)
library(nlstools)
library(nls2)
library(propagate)

###
data=read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_Inundation.csv")
###

Level_MTL = array(0,dim=c(length(data$Level),5))
for(i in 1:ncol(Level_MTL)){
	Level_MTL[,i] = data$Level + (i - 1)*0.2
	}
Level_MTL = as.data.frame(Level_MTL)
colnames(Level_MTL) = c("MTLpresent","MTL.2","MTL.4","MTL.6","MTL.8")

###
Sensor_Elevation = 0.4
Mean_Tidal_Level = 1.34
###
High_Water = max(Level_MTL, na.rm = TRUE)
Ints = as.integer((High_Water-Sensor_Elevation)/0.05)
	
ID_MTL = array(0,dim=c(Ints,6))
for (i in 1:Ints){
	for (j in 1:5){
		Threshold = (Sensor_Elevation - 0.05) + i*0.05
		ID_MTL[i,1] = Threshold
		ID_MTL[i,j+1] = sum(Level_MTL[,j] > ID_MTL[i,1], na.rm=TRUE)/sum(Level_MTL[,j] > 0, na.rm=TRUE)
		}
	}	
ID_MTL = as.data.frame(ID_MTL)
Colnames_MTL = array(0, ncol(ID_MTL))
Colnames_MTL[1] = "Elevation"
for (i in 1:5){
	Colnames_MTL[i+1] = Mean_Tidal_Level + (i - 1)*0.2
	}
colnames(ID_MTL) = Colnames_MTL
for(i in 2:ncol(ID_MTL)){
	ID_MTL[,i][ID_MTL[,i] > 0.8] = NA
	}
ID_MTL = as.data.frame(ID_MTL)
ID2_MTL = data.frame(ID_MTL[1], stack(ID_MTL[2:ncol(ID_MTL)]))	
plot(ID2_MTL$Elevation, ID2_MTL$values)	
ID2_MTL$ind = as.numeric(as.character(ID2_MTL$ind))
colnames(ID2_MTL) = c("El_m","Inund_perc","MTL")
ID2_MTL = na.omit(ID2_MTL)

NLS_MTL=nls(Inund_perc ~ 1/((1 + (mtl_a*exp(mtl_b*MTL))*exp(-B*El_m))^(1/v)), 
	data=ID2_MTL, 
	start = list(B=-0.9, v=0.1, mtl_a=0.058, mtl_b=-0.9),
	lower=c(-10, 	0.05, 	0.00001, 	-5),
	upper=c(-0.00001, 1, 		1, 		5),
	algorithm="port")
	summary(NLS_MTL)
NLS_MTL_boot = nlsBoot(NLS_MTL, niter=500)
	plot(NLS_MTL_boot, type="boxplot")
	summary(NLS_MTL_boot)

Inund = function(B, v, mtl_a, mtl_b, MTL) {
	I = 1/((1 + (mtl_a*exp(mtl_b*MTL))*exp(-B*Elev_pred))^(1/v))
	}
Elev_pred = seq(-1,5,0.01)
plot(Elev_pred, Inund(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], Mean_Tidal_Level), 
	type='l', lwd=5, col="steelblue4", xlab="Elevation (m NAVD88)", ylab="Inundation Duration", 
	xlim=c(0,5), ylim=c(0,1), cex.lab=1.35, cex.axis=1.25)
lines(Elev_pred, Inund(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], Mean_Tidal_Level+0.2), lwd=5, col="steelblue3")
lines(Elev_pred, Inund(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], Mean_Tidal_Level+0.4), lwd=5, col="steelblue2")
lines(Elev_pred, Inund(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], Mean_Tidal_Level+0.6), lwd=5, col="steelblue1")
lines(Elev_pred, Inund(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], Mean_Tidal_Level+0.8), lwd=5, col="lightsteelblue2")
points(ID2_MTL$El_m, ID2_MTL$Inund_perc, pch=20, col="black", cex = 1.1)
legend(3.2,0.99, c("Raw","MTL Present Day", "MTL + 0.2 m", "MTL + 0.4 m", "MTL + 0.6 m", "MTL + 0.8 m"),
	lty = c(3, 1, 1, 1, 1, 1), lwd = 4, col = c("black", "steelblue4", "steelblue3", "steelblue2", "steelblue1", "lightsteelblue2"))

NLS_Pred = predictNLS(NLS_MTL, newdata=data.frame(MTL = Mean_Tidal_Level, El_m = Elev_pred) ,interval="prediction", alpha=0.05)
Inund_pred = Inund(coef(NLS_MTL)[1], coef(NLS_MTL)[2], coef(NLS_MTL)[3], coef(NLS_MTL)[4], Mean_Tidal_Level)	
NLS_Pred2 = as.data.frame(NLS_Pred[1])
Plot_Theme = theme_classic() +
	theme(legend.justification = c(0,0),
	legend.position = c(3,0.75),
	legend.background = element_rect(fill="transparent"),
	legend.text = element_text(size=12, color="gray25",family="sans"),
	text=element_text(size=14),
	axis.title.x = element_text(size=14, color="gray25", face="bold",family="sans"),
	axis.title.y = element_text(size=14, color="gray25", face="bold",family="sans"),
	axis.text.x = element_text(size=12, vjust=0.5, color="gray25",family="sans"),
	axis.text.y = element_text(size=12, color="gray25",family="sans"),
	axis.ticks.x = element_blank(),
	axis.ticks.y=element_blank(),
	strip.text.x = element_text(size = 12, color="gray25",family="sans"),
	strip.background = element_rect(color="transparent"),
	line = element_blank(),
	panel.spacing = unit(0.1, "lines"))
Inundation_plot = ggplot() +
	geom_line(aes(x = Elev_pred, y = Inund_pred), linetype = 1, size = 1) +
	geom_ribbon(aes(x = Elev_pred, ymin = NLS_Pred2$summary.Prop.2.5., ymax = NLS_Pred2$summary.Prop.97.5.), fill="gray20", alpha = 0.15) +
	geom_point(aes(x = ID_MTL[,1], y = ID_MTL[,2]), size=2, color="red") +
	scale_color_manual(name = "Raw Data", values = "red") +
	scale_y_continuous(name="Inundation Duration", limit = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
	scale_x_continuous(name="Elevation (m NAVD88)", limit = c(-1,4), breaks = c(-1,0,1,2,3,4)) +
	Plot_Theme

