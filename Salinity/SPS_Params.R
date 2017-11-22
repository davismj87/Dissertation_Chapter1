###SPS_Params.R###
###Code used to parameterize soil pore salinity by inundation duration###

###Clearing the working memory and loading necessary packages###

rm(list = ls())
library(ggplot2)
library(nls2)
library(nlstools)
library(propagate)

###Importing the dataset###

	#Column 1 (ID) is sample ID including year, site, and transect
	#Column 2 (Date) was the date each sample was collected
	#Column 3 (Chan_Dist_m) is the Euclidean distance of each sample from the edge of the channel
	#Column 4 (Elev_m) is the surface elevation in meterns NAVD 88 of each sample
	#Column 5 (River_m) is the path distance of each sample to the edge of the river in meters
	#Column 6 (Delta_m) is the path distance of each sample to open water in meters
	#Column 7 (Salinity) is the measured soil pore salinity of each sample

	SPS=read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_SPS.csv")
	head(SPS)

###Preparing the data###

	#First we need to calculate inundation duration using the elevation data
	#Inundation duration is the proportion of time a cell of elevation x will be covered by the tide

		SPS$Inundation = 1/(1 + (0.033*exp(-0.9*1.34))*exp(0.9*SPS$Elev_m))^(1/0.06)
	
		#Checking output

			plot(SPS$Elev_m, SPS$Inundation, pch=19)

	#Now we look at the distribution of the salinity data and run some preliminary plots

		hist(SPS$Salinity)
		
			#The data are fairly normal

		plot(SPS$Inundation, SPS$Salinity)
			abline(lm(SPS$Salinity~SPS$Inundation))

			#Salinity appears to increase with increasing inundation duration

		plot(SPS$Chan_Dist_m, SPS$Salinity)
			abline(lm(SPS$Salinity~SPS$Chan_Dist_m))
			
			#The relationship between salinity and distance from channel edge is unclear

		plot(SPS$River_m, SPS$Salinity)
			abline(lm(SPS$Salinity~SPS$River_m))

			#The salinity increases with increasing distance to the river with an upper limit

		plot(SPS$Delta_m, SPS$Salinity)
			abline(lm(SPS$Salinity~SPS$Delta_m))

			#The salinity decreases with increasing distance to the delta (with a lower limit)

		SPS$Ratio=SPS$River_m/SPS$Delta_m
		plot(SPS$Ratio, SPS$Salinity)
			abline(lm(SPS$Salinity~SPS$Ratio))

			#The ratio between the distance to the river and distance to the delta may
			#also be an effective, simplified way to model soil pore salinity

###Running linear models###

	#First let's run some linear models to parse out significant driving forces

		SPSlm1 = lm(Salinity ~ Inundation, data=SPS)
		SPSlm2 = lm(Salinity ~ River_m, data=SPS)
		SPSlm3 = lm(Salinity ~ Delta_m, data=SPS)
		SPSlm4 = lm(Salinity ~ Chan_Dist_m, data=SPS)
		SPSlm5 = lm(Salinity ~ Ratio, data=SPS)

		AIC(SPSlm1,SPSlm2,SPSlm3,SPSlm4,SPSlm5)
		as.array(c(sigma(SPSlm1),sigma(SPSlm2),sigma(SPSlm3),sigma(SPSlm4),sigma(SPSlm5)))

		summary(SPSlm5)

			#The ratio of mouth to delta distance appears to be the best-fit model
			#This makes sense, because they are mixing sources of fresh and salt water
			#Channel distance is not significant

		SPSlm6 = lm(Salinity ~ Ratio + Inundation, data=SPS)	
		SPSlm7 = lm(Salinity ~ Ratio*Inundation, data=SPS)
		SPSlm8 = lm(Salinity ~ River_m + Inundation, data=SPS)
		SPSlm9 = lm(Salinity ~ River_m*Inundation, data=SPS)
		SPSlm10 = lm(Salinity ~ Delta_m + Inundation, data=SPS)
		SPSlm11 = lm(Salinity ~ Delta_m*Inundation, data=SPS)	
		SPSlm12 = lm(Salinity ~ Delta_m + River_m + Inundation, data=SPS)
		SPSlm13 = lm(Salinity ~ River_m*Delta_m*Inundation, data=SPS)

		AIC(SPSlm6,SPSlm7,SPSlm8,SPSlm9,SPSlm10,SPSlm11,SPSlm12,SPSlm13)
		as.array(c(sigma(SPSlm6),sigma(SPSlm7),sigma(SPSlm8),sigma(SPSlm9),sigma(SPSlm10)
			,sigma(SPSlm11),sigma(SPSlm12),sigma(SPSlm13)))

		summary(SPSlm7)
		summary(SPSlm13)

			#Inundation by dist from mouth by dist from delta is best fit model
			#Even though not all interactive relationships are significant
		
	###GG Plots###

		#River x Inundation

			ggplot(SPS, aes(x = Inundation, y = Salinity, color = River_m))+
			geom_point(size=3)+
			labs(x = 'Inundation', y = 'Salinity')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

			ggplot(SPS, aes(x = River_m, y = Salinity, color = Inundation))+
			geom_point(size=3)+
			labs(x = 'River', y = 'Salinity')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

		#Delta x Inundation

			ggplot(SPS, aes(x = Inundation, y = Salinity, color = Delta_m))+
			geom_point(size=3)+
			labs(x = 'Inundation', y = 'Salinity')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

			ggplot(SPS, aes(x = Delta_m, y = Salinity, color = Inundation))+
			geom_point(size=3)+
			labs(x = 'River', y = 'Salinity')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

		#Ratio x Inundation

			ggplot(SPS, aes(x = Inundation, y = Salinity, color = Ratio))+
			geom_point(size=3)+
			labs(x = 'Inundation', y = 'Salinity')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

			ggplot(SPS, aes(x = Ratio, y = Salinity, color = Inundation))+
			geom_point(size=3)+
			labs(x = 'River', y = 'Salinity')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

		#Model predictions using the "predict" function for the full model

			Pred.full=model.matrix(SPSlm13)
			head(Pred.full)
			Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(SPSlm13, Pred.full)
			plot(Pred.full$Inundation, SPS$Salinity)
				points(Pred.full$Inundation, Pred.full$Predicted, col="red", pch=19)
			plot(Pred.full$River_m, SPS$Salinity)
				points(Pred.full$River_m, Pred.full$Predicted, col="red", pch=19)
			plot(Pred.full$Delta_m, SPS$Salinity)
				points(Pred.full$Delta_m, Pred.full$Predicted, col="red", pch=19)

				#Not bad, but it overestimates salinities close to the river.

		#Model predictions using the "predict" function for the ratio model

			Pred.ratio=model.matrix(SPSlm7)
			head(Pred.ratio)
			Pred.ratio=as.data.frame(Pred.ratio)
			Pred.ratio$Predicted=predict(SPSlm7, Pred.ratio)
			plot(Pred.ratio$Inundation, SPS$Salinity)
				points(Pred.ratio$Inundation, Pred.ratio$Predicted, col="red", pch=19)
			plot(Pred.ratio$Ratio, SPS$Salinity)
				points(Pred.ratio$Ratio, Pred.ratio$Predicted, col="red", pch=19)

				#Not bad, but grossly overestimates salinities at low river/delta ratios.

###Maximum Likelihood Estimation###

	#Exponential Decay

		SPS1 = function(P_C, P_k){
			SPS_1 = P_C*(1 - exp(-P_k*SPS$River_m))
			}
			trySPS1 = SPS1(30, 0.1)

		SPSmle1 = nls(Salinity ~ P_C*(1-exp(-P_k*River_m)), 
			data=SPS, 
			start = list(P_C = 30, P_k = 0.1))
		summary(SPSmle1)
		plot(SPS$River_m, SPS$Salinity)
			points(SPS$River_m, SPS1(coef(SPSmle1)[1], coef(SPSmle1)[2]), pch=19, col="red")
		AIC(SPSlm2, SPSmle1)
		sigma(SPSmle1)

			#The exponential decal model performs considerably better!

	#Exponential Decay plus Inundation Frequency

		SPS2 = function(i1, i2, P_C, P_k){
			SPS_2 = (i1 + i2*SPS$Inundation) + P_C*(1 - exp(-P_k*SPS$River_m))
			}
			trySPS2 = SPS2(0, 5, 30, 0.1)

		SPSmle2 = nls(Salinity ~ i1 + i2*Inundation + P_C*(1-exp(-P_k*River_m)), 
			data=SPS, 
			start = list(i1 = 0, i2 = 5, P_C = 30, P_k = 0.1),
			lower=c(-200, -100, 0, -1),
			upper=c(200, 100, 60, 1),
			algorithm="port")
		summary(SPSmle2)
		plot(SPS$River_m, SPS$Salinity)
			points(SPS$River_m, SPS2(coef(SPSmle2)[1], coef(SPSmle2)[2], coef(SPSmle2)[3], coef(SPSmle2)[4]), pch=19, col="red")
		AIC(SPSlm8, SPSmle1, SPSmle2)
		sigma(SPSmle2)

			#It performs well, but inundation terms are non-significant

	#Exponential Decay times Inundation Frequency

		SPS3 = function(i1, P_C, P_k){
			SPS_3 = i1*SPS$Inundation*(1 - exp(-P_k*SPS$River_m)) + P_C*(1 - exp(-P_k*SPS$River_m))
			}
			trySPS3 = SPS3(50, 20, 0.006)

		SPSmle3 = nls(Salinity ~ i1*Inundation*(1 - exp(-P_k*River_m)) + P_C*(1 - exp(-P_k*River_m)), 
			data=SPS, 
			start = list(i1 = 0, P_C = 30, P_k = 0.1),
			lower=c(-150, 0, -1),
			upper=c(150, 60, 1),
			algorithm="port")
		summary(SPSmle3)
		AIC(SPSlm9, SPSmle1, SPSmle2, SPSmle3)
		sigma(SPSmle3)

		plot(SPS$River_m, SPS$Salinity)
			points(SPS$River_m, SPS3(coef(SPSmle3)[1], coef(SPSmle3)[2], coef(SPSmle3)[3]), pch=19, col="red")
			lm9Predicted=predict(SPSlm9, as.data.frame(model.matrix(SPSlm9)))
			points(SPS$River_m, lm9Predicted, pch=19, col="blue")

	#Exponential Decay times Inundation Frequency plus Delta

		SPS4 = function(d1, d2, i1, P_C, P_k){
			SPS_4 = d1 + d2*SPS$Delta_m + i1*SPS$Inundation*(1 - exp(-P_k*SPS$River_m)) + P_C*(1 - exp(-P_k*SPS$River_m))
			}
			trySPS4 = SPS4(1, 0.01, 3.9, 26, 0.006)

		SPSmle4 = nls(Salinity ~ d1 + d2*Delta_m + i1*Inundation*(1 - exp(-P_k*River_m)) + P_C*(1 - exp(-P_k*River_m)), 
			data=SPS, 
			start = list(d1 = 1, d2 = 0.01, i1 = 0, P_C = 30, P_k = 0.1),
			lower=c(0, -0.1, -150, 0, -1),
			upper=c(20, 0.1, 150, 60, 1),
			algorithm="port")
		summary(SPSmle4)
		AIC(SPSlm13, SPSmle1, SPSmle2, SPSmle3, SPSmle4)
		sigma(SPSmle4)

		plot(SPS$River_m, SPS$Salinity, xlab="River (m)",ylab="Salinity (PSU)")
			points(SPS$River_m, SPS4(coef(SPSmle4)[1], coef(SPSmle4)[2], coef(SPSmle4)[3], coef(SPSmle4)[4], coef(SPSmle4)[5]), pch=19, col="red")
			lm13Predicted=predict(SPSlm13, as.data.frame(model.matrix(SPSlm13)))
			points(SPS$River_m, lm13Predicted, pch=19, col="blue")
		
	#Exponential Decay times Inundation Frequency times Delta

		SPS5 = function(d1, i1, P_C, P_C2, P_k){
			SPS_5 = P_C*(1 - exp(-P_k*SPS$River_m)) + d1*SPS$Delta_m*(1 - exp(-P_k*SPS$River_m)) + i1*SPS$Inundation*(1 - exp(-P_k*SPS$River_m)) + P_C2*SPS$Delta_m*SPS$Inundation*(1 - exp(-P_k*SPS$River_m))
			}
			trySPS5 = SPS5(0.001, 2, 26, 0.01, 0.006)

		SPSmle5 = nls(Salinity ~ P_C*(1 - exp(-P_k*River_m)) + d1*Delta_m*(1 - exp(-P_k*River_m)) + i1*Inundation*(1 - exp(-P_k*River_m)) + P_C2*Delta_m*Inundation*(1 - exp(-P_k*River_m)), 
			data=SPS, 
			start = list(d1 = -0.001, i1 = -2, P_C = 26, P_C2 = 0.01, P_k = 0.006),
			lower=c(-1, -30, -50, -10, -10),
			upper=c(1, 20, 60, 10, 10),
			algorithm="port")
		summary(SPSmle5)
		AIC(SPSlm13, SPSmle1, SPSmle2, SPSmle3, SPSmle4, SPSmle5)
		sigma(SPSmle5)

		plot(SPS$Ratio, SPS$Salinity)
			points(SPS$River_m, SPS5(coef(SPSmle5)[1], coef(SPSmle5)[2], coef(SPSmle5)[3], coef(SPSmle5)[4], coef(SPSmle5)[5]), pch=19, col="red")
			points(SPS$River_m, lm13Predicted, pch=19, col="blue")

	#Exponential Decay Ratio

		SPS6 = function(P_C, P_k){
			SPS_6 = P_C*(1 - exp(-P_k*SPS$Ratio))
			}
			trySPS6 = SPS6(28, 4)

		SPSmle6 = nls(Salinity ~ P_C*(1-exp(-P_k*Ratio)), 
			data=SPS, 
			start = list(P_C = 30, P_k = 0.1))
		summary(SPSmle6)
		plot(SPS$Ratio, SPS$Salinity)
			points(SPS$Ratio, SPS6(coef(SPSmle6)[1], coef(SPSmle6)[2]), pch=19, col="red")
		AIC(SPSmle1,SPSmle6)
		sigma(SPSmle6)

	#Exponential Decay Ratio + Inundation

		SPS7 = function(i1, i2, P_C, P_k){
			SPS_7 = (i1 + i2*SPS$Inundation) + P_C*(1 - exp(-P_k*SPS$Ratio))
			}
			trySPS7 = SPS7(0, 5, 28, 4)

		SPSmle7 = nls(Salinity ~ (i1 + i2*Inundation) + P_C*(1-exp(-P_k*Ratio)), 
			data=SPS, 
			start = list(i1 = 0, i2 = 5, P_C = 30, P_k = 1),
			lower=c(-200, -100, 0, -100),
			upper=c(200, 100, 60, 100),
			algorithm="port")
		summary(SPSmle7)
		plot(SPS$Ratio, SPS$Salinity)
			points(SPS$Ratio, SPS7(coef(SPSmle7)[1], coef(SPSmle7)[2], coef(SPSmle7)[3], coef(SPSmle7)[4]), pch=19, col="red")
		AIC(SPSmle2,SPSmle7)
		sigma(SPSmle7)

	#Exponential Decay Ratio x Inundation

		SPS8 = function(i1, P_C, P_k){
			SPS_8 = i1*SPS$Inundation*(1 - exp(-P_k*SPS$Ratio)) + P_C*(1 - exp(-P_k*SPS$Ratio))
			}
			trySPS8 = SPS8(50, 28, 4)

		SPSmle8 = nls(Salinity ~ i1*Inundation*(1 - exp(-P_k*Ratio)) + P_C*(1 - exp(-P_k*Ratio)), 
			data=SPS, 
			start = list(i1 = 50, P_C = 28, P_k = 4),
			lower=c(-150, 0, -100),
			upper=c(150, 60, 100),
			algorithm="port")
		summary(SPSmle8)
		AIC(SPSmle3,SPSmle8)
		sigma(SPSmle8)

		plot(SPS$Ratio, SPS$Salinity)
			points(SPS$Ratio, SPS8(coef(SPSmle8)[1], coef(SPSmle8)[2], coef(SPSmle8)[3]), pch=19, col="red")

	#So the best-fit model according to AIC values is the exponential decay function of
	# River x Inundation + Delta

###Plotting predicted output###

	#Simulating data

		River_test = seq(0,1750,5)
		Imean = mean(SPS$Inundation, na.rm=TRUE)
		Inundation_test = array(Imean,length(River_test))
		#Inundation_test = rnorm(length(River_test), mean=mean(SPS$Inundation, na.rm=TRUE), sd=sd(SPS$Inundation, na.rm=TRUE))
		Dmean = mean(SPS$Delta, na.rm=TRUE)
		Dsd = sd(SPS$Delta, na.rm=TRUE)
		Delta_test = array(Dmean,length(River_test))
		#Delta_test = rgamma(length(River_test), shape = (Dmean/Dsd)^2, scale = (Dsd^2/Dmean))
		newdat=as.data.frame(cbind(Delta_test, Inundation_test, River_test))
		colnames(newdat) = c("Delta_m","Inundation","River_m")

	#Linear model

		LM_pred = predict(SPSlm13, newdata=newdat, type="response", se.fit=TRUE)
			plot(SPS$River_m, SPS$Salinity)
			plot(River_test, LM_pred$fit, type="l", lwd=2, ylim=c(0,60))
				lines(River_test, LM_pred$fit + LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
				lines(River_test, LM_pred$fit - LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
				points(SPS$River_m, SPS$Salinity, pch=19)

	#Propagate and nlstools#

		Saljack = nlsJack(SPSmle4)
		Saljack$jackCI
		plot(Saljack$jackCI)

	###GG plot##

		MLE_pred = predictNLS(SPSmle4, newdata=data.frame(Delta_m=1500, Inundation=0.1, River_m=River_test), interval="prediction", alpha=0.05)
		MLE_pred_out=as.data.frame(MLE_pred[1])
		head(MLE_pred_out)

		SubSPS = subset(SPS, !Delta_m > 2000)
		SubSPS = subset(SubSPS, !Delta_m < 1000)
		SubSPS = subset(SubSPS, !Inundation > 0.25)
		SubSPS = subset(SubSPS, !Inundation < 0.05)

		ggplot() +
			geom_area(aes(x=River_test, y=MLE_pred_out$summary.Sim.97.5.), fill="lightgray") +
			geom_area(aes(x=River_test, y=MLE_pred_out$summary.Sim.2.5.), fill="white") +
			geom_line(aes(x=River_test, y=MLE_pred_out$summary.Sim.Mean), lwd=1.5) +
			geom_point(aes(x=SPS$River_m, y=SPS$Salinity), pch=19, col="black") +
			geom_point(aes(x=SubSPS$River_m, y=SubSPS$Salinity), pch=19, col="red") +
			scale_x_continuous(limits=c(0,1750)) +
			scale_y_continuous(limits=c(0,50)) +
			xlab("River (m)") +
			ylab("Salinity (PSU)") +
			ggtitle("Inundation = 10%, Delta = 1500 m") +
			theme_classic()+
			theme(plot.title = element_text(hjust = 0.5))


#######################
###CONSOLIDATED CODE###
#######################
		
rm(list = ls())
library(ggplot2)
library(nls2)
library(nlstools)
library(propagate)

#
SPS=read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_SPS.csv")
#

SPS$Inundation = 1/(1 + (0.033*exp(-0.9*1.34))*exp(0.9*SPS$Elev_m))^(1/0.06)
SPS$Ratio=SPS$River_m/SPS$Delta_m

SPSlm1 = lm(Salinity ~ Inundation, data=SPS)
SPSlm2 = lm(Salinity ~ River_m, data=SPS)
SPSlm3 = lm(Salinity ~ Delta_m, data=SPS)
SPSlm4 = lm(Salinity ~ Chan_Dist_m, data=SPS)
SPSlm5 = lm(Salinity ~ Ratio, data=SPS)

AIC(SPSlm1,SPSlm2,SPSlm3,SPSlm4,SPSlm5)
as.array(c(sigma(SPSlm1),sigma(SPSlm2),sigma(SPSlm3),sigma(SPSlm4),sigma(SPSlm5)))
summary(SPSlm5)

SPSlm6 = lm(Salinity ~ Ratio + Inundation, data=SPS)	
SPSlm7 = lm(Salinity ~ Ratio*Inundation, data=SPS)
SPSlm8 = lm(Salinity ~ River_m + Inundation, data=SPS)
SPSlm9 = lm(Salinity ~ River_m*Inundation, data=SPS)
SPSlm10 = lm(Salinity ~ Delta_m + Inundation, data=SPS)
SPSlm11 = lm(Salinity ~ Delta_m*Inundation, data=SPS)	
SPSlm12 = lm(Salinity ~ Delta_m + River_m + Inundation, data=SPS)
SPSlm13 = lm(Salinity ~ River_m*Delta_m*Inundation, data=SPS)

AIC(SPSlm6,SPSlm7,SPSlm8,SPSlm9,SPSlm10,SPSlm11,SPSlm12,SPSlm13)
as.array(c(sigma(SPSlm6),sigma(SPSlm7),sigma(SPSlm8),sigma(SPSlm9),sigma(SPSlm10)
	,sigma(SPSlm11),sigma(SPSlm12),sigma(SPSlm13)))
summary(SPSlm7)
summary(SPSlm13)

ggplot(SPS, aes(x = River_m, y = Salinity, color = Inundation))+
	geom_point(size=3)+
	labs(x = 'River', y = 'Salinity')+
	scale_color_continuous(guide = FALSE)+
	theme_classic()

SPS1 = function(P_C, P_k){
	SPS_1 = P_C*(1 - exp(-P_k*SPS$River_m))
	}
	trySPS1 = SPS1(30, 0.1)

SPSmle1 = nls(Salinity ~ P_C*(1-exp(-P_k*River_m)), 
	data=SPS, 
	start = list(P_C = 30, P_k = 0.1))
	summary(SPSmle1)
	plot(SPS$River_m, SPS$Salinity)
		points(SPS$River_m, SPS1(coef(SPSmle1)[1], coef(SPSmle1)[2]), pch=19, col="red")
	AIC(SPSlm2, SPSmle1)
	sigma(SPSmle1)

SPS2 = function(i1, i2, P_C, P_k){
	SPS_2 = (i1 + i2*SPS$Inundation) + P_C*(1 - exp(-P_k*SPS$River_m))
	}
	trySPS2 = SPS2(0, 5, 30, 0.1)

SPSmle2 = nls(Salinity ~ i1 + i2*Inundation + P_C*(1-exp(-P_k*River_m)), 
	data=SPS, 
	start = list(i1 = 0, i2 = 5, P_C = 30, P_k = 0.1),
	lower=c(-200, -100, 0, -1),
	upper=c(200, 100, 60, 1),
	algorithm="port")
	summary(SPSmle2)
	plot(SPS$River_m, SPS$Salinity)
		points(SPS$River_m, SPS2(coef(SPSmle2)[1], coef(SPSmle2)[2], coef(SPSmle2)[3], coef(SPSmle2)[4]), pch=19, col="red")
	AIC(SPSlm8, SPSmle1, SPSmle2)
	sigma(SPSmle2)

SPS3 = function(i1, P_C, P_k){
	SPS_3 = i1*SPS$Inundation*(1 - exp(-P_k*SPS$River_m)) + P_C*(1 - exp(-P_k*SPS$River_m))
	}
	trySPS3 = SPS3(50, 20, 0.006)

SPSmle3 = nls(Salinity ~ i1*Inundation*(1 - exp(-P_k*River_m)) + P_C*(1 - exp(-P_k*River_m)), 
	data=SPS, 
	start = list(i1 = 0, P_C = 30, P_k = 0.1),
	lower=c(-150, 0, -1),
	upper=c(150, 60, 1),
	algorithm="port")
	summary(SPSmle3)
	AIC(SPSlm9, SPSmle1, SPSmle2, SPSmle3)
	sigma(SPSmle3)

	plot(SPS$River_m, SPS$Salinity)
		points(SPS$River_m, SPS3(coef(SPSmle3)[1], coef(SPSmle3)[2], coef(SPSmle3)[3]), pch=19, col="red")
		lm9Predicted=predict(SPSlm9, as.data.frame(model.matrix(SPSlm9)))
		points(SPS$River_m, lm9Predicted, pch=19, col="blue")

SPS4 = function(d1, d2, i1, P_C, P_k){
	SPS_4 = d1 + d2*SPS$Delta_m + i1*SPS$Inundation*(1 - exp(-P_k*SPS$River_m)) + P_C*(1 - exp(-P_k*SPS$River_m))
	}
	trySPS4 = SPS4(1, 0.01, 3.9, 26, 0.006)

SPSmle4 = nls(Salinity ~ d1 + d2*Delta_m + i1*Inundation*(1 - exp(-P_k*River_m)) + P_C*(1 - exp(-P_k*River_m)), 
	data=SPS, 
	start = list(d1 = 1, d2 = 0.01, i1 = 0, P_C = 30, P_k = 0.1),
	lower=c(0, -0.1, -150, 0, -1),
	upper=c(20, 0.1, 150, 60, 1),
	algorithm="port")
	summary(SPSmle4)
	AIC(SPSlm13, SPSmle1, SPSmle2, SPSmle3, SPSmle4)
	sigma(SPSmle4)

	plot(SPS$River_m, SPS$Salinity, xlab="River (m)",ylab="Salinity (PSU)")
		points(SPS$River_m, SPS4(coef(SPSmle4)[1], coef(SPSmle4)[2], coef(SPSmle4)[3], coef(SPSmle4)[4], coef(SPSmle4)[5]), pch=19, col="red")
		lm13Predicted=predict(SPSlm13, as.data.frame(model.matrix(SPSlm13)))
		points(SPS$River_m, lm13Predicted, pch=19, col="blue")
		
SPS5 = function(d1, i1, P_C, P_C2, P_k){
	SPS_5 = P_C*(1 - exp(-P_k*SPS$River_m)) + d1*SPS$Delta_m*(1 - exp(-P_k*SPS$River_m)) + i1*SPS$Inundation*(1 - exp(-P_k*SPS$River_m)) + P_C2*SPS$Delta_m*SPS$Inundation*(1 - exp(-P_k*SPS$River_m))
	}
	trySPS5 = SPS5(0.001, 2, 26, 0.01, 0.006)

SPSmle5 = nls(Salinity ~ P_C*(1 - exp(-P_k*River_m)) + d1*Delta_m*(1 - exp(-P_k*River_m)) + i1*Inundation*(1 - exp(-P_k*River_m)) + P_C2*Delta_m*Inundation*(1 - exp(-P_k*River_m)), 
	data=SPS, 
	start = list(d1 = -0.001, i1 = -2, P_C = 26, P_C2 = 0.01, P_k = 0.006),
	lower=c(-1, -30, -50, -10, -10),
	upper=c(1, 20, 60, 10, 10),
	algorithm="port")
	summary(SPSmle5)
	AIC(SPSlm13, SPSmle1, SPSmle2, SPSmle3, SPSmle4, SPSmle5)
	sigma(SPSmle5)

	plot(SPS$Ratio, SPS$Salinity)
		points(SPS$River_m, SPS5(coef(SPSmle5)[1], coef(SPSmle5)[2], coef(SPSmle5)[3], coef(SPSmle5)[4], coef(SPSmle5)[5]), pch=19, col="red")
		points(SPS$River_m, lm13Predicted, pch=19, col="blue")

SPS6 = function(P_C, P_k){
	SPS_6 = P_C*(1 - exp(-P_k*SPS$Ratio))
	}
	trySPS6 = SPS6(28, 4)

SPSmle6 = nls(Salinity ~ P_C*(1-exp(-P_k*Ratio)), 
	data=SPS, 
	start = list(P_C = 30, P_k = 0.1))
	summary(SPSmle6)
	plot(SPS$Ratio, SPS$Salinity)
		points(SPS$Ratio, SPS6(coef(SPSmle6)[1], coef(SPSmle6)[2]), pch=19, col="red")
	AIC(SPSmle1,SPSmle6)
	sigma(SPSmle6)

SPS7 = function(i1, i2, P_C, P_k){
	SPS_7 = (i1 + i2*SPS$Inundation) + P_C*(1 - exp(-P_k*SPS$Ratio))
	}
	trySPS7 = SPS7(0, 5, 28, 4)

SPSmle7 = nls(Salinity ~ (i1 + i2*Inundation) + P_C*(1-exp(-P_k*Ratio)), 
	data=SPS, 
	start = list(i1 = 0, i2 = 5, P_C = 30, P_k = 1),
	lower=c(-200, -100, 0, -100),
	upper=c(200, 100, 60, 100),
	algorithm="port")
	summary(SPSmle7)
	plot(SPS$Ratio, SPS$Salinity)
		points(SPS$Ratio, SPS7(coef(SPSmle7)[1], coef(SPSmle7)[2], coef(SPSmle7)[3], coef(SPSmle7)[4]), pch=19, col="red")
	AIC(SPSmle2,SPSmle7)
	sigma(SPSmle7)

SPS8 = function(i1, P_C, P_k){
	SPS_8 = i1*SPS$Inundation*(1 - exp(-P_k*SPS$Ratio)) + P_C*(1 - exp(-P_k*SPS$Ratio))
	}
	trySPS8 = SPS8(50, 28, 4)

SPSmle8 = nls(Salinity ~ i1*Inundation*(1 - exp(-P_k*Ratio)) + P_C*(1 - exp(-P_k*Ratio)), 
	data=SPS, 
	start = list(i1 = 50, P_C = 28, P_k = 4),
	lower=c(-150, 0, -100),
	upper=c(150, 60, 100),
	algorithm="port")
	summary(SPSmle8)
	AIC(SPSmle3,SPSmle8)
	sigma(SPSmle8)

	plot(SPS$Ratio, SPS$Salinity)
		points(SPS$Ratio, SPS8(coef(SPSmle8)[1], coef(SPSmle8)[2], coef(SPSmle8)[3]), pch=19, col="red")

River_test = seq(0,1750,5)
Imean = mean(SPS$Inundation, na.rm=TRUE)
Inundation_test = array(Imean,length(River_test))
Dmean = mean(SPS$Delta, na.rm=TRUE)
Dsd = sd(SPS$Delta, na.rm=TRUE)
Delta_test = array(Dmean,length(River_test))
newdat=as.data.frame(cbind(Delta_test, Inundation_test, River_test))
colnames(newdat) = c("Delta_m","Inundation","River_m")

LM_pred = predict(SPSlm13, newdata=newdat, type="response", se.fit=TRUE)
	plot(SPS$River_m, SPS$Salinity)
	plot(River_test, LM_pred$fit, type="l", lwd=2, ylim=c(0,60))
		lines(River_test, LM_pred$fit + LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
		lines(River_test, LM_pred$fit - LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
		points(SPS$River_m, SPS$Salinity, pch=19)

Saljack = nlsJack(SPSmle4)
Saljack$jackCI
plot(Saljack$jackCI)

MLE_pred = predictNLS(SPSmle4, newdata=data.frame(Delta_m=1500, Inundation=0.1, River_m=River_test), interval="prediction", alpha=0.05)
	MLE_pred_out=as.data.frame(MLE_pred[1])
	head(MLE_pred_out)

	SubSPS = subset(SPS, !Delta_m > 2000)
	SubSPS = subset(SubSPS, !Delta_m < 1000)
	SubSPS = subset(SubSPS, !Inundation > 0.25)
	SubSPS = subset(SubSPS, !Inundation < 0.05)

ggplot() +
	geom_area(aes(x=River_test, y=MLE_pred_out$summary.Sim.97.5.), fill="lightgray") +
	geom_area(aes(x=River_test, y=MLE_pred_out$summary.Sim.2.5.), fill="white") +
	geom_line(aes(x=River_test, y=MLE_pred_out$summary.Sim.Mean), lwd=1.5) +
	geom_point(aes(x=SPS$River_m, y=SPS$Salinity), pch=19, col="black") +
	geom_point(aes(x=SubSPS$River_m, y=SubSPS$Salinity), pch=19, col="red") +
	scale_x_continuous(limits=c(0,1750)) +
	scale_y_continuous(limits=c(0,50)) +
	xlab("River (m)") +
	ylab("Salinity (PSU)") +
	ggtitle("Inundation = 10%, Delta = 1500 m") +
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))





	