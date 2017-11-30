###Elev_Params.R###
###Code used to parameterize vegetation by inundation duration and soil pore salinity###

###Clearing the working memory and loading necessary packages###

rm(list=ls(all=TRUE))
library(ggplot2)
library(bbmle)
library(nlstools)
library(propagate)
library(MASS)

###Importing the dataset###

	#Column 1 (ID) is the SET ID, including pipe, direction, pin, and year interval
	#Column 2 (Elev_m) is the SET elevation in meters NAVD88
	#Column 3 (Nisq_m) is the distance from the river in meters
	#Column 4 (Delta_m) is the distance from open water (delta) in meters
	#Column 5 (Chan_Dist_m) is the distance from the edge of the nearest channel in meters
	#Column 6 (Sed_fall) is the average fall (Sept-Dec) suspended sediment input in metric tons
	#Column 7 (Sed_yr) is the average yearly suspended sediment input in metric tons
	#Column 8 (Stem_dens) is the mean stem density (stems/m2) at each SET
	#Column 9 (delE) is the yearly elevation change (m) standardized by the number of years between measurements
	
	Elev = read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_Elevation.csv")
	head(Elev)

###Preparing the data###

	#First we need to calculate inundation duration using the elevation data
	#Inundation duration is the proportion of time a cell of elevation x will be covered by the tide

		Elev$Inundation = 1/(1 + (0.033*exp(-0.9*1.34))*exp(0.9*Elev$Elev_m))^(1/0.06)
	
		#Checking output

			plot(Elev$Elev_m, Elev$Inundation, pch=19)

	#Now we look at the distribution of the elevation change data and run some preliminary plots

		hist(Elev$delE)
			
			#The data are quite normal
		
		hist(Elev$Inundation)
		plot(Elev$Inundation, Elev$delE, xlim=c(0,0.6))
			abline(lm(Elev$delE ~ Elev$Inundation))
		plot(log(Elev$Inundation+0.01), Elev$delE)
			abline(lm(Elev$delE ~ log(Elev$Inundation+0.01)))

			#May be slight negative relationship
			#Best represented as exponential or gaussian function with asymptote of 0

		hist(Elev$Nisq_m)
		plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
			abline(lm(Elev$delE ~ Elev$Nisq_m))
		plot(Elev$Nisq_m, Elev$Inundation)
			abline(lm(Elev$Inundation ~ Elev$Nisq_m))

			#Slight positive relationship (sediment is being distributed further on the delta)
			#Can be linear


		hist(Elev$Delta_m)
		plot(Elev$Delta_m, Elev$delE, xlim=c(0,2000))
			abline(lm(Elev$delE ~ Elev$Delta_m))
		plot(log(Elev$Delta_m), Elev$delE, xlim=c(5.5,7.5))
			abline(lm(Elev$delE ~ log(Elev$Delta_m)))
		plot(Elev$Delta_m, Elev$Inundation)
			abline(lm(Elev$Inundation ~ Elev$Delta_m))
			summary(lm(Elev$Inundation ~ Elev$Delta_m))

			#Slight positive relationship, but delta covaries significantly with inundation frequency
			#Can be linear

		hist(Elev$Chan_Dist_m)
		plot(Elev$Chan_Dist_m, Elev$delE)
			abline(lm(Elev$delE ~ Elev$Chan_Dist_m))

			#May be slight negative relationship
			#Probably exponential decline

		plot(Elev$Sed_fall, Elev$delE)
			abline(lm(Elev$delE ~ Elev$Sed_fall))

			#Positive relationship
			#Can be linear

		plot(Elev$Stem_dens, Elev$delE)
			abline(lm(Elev$delE ~ Elev$Stem_dens))
		plot(sqrt(Elev$Stem_dens), Elev$delE)
			abline(lm(Elev$delE ~ sqrt(Elev$Stem_dens)))
		
			#May be slight positive relationship
			#Can be linear, likely exponential with a positive asymptote

###Running linear models###

	#First let's run some linear models to parse out significant driving forces

		ELVlm1 = lm(delE ~ Inundation, data=Elev)
		ELVlm2 = lm(delE ~ Nisq_m, data=Elev)
		ELVlm3 = lm(delE ~ Delta_m, data=Elev)
		ELVlm4 = lm(delE ~ Chan_Dist_m, data=Elev)
		ELVlm5 = lm(delE ~ Sed_fall, data=Elev)
		ELVlm6 = lm(delE ~ Stem_dens, data=Elev)
		AIC(lm(delE ~ Sed_yr, data=Elev))

		AIC(ELVlm1,ELVlm2,ELVlm3,ELVlm4,ELVlm5,ELVlm6)
		sigma(ELVlm5)
		summary(ELVlm5)

			#Sediment is the biggest driving force, followed by dist from the river

		ELVlm7 = lm(delE ~ Sed_fall + Inundation, data=Elev)
		ELVlm8 = lm(delE ~ Sed_fall*Inundation, data=Elev)
		ELVlm9 = lm(delE ~ Sed_fall + Nisq_m, data=Elev)
		ELVlm10 = lm(delE ~ Sed_fall*Nisq_m, data=Elev)
		ELVlm11 = lm(delE ~ Sed_fall + Delta_m, data=Elev)
		ELVlm12 = lm(delE ~ Sed_fall*Delta_m, data=Elev)
		ELVlm13 = lm(delE ~ Sed_fall + Chan_Dist_m, data=Elev)
		ELVlm14 = lm(delE ~ Sed_fall*Chan_Dist_m, data=Elev)
		ELVlm15 = lm(delE ~ Sed_fall + Stem_dens, data=Elev)
		ELVlm16 = lm(delE ~ Sed_fall*Stem_dens, data=Elev)

		AIC(ELVlm7,ELVlm8,ELVlm9,ELVlm10,ELVlm11,ELVlm12,ELVlm13,ELVlm14,ELVlm15,ELVlm16)
		summary(ELVlm10)
		sigma(ELVlm10)

			Pred.full=model.matrix(ELVlm10)
				head(Pred.full)
				Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(ELVlm10, Pred.full)
			plot(Elev$Sed_fall, Elev$delE)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="red", pch=19)
			plot(Elev$Nisq_m, Elev$delE)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="red", pch=19)

				#Interactive sediment and distance from the nisqually river
				#This makes sense. More sediment = more accretion. The plume extends further from the river.

		ELVlm17 = lm(delE ~ Sed_fall*Nisq_m + Inundation, data=Elev)
		ELVlm18 = lm(delE ~ Sed_fall*Nisq_m*Inundation, data=Elev)
		ELVlm19 = lm(delE ~ Sed_fall*Nisq_m + Delta_m, data=Elev)
		ELVlm20 = lm(delE ~ Sed_fall*Nisq_m*Delta_m, data=Elev)
		ELVlm21 = lm(delE ~ Sed_fall*Nisq_m + Chan_Dist_m, data=Elev)
		ELVlm22 = lm(delE ~ Sed_fall*Nisq_m*Chan_Dist_m, data=Elev)
		ELVlm23 = lm(delE ~ Sed_fall*Nisq_m + Stem_dens, data=Elev)
		ELVlm24 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens, data=Elev)

		AIC(ELVlm17,ELVlm18,ELVlm19,ELVlm20,ELVlm21,ELVlm22,ELVlm23,ELVlm24)
		summary(ELVlm24)
		sigma(ELVlm18)

			Pred.full=model.matrix(ELVlm24)
				head(Pred.full)
				Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(ELVlm24, Pred.full)
			plot(Elev$Sed_fall, Elev$delE)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="red", pch=19)
			plot(Elev$Nisq_m, Elev$delE)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="red", pch=19)

			#Adding stem density (interactive) improves the model. Delta_m and Inundation also improve the model.

		ELVlm25 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Inundation, data=Elev)
		ELVlm26 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Inundation, data=Elev)
		ELVlm27 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Chan_Dist_m, data=Elev)
		ELVlm28 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Chan_Dist_m, data=Elev)
		ELVlm29 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Delta_m, data=Elev)
		ELVlm30 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Delta_m, data=Elev)
		ELVlm31 = lm(delE ~ Sed_fall*Nisq_m*Inundation + Stem_dens, data=Elev)
		ELVlm32 = lm(delE ~ Sed_fall*Nisq_m*Inundation + Chan_Dist_m, data=Elev)
		ELVlm33 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m, data=Elev)

		AIC(ELVlm25,ELVlm26,ELVlm27,ELVlm28,ELVlm29,ELVlm30,ELVlm31,ELVlm32,ELVlm33)
		sigma(ELVlm33)

			#Best-fit models = *26*, 27, 28, *30*, *33*

			Pred.full=model.matrix(ELVlm33)
				head(Pred.full)
				Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(ELVlm33, Pred.full)
			plot(Elev$Sed_fall, Elev$delE)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="red", pch=19)
			plot(Elev$Nisq_m, Elev$delE)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="red", pch=19)

		ELVlm34 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m + Stem_dens, data=Elev)
		ELVlm35 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m*Stem_dens, data=Elev)
		ELVlm36 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Inundation + Chan_Dist_m, data=Elev)
		ELVlm37 = lm(delE ~ Sed_fall*Nisq_m*Inundation + Chan_Dist_m + Stem_dens, data=Elev)
		ELVlm38 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Chan_Dist_m + Inundation, data=Elev)
		ELVlm39 = lm(delE ~ Sed_fall + Nisq_m + Stem_dens + Chan_Dist_m + Inundation, data=Elev)
		ELVlm40 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m*sqrt(Stem_dens), data=Elev)
		ELVlm41 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m*log(Stem_dens+1), data=Elev)
		ELVlm42 = lm(delE ~ Sed_fall*Nisq_m*Inundation*log(Chan_Dist_m+1)*log(Stem_dens+1), data=Elev)
		ELVlm43 = lm(delE ~ Sed_fall + Nisq_m*Inundation*log(Chan_Dist_m+1)*log(Stem_dens+1), data=Elev)

		summary(ELVlm35)
		AIC(ELVlm34,ELVlm35,ELVlm36,ELVlm37,ELVlm38,ELVlm39,ELVlm40,ELVlm41,ELVlm42,ELVlm43)
		sigma(ELVlm35)
		sigma(ELVlm42)

		ELVlm_step = stepAIC(ELVlm1, scope = list(upper = ~Sed_fall*Nisq_m*Stem_dens*Chan_Dist_m*Inundation, 
			lower = ~1), direction="both")
			ELVlm_step$anova
			AIC(ELVlm_step)
			sigma(ELVlm_step)

			#The best-fit model is the full model w/ some log transformed variables, 
			#but it has LOTS of parameters

			Pred.full=model.matrix(ELVlm35)
				head(Pred.full)
				Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(ELVlm35, Pred.full)
			plot(Elev$Sed_fall, Elev$delE)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="red", pch=19)
			plot(Elev$Nisq_m, Elev$delE)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="red", pch=19)
			plot(Elev$Stem_dens, Elev$delE)
				points(Pred.full$Stem_dens, Pred.full$Predicted, col="red", pch=19)
			plot(Elev$Inundation, Elev$delE)
				points(Pred.full$Inundation, Pred.full$Predicted, col="red", pch=19)

		#Plotting predicted output

			Stem_test = seq(0,3000,5)
			Inundation_test = seq(0,0.75,0.001)
			Nisq_test = seq(0,2000,4)
			Sed_test = seq(0,500,1)
			Chan_test = seq(0,100,0.5)

			newdat = data.frame(Stem_dens = Stem_test,
				Inundation = 0.2, 
				Nisq_m = 750, 
				Sed_fall = 200, 
				Chan_Dist_m = 30)
			FULL_pred_Stem = predict(ELVlm35, newdata=newdat,type="response", se.fit=TRUE)

			plot(Stem_test, FULL_pred_Stem$fit, type="l", lwd=2, ylim=c(-0.1,0.1),
				xlab="Stem dens (/m2)", ylab="Elevation change (m/yr)")
				lines(Stem_test, FULL_pred_Stem$fit + FULL_pred_Stem$se.fit, col="gray",pch=19,lwd=2,lty=2)
				lines(Stem_test, FULL_pred_Stem$fit - FULL_pred_Stem$se.fit, col="gray",pch=19,lwd=2,lty=2)
				points(Elev$Stem_dens, Elev$delE, pch=19)

				#It works all right, but I don't buy some of the relationships

	###GG plots###

			ggplot(Elev, aes(x = Sed_fall, y = delE, color = Nisq_m))+
			geom_point(size=3)+
			labs(x = 'Sediment (metric tons)', y = 'Elevation change (m/yr)')+
			scale_color_continuous()+
			theme_classic()	

			ggplot(Elev, aes(x = Stem_dens, y = delE, color = Nisq_m))+
			geom_point(size=3)+
			labs(x = 'Stem Density (/m2)', y = 'Elevation change (m/yr)')+
			scale_color_continuous()+
			theme_classic()	

			ggplot(Elev, aes(x = Stem_dens, y = delE, color = Inundation))+
			geom_point(size=3)+
			labs(x = 'Stem Density (/m2)', y = 'Inundation duration (%)')+
			scale_color_continuous()+
			theme_classic()

			ggplot(Elev, aes(x = Stem_dens, y = delE, color = Sed_fall))+
			geom_point(size=3)+
			labs(x = 'Stem Density (/m2)', y = 'Sediment (metric tons)')+
			scale_color_continuous()+
			theme_classic()			

###Maximum Likelihood Estimation###

	###Stem Density (m-2)###

		summary(ELVlm6)

		plot(Elev$Stem_dens, Elev$delE)
			points(Pred.full$Stem_dens, Pred.full$Predicted, col="red", pch=19)

		#Exponential Decay (increasing)

			ELV_func1 = function(V_C,V_k){
				El_ch = V_C - V_C*exp(-V_k*Elev$Stem_dens)
				}
			plot(Elev$Stem_dens, ELV_func1(0.02, 0.001))

			ELVmle1=nls(delE ~ V_C - V_C*exp(-V_k*Stem_dens), 
				data=Elev, 
				start = list(V_C = 0.02, V_k = 0.01),
				lower=c(-1, -0.49),
				upper=c(1, 0.49),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle1)
			plot(Elev$Stem_dens, Elev$delE)
				points(Elev$Stem_dens, ELV_func1(coef(ELVmle1)[1], coef(ELVmle1)[2]), pch=19, col="red")

		#Exponential Decay

			ELV_func2 = function(V_C,V_k){
				El_ch = V_C*exp(-V_k*Elev$Stem_dens)
				}
			plot(Elev$Stem_dens, ELV_func2(0.02, 0.001))

			ELVmle2=nls(delE ~ V_C*exp(-V_k*Stem_dens), 
				data=Elev, 
				start = list(V_C = 0.02, V_k = 0.001),
				lower=c(-1, -10),
				upper=c(1, 10),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle2)
			plot(Elev$Stem_dens, Elev$delE)
				points(Elev$Stem_dens, ELV_func2(coef(ELVmle2)[1], coef(ELVmle2)[2]), pch=19, col="red")

		#Logarithmic

			ELV_func3 = function(V_a,V_b){
				El_ch = V_a + V_b*log(Stem_test+1)
				}
			plot(Stem_test, ELV_func3(0.01, 0.001), type="l", lwd=2)

			ELVmle3=nls(delE ~ V_a + V_b*log(Stem_dens+1), 
				data=Elev, 
				start = list(V_a = 0.01, V_b = 0.001),
				lower=c(-1, -1),
				upper=c(1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle3)
			plot(Elev$Stem_dens, Elev$delE)
				lines(Stem_test, ELV_func3(coef(ELVmle3)[1], coef(ELVmle3)[2]), col="red", lwd=3)
				abline(lm(Elev$delE~Elev$Stem_dens), col="blue", lwd=2, lty=2)
			
			AIC(ELVlm6,ELVmle1,ELVmle2,ELVmle3)

			#Logarithmic growth works well as long as +1 is added to stem density to avoid zeros

	###Distance from the Nisqually River###

		summary(ELVlm2)

		plot(Elev$Nisq_m, Elev$delE)
			points(Pred.full$Nisq_m, Pred.full$Predicted, col="red", pch=19)

		#Logistic#

			ELV_func4 = function(N_up, N_steep, N_mid){
				El_ch = N_up/(1 + exp(-N_steep*(Nisq_test-N_mid)))
				}
			plot(Nisq_test, ELV_func4(0.1, 0.005, 500), type="l", lwd=2)

			ELVmle4=nls(delE ~ N_up/(1 + exp(-N_steep*(Nisq_m-N_mid))), 
				data=Elev, 
				start = list(N_up=0.01, N_steep=0.005, N_mid=500),
				lower=c(0, 0, 0),
				upper=c(0.2, 1, 2000),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle4)
			plot(Elev$Nisq_m, Elev$delE)
				lines(Nisq_test, ELV_func4(coef(ELVmle4)[1], coef(ELVmle4)[2], coef(ELVmle4)[3]), col="red", lwd=3)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)

		#Exponential Growth#

			ELV_func5 = function(N_C, N_k){
				El_ch = N_C*exp(N_k*Nisq_test)
				}
			plot(Nisq_test, ELV_func5(0.001,0.001), type="l", lwd=2)

			ELVmle5=nls(delE ~ N_C*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(N_C=0.01, N_k=0.005),
				lower=c(0, 0),
				upper=c(1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle5)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				lines(Nisq_test, ELV_func5(coef(ELVmle5)[1], coef(ELVmle5)[2]), col="red", lwd=3)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)

		#Exponential Decay#

			ELV_func6 = function(N_C, N_k){
				El_ch = N_C*(1 - exp(-N_k*Nisq_test))
				}
			plot(Nisq_test, ELV_func6(0.05,0.001), type="l", lwd=2)

			ELVmle6=nls(delE ~ N_C*(1 - exp(-N_k*Nisq_m)), 
				data=Elev, 
				start = list(N_C=0.01, N_k=0.000005),
				lower=c(0, 0),
				upper=c(1.5, 0.0001),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle6)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				lines(Nisq_test, ELV_func6(coef(ELVmle6)[1], coef(ELVmle6)[2]), col="red", lwd=3)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)

			AIC(ELVlm2,ELVmle4,ELVmle5,ELVmle6)
	
			#Exponential decay and logistic did not fully converge, but exponential growth is still more
			#suitable than linear

	###Suspended Sediment Input###

		summary(ELVlm5)

		plot(Elev$Sed_fall, Elev$delE)
			points(Pred.full$Sed_fall, Pred.full$Predicted, col="red", pch=19)

		#Exponential Decay#

			ELV_func7 = function(S_C, S_k){
				El_ch = S_C*(1 - exp(-S_k*Sed_test))
				}
			plot(Sed_test, ELV_func7(0.05,0.001), type="l", lwd=2)

			ELVmle7=nls(delE ~ S_C*(1 - exp(-S_k*Sed_fall)), 
				data=Elev, 
				start = list(S_C=0.05, S_k=0.001),
				lower=c(0, 0),
				upper=c(2, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle7)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,500))
				lines(Sed_test, ELV_func7(coef(ELVmle7)[1], coef(ELVmle7)[2]), col="red", lwd=3)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)

		#Logarithmic

			ELV_func8 = function(S_a,S_b){
				El_ch = S_a + S_b*log(Sed_test+1)
				}
			plot(Sed_test, ELV_func8(0.01, 0.00001), type="l", lwd=2)

			ELVmle8=nls(delE ~ S_a + S_b*log(Sed_fall+1), 
				data=Elev, 
				start = list(S_a = 0.01, S_b = 0.001),
				lower=c(-1, -1),
				upper=c(1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle8)
			plot(Elev$Sed_fall, Elev$delE)
				lines(Sed_test, ELV_func8(coef(ELVmle8)[1], coef(ELVmle8)[2]), col="red", lwd=3)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)

			AIC(ELVlm5,ELVmle7,ELVmle8)

			#When it comes down to it, they're really similar, linear works fine

	###Stem Density x Suspended Sediment x Distance to the River###

		#Sediment + River#

			summary(ELVlm9)

			ELV_func9 = function(Int, S_a, N_C, N_k){
				El_ch = Int + S_a*Elev$Sed_fall + N_C*exp(N_k*Elev$Nisq_m) 
				}
			plot(Elev$Sed_fall, ELV_func9(0,0.00016,0.0001,0.001))		
			plot(Elev$Nisq_m, ELV_func9(0,0.00016,0.0001,0.001))

			ELVmle9=nls(delE ~ Int + S_a*Sed_fall + N_C*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(Int=0, S_a=0.0001, N_C=0.01, N_k=0.005),
				lower=c(-0.1, 0, 0, 0),
				upper=c(0.1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle9)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func9(coef(ELVmle9)[1], coef(ELVmle9)[2], coef(ELVmle9)[3], 
					coef(ELVmle9)[4]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func9(coef(ELVmle9)[1], coef(ELVmle9)[2], coef(ELVmle9)[3], 
					coef(ELVmle9)[4]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)

			AIC(ELVlm9, ELVmle9)

				#Slightly better

		#Sediment x River#

			summary(ELVlm10)

			ELV_func10 = function(S_a, N_C, N_k){
				El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) 
				}
			plot(Elev$Sed_fall, ELV_func10(0.000016,-0.0001,0.0001))		
			plot(Elev$Nisq_m, ELV_func10(0.000016,-0.0001,0.0001))

			ELVmle10=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(S_a=0.0001, N_C=0.01, N_k=0.005),
				lower=c(-1, -1, -1),
				upper=c(1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle10)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func10(coef(ELVmle10)[1], coef(ELVmle10)[2], coef(ELVmle10)[3]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func10(coef(ELVmle10)[1], coef(ELVmle10)[2], coef(ELVmle10)[3]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)

			AIC(ELVlm10, ELVmle9, ELVmle10)
			sigma(ELVmle10)

			mle10_pred = predictNLS(ELVmle10, newdata=data.frame(Nisq_m = Nisq_test, Sed_fall = 400), interval="prediction", alpha=0.05)
			mle10_pred_out=as.data.frame(mle10_pred[1])
			head(mle10_pred_out)

			SubELV = subset(Elev, !Sed_fall < 300)
			SubELV = subset(SubELV, !Sed_fall > 500)

			ggplot() +
				geom_area(aes(x=Nisq_test, y=(mle10_pred_out$summary.Sim.97.5.)), fill="lightgray") +
				geom_area(aes(x=Nisq_test, y=(mle10_pred_out$summary.Sim.2.5.)), fill="white") +
				geom_line(aes(x=Nisq_test, y=(mle10_pred_out$summary.Sim.Mean)), lwd=1) +
				geom_point(aes(x=Elev$Nisq_m, y=Elev$delE), pch=19, col="black") +
				geom_point(aes(x=SubELV$Nisq_m, y=SubELV$delE), pch=19, col="red") +
				scale_x_continuous(limits=c(0,2000)) +
				scale_y_continuous(limits=c(-0.06,0.06)) +
				xlab("Distance from the River (m)") +
				ylab("delE (m/yr)") +
				ggtitle("Sediment = 400 metric tons") +
				theme_classic()+
				theme(plot.title = element_text(hjust = 0.5))

				#Much better

		#Sediment x River (exponential) + Stems (logarithmic)#

			summary(ELVlm23)

			ELV_func11 = function(Int, V_a, S_a, N_C, N_k){
				El_ch = Int + V_a*log(Elev$Stem_dens+1) + N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) 
				}
			plot(Elev$Sed_fall, ELV_func11(0, 0.0001, 0.0000039,-0.0004,0.0016))		
			plot(Elev$Nisq_m, ELV_func11(0, 0.0001, 0.0000039,-0.0004,0.0016))
			plot(Elev$Stem_dens, ELV_func11(0, 0.0001, 0.0000039,-0.0004,0.0016))

			ELVmle11=nls(delE ~ Int + V_a*log(Stem_dens+1) + N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(Int = 0, V_a = 0.0001, S_a=0.0001, N_C=0.01, N_k=0.005),
				lower=c(-1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle11)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func11(coef(ELVmle11)[1], coef(ELVmle11)[2], coef(ELVmle11)[3]
					, coef(ELVmle11)[4], coef(ELVmle11)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func11(coef(ELVmle11)[1], coef(ELVmle11)[2], coef(ELVmle11)[3]
					, coef(ELVmle11)[4], coef(ELVmle11)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)
			plot(Elev$Stem_dens, Elev$delE, xlim=c(0,2000))
				points(Elev$Stem_dens, ELV_func11(coef(ELVmle11)[1], coef(ELVmle11)[2], coef(ELVmle11)[3]
					, coef(ELVmle11)[4], coef(ELVmle11)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Stem_dens), col="blue", lwd=2, lty=2)

			AIC(ELVlm23, ELVmle11)
			sigma(ELVmle11)

		#Sediment x River (exponential) x Stems (logarithmic)#

			summary(ELVlm24)

			ELV_func12 = function(S_a, V_a, S_b, N_C, N_k){
				El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_a*log(Elev$Stem_dens+1)*exp(N_k*Elev$Nisq_m) + S_b*log(Elev$Stem_dens+1)*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}
			plot(Elev$Sed_fall, ELV_func11(0.01, 0.0001, 0.0000039,-0.0004,0.0016))		
			plot(Elev$Nisq_m, ELV_func11(0.01, 0.0001, 0.0000039,-0.0004,0.0016))
			plot(Elev$Stem_dens, ELV_func11(0.01, 0.0001, 0.0000039,-0.0004,0.0016))

			ELVmle12=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + V_a*log(Stem_dens+1)*exp(N_k*Nisq_m) + S_b*log(Stem_dens+1)*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(S_a = 0.01, V_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016),
				lower=c(-1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle12)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func12(coef(ELVmle12)[1], coef(ELVmle12)[2], coef(ELVmle12)[3]
					, coef(ELVmle12)[4], coef(ELVmle12)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="blue", pch=18, size=0.5)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func12(coef(ELVmle12)[1], coef(ELVmle12)[2], coef(ELVmle12)[3]
					, coef(ELVmle12)[4], coef(ELVmle12)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="blue", pch=18, size=0.5)
			plot(Elev$Stem_dens, Elev$delE, xlim=c(0,2000))
				points(Elev$Stem_dens, ELV_func12(coef(ELVmle12)[1], coef(ELVmle12)[2], coef(ELVmle12)[3]
					, coef(ELVmle12)[4], coef(ELVmle12)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Stem_dens), col="blue", lwd=2, lty=2)
				points(Pred.full$Stem_dens, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			AIC(ELVlm24, ELVmle11, ELVmle12)
			sigma(ELVmle12)

	###Inundation Duration###

		summary(ELVlm1)

		plot(Elev$Inundation, Elev$delE)
			points(Pred.full$Inundation, Pred.full$Predicted, col="red", pch=19)

		#Logistic Growth#

			ELV_func13 = function(I_up, I_steep, I_mid){
				El_ch = I_up/(1 + exp(-I_steep*(Elev$Inundation-I_mid)))
				}
			plot(Elev$Inundation, ELV_func13(0.02,10,0.2))

			ELVmle13=nls(delE ~ I_up/(1 + exp(-I_steep*(Inundation-I_mid))), 
				data=Elev, 
				start = list(I_up=0.02, I_steep=10, I_mid=0.2),
				lower=c(-10, -100, -2),
				upper=c(10, 100, 2),
				algorithm="port")
			summary(ELVmle13)
			plot(Elev$Inundation, Elev$delE)
				points(Elev$Inundation, ELV_func13(coef(ELVmle13)[1], coef(ELVmle13)[2], coef(ELVmle13)[3]), col="red", pch=19)
			AIC(ELVlm1, ELVmle13)

			#Not really useful

		#Gaussian#

			ELV_func14 = function(Ialpha, Ibeta, Igamma){
				El_ch = (Ialpha*exp((-(Elev$Inundation - Ibeta)^2)/(2*(Igamma^2))))
				}
			plot(Elev$Inundation, ELV_func14(-0.001, 0.2, 0.1))

			ELVmle14=nls(delE ~ (Ialpha*exp((-(Inundation - Ibeta)^2)/(2*(Igamma^2)))), 
				data=Elev, 
				start = list(Ialpha = -0.01, Ibeta = 0.2, Igamma = 0.1),
				lower=c(-1, 0, 0),
				upper=c(1, 1, 1),
				algorithm="port")
			summary(ELVmle14)
			plot(Elev$Inundation, Elev$delE)
				points(Elev$Inundation, ELV_func14(coef(ELVmle14)[1], coef(ELVmle14)[2], coef(ELVmle14)[3]), col="red", pch=19)
			AIC(ELVlm1, ELVmle13, ELVmle14)

			#Also doesn't really work...

		#Exponential#

			ELV_func15 = function(I_C, I_k){
				El_ch = I_C*exp(I_k*Elev$Inundation)
				}
			plot(Elev$Inundation, ELV_func15(0.1, 0.1))

			ELVmle15=nls(delE ~ I_C*exp(I_k*Inundation), 
				data=Elev, 
				start = list(I_C=0.1, I_k=0.1),
				lower=c(-1, -10),
				upper=c(1, 10),
				algorithm="port")
			summary(ELVmle15)
			plot(Elev$Inundation, Elev$delE)
				points(Elev$Inundation, ELV_func15(coef(ELVmle15)[1], coef(ELVmle15)[2]), col="red", pch=19)
			AIC(ELVlm1, ELVmle13, ELVmle14, ELVmle15)

			#NO

		#Logarithmic#

			ELV_func16 = function(I_C, I_k){
				El_ch = I_C - I_C*exp(-I_k*Elev$Inundation)
				}
			plot(Elev$Inundation, ELV_func16(0.1, 0.1))

			ELVmle16=nls(delE ~ I_C - I_C*exp(-I_k*Inundation), 
				data=Elev, 
				start = list(I_C=0.1, I_k=0.1),
				lower=c(-1, -10),
				upper=c(1, 10),
				algorithm="port")
			summary(ELVmle16)
			plot(Elev$Inundation, Elev$delE)
				points(Elev$Inundation, ELV_func16(coef(ELVmle16)[1], coef(ELVmle16)[2]), col="red", pch=19)
			AIC(ELVlm1, ELVmle13, ELVmle14, ELVmle15, ELVmle16)

		#Linear does fine. There's too much variation in the data to model it any other way.

	###Inundation x Suspended Sediment x Distance to the River###

		#Sediment x River (exponential) + Inundation#

			summary(ELVlm17)

			ELV_func17 = function(Int, I_a, S_a, N_C, N_k){
				El_ch = Int + I_a*Elev$Inundation + N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) 
				}
			plot(Elev$Sed_fall, ELV_func11(0, 0.0001, 0.0000039,-0.0004,0.0016))		
			plot(Elev$Nisq_m, ELV_func11(0, 0.0001, 0.0000039,-0.0004,0.0016))
			plot(Elev$Inundation, ELV_func11(0, 0.0001, 0.0000039,-0.0004,0.0016))

			ELVmle17=nls(delE ~ Int + I_a*Elev$Inundation + N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(Int = 0, I_a = 0.0001, S_a=0.0001, N_C=0.01, N_k=0.005),
				lower=c(-1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle17)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func17(coef(ELVmle17)[1], coef(ELVmle17)[2], coef(ELVmle17)[3]
					, coef(ELVmle17)[4], coef(ELVmle17)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func17(coef(ELVmle17)[1], coef(ELVmle17)[2], coef(ELVmle17)[3]
					, coef(ELVmle17)[4], coef(ELVmle17)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)
			plot(Elev$Inundation, Elev$delE, xlim=c(0,0.6))
				points(Elev$Inundation, ELV_func17(coef(ELVmle17)[1], coef(ELVmle17)[2], coef(ELVmle17)[3]
					, coef(ELVmle17)[4], coef(ELVmle17)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Inundation), col="blue", lwd=2, lty=2)

			AIC(ELVlm23, ELVmle11, ELVlm17, ELVmle17)
			sigma(ELVmle17)

		#Sediment x River (exponential) x Inundation#

			summary(ELVlm18)

			ELV_func18 = function(S_a, I_a, S_b, N_C, N_k){
				El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + S_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}
			plot(Elev$Sed_fall, ELV_func18(0.01, 0.00001, 0.00000039,-0.0004,0.0016))		
			plot(Elev$Nisq_m, ELV_func18(0.01, 0.00001, 0.00000039,-0.0004,0.0016))
			plot(Elev$Inundation, ELV_func18(0.01, 0.00001, 0.00000039,-0.0004,0.0016))

			ELVmle18=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + I_a*Inundation*exp(N_k*Nisq_m) + S_b*Inundation*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(S_a = 0.01, I_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016),
				lower=c(-1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle18)
			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func18(coef(ELVmle18)[1], coef(ELVmle18)[2], coef(ELVmle18)[3]
					, coef(ELVmle18)[4], coef(ELVmle18)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="blue", pch=18, size=0.5)
			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func18(coef(ELVmle18)[1], coef(ELVmle18)[2], coef(ELVmle18)[3]
					, coef(ELVmle18)[4], coef(ELVmle18)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="blue", pch=18, size=0.5)
			plot(Elev$Inundation, Elev$delE, xlim=c(0,0.6))
				points(Elev$Inundation, ELV_func18(coef(ELVmle18)[1], coef(ELVmle18)[2], coef(ELVmle18)[3]
					, coef(ELVmle18)[4], coef(ELVmle18)[5]), col="red", pch=19)
				abline(lm(Elev$delE~Elev$Inundation), col="blue", lwd=2, lty=2)
				points(Pred.full$Inundation, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			AIC(ELVlm23, ELVlm24, ELVmle11, ELVlm17, ELVmle17, ELVmle18)
			sigma(ELVmle18)

				#Sediment x River Distance (exponential) x Stem Density (logarithmic) is still BFM

		#Sediment x River (exponential) x Stems (logarithmic) + Inundation#

			summary(ELVlm25)

			ELV_func19 = function(Int, I_a, S_a, V_a, S_b, N_C, N_k){
				El_ch = Int + I_a*Elev$Inundation + N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_a*log(Elev$Stem_dens+1)*exp(N_k*Elev$Nisq_m) + S_b*log(Elev$Stem_dens+1)*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}
			plot(Elev$Sed_fall, ELV_func19(0, 0.001, 0.00001, 0.0005, -0.0000013,-0.002,0.0013))		
			plot(Elev$Nisq_m, ELV_func19(0, 0.001, 0.00001, 0.0005, -0.0000013,-0.002,0.0013))
			plot(Elev$Stem_dens, ELV_func19(0, 0.001, 0.00001, 0.0005, -0.0000013,-0.002,0.0013))

			ELVmle19=nls(delE ~ Int + I_a*Inundation + N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + V_a*log(Stem_dens+1)*exp(N_k*Nisq_m) + S_b*log(Stem_dens+1)*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(Int = 0, I_a = 0.001, S_a = 0.00001, V_a = 0.0005, S_b=-0.0000013, N_C=-0.002, N_k=0.0013),
				lower=c(-1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle19)

			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func19(coef(ELVmle19)[1], coef(ELVmle19)[2], coef(ELVmle19)[3],
					coef(ELVmle19)[4], coef(ELVmle19)[5], coef(ELVmle19)[6], coef(ELVmle19)[7]), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func19(coef(ELVmle19)[1], coef(ELVmle19)[2], coef(ELVmle19)[3],
					coef(ELVmle19)[4], coef(ELVmle19)[5], coef(ELVmle19)[6], coef(ELVmle19)[7]), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			plot(Elev$Stem_dens, Elev$delE, xlim=c(0,2000))
				points(Elev$Stem_dens, ELV_func19(coef(ELVmle19)[1], coef(ELVmle19)[2], coef(ELVmle19)[3],
					coef(ELVmle19)[4], coef(ELVmle19)[5], coef(ELVmle19)[6], coef(ELVmle19)[7]), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Stem_dens), col="blue", lwd=2, lty=2)
				points(Pred.full$Stem_dens, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			AIC(ELVlm25, ELVmle19)
			sigma(ELVmle19)

		#Sediment x River (exponential) x Stems (logarithmic) x Inundation#

			summary(ELVlm26)

			ELV_func20 = function(S_a, V_a, S_b, N_C, N_k, I_a, I_b, I_c, I_d){
				El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_a*log(Elev$Stem_dens+1)*exp(N_k*Elev$Nisq_m) + S_b*log(Elev$Stem_dens+1)*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + 
					I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + I_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_c*Elev$Inundation*log(Elev$Stem_dens+1)*exp(N_k*Elev$Nisq_m) + I_d*Elev$Inundation*log(Elev$Stem_dens+1)*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}
			plot(Elev$Sed_fall, ELV_func20(0.00001, 0.0005, -0.0000013,-0.002,0.0013, 0.00001, 0.00001, 0.00001, 0.00001))		
			plot(Elev$Nisq_m, ELV_func20(0.00001, 0.0005, -0.0000013,-0.002,0.0013, 0.00001, 0.00001, 0.00001, 0.00001))
			plot(Elev$Stem_dens, ELV_func20(0.00001, 0.0005, -0.0000013,-0.002,0.0013, 0.00001, 0.00001, 0.00001, 0.00001))

			ELVmle20=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + V_a*log(Stem_dens+1)*exp(N_k*Nisq_m) + S_b*log(Stem_dens+1)*Sed_fall*exp(N_k*Nisq_m) +
					I_a*Inundation*exp(N_k*Nisq_m) + I_b*Inundation*Sed_fall*exp(N_k*Nisq_m) + I_c*Inundation*log(Stem_dens+1)*exp(N_k*Nisq_m) + I_d*Inundation*log(Stem_dens+1)*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(S_a = 0.00001, V_a = 0.0005, S_b=-0.0000013, N_C=-0.002, N_k=0.0013, I_a=0.0001, I_b=0.0001, I_c=0.0001, I_d=0.0001),
				lower=c(-1, -1, -1, -1, -1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle20)

			plot(Elev$Sed_fall, Elev$delE, xlim=c(0,600))
				points(Elev$Sed_fall, ELV_func20(coef(ELVmle20)[1], coef(ELVmle20)[2], coef(ELVmle20)[3],
					coef(ELVmle20)[4], coef(ELVmle20)[5], coef(ELVmle20)[6], coef(ELVmle20)[7], coef(ELVmle20)[8], 
					coef(ELVmle20)[9]), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Sed_fall), col="blue", lwd=2, lty=2)
				points(Pred.full$Sed_fall, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
				points(Elev$Nisq_m, ELV_func20(coef(ELVmle20)[1], coef(ELVmle20)[2], coef(ELVmle20)[3],
					coef(ELVmle20)[4], coef(ELVmle20)[5], coef(ELVmle20)[6], coef(ELVmle20)[7], coef(ELVmle20)[8], 
					coef(ELVmle20)[9]), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Nisq_m), col="blue", lwd=2, lty=2)
				points(Pred.full$Nisq_m, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			plot(Elev$Stem_dens, Elev$delE, xlim=c(0,2000))
				points(Elev$Stem_dens, ELV_func20(coef(ELVmle20)[1], coef(ELVmle20)[2], coef(ELVmle20)[3],
					coef(ELVmle20)[4], coef(ELVmle20)[5], coef(ELVmle20)[6], coef(ELVmle20)[7], coef(ELVmle20)[8], 
					coef(ELVmle20)[9]), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Stem_dens), col="blue", lwd=2, lty=2)
				points(Pred.full$Stem_dens, Pred.full$Predicted, col="blue", pch=18, size=0.5)

			AIC(ELVlm26, ELVmle20)
			sigma(ELVmle20)

				#MLE version is better, but not by a ton

	###Distance from channel edge###

		summary(ELVlm4)

		plot(Elev$Chan_Dist_m, Elev$delE)
			points(Pred.full$Chan_Dist_m, Pred.full$Predicted, col="red", pch=19)

		#Exponential Decay#

			ELV_func21 = function(C_up, C_decay){
				El_ch = C_up*exp(-C_decay*Elev$Chan_Dist_m)
				}
			plot(Elev$Chan_Dist_m, ELV_func21(0.01,-0.001))

			ELVmle21=nls(delE ~ C_up*exp(-C_decay*Chan_Dist_m), 
				data=Elev, 
				start = list(C_up = 0.001, C_decay = 0.001),
				lower=c(-1, -1),
				upper=c(1, 1),
				algorithm="port")
			summary(ELVmle21)
			plot(Elev$Chan_Dist_m, Elev$delE)
				points(Elev$Chan_Dist_m, ELV_func21(coef(ELVmle21)[1], coef(ELVmle21)[2]), ylim=c(0, 0.04), xlim=c(0,100), 
					col="red", pch=19)
				abline(lm(Elev$delE~Elev$Chan_Dist_m), col="blue", lwd=2, lty=2)
			AIC(ELVlm4, ELVmle21)

			#Exponential decay is not any better

		#Sediment x River (exponential) x Inundation x Channel Distance#

			summary(ELVlm33)

			ELV_func22 = function(S_a, I_a, S_b, N_C, N_k, C_a, C_b, C_c, C_d){
				El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + S_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
					C_a*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + C_b*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + C_c*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + C_d*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}
			plot(Elev$Sed_fall, ELV_func22(-0.000002, -0.009, 0.000027,0.0018,0.00014,0.000001,0.000001,0.000001,0.000001))		
			plot(Elev$Nisq_m, ELV_func22(-0.000002, -0.009, 0.000027,0.0018,0.0014,0.000001,0.000001,0.000001,0.000001))
			plot(Elev$Inundation, ELV_func22(-0.000002, -0.009, 0.000027,0.0018,0.000001,0.000001,0.000001,0.000001))

			ELVmle22=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + I_a*Inundation*exp(N_k*Nisq_m) + S_b*Inundation*Sed_fall*exp(N_k*Nisq_m) + 
				C_a*Chan_Dist_m*exp(N_k*Nisq_m) + C_b*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) + C_c*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + C_d*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(S_a = 0.01, I_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016, C_a=0.00001, C_b=0.00001, C_c=0.00001, C_d=0.00001),
				lower=c(-1, -1, -1, -1, -1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle22)

			AIC(ELVlm23, ELVlm24, ELVlm33, ELVmle11, ELVlm17, ELVmle17, ELVmle18, ELVmle22)
			sigma(ELVmle22)

				#The linear model is actually better at this point

		#Sediment x River (exponential) x Inundation x Channel Distance + Stem density (logarithmic)#

			summary(ELVlm34)

			ELV_func23 = function(V_a, V_b, S_a, I_a, S_b, N_C, N_k, C_a, C_b, C_c, C_d){
				El_ch = V_a + V_b*log(Elev$Stem_dens+1) + N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + S_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
					C_a*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + C_b*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + C_c*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + C_d*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}
			plot(Elev$Sed_fall, ELV_func23(0, 0.00001, -0.00000002, -0.009, 0.000027,0.0018,0.00014,0.000001,0.000001,0.000001,0.000001))		
			plot(Elev$Nisq_m, ELV_func23(0, 0.00001, -0.00000002, -0.009, 0.000027,0.0018,0.0014,0.000001,0.000001,0.000001,0.000001))
			plot(Elev$Inundation, ELV_func23(0, 0.00001, -0.00000002, -0.009, 0.000027,0.0018,0.000001,0.000001,0.000001,0.000001))

			ELVmle23=nls(delE ~ V_a + V_b*log(Stem_dens+1) + N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + I_a*Inundation*exp(N_k*Nisq_m) + S_b*Inundation*Sed_fall*exp(N_k*Nisq_m) + 
				C_a*Chan_Dist_m*exp(N_k*Nisq_m) + C_b*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) + C_c*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + C_d*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(V_a = 0.0001, V_b = 0.00001, S_a = 0.01, I_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016, C_a=0.00001, C_b=0.00001, C_c=0.00001, C_d=0.00001),
				lower=c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle23)

			AIC(ELVlm33, ELVlm34, ELVmle22, ELVmle23)
			sigma(ELVmle23)

			#But this one works way better...

		#Sediment x River (exponential) x Inundation x Channel Distance x Stem density (logarithmic)#

			summary(ELVlm35)

			ELV_func24 = function(V_a, V_b, V_c, V_d, V_e, V_f, V_g, V_h, S_a, I_a, S_b, N_C, N_k, C_a, C_b, C_c, C_d){
				El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + S_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
					C_a*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + C_b*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + C_c*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + C_d*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
					V_a*log(Elev$Stem_dens+1)*exp(N_k*Elev$Nisq_m) +  V_b*log(Elev$Stem_dens+1)*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_c*log(Elev$Stem_dens+1)*Elev$Inundation*exp(N_k*Elev$Nisq_m) +
					V_d*log(Elev$Stem_dens+1)*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_e*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + V_f*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
					V_g*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + V_h*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
				}

			ELVmle24=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + I_a*Inundation*exp(N_k*Nisq_m) + S_b*Inundation*Sed_fall*exp(N_k*Nisq_m) +
					C_a*Chan_Dist_m*exp(N_k*Nisq_m) + C_b*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) + C_c*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + C_d*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m) +
					V_a*log(Stem_dens+1)*exp(N_k*Nisq_m) +  V_b*log(Stem_dens+1)*Sed_fall*exp(N_k*Nisq_m) + V_c*log(Stem_dens+1)*Inundation*exp(N_k*Nisq_m) +
					V_d*log(Stem_dens+1)*Inundation*Sed_fall*exp(N_k*Nisq_m) + V_e*log(Stem_dens+1)*Chan_Dist_m*exp(N_k*Nisq_m) + V_f*log(Stem_dens+1)*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) +
					V_g*log(Stem_dens+1)*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + V_h*log(Stem_dens+1)*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m), 
				data=Elev, 
				start = list(V_a = 0.1, V_b = 0.1, V_c = 0.1, V_d = 0.1, V_e = 0.1, V_f = 0.1, V_g = 0.1, V_h = 0.1, S_a = 0.01, I_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016, C_a=0.00001, C_b=0.00001, C_c=0.00001, C_d=0.00001),
				lower=c(-1, -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
				upper=c(1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
				nls.control(maxiter=100),
				algorithm="port")
			summary(ELVmle24)

			AIC(ELVlm33, ELVlm34, ELVlm35, ELVmle22, ELVmle23, ELVmle24)
			sigma(ELVmle24)

				#Not any better than linear!


###Predicted Plots###

	###Changing inundation frequency###

		Inund_test=seq(0,1,0.01)
		Nisq_test=runif(length(Inund_test), 550, 550)
		Chan_test=runif(length(Inund_test), 30, 30)
		Sed_test=rnorm(length(Inund_test), 300, 50)
		Veg_test=rnorm(length(Inund_test), 500, 100)

		Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
			Nisq_m=Nisq_test,  
			Sed_fall=Sed_test, 
			Stem_dens=Veg_test, 
			Chan_Dist_m=Chan_test),
			interval="prediction", se.fit=TRUE)
		Elev_pred_out=as.data.frame(Elev_pred)
		
		plot(delE ~ Inundation, data = Elev, col = "black", pch=20,
			xlab="Inundation Frequency", ylab="Annual Elev Change (m)",
			main="Sed = 300, Stems = 500, Nisq = 550 m, Chan = 30 m", 
			cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,0.6), ylim=c(-0.1, 0.1))
		lines(Inund_test, Elev_pred_out$fit.fit, lwd=3)
		lines(Inund_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
		lines(Inund_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
		points(delE ~ Inundation, data = Elev, col = "black", pch=20)

	###Changing sediment input###

		Sed_test=seq(0, 600, 1)
		Inund_test=runif(length(Sed_test), 0.22, 0.22)
		Nisq_test=runif(length(Sed_test), 550, 550)
		Chan_test=runif(length(Sed_test), 30, 30)
		Veg_test=rnorm(length(Sed_test), 500, 100)

		Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
			Nisq_m=Nisq_test,  
			Sed_fall=Sed_test, 
			Stem_dens=Veg_test, 
			Chan_Dist_m=Chan_test),
			interval="prediction", se.fit=TRUE)
		Elev_pred_out=as.data.frame(Elev_pred)
	
		plot(delE ~ Sed_fall, data = Elev, col = "black", pch=20,
			xlab="Suspended Sediment (metric tons)", ylab="Annual Elev Change (m)",
			main="Inund = 22%, Stems = 500, Nisq = 550 m, Chan = 30 m", 
			cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,600), ylim=c(-0.1, 0.1))
		lines(Sed_test, Elev_pred_out$fit.fit, lwd=3)
		lines(Sed_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
		lines(Sed_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
		points(delE ~ Sed_fall, data = Elev, col = "black", pch=20)

	###Changing Stem Density###

		Veg_test=seq(0, 2000, 8)
		Inund_test=runif(length(Veg_test), 0.22, 0.22)
		Nisq_test=runif(length(Veg_test), 550, 550)
		Chan_test=runif(length(Veg_test), 30, 30)
		Sed_test=rnorm(length(Veg_test), 300, 50)

		Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
			Nisq_m=Nisq_test,  
			Sed_fall=Sed_test, 
			Stem_dens=Veg_test, 
			Chan_Dist_m=Chan_test),
			interval="prediction", se.fit=TRUE)
		Elev_pred_out=as.data.frame(Elev_pred)
	
		plot(delE ~ Stem_dens, data = Elev, col = "black", pch=20,
			xlab="Stem Density (m-2)", ylab="Annual Elev Change (m)",
			main="Inund = 22%, Sediment = 300, Nisq = 550 m, Chan = 30 m", 
			cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,2000), ylim=c(-0.1, 0.1))
		lines(Veg_test, Elev_pred_out$fit.fit, lwd=3)
		lines(Veg_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
		lines(Veg_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
		points(delE ~ Stem_dens, data = Elev, col = "black", pch=20)

	###Changing Channel Distance###

		Chan_test=seq(0, 100, 1)
		Inund_test=runif(length(Chan_test), 0.22, 0.22)
		Nisq_test=runif(length(Chan_test), 550, 550)
		Veg_test=rnorm(length(Chan_test), 500, 100)
		Sed_test=rnorm(length(Chan_test), 300, 50)

		Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
			Nisq_m=Nisq_test,  
			Sed_fall=Sed_test, 
			Stem_dens=Veg_test, 
			Chan_Dist_m=Chan_test),
			interval="prediction", se.fit=TRUE)
		Elev_pred_out=as.data.frame(Elev_pred)
	
		plot(delE ~ Chan_Dist_m, data = Elev, col = "black", pch=20,
			xlab="Channel Dist (m)", ylab="Annual Elev Change (m)",
			main="Inund = 22%, Sediment = 300, Nisq = 550 m, Stems = 500 (m-2)", 
			cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,100), ylim=c(-0.1, 0.1))
		lines(Chan_test, Elev_pred_out$fit.fit, lwd=3)
		lines(Chan_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
		lines(Chan_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
		points(delE ~ Chan_Dist_m, data = Elev, col = "black", pch=20)

#######################
###CONSOLIDATED CODE###
#######################

rm(list=ls(all=TRUE))
library(ggplot2)
library(bbmle)
library(nlstools)
library(propagate)
library(MASS)

#
Elev = read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_Elevation.csv")
#

Elev$Inundation = 1/(1 + (0.033*exp(-0.9*1.34))*exp(0.9*Elev$Elev_m))^(1/0.06)
plot(Elev$Elev_m, Elev$Inundation, pch=19)

hist(Elev$delE)
			
hist(Elev$Inundation)
plot(Elev$Inundation, Elev$delE, xlim=c(0,0.6))
	abline(lm(Elev$delE ~ Elev$Inundation))
plot(log(Elev$Inundation+0.01), Elev$delE)
	abline(lm(Elev$delE ~ log(Elev$Inundation+0.01)))

hist(Elev$Nisq_m)
plot(Elev$Nisq_m, Elev$delE, xlim=c(0,2000))
	abline(lm(Elev$delE ~ Elev$Nisq_m))
plot(Elev$Nisq_m, Elev$Inundation)
	abline(lm(Elev$Inundation ~ Elev$Nisq_m))

hist(Elev$Delta_m)
plot(Elev$Delta_m, Elev$delE, xlim=c(0,2000))
	abline(lm(Elev$delE ~ Elev$Delta_m))
plot(log(Elev$Delta_m), Elev$delE, xlim=c(5.5,7.5))
	abline(lm(Elev$delE ~ log(Elev$Delta_m)))
plot(Elev$Delta_m, Elev$Inundation)
	abline(lm(Elev$Inundation ~ Elev$Delta_m))
	summary(lm(Elev$Inundation ~ Elev$Delta_m))

hist(Elev$Chan_Dist_m)
plot(Elev$Chan_Dist_m, Elev$delE)
	abline(lm(Elev$delE ~ Elev$Chan_Dist_m))

plot(Elev$Sed_fall, Elev$delE)
	abline(lm(Elev$delE ~ Elev$Sed_fall))

plot(Elev$Stem_dens, Elev$delE)
	abline(lm(Elev$delE ~ Elev$Stem_dens))
plot(sqrt(Elev$Stem_dens), Elev$delE)
	abline(lm(Elev$delE ~ sqrt(Elev$Stem_dens)))
		
ELVlm1 = lm(delE ~ Inundation, data=Elev)
ELVlm2 = lm(delE ~ Nisq_m, data=Elev)
ELVlm3 = lm(delE ~ Delta_m, data=Elev)
ELVlm4 = lm(delE ~ Chan_Dist_m, data=Elev)
ELVlm5 = lm(delE ~ Sed_fall, data=Elev)
ELVlm6 = lm(delE ~ Stem_dens, data=Elev)
AIC(lm(delE ~ Sed_yr, data=Elev))

AIC(ELVlm1,ELVlm2,ELVlm3,ELVlm4,ELVlm5,ELVlm6)
summary(ELVlm5)

ELVlm7 = lm(delE ~ Sed_fall + Inundation, data=Elev)
ELVlm8 = lm(delE ~ Sed_fall*Inundation, data=Elev)
ELVlm9 = lm(delE ~ Sed_fall + Nisq_m, data=Elev)
ELVlm10 = lm(delE ~ Sed_fall*Nisq_m, data=Elev)
ELVlm11 = lm(delE ~ Sed_fall + Delta_m, data=Elev)
ELVlm12 = lm(delE ~ Sed_fall*Delta_m, data=Elev)
ELVlm13 = lm(delE ~ Sed_fall + Chan_Dist_m, data=Elev)
ELVlm14 = lm(delE ~ Sed_fall*Chan_Dist_m, data=Elev)
ELVlm15 = lm(delE ~ Sed_fall + Stem_dens, data=Elev)
ELVlm16 = lm(delE ~ Sed_fall*Stem_dens, data=Elev)

AIC(ELVlm7,ELVlm8,ELVlm9,ELVlm10,ELVlm11,ELVlm12,ELVlm13,ELVlm14,ELVlm15,ELVlm16)
summary(ELVlm10)
sigma(ELVlm10)

ELVlm17 = lm(delE ~ Sed_fall*Nisq_m + Inundation, data=Elev)
ELVlm18 = lm(delE ~ Sed_fall*Nisq_m*Inundation, data=Elev)
ELVlm19 = lm(delE ~ Sed_fall*Nisq_m + Delta_m, data=Elev)
ELVlm20 = lm(delE ~ Sed_fall*Nisq_m*Delta_m, data=Elev)
ELVlm21 = lm(delE ~ Sed_fall*Nisq_m + Chan_Dist_m, data=Elev)
ELVlm22 = lm(delE ~ Sed_fall*Nisq_m*Chan_Dist_m, data=Elev)
ELVlm23 = lm(delE ~ Sed_fall*Nisq_m + Stem_dens, data=Elev)
ELVlm24 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens, data=Elev)

AIC(ELVlm17,ELVlm18,ELVlm19,ELVlm20,ELVlm21,ELVlm22,ELVlm23,ELVlm24)
summary(ELVlm24)
sigma(ELVlm24)

ELVlm25 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Inundation, data=Elev)
ELVlm26 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Inundation, data=Elev)
ELVlm27 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Chan_Dist_m, data=Elev)
ELVlm28 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Chan_Dist_m, data=Elev)
ELVlm29 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Delta_m, data=Elev)
ELVlm30 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Delta_m, data=Elev)
ELVlm31 = lm(delE ~ Sed_fall*Nisq_m*Inundation + Stem_dens, data=Elev)
ELVlm32 = lm(delE ~ Sed_fall*Nisq_m*Inundation + Chan_Dist_m, data=Elev)
ELVlm33 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m, data=Elev)

AIC(ELVlm25,ELVlm26,ELVlm27,ELVlm28,ELVlm29,ELVlm30,ELVlm31,ELVlm32,ELVlm33)

ELVlm34 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m + Stem_dens, data=Elev)
ELVlm35 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m*Stem_dens, data=Elev)
ELVlm36 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens*Inundation + Chan_Dist_m, data=Elev)
ELVlm37 = lm(delE ~ Sed_fall*Nisq_m*Inundation + Chan_Dist_m + Stem_dens, data=Elev)
ELVlm38 = lm(delE ~ Sed_fall*Nisq_m*Stem_dens + Chan_Dist_m + Inundation, data=Elev)
ELVlm39 = lm(delE ~ Sed_fall + Nisq_m + Stem_dens + Chan_Dist_m + Inundation, data=Elev)
ELVlm40 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m*sqrt(Stem_dens), data=Elev)
ELVlm41 = lm(delE ~ Sed_fall*Nisq_m*Inundation*Chan_Dist_m*log(Stem_dens+1), data=Elev)
ELVlm42 = lm(delE ~ Sed_fall*Nisq_m*Inundation*log(Chan_Dist_m+1)*log(Stem_dens+1), data=Elev)

summary(ELVlm35)
AIC(ELVlm34,ELVlm35,ELVlm36,ELVlm37,ELVlm38,ELVlm39,ELVlm40,ELVlm41,ELVlm42)
sigma(ELVlm35)
sigma(ELVlm40)

ELVlm_step = stepAIC(ELVlm1, scope = list(upper = ~Sed_fall*Nisq_m*Stem_dens*Chan_Dist_m*Inundation, 
	lower = ~1), direction="both")
	ELVlm_step$anova
	AIC(ELVlm_step)
	sigma(ELVlm_step)

Pred.full=model.matrix(ELVlm35)
	head(Pred.full)
	Pred.full=as.data.frame(Pred.full)
Pred.full$Predicted=predict(ELVlm35, Pred.full)
	plot(Elev$Sed_fall, Elev$delE)
		points(Pred.full$Sed_fall, Pred.full$Predicted, col="red", pch=19)
	plot(Elev$Nisq_m, Elev$delE)
		points(Pred.full$Nisq_m, Pred.full$Predicted, col="red", pch=19)
	plot(Elev$Stem_dens, Elev$delE)
		points(Pred.full$Stem_dens, Pred.full$Predicted, col="red", pch=19)
	plot(Elev$Inundation, Elev$delE)
		points(Pred.full$Inundation, Pred.full$Predicted, col="red", pch=19)

Stem_test = seq(0,3000,5)
Inundation_test = seq(0,0.75,0.001)
Nisq_test = seq(0,2000,4)
Sed_test = seq(0,500,1)
Chan_test = seq(0,100,0.5)

newdat = data.frame(Stem_dens = Stem_test,
	Inundation = 0.2, 
	Nisq_m = 750, 
	Sed_fall = 200, 
	Chan_Dist_m = 30)
FULL_pred_Stem = predict(ELVlm35, newdata=newdat,type="response", se.fit=TRUE)

plot(Stem_test, FULL_pred_Stem$fit, type="l", lwd=2, ylim=c(-0.1,0.1),
	xlab="Stem dens (/m2)", ylab="Elevation change (m/yr)")
	lines(Stem_test, FULL_pred_Stem$fit + FULL_pred_Stem$se.fit, col="gray",pch=19,lwd=2,lty=2)
	lines(Stem_test, FULL_pred_Stem$fit - FULL_pred_Stem$se.fit, col="gray",pch=19,lwd=2,lty=2)
	points(Elev$Stem_dens, Elev$delE, pch=19)

ggplot(Elev, aes(x = Sed_fall, y = delE, color = Nisq_m))+
	geom_point(size=3)+
	labs(x = 'Sediment (metric tons)', y = 'Elevation change (m/yr)')+
	scale_color_continuous()+
	theme_classic()	

ggplot(Elev, aes(x = Stem_dens, y = delE, color = Nisq_m))+
	geom_point(size=3)+
	labs(x = 'Stem Density (/m2)', y = 'Elevation change (m/yr)')+
	scale_color_continuous()+
	theme_classic()	

ggplot(Elev, aes(x = Stem_dens, y = delE, color = Inundation))+
	geom_point(size=3)+
	labs(x = 'Stem Density (/m2)', y = 'Inundation duration (%)')+
	scale_color_continuous()+
	theme_classic()

ggplot(Elev, aes(x = Stem_dens, y = delE, color = Sed_fall))+
	geom_point(size=3)+
	labs(x = 'Stem Density (/m2)', y = 'Sediment (metric tons)')+
	scale_color_continuous()+
	theme_classic()			


ELV_func23 = function(V_a, V_b, S_a, I_a, S_b, N_C, N_k, C_a, C_b, C_c, C_d){
	El_ch = V_a + V_b*log(Elev$Stem_dens+1) + N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + S_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
		C_a*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + C_b*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + C_c*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + C_d*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
	}

ELVmle23=nls(delE ~ V_a + V_b*log(Stem_dens+1) + N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + I_a*Inundation*exp(N_k*Nisq_m) + S_b*Inundation*Sed_fall*exp(N_k*Nisq_m) + 
	C_a*Chan_Dist_m*exp(N_k*Nisq_m) + C_b*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) + C_c*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + C_d*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m), 
	data=Elev, 
	start = list(V_a = 0.0001, V_b = 0.00001, S_a = 0.01, I_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016, C_a=0.00001, C_b=0.00001, C_c=0.00001, C_d=0.00001),
	lower=c(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
	upper=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
	nls.control(maxiter=100),
	algorithm="port")
	summary(ELVmle23)

AIC(ELVlm34, ELVmle23)
sigma(ELVmle23)

ELV_func24 = function(V_a, V_b, V_c, V_d, V_e, V_f, V_g, V_h, S_a, I_a, S_b, N_C, N_k, C_a, C_b, C_c, C_d){
	El_ch = N_C*exp(N_k*Elev$Nisq_m) + S_a*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + I_a*Elev$Inundation*exp(N_k*Elev$Nisq_m) + S_b*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
	C_a*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + C_b*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + C_c*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + C_d*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
	V_a*log(Elev$Stem_dens+1)*exp(N_k*Elev$Nisq_m) +  V_b*log(Elev$Stem_dens+1)*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_c*log(Elev$Stem_dens+1)*Elev$Inundation*exp(N_k*Elev$Nisq_m) +
	V_d*log(Elev$Stem_dens+1)*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) + V_e*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*exp(N_k*Elev$Nisq_m) + V_f*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*Elev$Sed_fall*exp(N_k*Elev$Nisq_m) +
	V_g*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*Elev$Inundation*exp(N_k*Elev$Nisq_m) + V_h*log(Elev$Stem_dens+1)*Elev$Chan_Dist_m*Elev$Inundation*Elev$Sed_fall*exp(N_k*Elev$Nisq_m)
	}

ELVmle24=nls(delE ~ N_C*exp(N_k*Nisq_m) + S_a*Sed_fall*exp(N_k*Nisq_m) + I_a*Inundation*exp(N_k*Nisq_m) + S_b*Inundation*Sed_fall*exp(N_k*Nisq_m) +
	C_a*Chan_Dist_m*exp(N_k*Nisq_m) + C_b*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) + C_c*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + C_d*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m) +
	V_a*log(Stem_dens+1)*exp(N_k*Nisq_m) +  V_b*log(Stem_dens+1)*Sed_fall*exp(N_k*Nisq_m) + V_c*log(Stem_dens+1)*Inundation*exp(N_k*Nisq_m) +
	V_d*log(Stem_dens+1)*Inundation*Sed_fall*exp(N_k*Nisq_m) + V_e*log(Stem_dens+1)*Chan_Dist_m*exp(N_k*Nisq_m) + V_f*log(Stem_dens+1)*Chan_Dist_m*Sed_fall*exp(N_k*Nisq_m) +
	V_g*log(Stem_dens+1)*Chan_Dist_m*Inundation*exp(N_k*Nisq_m) + V_h*log(Stem_dens+1)*Chan_Dist_m*Inundation*Sed_fall*exp(N_k*Nisq_m), 
	data=Elev, 
	start = list(V_a = 0.1, V_b = 0.1, V_c = 0.1, V_d = 0.1, V_e = 0.1, V_f = 0.1, V_g = 0.1, V_h = 0.1, S_a = 0.01, I_a = 0.0001, S_b=0.0000039, N_C=-0.0004, N_k=0.0016, C_a=0.00001, C_b=0.00001, C_c=0.00001, C_d=0.00001),
	lower=c(-1, -1, -1, -1, -1, -1,-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
	upper=c(1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
	nls.control(maxiter=100),
	algorithm="port")
	summary(ELVmle24)

AIC(ELVlm34, ELVlm35, ELVmle23, ELVmle24)
sigma(ELVmle24)

Inund_test=seq(0,1,0.01)
Nisq_test=runif(length(Inund_test), 550, 550)
Chan_test=runif(length(Inund_test), 30, 30)
Sed_test=rnorm(length(Inund_test), 300, 50)
Veg_test=rnorm(length(Inund_test), 500, 100)

Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
	Nisq_m=Nisq_test,  
	Sed_fall=Sed_test, 
	Stem_dens=Veg_test, 
	Chan_Dist_m=Chan_test),
	interval="prediction", se.fit=TRUE)
Elev_pred_out=as.data.frame(Elev_pred)
		
plot(delE ~ Inundation, data = Elev, col = "black", pch=20,
	xlab="Inundation Frequency", ylab="Annual Elev Change (m)",
	main="Sed = 300, Stems = 500, Nisq = 550 m, Chan = 30 m", 
	cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,0.6), ylim=c(-0.1, 0.1))
	lines(Inund_test, Elev_pred_out$fit.fit, lwd=3)
	lines(Inund_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
	lines(Inund_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
	points(delE ~ Inundation, data = Elev, col = "black", pch=20)

Sed_test=seq(0, 600, 1)
	Inund_test=runif(length(Sed_test), 0.22, 0.22)
	Nisq_test=runif(length(Sed_test), 550, 550)
	Chan_test=runif(length(Sed_test), 30, 30)
	Veg_test=rnorm(length(Sed_test), 500, 100)

Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
	Nisq_m=Nisq_test,  
	Sed_fall=Sed_test, 
	Stem_dens=Veg_test, 
	Chan_Dist_m=Chan_test),
	interval="prediction", se.fit=TRUE)
Elev_pred_out=as.data.frame(Elev_pred)
	
plot(delE ~ Sed_fall, data = Elev, col = "black", pch=20,
	xlab="Suspended Sediment (metric tons)", ylab="Annual Elev Change (m)",
	main="Inund = 22%, Stems = 500, Nisq = 550 m, Chan = 30 m", 
	cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,600), ylim=c(-0.1, 0.1))
lines(Sed_test, Elev_pred_out$fit.fit, lwd=3)
lines(Sed_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
lines(Sed_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
points(delE ~ Sed_fall, data = Elev, col = "black", pch=20)

Veg_test=seq(0, 2000, 8)
Inund_test=runif(length(Veg_test), 0.22, 0.22)
Nisq_test=runif(length(Veg_test), 550, 550)
Chan_test=runif(length(Veg_test), 30, 30)
Sed_test=rnorm(length(Veg_test), 300, 50)

Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
	Nisq_m=Nisq_test,  
	Sed_fall=Sed_test, 
	Stem_dens=Veg_test, 
	Chan_Dist_m=Chan_test),
	interval="prediction", se.fit=TRUE)
Elev_pred_out=as.data.frame(Elev_pred)
	
plot(delE ~ Stem_dens, data = Elev, col = "black", pch=20,
	xlab="Stem Density (m-2)", ylab="Annual Elev Change (m)",
	main="Inund = 22%, Sediment = 300, Nisq = 550 m, Chan = 30 m", 
	cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,2000), ylim=c(-0.1, 0.1))
lines(Veg_test, Elev_pred_out$fit.fit, lwd=3)
lines(Veg_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
lines(Veg_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
points(delE ~ Stem_dens, data = Elev, col = "black", pch=20)

Chan_test=seq(0, 100, 1)
Inund_test=runif(length(Chan_test), 0.22, 0.22)
Nisq_test=runif(length(Chan_test), 550, 550)
Veg_test=rnorm(length(Chan_test), 500, 100)
Sed_test=rnorm(length(Chan_test), 300, 50)

Elev_pred = predict(ELVlm42, newdata=list(Inundation=Inund_test, 
	Nisq_m=Nisq_test,  
	Sed_fall=Sed_test, 
	Stem_dens=Veg_test, 
	Chan_Dist_m=Chan_test),
	interval="prediction", se.fit=TRUE)
	Elev_pred_out=as.data.frame(Elev_pred)
	
plot(delE ~ Chan_Dist_m, data = Elev, col = "black", pch=20,
	xlab="Stem Density (m-2)", ylab="Annual Elev Change (m)",
	main="Inund = 22%, Sediment = 300, Nisq = 550 m, Stems = 500 (m-2)", 
	cex.main=1.25, cex.lab=1.25, cex.axis=1.2, xlim=c(0,100), ylim=c(-0.1, 0.1))
lines(Chan_test, Elev_pred_out$fit.fit, lwd=3)
lines(Chan_test, Elev_pred_out$fit.upr, lwd=2, lty=2, col="darkgray")			
lines(Chan_test, Elev_pred_out$fit.lwr, lwd=2, lty=2, col="darkgray")
points(delE ~ Chan_Dist_m, data = Elev, col = "black", pch=20)

