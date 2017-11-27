###Vegetation_Params.R###
###Code used to parameterize vegetation by inundation duration and soil pore salinity###

###Clearing the working memory and loading necessary packages###

rm(list = ls())
library(ggplot2)
library(bbmle)
library(nlstools)
library(propagate)

###Importing the dataset###

	#Column 1 (ID) is the quadrat ID, including year, month, and site
	#Column 2 (Density) is the stem density of the 0.25 x 0.25 m quadrat
	#Column 3 (Chan_Dist_m) is the distance from the edge of the nearest channel in meters
	#Column 4 (Elev_m) is the quadrat elevation in meters NAVD88
	#Column 5 (Salinity) is the measured soil pore salinity at the quadrat

	#file.choose()
	Veg=read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_Vegetation.csv")
	head(Veg)

###Preparing the data###

	#First we need to calculate inundation duration using the elevation data
	#Inundation duration is the proportion of time a cell of elevation x will be covered by the tide

		Veg$Inundation = 1/(1 + (0.033*exp(-0.9*1.34))*exp(0.9*Veg$Elev_m))^(1/0.06)
	
		#Checking output

			plot(Veg$Elev_m, Veg$Inundation, pch=19)

	#Now we look at the distribution of the stem density data and run some preliminary plots

		hist(Veg$Density)
			
			#The data are zero-skewed poisson (count) data

		hist(log(Veg$Density + 1))
		hist(sqrt(Veg$Density))

			#Log-transforming or square root transformation helps, except for the zero-skewedness

		plot(Veg$Inundation, Veg$Density)
			abline(lm(Veg$Density~Veg$Inundation))

			#Stem density appears to follow a gaussian function with respect to inundation

		plot(Veg$Salinity, Veg$Density)
			abline(lm(Veg$Density~Veg$Salinity))

			#Stem density appears to decrease with increasing salinity

		plot(Veg$Chan_Dist_m, Veg$Density)
			abline(lm(Veg$Density~Veg$Chan_Dist_m))	
			
			#Probably no relationship here	

###Running linear models###

	#First let's run some linear models to parse out significant driving forces

		VEGlm1 = lm(Density ~ Inundation, data=Veg)
		VEGlm2 = lm(Density ~ Salinity, data=Veg)
		VEGlm3 = lm(Density ~ Inundation + Salinity, data=Veg)
		VEGlm4 = lm(Density ~ Inundation*Salinity, data=Veg)
		VEGlm5 = lm(log(Density+1) ~ Inundation, data=Veg)
		VEGlm6 = lm(log(Density+1) ~ Salinity, data=Veg)
		VEGlm7 = lm(log(Density+1) ~ Inundation + Salinity, data=Veg)
		VEGlm8 = lm(log(Density+1) ~ Inundation*Salinity, data=Veg)
		VEGlm9 = lm(sqrt(Density) ~ Inundation, data=Veg)
		VEGlm10 = lm(sqrt(Density) ~ Salinity, data=Veg)
		VEGlm11 = lm(sqrt(Density) ~ Inundation + Salinity, data=Veg)
		VEGlm12 = lm(sqrt(Density) ~ Inundation*Salinity, data=Veg)

		AIC(VEGlm1,VEGlm2,VEGlm3,VEGlm4,VEGlm5,VEGlm6,VEGlm7,VEGlm8,VEGlm9,VEGlm10,VEGlm11,VEGlm12)
		as.array(c(sigma(VEGlm1),sigma(VEGlm2),sigma(VEGlm3),sigma(VEGlm4),sigma(VEGlm5),
			sigma(VEGlm6),sigma(VEGlm7),sigma(VEGlm8)))

		summary(VEGlm7)

		#The additive model of Inundation and Salinity appears to be the best-fit model,
		#although salinity is not significant.

	###GG Plots###

		#Salinity x Inundation

			ggplot(Veg, aes(x = Inundation, y = Density, color = Salinity))+
			geom_point(size=3)+
			labs(x = 'Inundation', y = 'Stem Density (m-2)')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

			ggplot(Veg, aes(x = Salinity, y = Density, color = Inundation))+
			geom_point(size=3)+
			labs(x = 'Salinity', y = 'Stem Density (m-2)')+
			scale_color_continuous(guide = FALSE)+
			theme_classic()

		#Model predictions using the "predict" function for the best-fit model

			Pred.full=model.matrix(VEGlm4)
				head(Pred.full)
				Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(VEGlm4, Pred.full)
			plot(Veg$Inundation, Veg$Density)
				points(Pred.full$Inundation, Pred.full$Predicted, col="red", pch=19)
			plot(Veg$Salinity, Veg$Density)
				points(Pred.full$Salinity, Pred.full$Predicted, col="red", pch=19)

			Pred.full=model.matrix(VEGlm7)
				head(Pred.full)
				Pred.full=as.data.frame(Pred.full)
			Pred.full$Predicted=predict(VEGlm7, Pred.full)
			plot(Veg$Inundation, Veg$Density)
				points(Pred.full$Inundation, exp(Pred.full$Predicted), col="red", pch=19)
			plot(Veg$Salinity, Veg$Density)
				points(Pred.full$Salinity, exp(Pred.full$Predicted), col="red", pch=19)

				#Not bad, but it doesn't capture the Gaussian relationship we're looking for
				#at the lower end of inundation values

###Maximum Likelihood Estimation###

	#Gaussian Inundation

		#Testing the model

		VEG_func1 = function(alpha, beta, gamma){
			Stems = alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}
		plot(Veg$Inundation, VEG_func1(100, 0.15, 0.2))

			#where alpha = 1/sd*sqrt(2pi), beta = mean, and gamma = sd
			#Looks pretty good

		VEGmle1=nls(Density ~ alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(alpha=100, beta=0.15, gamma=0.2),
			lower=c(0, 0, 0),
			upper=c(10000, 1, 1),
			algorithm="port")
		summary(VEGmle1)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func1(coef(VEGmle1)[1], coef(VEGmle1)[2], coef(VEGmle1)[3]), pch=19, cex=0.75, col="red")
		AIC(VEGlm1, VEGmle1)
		sigma(VEGmle1)

		#The Gaussian model performs slightly better

	#Gaussian Additive

		VEG_func2 = function(s1, s2, alpha, beta, gamma){
			Stems = s1 + s2*Veg$Salinity + alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}
		plot(Veg$Inundation, VEG_func2(0, 2, 100, 0.15, 0.2))

		VEGmle2=nls(Density ~ s1 + s2*Salinity + alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(s1 = 0, s2 = 2, alpha=100, beta=0.15, gamma=0.2),
			lower=c(-100, 0, 0, 0, 0),
			upper=c(100, 1000, 10000, 1, 1),
			algorithm="port")
		summary(VEGmle2)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func2(coef(VEGmle2)[1], coef(VEGmle2)[2], coef(VEGmle2)[3], coef(VEGmle2)[4], coef(VEGmle2)[5]), pch=19, cex=0.75, col="red")
		AIC(VEGlm3, VEGmle2)
		sigma(VEGmle2)

		#Salinity parameters don't really converge

	#Gaussian Interactive

		VEG_func3 = function(alpha, s1, beta, gamma){
			Stems = (alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2))))
			}
		plot(Veg$Inundation, VEG_func3(50, 100, 0.12, 0.11))

		VEGmle3=nls(Density ~ (alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2)))), 
			data=Veg, 
			start = list(alpha = 5, s1 = 1, beta=0.15, gamma=0.2),
			lower=c(0, 0, 0, 0),
			upper=c(10000, 100, 1, 1),
			algorithm="port")
		summary(VEGmle3)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func3(coef(VEGmle3)[1], coef(VEGmle3)[2], coef(VEGmle3)[3], coef(VEGmle3)[4]), pch=19, cex=0.75, col="red")
		AIC(VEGlm4, VEGmle3)
		sigma(VEGmle3)

		#Salinity parameters don't really converge

	#Gaussian Interactive (Simplified)

		VEG_func4 = function(alpha, beta, gamma){
			Stems = alpha*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}
		plot(Veg$Inundation, VEG_func4(10, 0.12, 0.11))

		VEGmle4=nls(Density ~ alpha*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(alpha = 10, beta=0.15, gamma=0.2),
			lower=c(0, 0, 0),
			upper=c(1000, 1, 1),
			algorithm="port")
		summary(VEGmle4)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func4(coef(VEGmle4)[1], coef(VEGmle4)[2], coef(VEGmle4)[3]), pch=19, cex=0.75, col="red")
		AIC(VEGlm4, VEGmle4)
		sigma(VEGmle4)
		sigma(VEGlm4)

		#Adding salinity parameter captures range of variation, but still not good enough.

	#Gaussian Inundation Sqrt-Transformed

		#Testing the model

		VEG_func5 = function(alpha, beta, gamma){
			Stems = alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}

		VEGmle5=nls(sqrt(Density) ~ alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(alpha=10, beta=0.15, gamma=0.2),
			lower=c(0, 0, 0),
			upper=c(10000, 1, 1),
			algorithm="port")
		summary(VEGmle5)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func5(coef(VEGmle5)[1], coef(VEGmle5)[2], coef(VEGmle5)[3])^2, pch=19, cex=0.75, col="red")
			points(Veg$Inundation, VEG_func1(coef(VEGmle1)[1], coef(VEGmle1)[2], coef(VEGmle1)[3]), pch=19, cex=0.75, col="blue")
		AIC(VEGmle1, VEGmle5)

		#Square rooting the stem density gives a better result

	#Gaussian Additive Sqrt-Transformed

		VEG_func6 = function(s1, s2, alpha, beta, gamma){
			Stems = s1 + s2*Veg$Salinity + alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}

		VEGmle6=nls(sqrt(Density) ~ s1 + s2*Salinity + alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(s1 = 0, s2 = 2, alpha=100, beta=0.15, gamma=0.2),
			lower=c(-100, 0, 0, 0, 0),
			upper=c(100, 1000, 10000, 1, 1),
			algorithm="port")
		summary(VEGmle6)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func6(coef(VEGmle6)[1], coef(VEGmle6)[2], coef(VEGmle6)[3], coef(VEGmle6)[4], coef(VEGmle6)[5])^2, pch=19, cex=0.75, col="red")
			points(Veg$Inundation, VEG_func2(coef(VEGmle2)[1], coef(VEGmle2)[2], coef(VEGmle2)[3], coef(VEGmle2)[4], coef(VEGmle2)[5]), pch=19, cex=0.75, col="blue")
		AIC(VEGmle2, VEGmle6)

		#Salinity parameters don't really converge

	#Gaussian Interactive Sqrt-Transformed

		VEG_func7 = function(alpha, s1, beta, gamma){
			Stems = (alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2))))
			}

		VEGmle7=nls(sqrt(Density) ~ (alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2)))), 
			data=Veg, 
			start = list(alpha = 5, s1 = 1, beta=0.15, gamma=0.2),
			lower=c(0, 0, 0, 0),
			upper=c(500, 100, 1, 1),
			algorithm="port")
		summary(VEGmle7)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func7(coef(VEGmle7)[1], coef(VEGmle7)[2], coef(VEGmle7)[3], coef(VEGmle7)[4])^2, pch=19, cex=0.75, col="red")
			points(Veg$Inundation, VEG_func3(coef(VEGmle3)[1], coef(VEGmle3)[2], coef(VEGmle3)[3], coef(VEGmle7)[4]), pch=19, cex=0.75, col="blue")
		AIC(VEGmle3,VEGmle7)

		#Salinity parameters don't really converge

	#Gaussian Interactive (Simplified) Sqrt-Transformed

		VEG_func8 = function(alpha, beta, gamma){
			Stems = alpha*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}

		VEGmle8=nls(sqrt(Density) ~ alpha*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(alpha = 10, beta=0.15, gamma=0.2),
			lower=c(0, 0, 0),
			upper=c(1000, 1, 1),
			algorithm="port")
		summary(VEGmle8)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func8(coef(VEGmle8)[1], coef(VEGmle8)[2], coef(VEGmle8)[3])^2, pch=19, cex=0.75, col="red")
			points(Veg$Inundation, VEG_func4(coef(VEGmle4)[1], coef(VEGmle4)[2], coef(VEGmle4)[3]), pch=19, cex=0.75, col="blue")
		AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8)

		#Adding salinity parameter captures range of variation, but still not good enough.

	#Salinity Exponential Decay

		VEG_func9 = function(N,lambda){
			Stems = N*exp(-lambda*Veg$Salinity)
			}

		VEGmle9=nls(Density ~ N*exp(-lambda*Salinity), 
			data=Veg, 
			start = list(N = 3000, lambda=0),
			lower=c(0, -5),
			upper=c(10000, 5),
			algorithm="port")
		summary(VEGmle9)
		plot(Veg$Salinity, Veg$Density, pch=19)
			points(Veg$Salinity, VEG_func9(coef(VEGmle9)[1], coef(VEGmle9)[2]), pch=19, cex=0.75, col="red")
		AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9)
		sigma(VEGmle9)

	#Salinity Exponential Decay + Inundation Gaussian

		VEG_func10 = function(N,lambda, alpha, beta, gamma){
			Stems = N*exp(-lambda*Veg$Salinity) + alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}

		VEGmle10=nls(Density ~ N*exp(-lambda*Salinity) + alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(N = 3000, lambda=0, alpha=100, beta=0.12, gamma=0.11),
			lower=c(0, -5, 0, 0, 0),
			upper=c(10000, 5, 10000, 100, 100),
			algorithm="port")
		summary(VEGmle10)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func10(coef(VEGmle10)[1], coef(VEGmle10)[2], coef(VEGmle10)[3], coef(VEGmle10)[4], coef(VEGmle10)[5]), pch=19, cex=0.75, col="red")
		AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10)
		sigma(VEGmle10)

	#Salinity Exponential Decay x Inundation Gaussian

		VEG_func11 = function(N, lambda, beta, gamma){
			Stems = N*exp(-lambda*Veg$Salinity)*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
			}
		plot(Veg$Inundation, VEG_func11(1000, 0.05, 0.12, 0.11))

		VEGmle11=nls(Density ~ N*exp(-lambda*Salinity)*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(N = 3000, lambda=0, beta=0.12, gamma=0.11),
			lower=c(0, -5, 0, 0),
			upper=c(10000, 5, 100, 100),
			algorithm="port")
		summary(VEGmle11)
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func11(coef(VEGmle11)[1], coef(VEGmle11)[2], coef(VEGmle11)[3], coef(VEGmle11)[4]), pch=19, cex=0.75, col="red")
			points(Veg$Inundation, VEG_func8(coef(VEGmle8)[1], coef(VEGmle8)[2], coef(VEGmle8)[3])^2, pch=19, cex=0.75, col="blue")
		plot(Veg$Salinity, Veg$Density, pch=19)
			points(Veg$Salinity, VEG_func11(coef(VEGmle11)[1], coef(VEGmle11)[2], coef(VEGmle11)[3], coef(VEGmle11)[4]), pch=19, cex=0.75, col="red")
			points(Veg$Salinity, VEG_func8(coef(VEGmle8)[1], coef(VEGmle8)[2], coef(VEGmle8)[3])^2, pch=19, cex=0.75, col="blue")
		AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10,VEGmle11)
		sigma(VEGmle11)

		#THIS IS THE CORRECT MODEL TO USE

	#Gaussian Salinity

		VEG_func12 = function(alpha, beta, gamma){
			Stems = alpha*exp((-(Veg$Salinity - beta)^2)/(2*(gamma^2)))
			}
		plot(Veg$Salinity, VEG_func12(100,3, 3))

		VEGmle12=nls(Density ~ alpha*exp((-(Salinity - beta)^2)/(2*(gamma^2))), 
			data=Veg, 
			start = list(alpha = 100, beta=3, gamma=3),
			lower=c(0, 0, 0),
			upper=c(10000, 1000, 1000),
			algorithm="port")
		summary(VEGmle12)
		plot(Veg$Salinity, Veg$Density, pch=19)
			points(Veg$Salinity, VEG_func12(coef(VEGmle12)[1], coef(VEGmle12)[2], coef(VEGmle12)[3]), pch=19, cex=0.75, col="red")
		AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10,VEGmle11,VEGmle12)

	#Gaussian Salinity x Gaussian Inundation

		VEG_func13 = function(alpha, beta_s, gamma_s, beta_i, gamma_i){
			Stems = alpha*exp((-(Veg$Salinity - beta_s)^2)/(2*(gamma_s^2)))*exp((-(Veg$Inundation - beta_i)^2)/(2*(gamma_i^2)))
			}
		plot(Veg$Inundation, VEG_func13(500,6, 4,0.1, 0.1))

		VEGmle13=nls(Density ~ alpha*exp((-(Salinity - beta_s)^2)/(2*(gamma_s^2)))*exp((-(Inundation - beta_i)^2)/(2*(gamma_i^2))), 
			data=Veg, 
			start = list(alpha = 100, beta_s=3, gamma_s=3, beta_i=0.1, gamma_i=0.1),
			lower=c(0, 0, 0, 0, 0),
			upper=c(10000, 1000, 1000, 100, 100),
			algorithm="port")
		summary(VEGmle13)
		plot(Veg$Salinity, Veg$Density, pch=19)
			points(Veg$Salinity, VEG_func13(coef(VEGmle13)[1], coef(VEGmle13)[2], coef(VEGmle13)[3], coef(VEGmle13)[4], coef(VEGmle13)[5]), pch=19, cex=0.75, col="red")
		plot(Veg$Inundation, Veg$Density, pch=19)
			points(Veg$Inundation, VEG_func13(coef(VEGmle13)[1], coef(VEGmle13)[2], coef(VEGmle13)[3], coef(VEGmle13)[4], coef(VEGmle13)[5]), pch=19, cex=0.75, col="red")
		AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10,VEGmle11,VEGmle12, VEGmle13)

###Plotting Predicted Output###

	#Simulating data

		Inundation_test = seq(0,0.75,0.005)
		Salinity_test = array(mean(Veg$Salinity, na.rm=TRUE),length(Inundation_test))
		newdat=as.data.frame(cbind(Inundation_test, Salinity_test))
		colnames(newdat) = c("Inundation","Salinity")

	#Linear model

		LM_pred = predict(VEGlm7, newdata=newdat, type="response", se.fit=TRUE)
		plot(Veg$Inundation, Veg$Density)
		plot(Inundation_test, LM_pred$fit, type="l", lwd=2, ylim=c(0,20))
			lines(Inundation_test, LM_pred$fit + LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
			lines(Inundation_test, LM_pred$fit - LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
			points(Veg$Inundation, log(Veg$Density+1), pch=19)

	#Propagate and nlstools#

		Saljack = nlsJack(VEGmle11)
		Saljack$jackCI
		boxplot(Saljack$jackCI)

	###GG plot##

		MLE_pred = predictNLS(VEGmle11, newdata=data.frame(Salinity=20, Inundation=Inundation_test), interval="prediction", alpha=0.05)
		MLE_pred_out=as.data.frame(MLE_pred[1])
		head(MLE_pred_out)

		SubVEG = subset(Veg, !Salinity > 30)
		SubVEG = subset(SubVEG, !Salinity < 10)

		ggplot() +
			geom_area(aes(x=Inundation_test, y=(MLE_pred_out$summary.Sim.97.5.)), fill="lightgray") +
			geom_area(aes(x=Inundation_test, y=(MLE_pred_out$summary.Sim.2.5.)), fill="white") +
			geom_line(aes(x=Inundation_test, y=(MLE_pred_out$summary.Sim.Mean)), lwd=1.5) +
			geom_point(aes(x=Veg$Inundation, y=Veg$Density), pch=19, col="black") +
			geom_point(aes(x=SubVEG$Inundation, y=SubVEG$Density), pch=19, col="red") +
			scale_x_continuous(limits=c(0,0.75)) +
			scale_y_continuous(limits=c(0,4000)) +
			xlab("Inundation Duration (%)") +
			ylab("Stem Density (m-2)") +
			ggtitle("Salinity = 20 PSU") +
			theme_classic()+
			theme(plot.title = element_text(hjust = 0.5))


#######################
###CONSOLIDATED CODE###
#######################

rm(list = ls())
library(ggplot2)
library(bbmle)
library(nlstools)
library(propagate)

#
Veg=read.csv("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_Data\\CH1_Vegetation.csv")
#

Veg$Inundation = 1/(1 + (0.033*exp(-0.9*1.34))*exp(0.9*Veg$Elev_m))^(1/0.06)
	plot(Veg$Elev_m, Veg$Inundation, pch=19)

hist(Veg$Density)
hist(log(Veg$Density + 1))
hist(sqrt(Veg$Density))

plot(Veg$Inundation, Veg$Density)
	abline(lm(Veg$Density~Veg$Inundation))

plot(Veg$Salinity, Veg$Density)
	abline(lm(Veg$Density~Veg$Salinity))

plot(Veg$Chan_Dist_m, Veg$Density)
	abline(lm(Veg$Density~Veg$Chan_Dist_m))	
			
VEGlm1 = lm(Density ~ Inundation, data=Veg)
VEGlm2 = lm(Density ~ Salinity, data=Veg)
VEGlm3 = lm(Density ~ Inundation + Salinity, data=Veg)
VEGlm4 = lm(Density ~ Inundation*Salinity, data=Veg)
VEGlm5 = lm(log(Density+1) ~ Inundation, data=Veg)
VEGlm6 = lm(log(Density+1) ~ Salinity, data=Veg)
VEGlm7 = lm(log(Density+1) ~ Inundation + Salinity, data=Veg)
VEGlm8 = lm(log(Density+1) ~ Inundation*Salinity, data=Veg)
VEGlm9 = lm(sqrt(Density) ~ Inundation, data=Veg)
VEGlm10 = lm(sqrt(Density) ~ Salinity, data=Veg)
VEGlm11 = lm(sqrt(Density) ~ Inundation + Salinity, data=Veg)
VEGlm12 = lm(sqrt(Density) ~ Inundation*Salinity, data=Veg)

AIC(VEGlm1,VEGlm2,VEGlm3,VEGlm4,VEGlm5,VEGlm6,VEGlm7,VEGlm8,VEGlm9,VEGlm10,VEGlm11,VEGlm12)
as.array(c(sigma(VEGlm1),sigma(VEGlm2),sigma(VEGlm3),sigma(VEGlm4),sigma(VEGlm5),
	sigma(VEGlm6),sigma(VEGlm7),sigma(VEGlm8)))
summary(VEGlm7)

ggplot(Veg, aes(x = Inundation, y = Density, color = Salinity))+
	geom_point(size=3)+
	labs(x = 'Inundation', y = 'Stem Density (m-2)')+
	scale_color_continuous(guide = FALSE)+
	theme_classic()

ggplot(Veg, aes(x = Salinity, y = Density, color = Inundation))+
	geom_point(size=3)+
	labs(x = 'Salinity', y = 'Stem Density (m-2)')+
	scale_color_continuous(guide = FALSE)+
	theme_classic()

Pred.full=model.matrix(VEGlm4)
	Pred.full=as.data.frame(Pred.full)
Pred.full$Predicted=predict(VEGlm4, Pred.full)
	plot(Veg$Inundation, Veg$Density)
	points(Pred.full$Inundation, Pred.full$Predicted, col="red", pch=19)
	plot(Veg$Salinity, Veg$Density)
	points(Pred.full$Salinity, Pred.full$Predicted, col="red", pch=19)

Pred.full=model.matrix(VEGlm7)
	Pred.full=as.data.frame(Pred.full)
Pred.full$Predicted=predict(VEGlm7, Pred.full)
	plot(Veg$Inundation, Veg$Density)
	points(Pred.full$Inundation, exp(Pred.full$Predicted), col="red", pch=19)
	plot(Veg$Salinity, Veg$Density)
	points(Pred.full$Salinity, exp(Pred.full$Predicted), col="red", pch=19)

VEG_func1 = function(alpha, beta, gamma){
	Stems = alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
	plot(Veg$Inundation, VEG_func1(100, 0.15, 0.2))
VEGmle1=nls(Density ~ alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(alpha=100, beta=0.15, gamma=0.2),
	lower=c(0, 0, 0),
	upper=c(10000, 1, 1),
	algorithm="port")
	summary(VEGmle1)
AIC(VEGlm1, VEGmle1)

VEG_func2 = function(s1, s2, alpha, beta, gamma){
	Stems = s1 + s2*Veg$Salinity + alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
	plot(Veg$Inundation, VEG_func2(0, 2, 100, 0.15, 0.2))
VEGmle2=nls(Density ~ s1 + s2*Salinity + alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(s1 = 0, s2 = 2, alpha=100, beta=0.15, gamma=0.2),
	lower=c(-100, 0, 0, 0, 0),
	upper=c(100, 1000, 10000, 1, 1),
	algorithm="port")
	summary(VEGmle2)
AIC(VEGlm3, VEGmle2)

VEG_func3 = function(alpha, s1, beta, gamma){
	Stems = (alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2))))
	}
	plot(Veg$Inundation, VEG_func3(50, 100, 0.12, 0.11))
VEGmle3=nls(Density ~ (alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2)))), 
	data=Veg, 
	start = list(alpha = 5, s1 = 1, beta=0.15, gamma=0.2),
	lower=c(0, 0, 0, 0),
	upper=c(10000, 100, 1, 1),
	algorithm="port")
	summary(VEGmle3)
AIC(VEGlm4, VEGmle3)

VEG_func4 = function(alpha, beta, gamma){
	Stems = alpha*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
	plot(Veg$Inundation, VEG_func4(10, 0.12, 0.11))
VEGmle4=nls(Density ~ alpha*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(alpha = 10, beta=0.15, gamma=0.2),
	lower=c(0, 0, 0),
	upper=c(1000, 1, 1),
	algorithm="port")
	summary(VEGmle4)
AIC(VEGlm4, VEGmle4)
sigma(VEGmle4)

VEG_func5 = function(alpha, beta, gamma){
	Stems = alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
VEGmle5=nls(sqrt(Density) ~ alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(alpha=10, beta=0.15, gamma=0.2),
	lower=c(0, 0, 0),
	upper=c(10000, 1, 1),
	algorithm="port")
	summary(VEGmle5)
	plot(Veg$Inundation, Veg$Density, pch=19)
		points(Veg$Inundation, VEG_func5(coef(VEGmle5)[1], coef(VEGmle5)[2], coef(VEGmle5)[3])^2, pch=19, cex=0.75, col="red")
		points(Veg$Inundation, VEG_func1(coef(VEGmle1)[1], coef(VEGmle1)[2], coef(VEGmle1)[3]), pch=19, cex=0.75, col="blue")
AIC(VEGmle1, VEGmle5)

VEG_func6 = function(s1, s2, alpha, beta, gamma){
	Stems = s1 + s2*Veg$Salinity + alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
VEGmle6=nls(sqrt(Density) ~ s1 + s2*Salinity + alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(s1 = 0, s2 = 2, alpha=100, beta=0.15, gamma=0.2),
	lower=c(-100, 0, 0, 0, 0),
	upper=c(100, 1000, 10000, 1, 1),
	algorithm="port")
	summary(VEGmle6)
	plot(Veg$Inundation, Veg$Density, pch=19)
		points(Veg$Inundation, VEG_func6(coef(VEGmle6)[1], coef(VEGmle6)[2], coef(VEGmle6)[3], coef(VEGmle6)[4], coef(VEGmle6)[5])^2, pch=19, cex=0.75, col="red")
		points(Veg$Inundation, VEG_func2(coef(VEGmle2)[1], coef(VEGmle2)[2], coef(VEGmle2)[3], coef(VEGmle2)[4], coef(VEGmle2)[5]), pch=19, cex=0.75, col="blue")
AIC(VEGmle2, VEGmle6)

VEG_func7 = function(alpha, s1, beta, gamma){
	Stems = (alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2))))
	}
VEGmle7=nls(sqrt(Density) ~ (alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2)))) + (s1*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2)))), 
	data=Veg, 
	start = list(alpha = 5, s1 = 1, beta=0.15, gamma=0.2),
	lower=c(0, 0, 0, 0),
	upper=c(500, 100, 1, 1),
	algorithm="port")
	summary(VEGmle7)
	plot(Veg$Inundation, Veg$Density, pch=19)
		points(Veg$Inundation, VEG_func7(coef(VEGmle7)[1], coef(VEGmle7)[2], coef(VEGmle7)[3], coef(VEGmle7)[4])^2, pch=19, cex=0.75, col="red")
		points(Veg$Inundation, VEG_func3(coef(VEGmle3)[1], coef(VEGmle3)[2], coef(VEGmle3)[3], coef(VEGmle7)[4]), pch=19, cex=0.75, col="blue")
AIC(VEGmle3,VEGmle7)

VEG_func8 = function(alpha, beta, gamma){
	Stems = alpha*Veg$Salinity*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
VEGmle8=nls(sqrt(Density) ~ alpha*Salinity*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(alpha = 10, beta=0.15, gamma=0.2),
	lower=c(0, 0, 0),
	upper=c(1000, 1, 1),
	algorithm="port")
	summary(VEGmle8)
	plot(Veg$Inundation, Veg$Density, pch=19)
		points(Veg$Inundation, VEG_func8(coef(VEGmle8)[1], coef(VEGmle8)[2], coef(VEGmle8)[3])^2, pch=19, cex=0.75, col="red")
		points(Veg$Inundation, VEG_func4(coef(VEGmle4)[1], coef(VEGmle4)[2], coef(VEGmle4)[3]), pch=19, cex=0.75, col="blue")
AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8)

VEG_func9 = function(N,lambda){
	Stems = N*exp(-lambda*Veg$Salinity)
	}
VEGmle9=nls(Density ~ N*exp(-lambda*Salinity), 
	data=Veg, 
	start = list(N = 3000, lambda=0),
	lower=c(0, -5),
	upper=c(10000, 5),
	algorithm="port")
	summary(VEGmle9)
	plot(Veg$Salinity, Veg$Density, pch=19)
		points(Veg$Salinity, VEG_func9(coef(VEGmle9)[1], coef(VEGmle9)[2]), pch=19, cex=0.75, col="red")
AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9)

VEG_func10 = function(N,lambda, alpha, beta, gamma){
	Stems = N*exp(-lambda*Veg$Salinity) + alpha*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
	}
VEGmle10=nls(Density ~ N*exp(-lambda*Salinity) + alpha*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(N = 3000, lambda=0, alpha=100, beta=0.12, gamma=0.11),
	lower=c(0, -5, 0, 0, 0),
	upper=c(10000, 5, 10000, 100, 100),
	algorithm="port")
	summary(VEGmle10)
	plot(Veg$Inundation, Veg$Density, pch=19)
		points(Veg$Inundation, VEG_func10(coef(VEGmle10)[1], coef(VEGmle10)[2], coef(VEGmle10)[3], coef(VEGmle10)[4], coef(VEGmle10)[5]), pch=19, cex=0.75, col="red")
AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10)

VEG_func11 = function(N, lambda, beta, gamma){
	Stems = N*exp(-lambda*Veg$Salinity)*exp((-(Veg$Inundation - beta)^2)/(2*(gamma^2)))
		}
VEGmle11=nls(Density ~ N*exp(-lambda*Salinity)*exp((-(Inundation - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(N = 3000, lambda=0, beta=0.12, gamma=0.11),
	lower=c(0, -5, 0, 0),
	upper=c(10000, 5, 100, 100),
	algorithm="port")
	summary(VEGmle11)
plot(Veg$Inundation, Veg$Density, pch=19)
	points(Veg$Inundation, VEG_func11(coef(VEGmle11)[1], coef(VEGmle11)[2], coef(VEGmle11)[3], coef(VEGmle11)[4]), pch=19, cex=0.75, col="red")
	points(Veg$Inundation, VEG_func8(coef(VEGmle8)[1], coef(VEGmle8)[2], coef(VEGmle8)[3])^2, pch=19, cex=0.75, col="blue")
plot(Veg$Salinity, Veg$Density, pch=19)
	points(Veg$Salinity, VEG_func11(coef(VEGmle11)[1], coef(VEGmle11)[2], coef(VEGmle11)[3], coef(VEGmle11)[4]), pch=19, cex=0.75, col="red")
	points(Veg$Salinity, VEG_func8(coef(VEGmle8)[1], coef(VEGmle8)[2], coef(VEGmle8)[3])^2, pch=19, cex=0.75, col="blue")
AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10,VEGmle11)
sigma(VEGmle11)

VEG_func12 = function(alpha, beta, gamma){
	Stems = alpha*exp((-(Veg$Salinity - beta)^2)/(2*(gamma^2)))
	}
VEGmle12=nls(Density ~ alpha*exp((-(Salinity - beta)^2)/(2*(gamma^2))), 
	data=Veg, 
	start = list(alpha = 100, beta=3, gamma=3),
	lower=c(0, 0, 0),
	upper=c(10000, 1000, 1000),
	algorithm="port")
	summary(VEGmle12)
AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10,VEGmle11,VEGmle12)

VEG_func13 = function(alpha, beta_s, gamma_s, beta_i, gamma_i){
	Stems = alpha*exp((-(Veg$Salinity - beta_s)^2)/(2*(gamma_s^2)))*exp((-(Veg$Inundation - beta_i)^2)/(2*(gamma_i^2)))
	}
VEGmle13=nls(Density ~ alpha*exp((-(Salinity - beta_s)^2)/(2*(gamma_s^2)))*exp((-(Inundation - beta_i)^2)/(2*(gamma_i^2))), 
	data=Veg, 
	start = list(alpha = 100, beta_s=3, gamma_s=3, beta_i=0.1, gamma_i=0.1),
	lower=c(0, 0, 0, 0, 0),
	upper=c(10000, 1000, 1000, 100, 100),
	algorithm="port")
	summary(VEGmle13)
plot(Veg$Salinity, Veg$Density, pch=19)
	points(Veg$Salinity, VEG_func13(coef(VEGmle13)[1], coef(VEGmle13)[2], coef(VEGmle13)[3], coef(VEGmle13)[4], coef(VEGmle13)[5]), pch=19, cex=0.75, col="red")
plot(Veg$Inundation, Veg$Density, pch=19)
	points(Veg$Inundation, VEG_func13(coef(VEGmle13)[1], coef(VEGmle13)[2], coef(VEGmle13)[3], coef(VEGmle13)[4], coef(VEGmle13)[5]), pch=19, cex=0.75, col="red")
AIC(VEGmle1,VEGmle2,VEGmle3,VEGmle4,VEGmle5,VEGmle6,VEGmle7,VEGmle8,VEGmle9,VEGmle10,VEGmle11,VEGmle12, VEGmle13)

Inundation_test = seq(0,0.75,0.005)
Salinity_test = array(mean(Veg$Salinity, na.rm=TRUE),length(Inundation_test))
newdat=as.data.frame(cbind(Inundation_test, Salinity_test))
colnames(newdat) = c("Inundation","Salinity")

LM_pred = predict(VEGlm7, newdata=newdat, type="response", se.fit=TRUE)
plot(Veg$Inundation, Veg$Density)
plot(Inundation_test, LM_pred$fit, type="l", lwd=2, ylim=c(0,20))
	lines(Inundation_test, LM_pred$fit + LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
	lines(Inundation_test, LM_pred$fit - LM_pred$se.fit, col="gray",pch=19,lwd=2,lty=2)
	points(Veg$Inundation, log(Veg$Density+1), pch=19)

Saljack = nlsJack(VEGmle11)
Saljack$jackCI
plot(Saljack$jackCI)

MLE_pred = predictNLS(VEGmle11, newdata=data.frame(Salinity=20, Inundation=Inundation_test), interval="prediction", alpha=0.05)
MLE_pred_out=as.data.frame(MLE_pred[1])
head(MLE_pred_out)

SubVEG = subset(Veg, !Salinity > 30)
SubVEG = subset(SubVEG, !Salinity < 10)

ggplot() +
	geom_area(aes(x=Inundation_test, y=(MLE_pred_out$summary.Sim.97.5.)), fill="lightgray") +
	geom_area(aes(x=Inundation_test, y=(MLE_pred_out$summary.Sim.2.5.)), fill="white") +
	geom_line(aes(x=Inundation_test, y=(MLE_pred_out$summary.Sim.Mean)), lwd=1.5) +
	geom_point(aes(x=Veg$Inundation, y=Veg$Density), pch=19, col="black") +
	geom_point(aes(x=SubVEG$Inundation, y=SubVEG$Density), pch=19, col="red") +
	scale_x_continuous(limits=c(0,0.75)) +
	scale_y_continuous(limits=c(0,4000)) +
	xlab("Inundation Duration (%)") +
	ylab("Stem Density (m-2)") +
	ggtitle("Salinity = 20 PSU") +
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))
	
	
