library(shiny)
setwd("E:\\UW PHD\\Dissertation\\Chapter_1\\CH1_R\\Inundation-app\\")
WLL=read.csv("data/CH1_Inundation.csv")

###USER INTERFACE###
	
ui = fluidPage(
	titlePanel("Nisqually Inundation Model"),
	
	sidebarLayout(
		sidebarPanel(
			
			helpText("Calculate inundation duration (%) by MTL and elevation (m)"),
     	
			sliderInput(inputId = "mtl",
          	label = "Mean Tidal Level (MTL)",
     		value = 1.34, min = 1, max = 2)
		),
	
		mainPanel(
     		
			plotOutput(outputId = "plot")
		)
	)
     )

server = function(input, output){
	library(bbmle)
	library(nlstools)
	library(ggplot2)
	
	WLL=read.csv("data/CH1_Inundation.csv")
	Year = as.factor(WLL$Year)
	Elev_pred = seq(-1,5,0.01)

	Plot_Theme = theme_classic() +
	theme(legend.position = c(0.8,0.8),
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
	
	output$plot = renderPlot({
		
		title = "Indundation Duration"
		
			Sensor_Elevation = 0.4 - 1.34 + input$mtl
			W_Level = WLL$Level - 1.34 + input$mtl
			High_Water = max(W_Level)
			Ints = as.integer((High_Water-Sensor_Elevation)/0.05)
			Level = split(W_Level, Year)
			ID = array(0,dim=c(Ints,1+nlevels(Year)))
		
			for (i in 1:Ints){
				for (j in 1:nlevels(Year)){
					Threshold = (Sensor_Elevation - 0.05) + i*0.05
					ID[i,1] = Threshold
					ID[i,j+1] = sum(as.data.frame(Level[j]) > ID[i,1])/sum(as.data.frame(Level[j]) > 0)
					}
			}
			ID = as.data.frame(ID)
			colnames(ID) = c("Elevation","y2010","y2011","y2012","y2013","y2014","y2015")			
			ID2 = data.frame(ID[1], stack(ID[2:ncol(ID)]))			
			colnames(ID2) = c("El_m","Inund_perc","Yr")

			NLS_I=nls(Inund_perc ~ 1/((1 + Q*exp(-B*El_m))^(1/v)), 
				data=ID2, 
				start = list(B=-0.9, v=0.1, Q =0.017),
				lower=c(-10, 0.05, 0.0001),
				upper=c(-0.0001, 1, 1),
				algorithm="port")
			
			Inund = function(B, v, Q) {
				I = 1/((1 + Q*exp(-B*ID2$El_m))^(1/v))
				}
			
			ggplot() +
			geom_point(aes(x = ID2$El_m, y = ID2$Inund_perc, fill=""), size=2, color="red") +
			geom_line(aes(x = ID2$El_m, y = Inund(coef(NLS_I)[1], coef(NLS_I)[2], coef(NLS_I)[3]), color = ""), linetype = 1, size = 1) +
			scale_color_manual(name = "Modeled Data", values = "black") +
			scale_fill_manual(name = "Raw Data", values = "red") +
			scale_y_continuous(name="Inundation Duration", limit = c(0,1), breaks = c(0,0.25,0.5,0.75,1)) +
			scale_x_continuous(name="Elevation (m NAVD88)", limit = c(0,5), breaks = c(0,1,2,3,4,5)) +
			Plot_Theme
			
		})
		}

shinyApp(ui = ui, server = server)
