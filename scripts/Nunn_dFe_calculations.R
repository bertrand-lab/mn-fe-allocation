## Script written by Loay Jabre (2019)

# Calculating free iron (Fe') concentrations:
# 
# Total iron, AKA total amount of iron you add to the media with EDTA is different than Fe', AKA free iron, AKA inorganic iron that is available for phytoplankton to grow. You can calculate Fe' from total iron, if you know the temperature, light level, pH, number of light hours etc.. This is a difficult calculation, and there are many assumptions to make, because the dissociation constants are not measured experimentally at the low temperatures temperatures that Fragilariopsis grows at. So, using extrapolations from published literature, we can do the following: 
# 
# These are data from Sunda and Huntsman, 2003 "Effect of pH, light, and temperature on Fe-EDTA chelation and Fe hydrolysis in seawater" paper. I used the data from the table they provided to re-draw their graphs, so I can do an extrapolation. 
# 
# just entering the values from the table in Sunda and Huntsman 2003


library(magrittr)
library(dplyr)
library(ggplot2)

Temperature <- c(20,20,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10
)

PAR<- c(0,0,0,0,0,0,0,500,500,500,500,500,0,0,0,0,0,0,0,0,0,0,0,0,500,500,500,500,500,500,500,500,500,500
)

pH <- c(8.95,8.43,7.99,8.19,8.14,8.15,8.81,7.83,7.98,8.38,8.93,8.15,7.85,8.06,8.29,7.78,8.11,8.33,7.65,7.99,8.25,7.68,7.98,8.21,7.89,8.07,8.28,7.85,8.02,8.22,8.02,8.17,8.02,8.15
)

EDTA <- c(50,10,3,3,35,8,25,20,3,10,50,35,1,3,3,1,3,3,1,3,3,1,3,3,10,30,30,10,30,30,30,30,30,30
)

Fe_total <- c(20,20,20,20,9,9,9,9,20,20,20,9,20,20,20,8,8,8,8,8,8,8,8,8,20,20,20,20,20,20,8,8,8,8
)

FeIII <- c(2.73,0.77,0.238,0.57,0.028,0.096,1.12,0.071,1.45,1.29,2.55,0.086,0.466,0.456,0.93,0.149,0.19,0.376,0.101,0.137,0.296,0.112,0.105,0.247,0.819,0.451,0.7,0.785,0.623,0.695,0.178,0.237,0.19,0.227
)

logKd <- c(-5.1,-6.4,-7.44,-7.05,-6.97,-7.07,-5.46,-6.8,-6.63,-6.16,-5.14,-6.47,-7.62,-7.16,-6.83,-7.72,-7.14,-6.83,-7.89,
           -7.28,-6.94,-7.85,-7.4,-7.02,-6.37,-6.16,-5.96,-6.39,-6.02,-5.97,-6.17,-6.04,-6.14,-6.06
)

#combining all the values into a dataframe
dissociation_table <- data.frame(Temperature, PAR, pH, EDTA, Fe_total, FeIII, logKd)

#changing some of the datatypes into characters for plotting purposes
dissociation_table$Temperature <- as.character(dissociation_table$Temperature)
dissociation_table$PAR <- as.character(dissociation_table$PAR)

#this plot is almost identical to the plot provided in the paper. I believe that they missed a data point, at the 500PAR, 20C, so that line looks slightly different. 
ggplot (dissociation_table, aes( x = pH, y = logKd, color = Temperature, shape = PAR))+
  geom_point(size = 5, alpha = 0.9)+
  theme_bw()+ 
  scale_color_manual (values = c("blue", "red"))+
  scale_shape_manual(values  = c(16,1))+
  geom_smooth(method = 'lm', se = FALSE)+
  ylab(expression(Log~K[d]~(dissociation~constant)))+
  xlab(expression(pH))+
  theme(axis.text.x=element_text(face = "bold", size = 20, color = "black"),
        axis.title.x=element_text(size=24,face="bold"), 
        axis.text.y=element_text(face = "bold", size = 20, color = "black"),
        axis.title.y=element_text(size=24,face="bold"))+
  theme (legend.title = element_text (size = 16), legend.text = element_text(size = 12), legend.position = c(0.83, 0.2))+
  labs(col="Temperature (°C)")+
  labs (shape = expression(PAR~(mu*mol~photons~m^-2~s^-1))) +
  theme(panel.grid = element_blank())

#Because pH has an effect on Kd (dissociation constant), I chose pH 8-8.3 (similar to my growth media) to extrapolate the effects of temperature and light. 
#The Kd in the dark will be Kd(dark), the Kd in the light (regardless of light level), can be refered to as Khv

dissociation_constants_pH <- dissociation_table %>% filter(between (pH, 8 ,8.3))

ggplot (dissociation_constants_pH, aes( x = Temperature, y = logKd, shape = PAR, group= PAR))+
  geom_point(size = 5)+
  geom_smooth(method = 'lm', color = "black")+
  scale_shape_manual (values = c(16,1))+
  ylab(expression(Log~K[dissociation]~(dissociation~constant)))+
  xlab(expression(Temperature~(degree*C)))+
  theme_bw()+
  theme(axis.text.x=element_text(face = "bold", size = 20, color = "black"),
        axis.title.x=element_text(size=24,face="bold"), 
        axis.text.y=element_text(face = "bold", size = 20, color = "black"),
        axis.title.y=element_text(size=24,face="bold"))+
  theme (legend.title = element_text (size = 14), legend.text = element_text(size = 12), legend.position = c(0.83, 0.9))+
  labs (shape = expression(PAR~(mu*mol~photons~m^-2~s^-1))) +
  theme(panel.grid = element_blank())

#From here we see that in the dark, Kdissociation(dark) is not affected by temperature, but Kd(light) or Khv is  temperature sensitive. I will assume that the temperature effect on Kd in the light (Khv) is linear (no other data available to suggest otherwise). so, I modelled it as such. 

dissociation_constants_pH$PAR <- as.numeric(dissociation_constants_pH$PAR)
dissociation_constants_pH$Temperature <- as.numeric(dissociation_constants_pH$Temperature)
#choosing only the light 
Kd_data <- dissociation_constants_pH %>% 
  filter(between(PAR, 100, 501))

#doing a linear model and plotting the Khv (Kd in the light) value against some temperatures
templinear <- lm(logKd~ Temperature, data=Kd_data)
summary(templinear)
dataframe <- data.frame (Temperature=c(1,3,6,10,20))

z<- predict.lm(templinear, dataframe)
xdata<- data.frame(z)
rownames(xdata)<- c(1,3,6,10,20)

xdata$Temperature <- c(1,3,6,10,20)


#Now we can use that linear model to calculate Fe' in media

#Enter the total Fe concentrations here (Mol/L). We are assuming that FeEDTA = totalFe (whatever concentrations they have in the paper)
FeEDTA <- c(42e-9, 1.37e-6)

#Enter the EDTA* (EDTA) concentrations here (Mol/L): This is the total concentration of EDTA minus all of the concentrations of metals (including iron, so at high iron concentrations, you'll have less EDTA)
EDTA <- c(100e-6 - FeEDTA[1], 100e-6 - FeEDTA[2])

#Enter the temperatures (C)that your cultures were growting at 
Temperature <- c(19)
temp <- data.frame(Temperature)

#Enter the light intensity  (umol photons/m2/s)
Ihv <- 150

#enter the number of light hours per day (i.e. 24 = constant light)
h <- 24

#Kd is the dissociation constant in the dark. This is independent of temperature and will not change. 
Kd <- 9.5499E-8

iron <- data.frame(EDTA, FeEDTA, Kd, Ihv)
iron_2 <- merge (iron, temp)

iron_2$logKhv <- predict.lm(templinear, iron_2)

#enter the value from the above dataframe and constants into this equation from 'algal culturing techniques', p 48. 
iron_2$Fe <- (iron_2$FeEDTA*(Kd+Ihv/500*10^(iron_2$logKhv)*h/24))/iron_2$EDTA


#the iron_2 dataframe now has the Fe' concentrations that correspond to different temperatures, we can plot those below , and we can see that as temperature increases, Fe' decreases. 

