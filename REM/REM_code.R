###############################
##                           ##
##      Estimating density   ##
##          through          ##
##            REM            ##
##                           ##
##       Pablo Palencia      ##
##        01/12/2025         ##
##                           ##
###############################


# R code for run Random Encounter Model (REM) analysis
# Working directory

# This will open an explorer window for easy browsing
setwd(choose.dir())
# To get the current working directory and verify that it is correct
getwd()
# Import the .txt table into R when the data are separated by ";"

# Load data frames
dataREM <- read.table("Data.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # parameters data frame
operat <- read.table("Operativity.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # operativity matrix (to estimate survey effort)
df_coord <- read.table("Coordinates.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # camera trap locations (plots, maps etc.)

# Load functions
source("REM_functions.R") # importing some key functions to run the analysis

# Packages required to run the analyses
library(activity) # to estimate activity pattern and day range (available on CRAN)
library(trappingmotion) # to estimate speed and day range (available on github https://github.com/PabloPalencia/trappingmotion)
library(Distance) # to estimate detection zone (available on CRAN)
library(dplyr) # to work with data frames (available on CRAN)
library(ggplot2) # to plot encounter rates (available on CRAN)

##########################################################################################################################################
## ACTIVITY
###################################################################################################################################
# Convert time of dataREM to a numeric vector of radian time-of-day
time_rad <- gettime(dataREM$time, format = "%H:%M:%S")

# fit activity model
actmod <- fitact(time_rad, sample="model") #sample=model: large sample size (greater than 100-200); sample=dataREM: small sample size (less than 100); o sample=none: no bootstrapping.

# Plot activity patterns
par(mfrow=c(1,1)); plot(actmod)

actmod@act[1] # mean activity level
actmod@act[2] # SE activity level

##################################################################################################################
## TRAVEL SPEED (Palencia et al. 2019, Palencia et al. 2020, Rowcliffe et al. 2016)
###################################################################################################################

# Explore speed distribution. Sometimes it is necessary to remove extreme-low/high values
par(mfrow=c(1,2))
boxplot(dataREM$speed); hist(dataREM$speed)
dataREM.speed <- subset(dataREM, speed < 8)
#dataREM.speed <- subset(dataREM.speed, speed > 0.04)

identbhvs(dataREM.speed$speed) # identify movement states
table(behav_class$behaviour)
meanspeed(behav_class) # average movement speed of each state

###################################################################################################################################
## Day Range
###################################################################################################################################

dayrange(act=actmod@act[1], act_se=actmod@act[2], speed_data) #day range (daily distance traveled)

###################################################################################################################################
### DETECTION ZONE
###################################################################################################################################

# explore and define the truncation distance
hist(dataREM$radius)

w_rad <- 10
# half-normal
hn <- ds(dataREM$radius, transect = "point", key="hn", adjustment = NULL, truncation=w_rad)
hn_cos <- ds(dataREM$radius, transect = "point", key="hn", adjustment = "cos", nadj = 1, truncation=w_rad)
hn_herm <- ds(dataREM$radius, transect = "point", key="hn", adjustment = "herm", nadj = 1, truncation=w_rad)
hn_poly <- ds(dataREM$radius, transect = "point", key="hn", adjustment = "poly", nadj = 1, truncation=w_rad)

# hazard-rate
hr <- ds(dataREM$radius, transect = "point", key="hr", adjustment = NULL, truncation=w_rad)
hr_cos <- ds(dataREM$radius, transect = "point", key="hr", adjustment = "cos", nadj = 1, truncation=w_rad)
hr_herm <- ds(dataREM$radius, transect = "point", key="hr", adjustment = "herm", nadj = 1, truncation=w_rad)
hr_poly <- ds(dataREM$radius, transect = "point", key="hr", adjustment = "poly", nadj = 1, truncation=w_rad)

# model selection
AIC(hn, hn_cos, hn_herm, hn_poly, hr, hr_cos, hr_herm, hr_poly)

# select best model
# (mind the fact that if your data is spiked at zero, you have to be careful with the hazard-rate model (details in Buckland et al. 2001))
best_modRad <- hn 

# Estimating effective detection radius and (SE)
EfecRad <- EDRtransform(best_modRad)

EfecRad$EDR # mean (m)
EfecRad$se.EDR # SE (m)

# Plots 
par(mfrow=c(1,2))
plot(best_modRad, main="Best model", xlab="Radius (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 10))
plot(best_modRad, main="Best model", xlab="Radius (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 10))

# Effective Detection Angle 

dataREM$angle <- abs(dataREM$angle)
hist(dataREM$angle)
w_ang <- 0.41 # truncation angle in radians (=23.5 grados)

# half-normal
hn_Ang <- ds(dataREM$angle, transect = "line", key="hn", adjustment = NULL, truncation=w_ang)
hn_cosAng <- ds(dataREM$angle, transect = "line", key="hn", adjustment = "cos", nadj = 1, truncation=w_ang)
hn_hermAng <- ds(dataREM$angle, transect = "line", key="hn", adjustment = "herm", nadj = 1, truncation=w_ang)
hn_polyAng <- ds(dataREM$angle, transect = "line", key="hn", adjustment = "poly", nadj = 1, truncation=w_ang)

# hazard-rate
hr_Ang <- ds(dataREM$angle, transect = "line", key="hr", adjustment = NULL, truncation=w_ang)
hr_cosAng <- ds(dataREM$angle, transect = "line", key="hr", adjustment = "cos", nadj = 1, truncation=w_ang)
hr_hermAng <- ds(dataREM$angle, transect = "line", key="hr", adjustment = "herm", nadj = 1, truncation=w_ang)
hr_polyAng <- ds(dataREM$angle, transect = "line", key="hr", adjustment = "poly", nadj = 1, truncation=w_ang)

# uniform
uni_cosAng <- ds(dataREM$angle, transect = "line", key="uni", adjustment = "cos", nadj = 1, truncation=w_ang)
uni_hermAng <- ds(dataREM$angle, transect = "line", key="uni", adjustment = "herm", nadj = 1, truncation=w_ang)
uni_polyAng <- ds(dataREM$angle, transect = "line", key="uni", adjustment = "poly", nadj = 1, truncation=w_ang)

# model comparison
AIC(hn_Ang, hn_cosAng, hn_hermAng, hn_polyAng, hr_Ang, hr_cosAng, hr_hermAng, hr_polyAng, uni_cosAng, uni_hermAng, uni_polyAng)

# select best model
best_modAng <- uni_cosAng  

# Estimating effective detection radius and (SE)
summary_ang<- summary(best_modAng$ddf)

EfecAng_mean <- summary_ang$average.p*w_ang # mean (radians)
EfecAng_SE <-  summary_ang$average.p.se*w_ang # SE (radians)

# Plots 
par(mfrow=c(1,1))
plot(best_modAng, main="Best model", xlab="Angle (rad)",
     showpoints=FALSE, lwd=3, xlim=c(0, w_ang))

###################################################################################################################################
### ENCOUNTER RATE
###################################################################################################################################

dataREM.dens <- as.data.frame(dataREM)

# estimate survey effort
tm <- sum(operat$days) # Survey effort (camera days)

dataREM.dens$angle <- abs(dataREM.dens$angle)
dataREM.dens <- subset(dataREM.dens, radius <= w_rad & angle <= w_ang) # Based on truncation distance to estimate effective radius and angle
seq <- length(dataREM.dens[, 1]) # Number of encounters

seq_point<-(table(dataREM.dens$camera))
operat <- merge(operat, seq_point, by.x= "camera", by.y = "Var1", all.x = TRUE); operat[is.na(operat)] <- 0
operat$tr <- operat$Freq/operat$days

enc_rate<-operat[,c("Freq","days")] # Selecting columns

# Plot encounter rate
df_coord$ER <- enc_rate$Freq # add encounter rate to coordinates dataREM frame

ggplot(data = df_coord) +
  aes(x = x, y = y, size = ER, colour = ER == 0) + 
  geom_point(alpha = 1, shape = 16) +
  scale_size(breaks = c(0, 1, 5, 10, 20, 50), range = c(1, 40)) +
  scale_colour_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  ylab("Latitude (m)") +
  xlab("Longitude (m)") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 20, face = "bold"),
    axis.title.x = element_text(color = "black", size = 20, face = "bold"),
    axis.text.y = element_text(size = rel(2)),
    axis.text.x = element_text(size = rel(2))
  )

###################################################################################################################################
### DENSITY
###################################################################################################################################
# We include in a list average values of day range and detection zone
param <- list(DR = DR,
              r = EfecRad$EDR / 1000,
              theta = EfecAng_mean*2)

# We include in a list standard error values of day range and detection zone
paramse <- list(DR = DR_se,
                r = EfecRad$se.EDR / 1000,
                theta = EfecAng_SE*2)

# Density estimation
density<-bootTRD(enc_rate$Freq, enc_rate$days, param, paramse); density

# CV
density[,2]/density[,1]

# Log-normal CIs
df_lnorm_confint <- lnorm_confint(density[,1], density[,2])

# Saving results
results <- data.frame(seq, tm, DR, DR_se, EfecRad$EDR, EfecRad$se.EDR, EfecAng_mean*2, EfecAng_SE*2, density[,1], density[,2])
dimnames(results) <- list("Value", c("y(seq.)","t (days)", "s(km/day)","s_se(km/day)", "r(m)", "r_se(m)", "ang(rad)", "ang_se(rad)", "d(ind/km2)", "d_se(ind/km2)")); View(results)

