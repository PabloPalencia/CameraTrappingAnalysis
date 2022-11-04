###############################
##                           ##
##      Estimating density   ##
##         througth          ##
##            REM            ##
##                           ##
##       Pablo Palencia      ##
##        03/11/2022         ##
##                           ##
###############################


# R code for implement Random Encounter Model (REM), a method to calculate populations densities from camera-trap 
# data for species which do no exhibit individually-indentifiable markings (Rowcliffe et al. 2008). The code also
# includes functions to estimate travel speed (Palencia et al. 2021, Rowcliffe et al. 2016), the activity value 
# (Rowcliffe et al 2014) and standard error of obtained densities values.

# Working directory

# This will open an explorer window for easy browsing
setwd(choose.dir())
# To get the current working directory and verify that it is correct
getwd()
# Import the .txt table into R when the data are separated by ";"

# Load dataframes
dataREM <- read.table("Data.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # parameters dataframe
operat <- read.table("Operativity.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # operatity matrix (to estimate survey effort)
df_coord <- read.table("Coordinates.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # camera trap locations (plots, maps etc.)

# The columns of the data must have the following information:
# col1: Point_ID: ID of the camera-trap location
# col2: Sp: Species observed in each sequence
# col3: G_size: the total number of individuals observed in the sequence
# col4: Date: date of the sequence
# col5: H_first: exactly hour of the first photo considered for estimate travel speed
# col6: H_last: exactly hour of the last photo considered for estimate travel speed
# col7: T: the difference between time stamps of the first and last images considered, in hours
# col8: T.s.: the difference between time stamps of the first and last images considered, in seconds
# col9: Dist.m.: distance traveled in each sequence. Expressed in m
# col10: Speed.m.s: travel speed estimated for each sequence. Expressed in m/s
# col11: Interval.min.: Shorter distance between the animal and the camera. 
# 5 intervals are considered: 1(0-2.5m to the camera-trap),2(2.5-5m), 3(5-7.5m), 4(7.5-10m), 5(>10m)
# col12: Dist_det: Distance of detection. Distance between the animal and the camera in the first photo of each sequence
# col13: Ang_det: Angle of detection. Angle of detection in the first photo of each sequence. Expresed in decimal grades
# col14: Behaviour: Observed behaviour for each sequence. See details in Palencia et al. 2019

### FUNCTIONS (not included in packages. Developed by M. Rowcliffe & St. Andrews University)
source("REM_functions.R") # importing some key functions to run the analysis


##########################################################################################################################################
## ACTIVITY
###################################################################################################################################
library(activity)
data.acti <- as.data.frame(dataREM)

# Estimating radian time of day
data.acti$T_sec2 <- as.numeric(strptime(data.acti$H_first, format="%H:%M:%S") - as.POSIXct(format(Sys.Date())), units="secs")
data.acti$T_0_1 <- data.acti$T_sec2/86400; data.acti <- subset(data.acti, G_size!="NA") # remove NAs

# Replicating each sequence based on the group size
activity.repli <- data.acti[rep(row.names(data.acti), data.acti$G_size), 1:16] 

# Discariding observations further than 5 meters (see Rowcliffe et al. 2014 Methods Ecol. Evol, 5(11): 1170-1179)
activity.repli <- subset(activity.repli, Interval.min =="1"  | Interval.min =="2")

# Estimating activity rate
activityRES <- 2*pi*activity.repli$T_0_1
mod1 <- fitact(activityRES, sample="data") 

# Plot activity patterns
par(mfrow=c(1,1)); plot(mod1); show(mod1)

##################################################################################################################
## TRAVEL SPEED (Palencia et al. 2019, Palencia et al. 2020, Rowcliffe et al. 2016)
###################################################################################################################
library(trappingmotion)

data.speed <- subset(dataREM, Behaviour != "Curiosity") # discard animals that react to the camera
data.speed <- subset(data.speed, Speed.m.s!="NA"); data.speed$Speed.m.s <- as.numeric(as.character(data.speed$Speed.m.s)) # Remove NAs

identbhvs(data.speed$Speed.m.s) # identify movement states (see Palencia et al. 2021 - Methods Ecol. Evol.)
meanspeed(behav_class) # average movement speed of each state

###################################################################################################################################
## Day Range
###################################################################################################################################

dayrange(act=mod1@act[1], act_se=mod1@act[2], speed_data) #day range (daily distance travelled)

###################################################################################################################################
### DETECTION ZONE
###################################################################################################################################

library(Distance)

data_dz_r<-subset(dataREM, Dist_det >= 0 )
data_dz_ang<-subset(dataREM, Ang_det != "NA" )


# Efective Detection Radius 
w_rad <- 10
# half-normal
hn <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = NULL, truncation=w_rad) 
hn_cos <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "cos", nadj = 1, truncation=w_rad)
hn_herm <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "herm", nadj = 1, truncation=w_rad)
hn_poly <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "poly", nadj = 1, truncation=w_rad)

#hazard-rate
hr <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = NULL, truncation=w_rad) 
hr_cos <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "cos", nadj = 1, truncation=w_rad)
hr_herm <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "herm", nadj = 1, truncation=w_rad)
hr_poly <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "poly", nadj = 1, truncation=w_rad)

#model comparison
AIC(hn, hn_cos, hn_herm, hn_poly, hr, hr_cos, hr_herm, hr_poly)

# select best model
# mind the fact that if your data is spiked at zero, you have to be carefoul with the hazard-rate model (details in Buckland et al. 2001)
best_modRad <- hr

# Estimating effective detection radius and (SE)
EfecRad <- EDRtransform(best_modRad)

EfecRad$EDR # mean
EfecRad$se.EDR # SE

# Plots... 
par(mfrow=c(1,2))
plot(best_modRad, main="Best model", xlab="Distance (m)",
     showpoints=FALSE, lwd=3, xlim=c(0, 10))
plot(best_modRad, main="Best model", xlab="Distance (m)", pdf=TRUE,
     showpoints=FALSE, lwd=3, xlim=c(0, 10))
par(mfrow=c(1,1))

# Efective Detection Angle 

data_dz_ang$Ang_rad <- abs(data_dz_ang$Ang_det*0.0174533)
FOV <- 42 # field of view of the cameras (degrees)
w_ang <- FOV/2*0.0174533

# half-normal
hn_Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = NULL, truncation=w_ang) 
hn_cosAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "cos", nadj = 1, truncation=w_ang)
hn_hermAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "herm", nadj = 1, truncation=w_ang)
hn_polyAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "poly", nadj = 1, truncation=w_ang)

#hazard-rate
hr_Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = NULL, truncation=w_ang) 
hr_cosAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "cos", nadj = 1, truncation=w_ang)
hr_hermAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "herm", nadj = 1, truncation=w_ang)
hr_polyAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "poly", nadj = 1, truncation=w_ang)

#uniform
uni_cosAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="uni", adjustment = "cos", nadj = 1, truncation=w_ang)
uni_hermAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="uni", adjustment = "herm", nadj = 1, truncation=w_ang)
uni_polyAng <- ds(data_dz_ang$Ang_rad, transect = "line", key="uni", adjustment = "poly", nadj = 1, truncation=w_ang)

#model comparison
AIC(hn_Ang, hn_cosAng, hn_hermAng, hn_polyAng, hr_Ang, hr_cosAng, hr_hermAng, hr_polyAng, uni_cosAng, uni_hermAng, uni_polyAng)

# select best model
# mind the fact that if your data is spiked at zero, you have to be carefoul with the hazard-rate model (details in Buckland et al. 2001)
best_modAng <- hr_Ang 

# Estimating effective detection radius and (SE)
summary_ang<- summary(best_modAng$ddf)

EfecAng_mean <- summary_ang$average.p*w_ang # mean
EfecAng_SE <-  summary_ang$average.p.se*w_ang# SE

# Plots... 
par(mfrow=c(1,1))
plot(best_modAng, main="Best model", xlab="Angle (rad)",
     showpoints=FALSE, lwd=3, xlim=c(0, w_ang))

###################################################################################################################################
### TRAPPING-RATE
###################################################################################################################################

data.dens <- as.data.frame(dataREM)
#data.dens <- subset(data.dens, Point_ID != 48 & Point_ID != 100)# to remove some points

tm <- sum(colSums(operat, na.rm = FALSE, dims = 1))- sum(operat$CAM) # Survey effort (camera days)

data.dens <- subset(data.dens, Dist_det < 10.1 | is.na(Dist_det)) # Based on truncation distance to estimate EDD. NAs are new sequences
seq <- length(data.dens[, 1]) # Number of sequences 

operat$oper_days<-rowSums(operat, dims = 1)-operat$CAM
seq_point<-data.frame(table(data.dens$Point_ID))
operat <- merge(operat, seq_point, by.x= "CAM", by.y = "Var1", all.x = TRUE); operat[is.na(operat)] <- 0
operat$tr <- operat$Freq/operat$oper_days

library(dplyr)
tr<-operat[,c("Freq","oper_days")] # Selecting columns
tr <- subset(tr, oper_days > 0) # remove CT that doesnt work properly

#Plot encounter rate
df_coord$ER <- tr$Freq # add encounter rate to coordinates data frame

library(ggplot2)
ggplot(data=df_coord) + 
  aes(x=Long, y=Lat, size=ER) + 
  geom_point(alpha=1, shape=16, col='blue') +
  scale_size(breaks = c(0, 1, 5, 10, 20, 50), range = c(1, 20))+
  ylab("Latitude (m)") +
  xlab("Longitude (m)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2)))

###################################################################################################################################
### DENSITY
###################################################################################################################################
# Average values
param <- list(DR = DR,
              r = EfecRad$EDR / 1000,
              theta = EfecAng_mean*2)

# Standar error values
paramse <- list(DR = DR_se,
                r = EfecRad$se.EDR / 1000,
                theta = EfecAng_SE*2)

density<-bootTRD(tr$Freq, tr$oper_days, param, paramse); density

# Saving results
results <- data.frame(seq, tm, DR, DR_se, EfecRad$EDR, EfecRad$se.EDR, EfecAng_mean*2, EfecAng_SE*2, density[,1], density[,2])
dimnames(results) <- list("Value", c("y(seq.)","t (days)", "s(km/day)","s_se(km/day)", "r(m)", "r_se(m)", "ang(rad)", "ang_se(rad)", "d(ind/km2)", "d_se(ind/km2)")); View(results)
#write.table(results, "ResultsREM_MARTESsp_FITO.txt", sep=";", row.names = FALSE)

