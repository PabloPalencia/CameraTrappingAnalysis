###############################
##                           ##
##      Estimating density   ##
##         througth          ##
##            REM            ##
##                           ##
##       Pablo Palencia      ##
##        26/03/2021         ##
##                           ##
###############################


# R code for implement Random Encounter Model (REM), a method to calculate populations densities from camera-trap 
# data for species which do no exhibit individually-indentifiable markings (Rowcliffe et al. 2008). The code also
# includes functions to estimate travel speed (Palencia et al. 2019, Rowcliffe et al. 2016), the activity value 
# (Rowcliffe et al 2014) and standard error of obtained densities values.

# Working directory

# This will open an explorer window for easy browsing
setwd(choose.dir())
# To get the current working directory and verify that it is correct
getwd()
# Import the .txt table into R when the data are separated by ";"

# Load data
data <- read.table("LiebresQM_REM.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE)

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
# col15: 1: Individuals cross the midline of the field of view of the camera-trap. 0: Individuals don't cross the midline of the field of view of the camera-trap.
# col16: Notes: Observations 

### FUNCTIONS (not included in packages. Developed by M. Rowcliffe & St. Andrews University)

source("REM_functions.R") # uploading some key functions to run the analysis

##########################################################################################################################################
## ACTIVITY
###################################################################################################################################
library(activity)
data.acti <- as.data.frame(data)

# Transforming clock time to sun time and estimating radian time of day
data.acti$T_sec2 <- as.numeric(strptime(data.acti$H_first, format="%H:%M:%S") - as.POSIXct(format(Sys.Date())), units="secs")

# Transforming (January to April & November to December: -1; May to October: -2)
# Select one line based on the period of study
#data.acti$T_sec2<-ifelse(data.acti$T_sec>3600, data.acti$T_sec-3600, data.acti$T_sec+82800) # January to April & November to December
#data.acti$T_sec2<-ifelse(data.acti$T_sec>7200, data.acti$T_sec-7200, data.acti$T_sec+79200) # May to October

data.acti$T_0_1 <- data.acti$T_sec2/86400; data.acti <- subset(data.acti, G_size!="NA") # remove NAs

# Replicating each sequence based on the group size
activity.repli <- data.acti[rep(row.names(data.acti), data.acti$G_size), 1:18] 

# Discariding observations further than 5 meters
activity.repli <- subset(activity.repli, Interval.min =="1"  | Interval.min =="2")#| Interval.min =="3") 

# Estimating activity value
activityRES <- 2*pi*activity.repli$T_0_1
mod1 <- fitact(activityRES, sample="data") #sample=model: large sample size (greater than 100-200); sample=data: small sample size (less than 100); o sample=none: no bootstrapping. 

# Plot activity patterns
par(mfrow=c(1,1))
plot(mod1)
show(mod1)

##################################################################################################################
## TRAVEL SPEED (Palencia et al. 2019, Palencia et al. 2020, Rowcliffe et al. 2016)
###################################################################################################################
library(trappingmotion)

data.speed <- subset(data, Behaviour != "Curiosity")

# Remove NAs
data.speed <- subset(data.speed, Speed.m.s!="NA"); data.speed$Speed.m.s <- as.numeric(as.character(data.speed$Speed.m.s))

# Explore speed distribution. Sometimes it is necessary to remove extreme-low values
#summary(data.speed$Speed.m.s); hist(data.speed$Speed.m.s, breaks=30)
boxplot(data.speed$Speed.m.s)
sort(data.speed$Speed.m.s)
data.speed <- subset(data.speed, Speed.m.s > 0.01)
data.speed <- subset(data.speed, Speed.m.s < 4)

identbhvs(data.speed$Speed.m.s) # identify movement states (see Palencia et al. 2021 - Methods Ecol. Evol.)
table(behav_class$behaviour)
meanspeed(behav_class) # average movement speed of each state

###################################################################################################################################
## Day Range
###################################################################################################################################

dayrange(act=mod1@act[1], act_se=mod1@act[2], speed_data) #day range (daily distance travelled)

###################################################################################################################################
### DETECTION ZONE
###################################################################################################################################

library(Distance)

data_dz_r<-subset(data, Dist_det >= 0 )
data_dz_ang<-subset(data, Ang_det != "NA" )


# Efective Detection Radius 
w_rad <- 10
# half-normal
hn_cos0 <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "cos", order = 0, truncation=w_rad) 
hn_cos2 <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "cos", order = 2, truncation=w_rad)
hn_herm0 <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "herm", order = 0, truncation=w_rad)
hn_herm2 <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "herm", order = 2, truncation=w_rad)
hn_poly0 <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "poly", order = 0, truncation=w_rad)
hn_poly2 <- ds(data_dz_r$Dist_det, transect = "point", key="hn", adjustment = "poly", order = 2, truncation=w_rad)

#hazard-rate
hr_cos0 <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "cos", order = 0, truncation=w_rad) 
hr_cos2 <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "cos", order = 2, truncation=w_rad)
hr_herm0 <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "herm", order = 0, truncation=w_rad)
hr_herm2 <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "herm", order = 2, truncation=w_rad)
hr_poly0 <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "poly", order = 0, truncation=w_rad)
hr_poly2 <- ds(data_dz_r$Dist_det, transect = "point", key="hr", adjustment = "poly", order = 2, truncation=w_rad)

#model comparison
AIC(hn_cos0, hn_cos2, hn_herm0, hn_herm2, hn_poly0, hn_poly2, hr_cos0, hr_cos2, hr_herm0, hr_herm2, hr_poly0, hr_poly2)

# select best model
# mind the fact that if your data is spiked at zero, you have to be carefoul with the hazard-rate model (details in Buckland et al. 2001)
best_modRad <- hn_herm2

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
FOV <- 55 # field of view of the cameras (degrees)
w_ang <- FOV/2*0.0174533

# half-normal
hn_cos0Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "cos", order = 0, truncation=w_ang) 
hn_cos2Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "cos", order = 2, truncation=w_ang)
hn_herm0Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "herm", order = 0, truncation=w_ang)
hn_herm2Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "herm", order = 2, truncation=w_ang)
hn_poly0Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "poly", order = 0, truncation=w_ang)
hn_poly2Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hn", adjustment = "poly", order = 2, truncation=w_ang)

#hazard-rate
hr_cos0Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "cos", order = 0, truncation=w_ang) 
hr_cos2Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "cos", order = 2, truncation=w_ang)
hr_herm0Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "herm", order = 0, truncation=w_ang)
hr_herm2Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "herm", order = 2, truncation=w_ang)
hr_poly0Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "poly", order = 0, truncation=w_ang)
hr_poly2Ang <- ds(data_dz_ang$Ang_rad, transect = "line", key="hr", adjustment = "poly", order = 2, truncation=w_ang)

#model comparison
AIC(hn_cos0Ang, hn_cos2Ang, hn_herm0Ang, hn_herm2Ang, hn_poly0Ang, hn_poly2Ang, hr_cos0Ang, hr_cos2Ang, hr_herm0Ang, hr_herm2Ang, hr_poly0Ang, hr_poly2Ang)
AIC(hn_cos0Ang,  hn_herm0Ang, hn_herm2Ang, hn_poly0Ang, hn_poly2Ang)

# select best model
# mind the fact that if your data is spiked at zero, you have to be carefoul with the hazard-rate model (details in Buckland et al. 2001)
best_modAng <- hn_herm0Ang  

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

data.dens <- as.data.frame(data)
#data.dens <- subset(data.dens, Point_ID != 48 & Point_ID != 100)# to remove some points
# Load the operatity matrix to estimate survey effort
operat <- read.table("OperativityLiebresQM_REM.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE)
tm <- sum(colSums(operat, na.rm = FALSE, dims = 1))- sum(operat$CAM) # Survey effort (camera days)

data.dens <- subset(data.dens, Dist_det < 10.1 | is.na(Dist_det)) # Based on truncation distance to estime EDD. NAs are new sequences
seq <- length(data.dens[, 1]) # Number of sequences 

operat$oper_days<-rowSums(operat, dims = 1)-operat$CAM
seq_point<-data.frame(table(data.dens$Point_ID))
operat <- merge(operat, seq_point, by.x= "CAM", by.y = "Var1", all.x = TRUE); operat[is.na(operat)] <- 0
operat$tr <- operat$Freq/operat$oper_days

library(dplyr)
tr<-operat[,c("Freq","oper_days")] # Selecting columns
tr <- subset(tr, oper_days > 0) # remove CT that doesnt work propoerly

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
dimnames(results) <- list("Value", c("seq","tm", "dr(km/day)","dr_se(km/day)", "r(m)", "r_se(m)", "ang(rad)", "ang_se(rad)", "d(ind/km2)", "d_se(ind/km2)")); View(results)
#write.table(results, "ResultsREM_MARTESsp_FITO.txt", sep=";", row.names = FALSE)
