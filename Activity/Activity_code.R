###############################
##                           ##
##      Activity pattern     ##
##            &              ##
##       activity level      ##
##                           ##
##   from camera trap data   ##
##       Pablo Palencia      ##
##        04/11/2022         ##
##                           ##
###############################


# R code to estimate wildlife activity patterns and activity levels from camera
# trap data

# Working directory
# This will open an explorer window for easy browsing
setwd(choose.dir())

# Import the .txt table into R when the data are separated by ";"
# Load data frame
data_activity <- read.table("Data.txt", sep = ";", dec=".", header=TRUE, as.is=TRUE) # data frame

# Data frame structure:
# N_ind: Number of individuals recorded
# Time: Time (HH:MM:SS format) when animals were recorded

#install.packages("activity") # just in case you didnt install the package previously
library(activity) # load package

# Replicating each sequence based on the number of animals recorded
data_activity.r <- data_activity[rep(row.names(data_activity), data_activity$N_ind), 1:2] 

# Convert time of data to a numeric vector of radian time-of-day
data_activity.r$radtime <- gettime(data_activity.r$Time, "%H:%M:%S", "proportion")
radtime <- 2*pi*gettime(data_activity.r$Time, "%H:%M:%S", "proportion")


# fit activity model
actmod <- fitact(radtime, sample="data") 

# plot activity pattern
plot(actmod)

# activity level
actmod@act
