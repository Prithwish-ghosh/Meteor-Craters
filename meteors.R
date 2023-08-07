#
#
# The packages which are needed for our analysis
#
#
#install.packages("gofest" , "Directional" , "circular", "protr","tidyverse",
#                 "leaflet" , "mixtools" , "spam", "ggplot2" , "SpatialVx" , 
#                 "skmeans" , "NPCirc" , "rgl", "flexmix", "slam", "mclust",
#                 "lattice", "mvdalab", "Amelia", "fields", "maps", "maptools",
#                 "tm")
library(Directional)
library(circular)
library(tidyverse)
library(Matrix)
library(mixtools)
library(ggplot2)
library(skmeans)
library(NPCirc)
library(rgl)
library(mclust)
library(maps)
library(maptools)
#  Meteor lalnding(1) dataset 

library(readr)

cneos_fireball_data_final <- read_csv("cneos_fireball_data (2).csv")
fireball = data.frame(cneos_fireball_data_final)
head(fireball)
View(fireball)
fireball$Longitude = -fireball$b * (fireball$Longitud..dg..)
head(fireball)
fireball = data.frame(fireball$Peak.Brightness.Date.Time..UT. , fireball$Altitude..km. ,
                      fireball$Velocity..km.s. , fireball$vx , fireball$vy , fireball$vz , 
                      fireball$Total.Radiated.Energy..J. , fireball$Calculated.Total.Impact.Energy..kt.,
                      fireball$Latitude , fireball$Longitude)
tail(fireball)
View(fireball)
fireball = fireball[1:769,]
head(fireball)
sum(is.na(fireball$Latitude..deg..))
dim(fireball)
#View(fireball)
final_data = fireball
dim(final_data)
tail(final_data)
watson.test(fireball$fireball.Latitude , alpha = 0.05 , dist = "vonmises")
watson.test(fireball$fireball.Longitude , alpha = 0.05 , dist = "vonmises")


library(maps)
library(mice)

imputed_data <- mice(final_data, method = "lasso.norm", m = 5, maxit = 50)
completed_data <- complete(imputed_data)
#is.numeric(completed_data[, -1]$long)
head(completed_data)
sum(is.na(completed_data$final_data.vx))
completed_data$fireball.long = circular(completed_data$fireball.long)
completed_data$fireball.Latitude = circular(completed_data$fireball.Latitude)
watson.test(completed_data$fireball.long , alpha = 0.05 , dist = "vonmises")
watson.test(completed_data$fireball.Latitude , alpha = 0.05 , dist = "vonmises")

plot(density.circular(completed_data$fireball.long , bw = 50) , pch = 19)
plot(density.circular(completed_data$fireball.Latitude , bw = 50) , pch = 19)

Impact_dataset = read_csv("Impact.csv")
Impact_dataset = data.frame(Impact_dataset)
library(Directional)
head(Impact_dataset)
dataset_imp = data.frame(Impact_dataset$Name , Impact_dataset$Location , Impact_dataset$Country ,
                         Impact_dataset$Diameter, Impact_dataset$Age , Impact_dataset$Coordinates, 
                         Impact_dataset$LAT...9 , Impact_dataset$LON...13)
dataset_imp = na.omit(dataset_imp)
head(dataset_imp)
dataset_imp$Impact_dataset.LON...13 = circular(dataset_imp$Impact_dataset.LON...13)
dataset_imp$Impact_dataset.LAT...9 = circular(dataset_imp$Impact_dataset.LAT...9)

watson.test(dataset_imp$Impact_dataset.LON...13 , alpha = 0.05 , dist = "vonmises")
watson.test(dataset_imp$Impact_dataset.LAT...9 , alpha = 0.05 , dist = "vonmises")

plot(density.circular(dataset_imp$Impact_dataset.LON...13 , bw = 50))
plot(density.circular(dataset_imp$Impact_dataset.LAT...9 , bw = 50))

