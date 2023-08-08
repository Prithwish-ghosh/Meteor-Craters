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
library(CircStats)
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
fireball$Longitude = -1 * (fireball$long)
head(fireball)
tail(fireball)
fireball = fireball[1:769,]
tail(fireball, 20)
fireball = data.frame(fireball$Peak.Brightness.Date.Time..UT. , fireball$Altitude..km. ,
                      fireball$Velocity..km.s. , fireball$vx , fireball$vy , fireball$vz , 
                      fireball$Total.Radiated.Energy..J. , fireball$Calculated.Total.Impact.Energy..kt.,
                      fireball$Latitude , fireball$Longitude)
tail(fireball)
fireball$fireball.Longitude = circular(fireball$fireball.Longitude)
fireball$fireball.Latitude = circular(fireball$fireball.Latitude)
head(fireball)
sum(is.na(fireball$fireball.Latitude))
sum(is.na(fireball$fireball.Longitude))
dim(fireball)


#View(fireball)
final_data = fireball
dim(final_data)
tail(final_data)
watson.test(fireball$fireball.Latitude , alpha = 0.05 , dist = "vonmises")
watson.test(fireball$fireball.Longitude , alpha = 0.05 , dist = "vonmises")
rao.spacing.test(fireball$fireball.Longitude , alpha = 0.05)
rao.spacing.test(fireball$fireball.Latitude , alpha = 0.05)

library(maps)
world_coordinates = map_data("world")
ggplot() +
  geom_map(
    data = world_coordinates,map = world_coordinates,
    aes(long, lat , map_id = region),
    color = "black" , fill = "lightyellow"
  ) +
  geom_point(
    data = fireball , 
    aes(fireball$fireball.Longitude , fireball$fireball.Latitude , color = "red",),
    alpha = 1
  ) +
  theme(legend.position = "top")

library(mice)

imputed_data <- mice(final_data, method = "lasso.norm", m = 5, maxit = 50)
completed_data <- complete(imputed_data)
#is.numeric(completed_data[, -1]$long)
head(completed_data)
sum(is.na(completed_data))
completed_data$fireball.long = circular(completed_data$fireball.long)
completed_data$fireball.Latitude = circular(completed_data$fireball.Latitude)
watson.test(completed_data$fireball.Longitude , alpha = 0.05 , dist = "vonmises")
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

world_coordinates = map_data("world")
ggplot() +
  geom_map(
    data = world_coordinates,map = world_coordinates,
    aes(long, lat , map_id = region),
    color = "black" , fill = "lightyellow"
  ) +
  geom_point(
    data = dataset_imp , 
    aes(dataset_imp$Impact_dataset.LON...13 , dataset_imp$Impact_dataset.LAT...9 , color = "red",),
    alpha = 1
  ) +
  theme(legend.position = "top")

head(dataset_imp)
dataset_imp$Impact_dataset.LON...13 = circular(dataset_imp$Impact_dataset.LON...13)
dataset_imp$Impact_dataset.LAT...9 = circular(dataset_imp$Impact_dataset.LAT...9)

watson.test(dataset_imp$Impact_dataset.LON...13 , alpha = 0.1 , dist = "vonmises")
watson.test(dataset_imp$Impact_dataset.LAT...9 , alpha = 0.1 , dist = "vonmises")

plot(density.circular(dataset_imp$Impact_dataset.LON...13 , bw = 50))
plot(density.circular(dataset_imp$Impact_dataset.LAT...9 , bw = 50))

rao.spacing.test(dataset_imp$Impact_dataset.LON...13 , alpha = 0.05)
rao.spacing.test(dataset_imp$Impact_dataset.LAT...9 , alpha = 0.05)

watson.two.test(dataset_imp$Impact_dataset.LON...13 , dataset_imp$Impact_dataset.LAT...9 , alpha = 0.05)
watson.two.test(completed_data$fireball.Latitude , completed_data$fireball.long , alpha = 0.05)
rao.test(dataset_imp$Impact_dataset.LON...13 , dataset_imp$Impact_dataset.LAT...9 , alpha = 0.05)

pp1_vm = pp.plot(dataset_imp$Impact_dataset.LON...13 , ref.line = TRUE )
pp1_un = pp.unif.plot(dataset_imp$Impact_dataset.LON...13 , ref.line = TRUE)
pp2_vm = pp.plot(dataset_imp$Impact_dataset.LAT...9 , ref.line = TRUE)
pp2_un = pp.unif.plot(dataset_imp$Impact_dataset.LAT...9, ref.line = TRUE)

ff1_vm = pp.plot(completed_data$fireball.Longitude , ref.line = TRUE)
ff1_un = pp.unif.plot(completed_data$fireball.Longitude , ref.line = TRUE)
ff2_vm = pp.plot(completed_data$fireball.Latitude , ref.line = TRUE)
ff2_un = pp.unif.plot(completed_data$fireball.Latitude , ref.line = TRUE)

sum(is.na(completed_data))
sum(is.na(completed_data$fireball.Latitude))
View(completed_data)

x= completed_data$fireball.Latitude
y <- rvonmises(n=1000, mu=mean(completed_data$fireball.Latitude), kappa=est.kappa(completed_data$fireball.Latitude))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(completed_data$fireball.Latitude , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for fireball.Latitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

x= completed_data$fireball.Longitude
y <- rvonmises(n=1000, mu=mean(completed_data$fireball.Longitude), kappa=est.kappa(completed_data$fireball.Longitude))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(completed_data$fireball.Latitude , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for fireball.Longitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)


x= dataset_imp$Impact_dataset.LON...13
y <- rvonmises(n=1000, mu=mean(dataset_imp$Impact_dataset.LON...13), kappa=est.kappa(dataset_imp$Impact_dataset.LON...13))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(dataset_imp$Impact_dataset.LON...13 , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Crater longitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

x= dataset_imp$Impact_dataset.LAT...9
y <- rvonmises(n=1000, mu=mean(dataset_imp$Impact_dataset.LAT...9), kappa=est.kappa(dataset_imp$Impact_dataset.LAT...9))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(dataset_imp$Impact_dataset.LAT...9 , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Crater Latitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

