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
library(sf)
# options(rgl.useNULL = TRUE)  # For Macbook
library(rgl)
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

fire = cbind(fireball$fireball.Latitude , fireball$fireball.Longitude)
set.seed(2022)
EvMFs <- 
  function(K){
    movMF(fire, k = K, control= list(nruns = 20))
  }

Esd = lapply(1:10, EvMFs)
gt = sapply(Esd, BIC)
gt

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
    color = "black" , fill = "lightblue"
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

cra = cbind(dataset_imp$Impact_dataset.LAT...9 , dataset_imp$Impact_dataset.LON...13)
set.seed(2023)
vMFs <- 
  function(K){
    movMF(cra, k = K, control= list(nruns = 20))
  }

sd = lapply(1:10, vMFs)
gtr = sapply(sd, BIC)
gtr

#library(Directional)
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


pp = function (x,col,col1, ref.line = TRUE )
{
  n <- length(x)
  x <- sort(x%%(2 * pi))
  z <- c(1:n)/(n + 1)
  mu <- circ.mean(x)%%(2 * pi)
  kappa <- est.kappa(x)
  y <- c(1:n)
  for (i in 1:n) {
    y[i] <- pvm(x[i], mu, kappa)
  }
  plot(z, y, xlab = "Von Mises Distribution", ylab = "Empirical Distribution" , col = col , lwd = 0.75)
  if (ref.line) 
    abline(0, 1 , col = col1 , lwd = 3.0)
  data.frame(mu, kappa)
}

pp(completed_data$fireball.Longitude , col = "red" , col1 = "blue" , ref.line = TRUE)
pp(completed_data$fireball.Latitude , col = "red" , col1 = "blue", ref.line = TRUE)

ppunif = function (x, ref.line = TRUE, frac = NULL, xlab = "Uniform Distribution", 
                   ylab = "Empirical Distribution", col = NULL, col.inf = NULL, 
                   col.sup = NULL, col1,...) 
{
  x <- na.omit(x)
  if (length(x) == 0) {
    warning("No observations (at least after removing missing values)")
    return(NULL)
  }
  x <- conversion.circular(x, units = "radians", zero = 0, 
                           rotation = "counter", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  y <- sort(x%%(2 * pi))/(2 * pi)
  n <- length(y)
  z <- (1:n)/(n + 1)
  if (is.null(col)) 
    col <- rep(1, n)
  else col <- rep(col, length.out = n)
  if (!is.null(frac)) {
    if (!is.numeric(frac) || (frac < 0 | frac > 1)) {
      stop("'frac' must be in the interval [0,1]")
    }
    f <- round(frac * n)
    if (f) {
      zm <- -1 + ((n - f + 1):n)/(n + 1)
      zp <- 1 + (1:f)/(n + 1)
      ym <- -1 + y[(n - f + 1):n]
      yp <- 1 + y[1:f]
      y <- c(ym, y, yp)
      z <- c(zm, z, zp)
      if (is.null(col.inf)) 
        col.inf <- rep(2, f)
      else col.inf <- rep(col.inf, length.out = f)
      if (is.null(col.sup)) 
        col.sup <- rep(2, f)
      else col.sup <- rep(col.sup, length.out = f)
      col <- c(col.inf, col, col.sup)
    }
  }
  plot.default(z, y, xlab = xlab, ylab = ylab, col = col, ...)
  if (ref.line) {
    abline(0, 1 , col = col1, lwd = 3.0)
    if (!is.null(frac)) {
      abline(h = c(0, 1), lty = 3)
      abline(v = c(0, 1), lty = 3)
    }
  }
}

ppunif(completed_data$fireball.Longitude , col = "red" , col1 = "blue")
ppunif(completed_data$fireball.Latitude , col = "red" , col1 = "blue")

pp(dataset_imp$Impact_dataset.LON...13 , col = "red" , col1 = "blue")
pp(dataset_imp$Impact_dataset.LAT...9, col = "red", col1 = "blue")

ppunif(dataset_imp$Impact_dataset.LON...13 , col = "red" , col1 = "blue")
ppunif(dataset_imp$Impact_dataset.LAT...9 , col = "red" , col1 = "blue")
