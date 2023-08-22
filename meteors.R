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
library(hrbrthemes)  # hrbrmstr themes
library(magick)      # For animation
library(mapproj)     # Needed for projection
theme_set(theme_ipsum())
library(NPCirc)
library(rgl)
library(colorspace)
library(mclust)
library(maps)
library(maptools)
#  Meteor lalnding(1) dataset 

library(readr)

lat_lon_to_xyz <- function(lat, lon, radius = 1) {
  lat_rad <- lat * pi / 180  # Convert latitude to radians
  lon_rad <- lon * pi / 180  # Convert longitude to radians
  
  x <- radius * cos(lat_rad) * cos(lon_rad)
  y <- radius * cos(lat_rad) * sin(lon_rad)
  z <- radius * sin(lat_rad)
  
  return(cbind(x, y, z))
}

vmf_density_grid = 
function(u, ngrid = 100) {
  # Translate to (0,180) and (0,360)
  u[,1] <- u[,1] + 90
  u[,2] <- u[,2] + 180
  res <- vmf.kerncontour(u, thumb = "none", den.ret = T, full = T,
                         ngrid = ngrid)
  
  # Translate back to (-90, 90) and (-180, 180) and create a grid of
  # coordinates
  ret <- expand.grid(Lat = res$lat - 90, Long = res$long - 180)
  ret$Density <- c(res$den)
  ret
}

########################################## Fire Ball Dataset #######################################################



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

latitude <- fireball$fireball.Latitude
longitude <- fireball$fireball.Longitude

xyz_cord = lat_lon_to_xyz(latitude, longitude)
xyz_cord

fishkent(xyz_cord)

fireball.densities <- vmf_density_grid(fireball[,c("fireball.Latitude",
                                                   "fireball.Longitude")],
                                       ngrid = 300);

world <- map_data("world")
g.fireball <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "grey90", fill = "grey80") +
  geom_point(data = fireball,
             mapping = aes(x = fireball.Longitude, y = fireball.Latitude),
             color = "red", alpha = .5, size = .5, stroke = 0.1) +
  geom_density_2d(data = fireball,
                  aes(x = fireball.Longitude, y = fireball.Latitude),
                  color = "#B266FF", alpha = 1) +
  geom_contour(data = fireball.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  scale_y_continuous(breaks = (-2:2) * 30, limits = c(-90, 90)) +
  scale_x_continuous(breaks = (-4:4) * 45, limits = c(-180, 180)) +
  coord_map("mercator")

g.fireball

g.fireball <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "grey90", fill = "grey80") +
  geom_point(data = fireball,
             mapping = aes(x = fireball.Longitude, y = fireball.Latitude),
             color = "red", alpha = .5, size = .5, stroke = 0.1) +
  geom_density_2d(data = fireball,
                  aes(x = fireball.Longitude, y = fireball.Latitude),
                  color = "#B266FF", alpha = 1) +
  geom_contour(data = fireball.densities, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  scale_y_continuous(breaks = (-2:2) * 30, limits = c(-90, 90)) +
  scale_x_continuous(breaks = (-4:4) * 45, limits = c(-180, 180)) +
  coord_map("orthographic", orientation = c(-10, 0, 0)) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) +
  scale_y_continuous(breaks = seq(-90, 90, 45)) +
  ggtitle("Orthographic Projection of Spherical Density", "Top / Front View") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.ontop = TRUE,
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_line(color = "black" ),
        panel.background = element_rect(fill = NA))
g.fireball


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




############################################### Impact Crater Dataset ####################################################




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

Crater = dataset_imp[, c(7 , 8)]

Crater.Density = vmf_density_grid(Crater, ngrid = 300)

g.Crater <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "grey90", fill = "grey80") +
  geom_point(data = Crater,
             mapping = aes(x = Crater[,1], y = Crater[,2]),
             color = "red", alpha = .5, size = .5, stroke = 0.1) +
  geom_density_2d(data = Crater,
                  aes(x = Crater[,1], y = Crater[,2]),
                  color = "#B266FF", alpha = 1) +
  geom_contour(data = Crater.Density, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  scale_y_continuous(breaks = (-2:2) * 30, limits = c(-90, 90)) +
  scale_x_continuous(breaks = (-4:4) * 45, limits = c(-180, 180)) +
  coord_map("mercator")

g.Crater

g.Crater <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "grey90", fill = "grey80") +
  geom_point(data = Crater,
             mapping = aes(x = Crater[,1], y = Crater[,2]),
             color = "red", alpha = .5, size = .5, stroke = 0.1) +
  geom_density_2d(data = Crater,
                  aes(x = Crater[,1], y = Crater[,2]),
                  color = "#B266FF", alpha = 1) +
  geom_contour(data = Crater.Density, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  scale_y_continuous(breaks = (-2:2) * 30, limits = c(-90, 90)) +
  scale_x_continuous(breaks = (-4:4) * 45, limits = c(-180, 180)) +
  #coord_map("mercator")+
  coord_map("orthographic", orientation = c(-20, 60, 0)) +
  ggtitle("Orthographic Projection of Spherical Density", "Top / Front View") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.ontop = TRUE,
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_line(color = "black"),
        panel.background = element_rect(fill = NA)) 


g.Crater


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



############################################### Meteor Landing Dataset ###################################################



dataset_landing = read_csv("meteorite-landings.csv")
dataset_landing = data.frame(dataset_landing)
dataset_landing = na.omit(dataset_landing)
library(maps)
library(ggplot2)
world_coordinates = map_data("world")
ggplot() +
  geom_map(
    data = world_coordinates,map = world_coordinates,
    aes(long, lat , map_id = region),
    color = "black" , fill = "lightyellow"
  ) +
  geom_point(
    data = dataset_landing ,
    aes(dataset_landing$reclong , dataset_landing$reclat , color = "red",),
    alpha = 0.4
  ) +
  theme(legend.position = "top")

library(circular)
watson.test(dataset_landing$reclong , alpha = 0.01 , dist = "vonmises")
watson.test(dataset_landing$reclat , alpha = 0.01 , dist = "vonmises")
watson.test(dataset_landing$reclong , alpha = 0.01 , dist = "uniform")
watson.test(dataset_landing$reclat , alpha = 0.01 , dist = "uniform")
watson.two(dataset_landing$reclong , dataset_landing$reclat , alpha = 0.01 , plot = TRUE)

dataset_landing = na.omit(dataset_landing)
library(movMF)
library(CircStats)

x= dataset_landing$reclong
y <- rvonmises(n=1000, mu=mean(dataset_landing$reclong), kappa=est.kappa(dataset_landing$reclong))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(dataset_landing$reclong , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1.5), ylim=c(-1.3, 1.5), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Meteor Landing.Longitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

x= dataset_landing$reclat
y <- rvonmises(n=1000, mu=mean(dataset_landing$reclat), kappa=est.kappa(dataset_landing$reclat))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(dataset_landing$reclat , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1.5), ylim=c(-1.3, 1.5), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Meteor Landing.Latitude data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)
ppunif(dataset_landing$reclong , col = "red" , col1 = "blue" , ref.line = TRUE)
ppunif(dataset_landing$reclat , col = "red" , col1 = "blue" , ref.line = TRUE)

pp(dataset_landing$reclong , col = "red" , col1 = "blue" , ref.line = TRUE)
pp(dataset_landing$reclat , col = "red" , col1 = "blue" , ref.line = TRUE)

dataset_landing = na.omit(dataset_landing)
lq = cbind(circular(dataset_landing$reclat) , circular(dataset_landing$reclong))
lq
sum(is.na(lan))
sum(is.null(lan))
sum(is.finite(lq))
dim(lq)

set.seed(2023)

rows_with_zeros <- apply(lq, 1, function(row) any(row == 0))
sum(rows_with_zeros)
dim(lq)
mat_no_zeros <- lq[!rows_with_zeros, ]

# Print the resulting matrix
print(mat_no_zeros)

vMFs_Landing <- 
  function(K){
    movMF(na.omit(mat_no_zeros), k = K, control= list(nruns = 20))
  }
sdl = lapply(1:20, vMFs_Landing)
sdl
sapply(sdl, BIC)

dataset_landing <- lq[!rows_with_zeros, ]
dataset_landing = as.matrix(dataset_landing) 
dataset_landing


num_samples <- 15000

# Randomly sample rows from the matrix
sampled_rows <- sample(1:nrow(dataset_landing), size = num_samples, replace = FALSE)

# Extract the sampled rows
Landing <- dataset_landing[sampled_rows, ]
dim(Landing)
Landing.dens = vmf_density_grid(Landing , ngrid = 300)

Landing = data.frame(Landing)

world <- map_data("world")

library(dplyr)

Landing_data_partitioned <- earthquake_data %>%
  mutate(partition = case_when(
    Landing$dataset_landing.reclat >= -65  ~ "Western Hemisphere",
    Landing$dataset_landing.reclat < 90
  ))

earthquake_data_partitioned

G_Asian_data <- Landing_data_partitioned %>%
  filter(partition == "Western Hemisphere")
                         
g.am_landing <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "grey90", fill = "grey90") +
  geom_point(data = Landing,
             mapping = aes(x = dataset_landing.reclong, y = dataset_landing.reclat),
             color = "red", alpha = .5, size = 1, stroke = 0.1) +
  geom_density_2d(data = G_Asian_data,
                  aes(x = dataset_landing.reclong, y = dataset_landing.reclat),
                  color = "#b266ff", alpha = 2) +
  geom_contour(data = south_a.dens, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  geom_density_2d(data = Landing , aes(x = dataset_landing.reclong, y = dataset_landing.reclat),
                  color = "darkgreen")+
  geom_contour(data = Landing.dens , aes(x=Long, y=Lat, z=Density),
               color = "deeppink4")+
  coord_map("mercator")

g.am_landing

                         
g.am_landing <- ggplot() +
  geom_map(data = world, map = world,
           mapping = aes(map_id = region),
           color = "grey90", fill = "grey90") +
  geom_point(data = Landing,
             mapping = aes(x = dataset_landing.reclong, y = dataset_landing.reclat),
             color = "red", alpha = .5, size = 1, stroke = 0.1) +
  geom_density_2d(data = G_Asian_data,
                  aes(x = dataset_landing.reclong, y = dataset_landing.reclat),
                  color = "#b266ff", alpha = 2) +
  geom_contour(data = south_a.dens, aes(x=Long, y=Lat, z=Density),
               color = "blue") +
  geom_density_2d(data = Landing , aes(x = dataset_landing.reclong, y = dataset_landing.reclat),
                  color = "darkgreen")+
  geom_contour(data = Landing.dens , aes(x=Long, y=Lat, z=Density),
               color = "deeppink4")+
  coord_map("orthographic", orientation = c(-10, 0, 0)) +
  scale_x_continuous(breaks = seq(-180, 180, 20)) +
  scale_y_continuous(breaks = seq(-90, 90, 45)) +
  ggtitle("Orthographic Projection of Spherical Density", "Top / Front View") +
  xlab("") +
  ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.ontop = TRUE,
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_line(color = "black" ),
        panel.background = element_rect(fill = NA))

g.am_landing             
