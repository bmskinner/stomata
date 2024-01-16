# analyse stomata CoMs to make a network

library(tidyverse)
library(rspatial)
library(raster)

in.file = "D:/git/stomata/stomatal_image/September 11 (part 1)- 55 samples/output/CoMs.txt"
data <- read.csv(in.file, header = T, sep = "\t")

test.img.data <- data %>% dplyr::filter(image=="1047-2-1.jpg")
pts <- test.img.data[,2:3]

ggplot(test.img.data, aes(x=x, y=y))+
  geom_point()


#### Homebrew method via https://rspatial.org/raster/analysis/8-pointpat.html ####

# for an image, calculate the density of points with a resolution 
# of 25 pixels ish
count.freq <- function(image.data){
  
  # Create a raster grid covering the image extent
  image.bounds <- raster(nrows=10, ncols=10, xmn=0, xmx=2592, ymn=0, ymx=1944, vals=0)
  # image.bounds <- rasterize(as.matrix(pts), image.bounds)
  # plot(image.bounds)
  # quads <- as(image.bounds, 'SpatialPolygons')
  # plot(quads, add=TRUE)
  # points(image.data, col='red', cex=.5)
  
  nc <- rasterize(image.data, image.bounds, fun='count', background=0)
  # plot(nc)
  
  nstom <- mask(nc, image.bounds)
  freq(nstom, useNA='no')
}

# For each image, calculate the density of points
results <- do.call(rbind, lapply(unique(data$image), function(i) count.freq(as.matrix(data[data$image==i,2:3]))))

# Aggregate and display
summ <- as.data.frame(results) %>% dplyr::group_by(value) %>%
  summarise(count = sum(count))
plot(summ)
lines(summ)


# The G-function - the cumulative distribution of the distances from events to their nearest neighboring event.
g.function <- function(image.data, image){
  # Calculate neighbour distances
  d <- dist(image.data)
  dm <- as.matrix(d)
  diag(dm) <- NA
  dmin <- apply(dm, 1, min, na.rm=TRUE)

  # get the unique distances (for the x-axis)
  distance <- sort(unique(round(dmin)))
  # compute how many cases there with distances smaller that each x
  Gd <- sapply(distance, function(x) sum(dmin < x))
  # normalize to get values between 0 and 1
  Gd <- Gd / length(dmin)

  data.frame(distance, Gd, image)
}

# K-function is designed to detect patterns at multiple scales (see Ripley 1976; and Haase 1995).
k.function <- function(image.data, image){
  d <- dist(image.data)
  distance <- seq(1, 3000, 100)
  Kd <- sapply(distance, function(x) sum(d < x)) # takes a while
  Kd <- Kd / (length(Kd))
  
  data.frame(distance, Kd, image)
}

g.results <- do.call(rbind, lapply(unique(data$image), function(i) g.function(as.matrix(data[data$image==i,2:3]), i)))
k.results <- do.call(rbind, lapply(unique(data$image), function(i) k.function(as.matrix(data[data$image==i,2:3]), i)))

# Plot G-function values across all images
ggplot(g.results, aes(x=distance, y=Gd, group=image))+
  geom_line(col="grey")+
  labs(x = "Distance between stomata", y = "Cumulative fraction")+
  theme_bw()

ggsave("g_function_stomata_distances.png", dpi=300, width = 170, height = 170, units = "mm")

# Plot K-function values across all images
ggplot(k.results, aes(x=distance, y=Kd, group=image))+
  geom_line(col="grey")+
  theme_bw()



# plot(distance, Gd, type='l', lwd=2, col='red', las=1,
#      ylab='F(d) or G(d)', xlab='Distance', yaxs="i", xaxs="i")
# lines(Fdistance, Fd, lwd=2, col='blue')
# lines(0:2000, expected, lwd=2)
# legend(1200, .3,
#        c(expression(italic("G")["d"]), expression(italic("F")["d"]), 'expected'),
#        lty=1, col=c('red', 'blue', 'black'), lwd=2, bty="n")
# 

#### Spatstat method ########