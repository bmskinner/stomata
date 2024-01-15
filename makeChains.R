# Test ideas for making chains of stomata
library(tidyverse)
library(magrittr)

MAX.DISTANCE <- 200

files <- list.files(path = "stomatal_images", pattern = "*.csv",  recursive = T, full.names = T)

data <- read.csv(files[1], header = T)

# Create the distance matrix for an image
euclidean <- function(x1, y1, x2, y2) sqrt( (x1-x2)^2 + (y1-y2)^2)

for(i in data$X.1){
  x <- data$X[i]
  y <- data$Y[i]
  
  calc.dists <- function(j) euclidean(x, y, data$X[j], data$Y[j])
  data[[paste0("d", i)]] <- sapply(data$X.1, calc.dists)
}

# Reformat long
filt <- tidyr::pivot_longer(data, cols = starts_with("d"), names_to = "Comparator") %>%
 dplyr::filter(value>0)
filt$stomata <- paste0("d", filt$X.1)


# Start making chains

# Calculate angle between three points
# point 2 is the centre
absolute.angle <- function(x1, y1, x2, y2, x3, y3){
  ab.x <- x2 - x1
  ab.y <- y2 - x1
  
  cb.x <- x2 - x3
  cb.y <- y2 - y3
  
  dot <-  ab.x * cb.x + ab.y * cb.y
  cross <- ab.x * cb.y - ab.y * cb.x
  
  alpha <- atan2(cross, dot)
  angle <- alpha * 180 / pi

  return(angle)
}

chains <- list()
for(i in unique(filt$stomata)){
  # Get the viable nearest stomata to link to
  sub.data <- filt %>% filter(stomata==i, value <= MAX.DISTANCE)

  # Get two levels of chain from current stomata
  mge <- merge(sub.data, filt, by.x = "Comparator", by.y = "stomata", all.y = F) %>%
    dplyr::filter(value.y <= MAX.DISTANCE) %>%
    merge(., filt, by.x = "Comparator.y", by.y = "stomata", all.y = F) %>%
    dplyr::select(-Comparator.y.y, -value) %>%
    dplyr::distinct()
  
  # Check angle is valid
  mge$theta <- mapply(absolute.angle, mge$X.x, mge$Y.x, mge$X.y, mge$Y.y, mge$X, mge$Y)
  
  mge %<>% dplyr::filter( (theta > 170 & theta < 190) |  (theta < -170 & theta > -190))
  
  if(nrow(mge)>0) {
    print("Found a chain")
    chains[[i]] <- c(mge$stomata, mge$Comparator.x, mge$Comparator.y)
  }
}


# Use a set of merges to get the three coordinates for angle calculation

# From a stomata CoM coordinate:
# Find any object within x radius
# Calculate line between the two points
# Project onwards with e.g. 5 degree arc a further x distance.
# Keep any objects that intersected.
# Calculate line between the two points
# Project onwards with e.g. 5 degree arc a further x distance.
# Continue to edge of image