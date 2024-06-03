# common functions
library(jsonlite)
library(tidyverse)
library(magrittr)
library(LearnGeom)
library(text.alignment)
library(igraph)
library(grid)
library(patchwork)
library(OpenImageR)
library(ggbeeswarm)
library(autoimage)
library(sf)
library(patchwork)
library(fs)
library(lwgeom)
library(parallel)
library(parallelsugar) # github 'nathanvan/parallelsugar', provides windows mclapply

# Save the given plot to the given path
save.ggplot <- function(plot, out.file, width = 170, height = 85){
  ggsave(out.file, plot = plot, dpi = 300, units = "mm", width = width, height = height)
}

# Save the given plot to the given path and return the plot
save.plot <- function(plot, image.file, suffix, width = 170, height = 85){
  out.path <- str_replace(image.file, ".jpg", suffix)
  out.path <- str_replace(out.path, "/", "_")
  out.path <- paste0("figure/", out.path)
  save.ggplot(out.file = out.path, plot = plot, width = width, height = height)
  plot
}


# Given two points in a list, return the leftmost (lowest x)
left <- function(p1, p2){
  if(p1$X<p2$X) return(p1)
  else(return(p2))
}
# Given two in a list, return the rightmost (highest x)
right <- function(p1, p2){
  if(p1$X>=p2$X) return(p1)
  else(return(p2))
}
# Given two in a list, return the upper (highest y)
upper <- function(p1, p2){
  if(p1$Y>p2$Y) return(p1)
  else(return(p2))
}
# Given two in a list, return the lower (lowest y)
lower <- function(p1, p2){
  if(p1$Y<=p2$Y) return(p1)
  else(return(p2))
}

# convert degrees to radians
deg2rad <- function(deg)(deg * pi) / (180)

# convert radians to degrees
rad2deg <- function(rad) (180 * rad) / pi

# Calculate the angle of line described by the given points to vertical
# Returns angle in degrees anticlockwise. The line will start from the upper point
angle.to.vertical <- function(x1, y1, x2, y2){
  a = 0; # same x coords from vertical
  b = 1; # arbitrary y offset
  c = x2 - x1;
  d = y2 - y1;
  
  atanA = atan2(a, b);
  atanB = atan2(c, d);
  
  rad2deg(atanA-atanB)
}

# Calculate the angle of line described by the given points to the horizontal
# Returns angle in degrees anticlockwise.
angle.to.horizontal <- function(x1, y1, x2, y2){

  a = 1; # same x coords from vertical
  b = 0; # arbitrary y offset
  c = x2 - x1;
  d = y2 - y1;
  
  atanA = atan2(a, b);
  atanB = atan2(c, d);
  
  rad2deg(atanA-atanB)
}


# Given a matrix containing X and Y columns constituting an object border, 
# calculate the max feret diameter and return the points in the border that
# lie at on this diameter
calculate.stomata.orientation <- function(points){
  
  # calc pairwise distances for half the points
  max.d <- 0
  p1 <- NA
  p2 <- NA
  
  half.border <- ceiling(nrow(points)/2)-1
  # Assume halfway round the border indexes is halfway around the border
  for(i in seq(1, half.border, 2)){
    
    # vary the check for a few points around the half border
    for(j in seq(i+half.border-5, i+half.border+5, 2)){
      
      if(j==i) next
      
      j <- max(1, j %% (nrow(points)))
      
      d <- euclidean(points[i, "X"], points[i, "Y"], points[j, "X"], points[j, "Y"])
      if(d>max.d){
        max.d <- d
        p1 <- points[i,]
        p2 <- points[j,]
      }
    }
  }

  list(max.feret = max.d, 
       p1    = list(X = p1["X"], "Y"=p1["Y"]),
       p2    = list(X = p2["X"], "Y"=p2["Y"]),
       dx    = abs(p2["X"] - p1["X"]),
       dy    = abs(p2["Y"] - p1["Y"]),
       x.com = mean(points[,"X"]),
       y.com = mean(points[,"Y"]))
}

# Find the maximum value in the density plot of the given vector
find.mode <-  function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# Calculate distance between two points
euclidean <- function(x1, y1, x2, y2) sqrt( (x1-x2)^2 + (y1-y2)^2)

# Read object borders from JSON
read.border.from.json <- function(file){
  data <- jsonlite::read_json(file)
  shape.number = 1
  process.shape <- function(shape){
    if(shape$shape_type!="polygon")  return(NULL)
    
    bounds <- split(unlist(shape$points, recursive = T), 1:2) %>% as.data.frame
    colnames(bounds) <- c("X", "Y")
    result <- data.frame("shape" = shape.number,
                         "x" = bounds$X,
                         "y" = bounds$Y,
                         "file" = file,
                         "folder" = dirname(file),
                         "imageName" = str_replace(basename(file), ".json", ".jpg"))
    shape.number <<- shape.number+1
    return(result)
  }
  
  
  do.call(rbind, lapply(data$shapes, process.shape)) %>% 
    as.data.frame %>%
    dplyr::mutate("imageWidth" = data$imageWidth,
                  "imageHeight" = data$imageHeight)
}

# Read object borders from YOLO predictions
# Convert outline to sf st_polygon
read.border.from.yolo <- function(file){
  
  border.data <- readr::read_tsv(file, col_names = TRUE, progress = FALSE, 
                                 col_types = cols())
  
  # Is this a bbox output, or a segment output?
  # bbox will have 6 columns, seg will have 4
  
  if(ncol(border.data)==6){
    cat("Reading YOLO bbox\n")
    data <- border.data %>%
      dplyr::group_by(Image) %>%
      dplyr::mutate(stomata  = paste0("s", sprintf("%02d", row_number()))) %>%
      dplyr::ungroup() %>%
      dplyr::select(stomata, x.com=x, y.com=y, w, h) %>%
      dplyr::mutate(xmin=x.com-w/2, xmax=x.com+w/2,
                    ymin=y.com-h/2, ymax=y.com+h/2)
    as.data.frame
    
    # Create polygons from border
    polys <- data %>%
      dplyr::group_by(stomata) %>%
      dplyr::rowwise() %>%
      dplyr::summarise(poly.matrix = list(matrix(c(xmin, xmin, xmax, xmin, xmin,
                                                   ymin, ymax, ymax, ymax, ymin), 
                                                 ncol=2, byrow=F)), .groups = 'drop')
    data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))
    data$feret <- lapply(data$polygons, \(x) find.max.feret.points(sf::st_coordinates(x)))
    data <- data %>%
      tidyr::unnest_wider(feret)

    return(data %>% dplyr::select(-w, -h, -(xmin:ymax)))
    
  } else {
    # Not 6 columns - interpret as a segment mask
    cat("Reading YOLO segment mask\n")
    data <- border.data %>%
      dplyr::group_by(Image, Object) %>%
      dplyr::mutate(stomata  = paste0("s", sprintf("%03d", Object))) %>%
      dplyr::ungroup() %>%
      dplyr::select(Image, stomata) %>%
      dplyr::distinct()
    
    # Create polygons from border
    cat("Creating polygons from segment mask\n")
    polys <- border.data %>% 
      dplyr::group_by(Image, Object) %>%
      dplyr::reframe(poly.matrix = list(matrix(c(x, x[1], y, y[1]), ncol=2, byrow=F)))
    data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))
    
    # Find the max diameter points from the polygon and angles to vertical/horizontal
    cat("Finding stomata orientations\n")
    data$feret <- lapply(data$polygons, \(x) calculate.stomata.orientation(sf::st_coordinates(x)))
    data <- data %>%
      tidyr::unnest_wider(feret) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(area      = sf::st_area(polygons),
                    perimeter = sf::st_perimeter(polygons),
                    circularity = 4 * pi * area / perimeter^2
                    )
      
    # Remove objects that are too small or too irregular to be stomata
    min.area <- median(data$area)/2
    max.area <- median(data$area)*2
    min.circ <- min(0.5, median(data$circularity)/2)
    data <- data %>%
      dplyr::filter( between(area, min.area, max.area) & circularity > min.circ) %>%
      dplyr::mutate(Folder = dirname(Image),
                    File   = basename(Image))

    return(data)
  }
  
  
}
