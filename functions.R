# common functions
# Given a matrix containing X and Y columns constituting an object border, 
# calculate the max feret diameter and return the points in the border that
# lie at on this diameter
find.max.feret.points <- function(points){
  
  # calc pairwise distances for half the points
  
  max.d <- 0
  p1 <- NA
  p2 <- NA
  for(i in seq(1, ceiling(nrow(points)/2), 3)){
    
    for(j in seq(floor(nrow(points)/2),  nrow(points)-1, 3)){
      d <- euclidean(points[i, "X"], points[i, "Y"], points[j, "X"], points[j, "Y"])
      if(d>max.d){
        max.d <- d
        p1 <- points[i,]
        p2 <- points[j,]
      }
    }
  }
  
  
  if(p1["X"] < p2["X"]){ p.min <- p1 } else { p.min <- p2}
  p.max <- ifelse(p2 == p.min, p1, p2)
  
  # angle of max Feret relative to vertical in image
  feretAngle <- LearnGeom::Angle(c(p.min["X"], p.min["Y"]+10), 
                                 c(p.min["X"], p.min["Y"]),  
                                 c(p.max["X"], p.max["Y"])) 
  
  list(max.d = max.d, 
       p1 = sf::st_point(x = c(p.min["X"], p.min["Y"])),
       p2 = sf::st_point(x = c(p.max["X"], p.max["Y"])),
       angle = feretAngle)
}

# convert degrees to radians
deg2rad <- function(deg)(deg * pi) / (180)

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
    return(data %>% dplyr::select(-w, -h, -(xmin:ymax)))
    
  } else {
    # Not 6 columns - interpret as a segment mask
    cat("Reading YOLO segment mask\n")
    data <- border.data %>%
      dplyr::group_by(Image, Object) %>%
      dplyr::mutate(stomata  = paste0("s", sprintf("%02d", Object))) %>%
      dplyr::ungroup() %>%
      dplyr::select(Image, stomata) %>%
      dplyr::distinct()
    
    # Create polygons from border
    polys <- border.data %>% 
      dplyr::group_by(Image, Object) %>%
      dplyr::reframe(poly.matrix = list(matrix(c(x, x[1], y, y[1]), ncol=2, byrow=F)))
    data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))
    
    # Find the max diameter points from the polygon
    data$feret <- lapply(data$polygons, \(x) find.max.feret.points(sf::st_coordinates(x)))
    data <- tidyr::unnest_wider(data, feret)
    
    return(data)
  }
  
  
}
