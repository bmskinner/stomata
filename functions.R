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
library(IRanges)

# Save the given plot to the given path
save.ggplot <- function(plot, out.file, width = 170, height = 85) {
  ggsave(out.file, plot = plot, dpi = 300, units = "mm", width = width, height = height)
}

# Save the given plot to the given path and return the plot
# Size for powerpoint slide: width= 253.8, height = 190.05
save.plot <- function(plot, image.file, suffix, width = 253.8, height = 190.05) { #
  out.path <- str_replace(image.file, ".jpg", suffix)
  out.path <- str_replace(out.path, "/", "_")
  out.path <- paste0("figure/", out.path)
  save.ggplot(out.file = out.path, plot = plot, width = width, height = height)
  plot
}


# Given two points in a list, return the leftmost (lowest x)
left <- function(p1, p2) {
  if (p1$X < p2$X) {
    return(p1)
  } else {
    (return(p2))
  }
}
# Given two points in a list, return the rightmost (highest x)
#
# p1 - a list with named elements x and y
# p2 - a list with named elements x and y
right <- function(p1, p2) {
  if (p1$X >= p2$X) {
    return(p1)
  } else {
    (return(p2))
  }
}
# Given two points in a list, return the upper (highest y)
#
# p1 - a list with named elements x and y
# p2 - a list with named elements x and y
upper <- function(p1, p2) {
  if (p1$Y > p2$Y) {
    return(p1)
  } else {
    (return(p2))
  }
}
# Given two in a list, return the lower (lowest y)
#
# p1 - a list with named elements x and y
# p2 - a list with named elements x and y
lower <- function(p1, p2) {
  if (p1$Y <= p2$Y) {
    return(p1)
  } else {
    (return(p2))
  }
}

# convert degrees to radians
deg2rad <- function(deg) (deg * pi) / (180)

# convert radians to degrees
rad2deg <- function(rad) (180 * rad) / pi

# Calculate the angle of line described by the given points to vertical
# Returns angle in degrees anticlockwise. The line will start from the upper point
angle.to.vertical <- function(x1, y1, x2, y2) {
  a <- 0 # same x coords from vertical
  b <- 1 # arbitrary y offset
  c <- x2 - x1
  d <- y2 - y1

  atanA <- atan2(a, b)
  atanB <- atan2(c, d)

  rad2deg(atanA - atanB)
}

# Calculate the angle of line described by the given points to the horizontal
# Returns angle in degrees anticlockwise.
angle.to.horizontal <- function(x1, y1, x2, y2) {
  a <- 1 # same x coords from vertical
  b <- 0 # arbitrary y offset
  c <- x2 - x1
  d <- y2 - y1

  atanA <- atan2(a, b)
  atanB <- atan2(c, d)

  rad2deg(atanA - atanB)
}

# Given a matrix containing X and Y columns constituting an bounding box border,
# identify the longest axis and return points that lie run on this axis through
# the centre of mass
calculate.bounding.box.orientation.points <- function(points) {
  # Get the CoM of the object. Marker points should align to this point.
  x.com <- mean(points[, "X"])
  y.com <- mean(points[, "Y"])

  p1 <- points[1, ]
  p2 <- points[2, ]
  p3 <- points[3, ]
  p4 <- points[4, ]

  # What are the lengths of the two axes in the rectangle?

  axis1 <- euclidean(p1["X"], p1["Y"], p2["X"], p2["Y"])
  axis2 <- euclidean(p1["X"], p1["Y"], p4["X"], p4["Y"])

  # Select the points lying on the longer axis
  pA <- p1

  if (axis1 > axis2) {
    pB <- p2
  } else {
    pB <- p4
  }
  long.axis <- ifelse(axis1 > axis2, axis1, axis2)

  cat(pB, "\n")

  list(
    max.feret = long.axis,
    p1 = list(X = pA["X"], "Y" = pA["Y"]),
    p2 = list(X = pB["X"], "Y" = pB["Y"]),
    dx = abs(pB["X"] - pA["X"]),
    dy = abs(pB["Y"] - pA["Y"]),
    x.com = x.com,
    y.com = y.com
  )
}


# Given a matrix containing X and Y columns constituting an object border,
# calculate the max feret diameter and return the points in the border that
# lie at on this diameter
calculate.stomata.orientation <- function(points) {
  # calc pairwise distances for half the points
  max.d <- 0
  p1 <- NA
  p2 <- NA

  half.border <- ceiling(nrow(points) / 2) - 1
  # Assume halfway round the border indexes is halfway around the border
  for (i in seq(1, half.border, 2)) {
    # vary the check for a few points around the half border
    for (j in seq(i + half.border - 5, i + half.border + 5, 2)) {
      if (j == i) next

      j <- max(1, j %% (nrow(points)))

      d <- euclidean(points[i, "X"], points[i, "Y"], points[j, "X"], points[j, "Y"])
      if (d > max.d) {
        max.d <- d
        p1 <- points[i, ]
        p2 <- points[j, ]
      }
    }
  }

  list(
    max.feret = max.d,
    p1 = list(X = p1["X"], "Y" = p1["Y"]),
    p2 = list(X = p2["X"], "Y" = p2["Y"]),
    dx = abs(p2["X"] - p1["X"]),
    dy = abs(p2["Y"] - p1["Y"]),
    x.com = mean(points[, "X"]),
    y.com = mean(points[, "Y"])
  )
}

# Find the maximum value in the density plot of the given vector
find.mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# Calculate distance between two points
euclidean <- function(x1, y1, x2, y2) sqrt((x1 - x2)^2 + (y1 - y2)^2)

# Read object borders from JSON
read.border.from.json <- function(file) {
  data <- jsonlite::read_json(file)
  shape.number <- 1
  process.shape <- function(shape) {
    if (shape$shape_type != "polygon") {
      return(NULL)
    }

    bounds <- split(unlist(shape$points, recursive = T), 1:2) %>% as.data.frame()
    colnames(bounds) <- c("X", "Y")
    result <- data.frame(
      "shape" = shape.number,
      "x" = bounds$X,
      "y" = bounds$Y,
      "file" = file,
      "folder" = dirname(file),
      "imageName" = str_replace(basename(file), ".json", ".jpg")
    )
    shape.number <<- shape.number + 1
    return(result)
  }

  do.call(rbind, lapply(data$shapes, process.shape)) %>%
    as.data.frame() %>%
    dplyr::mutate(
      "imageWidth" = data$imageWidth,
      "imageHeight" = data$imageHeight
    )
}



# Read object borders from YOLO predictions
# Convert outlines to sf::st_polygon objects
#
# file - the YOLO inferencing results text file
read.border.from.yolo <- function(file) {
  border.data <- readr::read_tsv(file,
    col_names = TRUE, progress = FALSE,
    col_types = cols()
  )

  # There are three possible input formats we need to handle:
  # bounding box (bbox) - a rectangle aligned with horizontal and vertical image axis
  # oriented bounding box (OBB) - a rectangle at an angle to best fit the object
  # segmentation - arbitrary number of points making an outline
  # The functions here handle each object type

  # bbox outputs have 6 columns with 1 row per object: Image, Object, x, y, w, h
  read.bbox.yolo.data <- function() {
    cat("Reading YOLO bbox\n")
    data <- border.data %>%
      dplyr::group_by(Image) %>%
      dplyr::mutate(stomata = paste0("s", sprintf("%02d", row_number()))) %>%
      dplyr::ungroup() %>%
      dplyr::select(stomata, x.com = x, y.com = y, w, h) %>%
      dplyr::mutate(
        xmin = x.com - w / 2, xmax = x.com + w / 2,
        ymin = y.com - h / 2, ymax = y.com + h / 2
      )
    as.data.frame

    # Create polygons from border
    polys <- data %>%
      dplyr::group_by(stomata) %>%
      dplyr::rowwise() %>%
      dplyr::summarise(poly.matrix = list(matrix(
        c(
          xmin, xmin, xmax, xmin, xmin,
          ymin, ymax, ymax, ymax, ymin
        ),
        ncol = 2, byrow = F
      )), .groups = "drop")
    data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))
    data$feret <- lapply(data$polygons, \(x) calculate.stomata.orientation(sf::st_coordinates(x)))
    data <- data %>%
      tidyr::unnest_wider(feret)

    data %>% dplyr::select(-w, -h, -(xmin:ymax))
  }

  # OBB outputs have one row per point in the OBB (i.e. 4 per object) with the
  # following columns:
  # side_a	side_b	obb_width	obb_height	aspect_ratio	image_name	Image	Object
  # max.feret	circularity	x_center	y_center	w	h	corner	x	y
  read.obb.yolo.data <- function() {
    cat("Reading YOLO OBB\n")

    data <- border.data %>%
      dplyr::select(Image, Object, x, y, corner) |>
      dplyr::group_by(Image, Object) |>
      dplyr::mutate(
        stomata = paste0("s", sprintf("%03d", cur_group_id())),
        xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y)
      ) |>
      tidyr::pivot_wider(names_from = corner, values_from = c(x, y)) |>
      dplyr::rowwise() |>
      dplyr::group_by(Image, stomata) |>
      dplyr::summarise(

        # Define the outline of the rectangle
        poly.matrix = list(matrix(
          c(
            x_1, x_2, x_3, x_4, x_1,
            y_1, y_2, y_3, y_4, y_1
          ),
          ncol = 2, byrow = F
        )),
        # Define a rectangle aligned to the vertical/horizontal axis of the image
        # that contains the OBB
        expanded.bbox = list(matrix(
          c(
            xmin, xmax, xmax, xmin, xmin,
            ymin, ymin, ymax, ymax, ymin
          ),
          ncol = 2, byrow = F
        )),
        .groups = "drop"
      )

    data$polygons <- lapply(data$expanded.bbox, function(x) sf::st_polygon(list(x)))

    # We don't want the max Feret diameter from the expanded bboxes; the diagonal
    # is entirely different. Instead, we want the marker points to be placed along
    # the longest axis of the expanded bbox.
    data$feret <- lapply(data$polygons, \(x) calculate.bounding.box.orientation.points(sf::st_coordinates(x)))
    data <- data |>
      tidyr::unnest_wider(feret) |>
      dplyr::rowwise() |>
      dplyr::mutate(
        area = sf::st_area(polygons),
        perimeter = sf::st_perimeter(polygons),
        circularity = 4 * pi * area / perimeter^2
      )

    data
  }

  # Segmentation outputs have 4 columns: Image, Object, x, y
  # Each point in an outline is a separate row
  read.segmentation.yolo.data <- function() {
    cat("Reading YOLO segment mask\n")
    data <- border.data %>%
      dplyr::group_by(Image, Object) %>%
      dplyr::mutate(stomata = paste0("s", sprintf("%03d", Object))) %>%
      dplyr::ungroup() %>%
      dplyr::select(Image, stomata) %>%
      dplyr::distinct()

    # Create polygons from border
    cat("Creating polygons from segment mask\n")
    polys <- border.data %>%
      dplyr::group_by(Image, Object) %>%
      dplyr::reframe(poly.matrix = list(matrix(c(x, x[1], y, y[1]), ncol = 2, byrow = F)))
    data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))

    # Find the max diameter points from the polygon and angles to vertical/horizontal
    cat("Finding stomata orientations\n")
    data$feret <- lapply(data$polygons, \(x) calculate.stomata.orientation(sf::st_coordinates(x)))
    data <- data %>%
      tidyr::unnest_wider(feret) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        area = sf::st_area(polygons),
        perimeter = sf::st_perimeter(polygons),
        circularity = 4 * pi * area / perimeter^2
      )

    # Remove objects that are too small or too irregular to be stomata
    min.area <- median(data$area) / 2
    max.area <- median(data$area) * 2
    min.circ <- min(0.5, median(data$circularity) / 2)
    data <- data %>%
      dplyr::filter(between(area, min.area, max.area) & circularity > min.circ) %>%
      dplyr::mutate(
        Folder = dirname(Image),
        File = basename(Image)
      )

    data
  }


  # Is this a bbox output, OBB output, or a segmentation output?
  if (ncol(border.data) == 4) {
    return(read.segmentation.yolo.data())
  }

  if (ncol(border.data) == 6) {
    return(read.bbox.yolo.data())
  }

  return(read.obb.yolo.data())
}


# Calculate the g function for an image
# image.data - a data.frame with columns for x and y coordinates of the CoM of objects to measure
# image the image file
g.function <- function(com.data, image.file.name) {
  # Calculate neighbour distances
  d <- dist(com.data)
  dm <- as.matrix(d)
  diag(dm) <- NA
  dmin <- apply(dm, 1, min, na.rm = TRUE)

  # get the unique distances (for the x-axis)
  distance <- sort(unique(round(dmin)))
  # compute how many cases there with distances smaller that each x
  Gd <- sapply(distance, function(x) sum(dmin < x))
  # normalize to get values between 0 and 1
  Gd <- Gd / length(dmin)

  data.frame(distance, Gd, Image = image.file.name)
}
