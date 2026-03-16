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
library(hexbin)

# Save the given plot to the given path
save.ggplot <- function(plot, out.file, width = 170, height = 85) {
  ggsave(out.file, plot = plot, dpi = 300, units = "mm", width = width, height = height)
}

# Save the given plot to the given path and return the plot
# Size for powerpoint slide: width= 253.8, height = 190.05
#
# plot - the ggplot to save
# image.file - the name of the image file to be saved
# suffix - a file name suffix
# width, height - dimensions in mm
save.plot <- function(plot, image.file, suffix, width = 253.8, height = 190.05) { #
  out.path <- str_replace(basename(image.file), ".jpg", suffix)
  out.path <- str_replace(out.path, "/", "_")
  out.path <- paste0("figure/", out.path)
  save.ggplot(out.file = out.path, plot = plot, width = width, height = height)
  plot
}


# Given two points in a list, return the leftmost (lowest x)
#
# p1 - a list with named elements X and Y
# p2 - a list with named elements X and Y
#
# Returns: the point with the lower X value
left <- function(p1, p2) {
  if (p1[["X"]] < p2[["X"]]) {
    return(p1)
  } else {
    return(p2)
  }
}
# Given two points in a list, return the rightmost (highest x)
#
# p1 - a list with named elements X and Y
# p2 - a list with named elements X and Y
#
# Returns: the point with the higher X value
right <- function(p1, p2) {
  if (p1[["X"]] >= p2[["X"]]) {
    return(p1)
  } else {
    return(p2)
  }
}
# Given two points in a list, return the upper (highest y)
#
# p1 - a list with named elements X and Y
# p2 - a list with named elements X and Y
#
# Returns: the point with the higher Y value
upper <- function(p1, p2) {
  if (p1[["Y"]] > p2[["Y"]]) {
    return(p1)
  } else {
    return(p2)
  }
}
# Given two in a list, return the lower (lowest y)
#
# p1 - a list with named elements X and Y
# p2 - a list with named elements X and Y
#
# Returns: the point with the lower Y value
lower <- function(p1, p2) {
  if (p1[["Y"]] <= p2[["Y"]]) {
    return(p1)
  } else {
    return(p2)
  }
}

# convert degrees to radians
#
# deg - the angle in degrees
#
# Returns: the angle in radians
deg2rad <- function(deg) (deg * pi) / (180)

# convert radians to degrees
#
# rad - the angle in radians
#
# Returns: the angle in degrees
rad2deg <- function(rad) (180 * rad) / pi

# Calculate the angle between two lines. Each line is defined by two points.
#
# l1.p1 - a point on line 1; a list with 'X' and 'Y' elements
# l1.p2- a point on line 1; a list with 'X' and 'Y' elements
# l2.p1 - a point on line 2; a list with 'X' and 'Y' elements
# l2.p2 - a point on line 2; a list with 'X' and 'Y' elements
#
# Returns: the absolute angle between the lines in degrees
angle.between.lines <- function(l1.p1, l1.p2, l2.p1, l2.p2) {
  a <- l1.p2[["X"]] - l1.p1[["X"]] # l1 x diff
  b <- l1.p2[["Y"]] - l1.p1[["Y"]] # l1 y diff

  c <- l2.p2[["X"]] - l2.p1[["X"]] # l2 x difference
  d <- l2.p2[["Y"]] - l2.p1[["Y"]] # l2 y difference

  atanA <- atan2(a, b)
  atanB <- atan2(c, d)

  return(rad2deg(atanA - atanB))
}


# Calculate the angle of line described by the given points to vertical
# Returns angle in degrees anticlockwise. The line will start from the upper point
angle.to.vertical <- function(x1, y1, x2, y2) {
  # L1 is the vertical axis
  a <- 0 # l1 x diff - same x coords from vertical
  b <- 1 # l1 y diff -arbitrary y offset

  # L2 is the line we input
  c <- x2 - x1 # l2 x difference
  d <- y2 - y1 # l2 y difference

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

angle.is.within.range <- function(angle, reference.angle, max.angle.delta) {
  if (is.null(reference.angle)) stop("Reference angle is missing")
  min.angle <- (reference.angle - max.angle.delta) %% 180
  max.angle <- (reference.angle + max.angle.delta) %% 180

  # cat(
  #   "Angle filter params:\nRef", reference.angle, "Max delta", max.angle.delta, "\n",
  #   "Min:", min.angle, "Max:", max.angle, "\n",
  #   "Input angles:", angle, "\n"
  # )

  # If we are close to 0 or 180, these will wrap; we must invert the range if so
  if (min.angle > max.angle) {
    # Handle wrapping of angles around 180
    return(angle >= min.angle | angle <= max.angle)
  }
  # otherwise normal filter
  return(angle >= min.angle & angle <= max.angle)
}

# Given a matrix containing X and Y columns constituting an bounding box border,
# identify the longest axis and return points that lie run on this axis through
# the centre of mass
calculate.bounding.box.orientation.points <- function(points) {
  # Get the CoM of the object.
  x.com <- mean(points[, "X"])
  y.com <- mean(points[, "Y"])

  #   axis12
  #  1 ----- 2
  #  |       |  axis14
  #  |       |
  #  4-------3

  p1 <- points[1, ]
  p2 <- points[2, ]
  p3 <- points[3, ]
  p4 <- points[4, ]

  # Create axes for the OBB border and the diagonals
  axis.12 <- sf::st_linestring(matrix(
    data = c(p1["X"], p1["Y"], p2["X"], p2["Y"]),
    byrow = TRUE, nrow = 2, ncol = 2
  ))
  axis.14 <- sf::st_linestring(matrix(
    data = c(p1["X"], p1["Y"], p4["X"], p4["Y"]),
    byrow = TRUE, nrow = 2, ncol = 2
  ))

  diagonal.13 <- sf::st_linestring(matrix(
    data = c(p1["X"], p1["Y"], p3["X"], p3["Y"]),
    byrow = TRUE, nrow = 2, ncol = 2
  ))
  diagonal.24 <- sf::st_linestring(matrix(
    data = c(p2["X"], p2["Y"], p4["X"], p4["Y"]),
    byrow = TRUE, nrow = 2, ncol = 2
  ))

  # What are the lengths of the axes in the rectangle?
  axis.12.length <- sf::st_length(axis.12)
  axis.14.length <- sf::st_length(axis.14)

  axis.13.length <- sf::st_length(diagonal.13)
  axis.24.length <- sf::st_length(diagonal.24)

  # The longest axis is mostly likely the general orientation of the stomata
  if (axis.12.length > axis.14.length) {
    long.axis <- axis.12
    short.axis <- axis.14
  } else {
    long.axis <- axis.14
    short.axis <- axis.12
  }

  # Calculate the angle between the longest axis and the two possible
  # diagonals. We assume that one of the diagonals is the desired orientation.

  # Which diagonal has the lowest angle to the long OBB axis?
  angle.to.diagonal.1 <- angle.between.lines(
    c("X" = long.axis[1], "Y" = long.axis[3]),
    c("X" = long.axis[2], "Y" = long.axis[4]),
    c("X" = diagonal.13[1], "Y" = diagonal.13[3]),
    c("X" = diagonal.13[2], "Y" = diagonal.13[4])
  )

  angle.to.diagonal.2 <- angle.between.lines(
    c("X" = long.axis[1], "Y" = long.axis[3]),
    c("X" = long.axis[2], "Y" = long.axis[4]),
    c("X" = diagonal.24[1], "Y" = diagonal.24[3]),
    c("X" = diagonal.24[2], "Y" = diagonal.24[4])
  )

  # Select the diagonal axis closest to the long axis
  if (abs(angle.to.diagonal.1) > abs(angle.to.diagonal.2)) {
    pA <- p1
    pB <- p3
  } else {
    pA <- p2
    pB <- p4
  }

  # # Select the points lying on the longer axis
  # pA <- p1
  #
  # if (axis.12.length > axis.14.length) {
  #   pB <- p2
  # } else {
  #   pB <- p4
  # }


  list(
    max.feret = sf::st_length(long.axis),
    p1 = list(X = pA["X"], "Y" = pA["Y"]),
    p2 = list(X = pB["X"], "Y" = pB["Y"]),
    dx = abs(pB["X"] - pA["X"]),
    dy = abs(pB["Y"] - pA["Y"]),
    d.long = sf::st_length(long.axis),
    d.short = sf::st_length(short.axis),
    angle.to.diagonal.1 = angle.to.diagonal.1,
    angle.to.diagonal.2 = angle.to.diagonal.2,
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


# Extend the given sf_linestring by the given amount
#
# line - the sf_linestring
# distance - the amount to extend
# dir - the direction to extend
extend.linestring <- function(line, distance, dir = "both") {
  if (!dir %in% c("head", "tail", "both")) {
    stop("Invalid input 'dir' must either be 'head', 'tail', or 'both'")
  }

  # convert to WK coords
  coords <- line
  # coords <- as.matrix(coords[c("x", "y")])

  # which index to keep
  to_keep <- dir != c("tail", "head")

  # dir coords index we want to keep
  dirs <- c(1, nrow(coords))[to_keep]

  # DETERMINE THE DIRECTION OF EACH COORDINATE
  x_coords <- unname(coords[, 1])
  y_coords <- unname(coords[, 2])

  # X/Y coordinate pairs
  x1 <- x_coords[1]
  y1 <- y_coords[1]
  x2 <- x_coords[2]
  y2 <- y_coords[2]

  directions <- c(atan2(y1 - y2, x1 - x2), atan2(y2 - y1, x2 - x1))

  # get directions of the direction of interest
  directions <- directions[to_keep]

  # if only a single distance, duplicate it, otherwise reverse the first 2 distances
  distances <- if (length(distance) == 1) {
    rep(distance, 2)
  } else {
    rev(distance[1:2])
  }

  # adjust dir point coordinates
  coords[dirs, ] <- coords[dirs, ] + distances[to_keep] * c(cos(directions), sin(directions))

  # # make a new linestring
  line <- sf::st_linestring(
    matrix(c(coords[, 1], coords[, 2]), byrow = FALSE, nrow = 2)
  )

  return(line)
}

# Find the maximum value in the density plot of the given vector
find.mode <- function(x) {
  # We are dealing with angles that may wrap around 180-0
  # Duplicate the counts so we get a proper peak at 180
  y <- c(x, x + 180)

  # Get histogram, find peak with max
  h <- hist(y, breaks = 72, plot = FALSE)
  h$breaks[which.max(h$counts)] %% 180
  #
  # This approach will not work if there is a very sharp peak in the histogram
  # beacuse a wide low peak may have overall higher density
  #   d <- density(y)
  #   d$x[which.max(d$y)] %% 180
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
        .groups = "drop"
      )

    # Create a polygon object with the OBB
    data$obbs <- lapply(data$poly.matrix, function(x) sf::st_polygon(list(x)))

    # Calculate the largest circle fitting inside the bounding box
    # This will help prevent spurious intersections between stomata being detected
    data$polygons <- lapply(data$obbs, st_inscribed_circle, dTolerance = 0.0001)

    # Calculate the Feret diameter of the polygons and other useful measures.
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
# com.data - a data.frame with columns for x and y coordinates of the CoM of objects to measure
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
