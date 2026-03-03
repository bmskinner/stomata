# Create chains from point coordinates

#### Imports #####
source("functions.R")

#### Constants #####
# min and maximum distance between stomata in pixels
MIN.DISTANCE <- 50
MAX.DISTANCE <- 500

# 2mer angle variation from the modal value in degrees
ANGLE.DELTA.2MER <- 15

# How far can a 3mer deviate from a 180 line in degrees?
ANGLE.DELTA.INTERNAL.3MER <- 10
# How much can the absolute angle of a 3mer deviate from the modal value?
ANGLE.DELTA.MODAL.3MER <- 10

# If we have remaining 2mers that might be part of a contig, how stringent an angle
# must they have to the modal 2mer angle in degrees?
ANGLE.DELTA.2MER.STRINGENT <- 5

#### Plotting Functions #####

plot.1mers <- function(mer.data, include.stomata = TRUE) {
  img.grob <- rasterGrob(mer.data$img, interpolate = TRUE)

  # Don't show message about new coord system by creating it, then setting it as
  # default before use
  cf <- coord_fixed(xlim = c(0, dim(mer.data$img)[2]), ylim = c(0, dim(mer.data$img)[1]))
  cf$default <- TRUE


  stomata <- data.frame(
    "x1" = sapply(mer.data$mer1$p1, \(p) p$X),
    "y1" = sapply(mer.data$mer1$p1, \(p) p$Y),
    "x2" = sapply(mer.data$mer1$p2, \(p) p$X),
    "y2" = sapply(mer.data$mer1$p2, \(p) p$Y)
  )


  p <- ggplot() +
    # Draw the input image as a background
    annotation_custom(img.grob,
      xmin = 0, xmax = dim(mer.data$img)[2],
      ymin = 0, ymax = dim(mer.data$img)[1]
    ) +
    cf +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Write the angles of the 1mer
    # geom_text(
    #   data = mer.data$mer1, aes(x = x.com, y = y.com + 20, label = sprintf("%.0f", abs.angle)),
    #   col = "red", size = 1.5
    # ) +
    # geom_text(data = mer.data$mer1, aes(x = x.com, y = y.com-20, label = sprintf("%.0f", horzAngle)),
    #           col = "red", size = 1.5)+

    # Draw the name of each stomata
    geom_text(data = mer.data$mer1, aes(x = x.com, y = y.com - 10, label = stomata), col = "pink1", size = 2) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = -20, r = -30, b = -20, l = -30) # remove whitespace outside image
    )

  if (include.stomata) {
    p <- p + # Draw the stomata points
      geom_point(data = stomata, aes(x = x1, y = y1), col = "blue", size = 2) +
      geom_point(data = stomata, aes(x = x2, y = y2), col = "orange", size = 2)
  }
  p
}


plot.2mers <- function(mer.data) {
  # Create image raster
  img.grob <- rasterGrob(mer.data$img, interpolate = TRUE)

  # Don't show message about new coord system by creating it, then setting it as
  # default before use
  cf <- coord_fixed(xlim = c(0, dim(mer.data$img)[2]), ylim = c(0, dim(mer.data$img)[1]))
  cf$default <- TRUE

  ggplot() +
    # Draw the input image as a background
    annotation_custom(img.grob,
      xmin = 0, xmax = dim(mer.data$img)[2],
      ymin = 0, ymax = dim(mer.data$img)[1]
    ) +
    cf +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Draw the 2mer
    geom_segment(
      data = mer.data$mer2, aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y),
      col = "blue", linewidth = 0.5
    ) +

    # Write the angle of the 2mer
    geom_text(
      data = mer.data$mer2, aes(
        x = (S2.x + S1.x) / 2, y = (S2.y + S1.y) / 2,
        label = paste(sprintf("%.1f°\n%.0fpx", angle.of.2mer, length))
      ),
      col = "red", size = 3
    ) +

    # Draw the centroid of each stomata in the 2mer
    geom_point(data = mer.data$mer2, aes(x = S1.x, y = S1.y), col = "blue", size = 2) +
    geom_point(data = mer.data$mer2, aes(x = S2.x, y = S2.y), col = "orange", size = 1) +

    # Draw the name of each stomata in the 2mer
    geom_text(data = mer.data$mer2, aes(x = S1.x, y = S1.y - 10, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer2, aes(x = S2.x, y = S2.y - 10, label = S2), col = "pink1", size = 2) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = -20, r = -30, b = -20, l = -30) # remove whitespace outside image
    )
}

plot.3mers <- function(mer.data) {
  img.grob <- rasterGrob(mer.data$img, interpolate = TRUE)

  cf <- coord_fixed(xlim = c(0, dim(mer.data$img)[2]), ylim = c(0, dim(mer.data$img)[1]))
  cf$default <- TRUE


  ggplot() +
    # Draw the input image as a background
    annotation_custom(img.grob,
      xmin = 0, xmax = dim(mer.data$img)[2],
      ymin = 0, ymax = dim(mer.data$img)[1]
    ) +
    cf +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Draw the centroid of each stomata
    geom_point(data = mer.data$mer1, aes(x = x.com, y = y.com), col = "black") +

    # Draw the 3mer
    geom_segment(
      data = mer.data$mer3, aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y),
      col = "blue", linewidth = 1.5, alpha = 0.8
    ) +
    geom_segment(
      data = mer.data$mer3, aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y),
      col = "orange", linewidth = 0.5
    ) +
    geom_text(
      data = mer.data$mer3, aes(
        x = (S1.x + S2.x) / 2, y = (S2.y + S1.y) / 2,
        label = sprintf("%.0f°\n|%.0f°|", mer3.internal.angle, mer3.angle.to.vertical)
      ),
      col = "red", size = 3
    ) +
    geom_text(data = mer.data$mer3, aes(x = S1.x, y = S1.y - 10, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer3, aes(x = S2.x, y = S2.y - 10, label = S2), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer3, aes(x = S3.x, y = S3.y - 10, label = S3), col = "pink1", size = 2) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = -20, r = -30, b = -20, l = -30) # remove whitespace outside image
    )
}

# plot contig data
plot.contigs <- function(mer.data) {
  img.grob <- rasterGrob(mer.data$img, interpolate = TRUE)

  cf <- coord_fixed(xlim = c(0, dim(mer.data$img)[2]), ylim = c(0, dim(mer.data$img)[1]))
  cf$default <- TRUE

  # visualise the contigs
  ggplot() +
    # Draw the input image as a background
    annotation_custom(img.grob,
      xmin = 0, xmax = dim(mer.data$img)[2],
      ymin = 0, ymax = dim(mer.data$img)[1]
    ) +
    cf +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Draw the 2mers in the chain
    geom_point(data = mer.data$contigs, aes(x = S1.x, y = S1.y, col = as.factor(Contig)), size = 2) +
    geom_point(data = mer.data$contigs, aes(x = S2.x, y = S2.y, col = as.factor(Contig)), size = 2) +
    geom_segment(data = mer.data$contigs, aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y, col = as.factor(Contig)), linewidth = 1) +
    labs(col = "Chain") +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = -20, r = -30, b = -20, l = -30) # remove whitespace outside image
    )
}

plot.rotated.contigs <- function(mer.data) {
  ggplot(mer.data$measured.contigs) +
    geom_segment(
      aes(
        x = Contig.start.x,
        y = Contig.start.y,
        xend = RotatedContigEnd[, 1],
        yend = RotatedContigEnd[, 2]
      ),
      col = "grey20"
    ) +
    # geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y), col = "grey")+
    geom_segment(aes(x = RotatedS1[, 1], y = RotatedS1[, 2], xend = RotatedS2[, 1], yend = RotatedS2[, 2]), col = "blue") +
    geom_point(aes(x = RotatedS1[, 1], y = RotatedS1[, 2]), col = "blue", size = 2) +
    geom_point(aes(x = RotatedS2[, 1], y = RotatedS2[, 2]), col = "blue", size = 2) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    )
}

#### Analysis Functions #####

# Give a file of stomata locations, create 1mers.
# the input file can be in .json (AnyLabelling) format
# or txt (YOLO bounding box format)
read.1mers <- function(file) {
  read.1mer.from.json <- function(file) {
    # Create 1mers from stomata outlines
    create.1mers <- function(border.data) {
      data <- border.data %>%
        dplyr::group_by(shape, file) %>%
        dplyr::summarise(
          x.com = mean(x),
          y.com = mean(y),
          .groups = "drop"
        ) %>%
        dplyr::mutate(stomata = paste0("s", sprintf("%02d", shape))) %>%
        dplyr::ungroup() %>%
        dplyr::select(stomata, x.com, y.com, file) %>%
        as.data.frame()

      # Create polygons from border
      polys <- border.data %>%
        dplyr::group_by(shape, file) %>%
        dplyr::summarise(poly.matrix = list(matrix(c(x, x[1], y, y[1]), ncol = 2, byrow = F)), .groups = "drop")
      data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))
      data %>% dplyr::select(-file)
    }

    border.data <- read.border.from.json(file)

    if (nrow(border.data) <= 1) {
      # return empty dist.contigs dataframe
      cat("  No objects in json file, returning\n")
      return(data.frame())
    }

    create.1mers(border.data)
  }

  if (str_ends(file, "json")) {
    mer1 <- read.1mer.from.json(file)
  } else {
    mer1 <- read.border.from.yolo(file)
  }

  mer1
}

# Check for overlapping objects and keep only the object with highest
# area
#
# mer.data - the complete data
# min.distance - only check for overlaps between stomata closer than this (pixels)
filter.1mers.for.overlaps <- function(mer.data, min.distance) {
  if (mer.data$write.debug.images) {
    save.plot(
      plot.1mers(mer.data, include.stomata = FALSE),
      mer.data$image.file, ".mer1.1.raw.png"
    )
  }

  # Pairwise comparison of stomata. Find any with overlapping bounds
  # Keep the larger
  distances <- expand.grid(
    s1 = mer.data$mer1$stomata, s2 = mer.data$mer1$stomata,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(s1 != s2) |>
    merge(mer.data$mer1, by.x = "s1", by.y = "stomata") |>
    merge(mer.data$mer1, by.x = "s2", by.y = "stomata") |>
    dplyr::rowwise() |>
    dplyr::mutate(
      distance.between.stomata = euclidean(x.com.x, y.com.x, x.com.y, y.com.y)
    ) |>
    dplyr::filter(distance.between.stomata < min.distance)

  if (nrow(distances) > 1) {
    distances <- distances |>
      dplyr::mutate(
        has.overlap = sf::st_overlaps(polygons.x, polygons.y),
        toRemove = ifelse(area.x < area.y, s1, s2)
      ) |>
      dplyr::filter(has.overlap == 1) |>
      dplyr::arrange(toRemove)

    stomata.to.remove <- unique(distances$toRemove)

    cat("Removing", length(stomata.to.remove), "1mers due to overlapping bounding boxes\n")

    filt <- mer.data$mer1[!(mer.data$mer1$stomata %in% stomata.to.remove), ]

    mer.data$distances <- distances
    mer.data$mer1 <- filt
  } else {
    cat("No overlapping bounding boxes detected\n")
  }

  if (mer.data$write.debug.images) {
    save.plot(
      plot.1mers(mer.data, include.stomata = FALSE),
      mer.data$image.file, ".mer1.2.overlaps.png"
    )
  }
  mer.data
}

# From given 1-mers, create 2mers. filter to those within a given distance of
# each other and annotate with absolute angles on image
#
# mer.data - the complete data object
# min.distance, max.distance - how close can stomata be to each other
create.2mers <- function(mer.data, min.distance, max.distance) {
  # Create all pairwise combinations of 1mers
  result <- expand.grid(mer.data$mer1$stomata, mer.data$mer1$stomata, stringsAsFactors = F) |>
    dplyr::distinct() |> # remove duplicates
    dplyr::filter(Var1 != Var2) |>
    merge(mer.data$mer1[, c("stomata", "x.com", "y.com")], by.x = "Var1", by.y = "stomata", all.y = F) |> # add coordinates for S1
    dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>% # rename for clarity
    merge(mer.data$mer1[, c("stomata", "x.com", "y.com")], by.x = "S2", by.y = "stomata", all.y = F) |> # add coordinates for S2
    dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) |> # rename for clarity

    # Apply a length filter to exclude 2mers that are too distant
    # The plurality of the remaining edges will be orientated along the chains
    dplyr::mutate(length = euclidean(S1.x, S1.y, S2.x, S2.y)) |>
    dplyr::filter(between(length, min.distance, max.distance)) |>
    # Calculate the angle of a 2mer based on the locations of the points at either end
    # Wrap at 180 degrees
    dplyr::mutate(angle.of.2mer = angle.to.vertical(S1.x, S1.y, S2.x, S2.y) %% 180)


  # What is the distribution of 2mer angles in the data?
  mer.data$modal.2mer.angle <- find.mode(result$angle.of.2mer)

  # Choose the first and second points of the 2mer depending on overall orientation
  # of the 2mers. If generally horizontal, go left to right. If generally vertical,
  # go top to bottom
  if (mer.data$modal.2mer.angle > 45 & mer.data$modal.2mer.angle < 135) {
    # the first stomata is always on the left in the image

    result <- result |>
      dplyr::mutate(
        s1.is.first = S1.x < S2.x,
        first = ifelse(s1.is.first, S1, S2),
        first.x = ifelse(s1.is.first, S1.x, S2.x),
        first.y = ifelse(s1.is.first, S1.y, S2.y),
        second = ifelse(s1.is.first, S2, S1),
        second.x = ifelse(s1.is.first, S2.x, S1.x),
        second.y = ifelse(s1.is.first, S2.y, S1.y),
        S1 = first,
        S1.x = first.x,
        S1.y = first.y,
        S2 = second,
        S2.x = second.x,
        S2.y = second.y
      )
  } else {
    # the first stomata is always at the top in the image

    result <- result |>
      dplyr::mutate(
        s1.is.first = S1.y < S2.y,
        first = ifelse(s1.is.first, S1, S2),
        first.x = ifelse(s1.is.first, S1.x, S2.x),
        first.y = ifelse(s1.is.first, S1.y, S2.y),
        second = ifelse(s1.is.first, S2, S1),
        second.x = ifelse(s1.is.first, S2.x, S1.x),
        second.y = ifelse(s1.is.first, S2.y, S1.y),
        S1 = first,
        S1.x = first.x,
        S1.y = first.y,
        S2 = second,
        S2.x = second.x,
        S2.y = second.y
      )
  }

  result <- result |>
    dplyr::select(S1, S1.x, S1.y, S2, S2.x, S2.y, length, angle.of.2mer) |>
    dplyr::mutate(kmer = paste0(S1, S2))


  # Convert stomata CoMs to list of coordinates
  # result$S1.point <- mapply(\(x, y, n) list(X = x, Y = y, name = n), result$S1.x, result$S1.y, result$S1, SIMPLIFY = FALSE)
  # result$S2.point <- mapply(\(x, y, n) list(X = x, Y = y, name = n), result$S2.x, result$S2.y, result$S2, SIMPLIFY = FALSE)

  # Make a name for the 2mer that can be used as an ID
  # result$kmer <- mapply(\(p1, p2) paste0(p1$name, p2$name), result$S1.point, result$S2.point)

  if (mer.data$write.debug.images) {
    dup.data <- data.frame("angle.of.2mer" = c(result$angle.of.2mer, result$angle.of.2mer + 180))

    histo.plot <- ggplot(dup.data, aes(x = angle.of.2mer)) +
      geom_histogram(aes(y = ..density..), binwidth = 5) +
      geom_density() +
      geom_vline(xintercept = mer.data$modal.2mer.angle) +
      geom_vline(xintercept = mer.data$modal.2mer.angle + 180) +
      labs(title = sprintf("Modal angle %.2f°", mer.data$modal.2mer.angle)) +
      scale_x_continuous(breaks = seq(0, 360, 45)) +
      theme_bw()
    save.plot(histo.plot, mer.data$image.file, ".mer2.histo.png")
  }


  if (nrow(result) > 0) {
    result$mer2id <- 1:nrow(result)
  }
  cat(
    "Created", nrow(result), "2mers from", nrow(mer.data$mer1), "1mers within length bounds",
    min.distance, "-", max.distance, "pixels\n"
  )

  # cat("Created", nrow(result), "2mers from", nrow(mer.data$mer1), "1mers\n")
  mer.data$mer2 <- result

  if (mer.data$write.debug.images) save.plot(plot.2mers(mer.data), mer.data$image.file, ".mer2.1.raw.png")

  return(mer.data)
}


# Filter 2mers to those lying within a given delta to the modal angle of edges
# in the image
#
# mer.data - the complete data
# max.angle.delta - the maximum deviation from the modal angle in degrees
filter.2mers.by.angle <- function(mer.data, max.angle.delta) {
  filt <- mer.data$mer2 |>
    dplyr::filter(angle.is.within.range(angle.of.2mer, mer.data$modal.2mer.angle, max.angle.delta))

  cat(
    sprintf(
      "Filtered from %i to %i 2mers within %.2f° of the modal 2mer angle (%.2f°)\n",
      nrow(mer.data$mer2), nrow(filt), max.angle.delta, mer.data$modal.2mer.angle
    )
  )
  mer.data$mer2 <- filt

  # Update the modal angle
  mer.data$modal.2mer.angle <- find.mode(mer.data$mer2$angle.of.2mer)

  if (mer.data$write.debug.images) save.plot(plot.2mers(mer.data), mer.data$image.file, ".mer2.2.filt.modal.angle.png")
  mer.data
}


# Remove 2mers intersecting a third stomata
#
# mer.data - the complete data
filter.2mers.by.intersections <- function(mer.data) {
  # Remove edges in the graph that intersect a third stomata
  # can use sf: https://stackoverflow.com/questions/61703791/determine-lines-that-intersect-a-polygon-in-r
  # Create line objects
  lines <- lapply(1:nrow(mer.data$mer2), function(i) {
    sf::st_linestring(rbind(
      c(mer.data$mer2$S1.x[i], mer.data$mer2$S1.y[i]),
      c(mer.data$mer2$S2.x[i], mer.data$mer2$S2.y[i])
    ))
  })

  # Check each 2mer for intersections with a stomata
  # If we intersect, then this is not a valid 2mer
  mer.data$mer2$intersects <- lapply(lines, function(l) sum(unlist(sapply(mer.data$mer1$polygons, sf::st_intersects, y = l))))

  filt <- mer.data$mer2 %>% dplyr::filter(intersects <= 2)

  cat("Filtered from", nrow(mer.data$mer2), "to", nrow(filt), "2mers with no internally intersecting stomata\n")

  mer.data$mer2 <- filt

  if (mer.data$write.debug.images) save.plot(plot.2mers(mer.data), mer.data$image.file, ".mer2.3.filt.intersection.png")

  return(mer.data)
}


# From given 2mers, create 3mers. Filter to those within an angle delta of 180
# degrees
#
# mer.data - the complete data
create.3mers <- function(mer.data) {
  # Merge the 2mers
  result <- merge(mer.data$mer2, mer.data$mer2,
    by.x = c("S2", "S2.x", "S2.y"),
    by.y = c("S1", "S1.x", "S1.y"),
    all.y = F, suffixes = c("A", "B")
  ) %>%
    dplyr::select(S1, S1.x, S1.y,
      S2, S2.x, S2.y,
      S3 = S2B, S3.x = S2.xB, S3.y = S2.yB,
      S1S2 = lengthA, S2S3 = lengthB,
      S1S2.angle = angle.of.2merA, S2S3.angle = angle.of.2merB
    ) %>% # rename for clarity

    # Calculate angle of 3mer. Check with stomata are top and bottom
    dplyr::rowwise() %>%
    dplyr::mutate(

      # Internal angle of the 3mer (closeness to straight line)
      mer3.internal.angle = LearnGeom::Angle(c(S1.x, S1.y), c(S2.x, S2.y), c(S3.x, S3.y)),

      # From create.2mers we know that S1-S2-S3 is already ordered L-R or T-B
      mer3.angle.to.vertical = angle.to.vertical(S1.x, S1.y, S3.x, S3.y) %% 180
    ) %>%
    dplyr::mutate(
      merFirst = paste0(S1, S2),
      merLast = paste0(S2, S3),
      mer3Id = paste0(S1, S2, S3)
    ) %>%
    # Remove the unused columns
    dplyr::distinct()

  cat("Created", nrow(result), "3mers from", nrow(mer.data$mer2), "2mers\n")

  mer.data$mer3 <- result

  if (mer.data$write.debug.images) save.plot(plot.3mers(mer.data), mer.data$image.file, ".mer3.1.raw.png")

  mer.data
}

# Filter 3mers to those near 180 degrees
#
# mer.data - the complete data
filter.3mers.by.straightness <- function(mer.data, max.angle.delta) {
  if (mer.data$write.debug.images) {
    histo.plot <- ggplot(mer.data$mer3, aes(x = mer3.internal.angle)) +
      geom_histogram(aes(y = ..density..), binwidth = 5) +
      geom_density() +
      scale_x_continuous(breaks = seq(0, 360, 5)) +
      theme_bw()
    save.plot(histo.plot, mer.data$image.file, ".mer3.angle.histo.png")
  }

  filt <- mer.data$mer3 %>%
    # First keep only the 3mers in straight lines
    dplyr::filter(mer3.internal.angle > 180 - max.angle.delta)

  mer.data$modal.mer3.angle <- find.mode(filt$mer3.angle.to.vertical)

  filt <- filt %>%
    # Look for 2mers present more than once
    # Drop the 3mers with the lowest angle.
    dplyr::group_by(merFirst) %>%
    dplyr::arrange(merFirst, desc(mer3.internal.angle)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::group_by(merLast) %>%
    dplyr::arrange(merLast, desc(mer3.internal.angle)) %>%
    dplyr::slice_head(n = 1)

  cat(
    sprintf(
      "Filtered from %i to %i 3mers within %.2f° of a straight line\n",
      nrow(mer.data$mer3), nrow(filt), max.angle.delta
    )
  )

  mer.data$mer3 <- filt
  if (mer.data$write.debug.images) save.plot(plot.3mers(mer.data), mer.data$image.file, ".mer3.2.straight.png")
  mer.data
}

# Filter 3mers to those near 180 degrees
#
# mer.data - the complete data
filter.3mers.by.orientation <- function(mer.data, max.angle.delta) {
  if (mer.data$write.debug.images) {
    histo.plot <- ggplot(mer.data$mer3, aes(x = mer3.angle.to.vertical)) +
      geom_histogram(aes(y = ..density..), binwidth = 5) +
      geom_density() +
      geom_vline(xintercept = mer.data$modal.mer3.angle) +
      labs(title = sprintf("Modal angle %.2f°", mer.data$modal.mer3.angle)) +
      scale_x_continuous(breaks = seq(0, 360, 5)) +
      theme_bw()
    save.plot(histo.plot, mer.data$image.file, ".mer3.angle.vertical.histo.png")
  }

  filt <- mer.data$mer3 |>
    # Keep only those 3mers within tolerance for modal absolute angle
    dplyr::filter(angle.is.within.range(
      mer3.angle.to.vertical,
      mer.data$modal.mer3.angle,
      max.angle.delta
    ))

  cat(
    sprintf(
      "Filtered from %s to %s 3mers within %.2f° of the modal 3mer angle (%.2f°)\n",
      nrow(mer.data$mer3),
      nrow(filt),
      max.angle.delta,
      mer.data$modal.mer3.angle
    )
  )

  mer.data$mer3 <- filt
  if (mer.data$write.debug.images) save.plot(plot.3mers(mer.data), mer.data$image.file, ".mer3.3.orientation.png")
  mer.data
}


# Remove any 3mers where the endpoint is in the middle of a chain
# to prevent branches
#
# mer.data - the complete data
filter.3mers.by.branch <- function(mer.data) {
  # Identify 3mers in which a K1S3 is also K2S2, but K2S1 is not K1S2.
  # In this case, K1 is terminating in the middle of another chain. Remove K1.
  # mismatches <- mer.data$mer3 |>
  #   merge(mer.data$mer3, by.x = "S3", by.y = "S2", all.y = F, suffixes = c(".K1", ".K2")) |>
  #   dplyr::filter(S1.K2 != S2)
  #
  # print(head(mismatches))
  #
  # filt <- mer.data$mer3 |>
  #   dplyr::filter(!(mer3Id %in% mismatches$mer3Id.K1))

  filt <- mer.data$mer3 %>%
    tidyr::pivot_longer(c(S1, S2, S3), names_to = "StomataPosition", values_to = "Stomata") %>%
    dplyr::group_by(StomataPosition, Stomata) %>%
    # Calculate how far the 3mer orientation differs from the modal orientation
    # We assume that the best 3mer to keep will be closest to the modal value
    dplyr::mutate(abs.mer3.angle.diff = abs(mer.data$modal.mer3.angle - mer3.angle.to.vertical)) %>%
    dplyr::arrange(StomataPosition, Stomata, abs.mer3.angle.diff) %>%
    dplyr::slice_head(n = 1) %>% # take only the 3mer closest to modal angle
    tidyr::pivot_wider(names_from = StomataPosition, values_from = Stomata) %>%
    na.omit() # remove the rows with NAs due to our slice

  cat("Filtered from", nrow(mer.data$mer3), "to", nrow(filt), "3mers removing branches\n")
  mer.data$mer3 <- filt

  if (mer.data$write.debug.images) save.plot(plot.3mers(mer.data), mer.data$image.file, ".mer3.4.prune.png")

  return(mer.data)
}
# There may be 3mers in which a stomata is a terminus of one 3mer and the centre
# of another 3mer. We need to remove the 3mers for which a shared stomata is
# terminal.
# ___                ___
#  \           ->
# ______             _____
#
# mer.data - the complete data
filter.3mers.by.terminal <- function(mer.data) {
  # Look for 1mers that are are the S2 of 2 different 2mers

  filt <- mer.data$mer3 %>%
    tidyr::pivot_longer(c(merFirst, merLast), names_to = "mer2Type", values_to = "mer2") %>%
    dplyr::mutate(
      StomataFirst = ifelse(mer2Type == "merFirst", S1, S2),
      StomataLast = ifelse(mer2Type == "merFirst", S2, S3)
    ) %>%
    dplyr::group_by(StomataFirst) %>% # Take all the starting stomata of the 3mer
    dplyr::mutate(
      CharactersPerStomata = nchar(StomataLast), # how many letters in a stomata id
      AllConnectedStomata = paste(unique(StomataLast), collapse = ""), # make a string of all stomata in a 2mer with the last stomata
      NumberOfConnectedStomata = nchar(AllConnectedStomata) / CharactersPerStomata
    ) %>% # how many stomata join this
    dplyr::group_by(AllConnectedStomata) %>%
    dplyr::arrange(AllConnectedStomata, desc(mer3.internal.angle)) %>% # order so the one to keep is top
    dplyr::mutate(rownum = row_number()) %>% # scruffy; figure out which row in group
    dplyr::filter(NumberOfConnectedStomata == 1 | rownum == 1) %>% # keep the first from multi2mers, or every row from unique 2mers
    dplyr::ungroup() %>%
    dplyr::select(-c(StomataFirst, StomataLast, AllConnectedStomata, CharactersPerStomata, NumberOfConnectedStomata, rownum)) %>%
    tidyr::pivot_wider(names_from = "mer2Type", values_from = "mer2") %>%
    na.omit()

  # Now do the same for when the branch is at the start, not the end

  filt <- filt %>%
    tidyr::pivot_longer(c(merFirst, merLast), names_to = "mer2Type", values_to = "mer2") %>%
    dplyr::mutate(
      StomataFirst = ifelse(mer2Type == "merFirst", S1, S2),
      StomataLast = ifelse(mer2Type == "merFirst", S2, S3)
    ) %>%
    dplyr::group_by(StomataLast) %>%
    dplyr::mutate(
      CharactersPerStomata = nchar(StomataFirst), # how many letters in a stomata id
      AllConnectedStomata = paste(unique(StomataFirst), collapse = ""), # make a string of all stomata in a 2mer with the last stomata
      NumberOfConnectedStomata = nchar(AllConnectedStomata) / CharactersPerStomata
    ) %>% # how many stomata join this
    dplyr::group_by(AllConnectedStomata) %>%
    dplyr::arrange(AllConnectedStomata, desc(mer3.internal.angle)) %>% # order so the one to keep is top
    dplyr::mutate(rownum = row_number()) %>% # scruffy; figure out which row in group
    dplyr::filter(NumberOfConnectedStomata == 1 | rownum == 1) %>% # keep the first from multi2mers, or every row from unique 2mers
    dplyr::ungroup() %>%
    dplyr::select(-c(StomataFirst, StomataLast, AllConnectedStomata, CharactersPerStomata, NumberOfConnectedStomata, rownum)) %>%
    tidyr::pivot_wider(names_from = "mer2Type", values_from = "mer2") %>%
    na.omit()

  cat("Filtered from", nrow(mer.data$mer3), "to", nrow(filt), "3mers removing kmers with shared terminus\n")

  mer.data$mer3 <- filt

  if (mer.data$write.debug.images) save.plot(plot.3mers(mer.data), mer.data$image.file, ".mer3.5.terminal.png")

  return(mer.data)
}


# create a deBruijn graph from kmers
#
# mer.data - the complete data
create.debruijn.graph <- function(mer.data) {
  # Find the unique 2mers. This assumes the 2mers are consistently ordered (i.e.
  # we cannot find both s1-s2 and s2-s1 in the data)
  uks <- unique(c(mer.data$merFirst, mer.data$merLast)) # unique 2mers
  if (length(uks) == 0) {
    return(make_empty_graph())
  }
  #
  # # Assign an id to each kmer
  k.index <- function(kmer) which(uks == kmer)
  mer.data$merLid <- sapply(mer.data$merFirst, k.index)
  mer.data$merRid <- sapply(mer.data$merLast, k.index)

  # Make a deBruijn graph from kmers
  g <- make_empty_graph()
  g <- add_vertices(g, length(uks))
  V(g)$kmer <- uks

  # Link kmers by edges
  for (i in 1:nrow(mer.data)) {
    g <- g + edges(c(mer.data$merLid[i], mer.data$merRid[i]))
  }
  g
}

# Create contigs from kmers. Each unlinked group of 3mers is separated to a
# separate contig and the stomata within the contig are ordered. Any unlinked
# 2mers are combined into a contig where possible
#
# mer.data - the complete data
create.contigs <- function(mer.data) {
  if (nrow(mer.data$mer3) == 0) {
    return(mer.data$mer3 %>% dplyr::mutate(Contig = list()))
  }

  g <- create.debruijn.graph(mer.data$mer3)

  # Decompose unlinked contigs
  gphs <- decompose.graph(g)

  # Get the number of the contig each 2mer belongs to
  get.contig.number <- function(i) {
    gph <- gphs[[i]]
    contig <- vertex_attr(gph, "kmer")
    contig.num <- rep(i, length(contig))
    names(contig.num) <- contig
    contig.num
  }

  # Assign each contig a number
  contig.numbers <- do.call(c, lapply(1:length(gphs), get.contig.number))
  result <- mer.data$mer3 %>%
    dplyr::mutate(Contig = map_int(merFirst, function(x) contig.numbers[names(contig.numbers) == x]))

  # Keep only the list of 2mers belonging to each contig
  result <- result %>%
    dplyr::select(Contig, merFirst, merLast, Contig) %>%
    tidyr::pivot_longer(c(merFirst, merLast), names_to = "kmer_type", values_to = "kmer") %>%
    dplyr::select(Contig, kmer)

  # What about leftover 2mers that should be in chains? Even more stringent
  # angle filter to get only those that are on the correct orientation
  remaining.2mers <- mer.data$mer2 %>%
    dplyr::filter(!(kmer %in% result$kmer))

  # If there are any valid 2mers left, try to bind them in
  if (nrow(remaining.2mers) > 0) {
    remaining.2mers <- remaining.2mers %>%
      dplyr::filter(angle.is.within.range(
        angle.of.2mer,
        mer.data$modal.2mer.angle,
        ANGLE.DELTA.2MER.STRINGENT
      )) %>%
      dplyr::mutate(Contig = max(result$Contig) + row_number()) %>%
      dplyr::select(Contig, kmer)

    result <- rbind(result, remaining.2mers)
  }

  mer.data$contigs <- merge(result, mer.data$mer2, all.y = FALSE, by = "kmer") %>%
    dplyr::mutate(
      Image = mer.data$image.file,
      Folder = basename(dirname(Image)),
      File = basename(Image)
    )


  cat("Created", length(unique(mer.data$contigs$Contig)), "contigs from 3mers\n")

  if (mer.data$write.chain.image) save.plot(plot.contigs(mer.data), mer.data$image.file, ".result.chains.png")

  mer.data
}

# Calculate straightness, angles and lengths of contigs
#
# mer.data - the complete data
measure.contigs <- function(mer.data) {
  # Calculate average distances per contig
  # and angle variation within the contig
  dist.contigs <- mer.data$contigs %>%
    dplyr::mutate(File = mer.data$image.file) %>%
    dplyr::group_by(Contig) %>%
    # Find the first and last stomata in the contig
    dplyr::mutate(
      Contig.start.index = which.min(S1.x),
      Contig.start.x = S1.x[Contig.start.index],
      Contig.start.y = S1.y[Contig.start.index],
      Contig.end.index = which.max(S2.x),
      Contig.end.x = S2.x[Contig.end.index],
      Contig.end.y = S2.y[Contig.end.index],
      n.stomata.in.contig = n()
    ) %>%
    # Calculate the absolute angle of the contig in the image against the vertical
    # Use this to calculate the angle against the horizontal
    dplyr::rowwise() %>%
    dplyr::mutate(
      Contig.abs.angle = angle.to.horizontal(
        Contig.start.x, Contig.start.y,
        Contig.end.x, Contig.end.y
      ),
      Contig.abs.length = euclidean(Contig.start.x, Contig.start.y, Contig.end.x, Contig.end.y),
      Contig.abs.angle.radians = -deg2rad(Contig.abs.angle),

      # Rotate the end of the contig about the contig start point
      RotatedContigEnd = autoimage::rotate(matrix(c(Contig.end.x, Contig.end.y), nrow = 1),
        Contig.abs.angle.radians,
        pivot = c(Contig.start.x, Contig.start.y)
      ),

      # Rotate the kmer coordinates about the contig start point
      RotatedS1 = autoimage::rotate(matrix(c(S1.x, S1.y), nrow = 1),
        Contig.abs.angle.radians,
        pivot = c(Contig.start.x, Contig.start.y)
      ),
      RotatedS2 = autoimage::rotate(matrix(c(S2.x, S2.y), nrow = 1),
        Contig.abs.angle.radians,
        pivot = c(Contig.start.x, Contig.start.y)
      )
    ) %>%
    # Calculate the deviation between the contig line and the individual points
    dplyr::mutate(
      S1.deviance = RotatedS1[, 2] - Contig.start.y,
      S2.deviance = RotatedS2[, 2] - Contig.start.y
    ) %>%
    # Ensure no duplicate kmers
    dplyr::select(-mer2id) %>%
    dplyr::distinct()

  # How many stomata are not in a chain?
  mer.data$stomata.in.contigs <- unique(c(mer.data$contigs$S1, mer.data$contigs$S2))
  mer.data$stomata.unassigned <- mer.data$mer1$stomata[!(mer.data$mer1$stomata %in% mer.data$stomata.in.contigs)]

  mer.data$nStomata <- length(unique(mer.data$mer1$stomata))
  mer.data$nUnassignedStomata <- length(unique(mer.data$stomata.unassigned))
  mer.data$fUnassignedStomata <- length(unique(mer.data$stomata.unassigned)) / mer.data$nStomata

  mer.data$measured.contigs <- dist.contigs


  if (mer.data$write.debug.images) save.plot(plot.rotated.contigs(mer.data), mer.data$image.file, ".result.rotated.png")

  mer.data
}



#### Main function #####

# Given YOLO 1mers read by read.1mers() for a single image, create contigs
# The 1mers object should contain image, and bounding polygon
process.yolo.predictions <- function(yolo.output.file, write.chain.image = FALSE, write.debug.images = FALSE) {
  # Assuming the yolo inferences are in a tsv file, check if there is already
  # a serialised output and skip any already done
  rds.output <- gsub("tsv", "Rds", yolo.output.file)
  if (file.exists(rds.output)) {
    cat("Skipping existing file", yolo.output.file, "\n")
    return()
  }

  # Store all data in one object for sub-functions to access
  data <- list()
  data$write.chain.image <- write.chain.image
  data$write.debug.images <- write.debug.images

  cat("Detecting stomata inferences in ", yolo.output.file, "\n")
  data$yolo.output.file <- yolo.output.file

  data$mer1 <- read.1mers(yolo.output.file) %>%
    dplyr::mutate(
      Image = str_remove(yolo.output.file, ".tsv"),
      File = basename(yolo.output.file),
      Folder = basename(dirname(yolo.output.file)),
    )
  cat("Detected", nrow(data$mer1), "stomata inferences in ", yolo.output.file, "\n")

  # Read the image for making annotations
  data$image.file <- unique(data$mer1$Image)
  if (!file.exists(data$image.file)) stop(paste("The image file", data$image.file, "was not found"))
  data$img <- OpenImageR::readImage(data$image.file)
  data$img <- OpenImageR::flipImage(data$img, mode = "vertical") # to draw Y-axis as expected

  cat("Analysing stomata patterning in", data$image.file, "\n")

  data <- filter.1mers.for.overlaps(data, MIN.DISTANCE)

  # Filter down the 2mers to those plausibly in chains
  data <- create.2mers(data, min.distance = MIN.DISTANCE, max.distance = MAX.DISTANCE)
  data <- filter.2mers.by.angle(data, ANGLE.DELTA.2MER)
  data <- filter.2mers.by.intersections(data)

  # Join the 2mers to create a 3mer chain
  data <- create.3mers(data)

  # Prune the 3mers. If two 3mers share a 2mer, keep the 3mer that is straightest
  data <- filter.3mers.by.straightness(data, ANGLE.DELTA.INTERNAL.3MER)
  data <- filter.3mers.by.orientation(data, ANGLE.DELTA.MODAL.3MER)
  data <- filter.3mers.by.branch(data)
  data <- filter.3mers.by.terminal(data)


  # Create contigs from the 3mers based on overlapping 2mers
  data <- create.contigs(data)

  # Calculate summary values of contigs
  data <- measure.contigs(data)

  # Summarise metadata for serialisation
  data$metadata <- data.frame(
    Image = data$image.file,
    File = basename(data$image.file),
    Folder = basename(dirname(data$image.file)),
    nStomata = data$nStomata,
    nUnassignedStomata = data$nUnassignedStomata,
    fUnassignedStomata = data$fUnassignedStomata,
    nContigs = length(unique(data$measured.contigs$Contig))
  )

  # The entire dataset is too large to save at scale and mostly not needed.
  # Keep the relevant contig information
  output.data <- list(
    "contigs" = data$contigs,
    "measurments" = data$measured.contigs,
    "metadata" = data$metadata
  )

  saveRDS(output.data, rds.output)

  data
}

#### Run the analysis ####
fs::dir_create("analysis")

# (Use the code below if you have single yolo file containing inferencing of multiple images)
# The complete YOLO output is too large to process in one go - split to a single
# output file per input image. Only needs to be done once for the YOLO output
# split.yolo.output.to.single.image.files <- function(file) {
#   border.data <- readr::read_tsv(file,
#     col_names = TRUE, progress = FALSE,
#     col_types = cols()
#   ) %>%
#     dplyr::group_by(Image) %>%
#     dplyr::mutate(
#       Folder = basename(dirname(Image)),
#       File = basename(Image)
#     ) %>%
#     dplyr::group_by(Folder) %>%
#     dplyr::group_walk(~ fs::dir_create(path = fs::path("analysis", .y$Folder))) %>%
#     dplyr::group_by(Folder, File) %>%
#     dplyr::group_walk(~ write_tsv(.x, file = fs::path("analysis", .y$Folder, .y$File, ext = "tsv")))
# }
# split.yolo.output.to.single.image.files("output.txt")


# Read the split YOLO files that contain border info for a single image each
# Serialise the detected chains for later use
# input.yolo.files <- list.files(path = ".", pattern = "*.tsv", include.dirs = TRUE, full.names = TRUE, recursive = FALSE)
input.yolo.files <- list.files(path = "analysis", pattern = "*.tsv", include.dirs = TRUE, full.names = TRUE, recursive = TRUE)


# Test on just the one
data <- process.yolo.predictions(input.yolo.files[1], write.chain.image = TRUE, write.debug.images = TRUE)

#### Run parallel ####

# Parallelsugar gives parallel syntax for Windows. Behaves like parallel on
# other platforms
parallelsugar::mclapply(input.yolo.files, process.yolo.predictions,
  mc.cores = 4,
  write.chain.image = TRUE, write.debug.images = TRUE
)
