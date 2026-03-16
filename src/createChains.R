# Create chains from point coordinates

#### Imports #####
source("src/functions.R")

#### Constants #####
# min and maximum distance between stomata in pixels
MIN.DISTANCE <- 50
MAX.DISTANCE <- 500

# 2mer angle variation from the modal value in degrees
ANGLE.DELTA.2MER <- 20

# How far can a 3mer deviate from a 180 line in degrees?
ANGLE.DELTA.INTERNAL.3MER <- 10
# How much can the absolute angle of a 3mer deviate from the modal value?
ANGLE.DELTA.MODAL.3MER <- 10

# If we have remaining 2mers that might be part of a contig, how stringent an angle
# must they have to the modal 2mer angle in degrees?
ANGLE.DELTA.2MER.STRINGENT <- 8

#### Plotting Functions #####

# Create a ggplot with the image drawn as a background layer
create.background.image.plot <- function(mer.data) {
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
    cf
}

plot.1mers <- function(mer.data, include.stomata = TRUE) {
  p <- create.background.image.plot(mer.data)

  stomata <- data.frame(
    "x1" = sapply(mer.data$mer1$p1, \(p) p$X),
    "y1" = sapply(mer.data$mer1$p1, \(p) p$Y),
    "x2" = sapply(mer.data$mer1$p2, \(p) p$X),
    "y2" = sapply(mer.data$mer1$p2, \(p) p$Y)
  )

  p <- p +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Draw the name of each stomata
    geom_text(data = mer.data$mer1, aes(x = x.com, y = y.com - 15, label = stomata), col = "pink1", size = 2) +
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
  create.background.image.plot(mer.data) +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Draw the 2mer
    geom_segment(
      data = mer.data$mer2, aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y),
      col = "blue", linewidth = 0.5, arrow = arrow(type = "closed", angle = 20, length = unit(3, "mm"))
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
    geom_text(data = mer.data$mer2, aes(x = S1.x, y = S1.y - 15, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer2, aes(x = S2.x, y = S2.y - 15, label = S2), col = "pink1", size = 2) +
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
  create.background.image.plot(mer.data) +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +

    # Draw the centroid of each stomata
    geom_point(data = mer.data$mer1, aes(x = x.com, y = y.com), col = "black") +

    # Draw the 3mer
    geom_segment(
      data = mer.data$mer3, aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y),
      col = "blue", linewidth = 1.5, alpha = 0.8, arrow = arrow(type = "closed", angle = 20, length = unit(3, "mm"))
    ) +
    geom_segment(
      data = mer.data$mer3, aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y),
      col = "orange", linewidth = 0.5, arrow = arrow(type = "closed", angle = 20, length = unit(3, "mm"))
    ) +
    geom_text(
      data = mer.data$mer3, aes(
        x = S2.x, y = S2.y,
        label = sprintf("%.0f°\n\n|%.0f°|", mer3.internal.angle, mer3.angle.to.vertical)
      ),
      col = "red", size = 3
    ) +
    geom_text(data = mer.data$mer3, aes(x = S1.x, y = S1.y - 15, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer3, aes(x = S2.x, y = S2.y - 15, label = S2), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer3, aes(x = S3.x, y = S3.y - 15, label = S3), col = "pink1", size = 2) +
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

plot.graph.2mers <- function(mer.data) {
  create.background.image.plot(mer.data) +

    # Draw the border polygons
    geom_sf(data = st_multipolygon(mer.data$mer1$polygons), fill = "darkgreen", alpha = 0.6) +
    # Draw the 2mer
    geom_segment(
      data = mer.data$contigable.2mers, aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y, colour = source),
      linewidth = 0.5, arrow = arrow(type = "closed", angle = 20, length = unit(3, "mm"))
    ) +
    scale_colour_manual(values = c("In 3mer" = "black", "Standalone" = "blue")) +

    # Write the angle of the 2mer
    geom_text(
      data = mer.data$contigable.2mers, aes(
        x = (S2.x + S1.x) / 2, y = (S2.y + S1.y) / 2,
        label = paste(sprintf("%.1f°\n%.0fpx", angle.of.2mer, length))
      ),
      col = "red", size = 3
    ) +

    # Draw the centroid of each stomata in the 2mer
    geom_point(data = mer.data$contigable.2mers, aes(x = S1.x, y = S1.y), col = "blue", size = 2) +
    geom_point(data = mer.data$contigable.2mers, aes(x = S2.x, y = S2.y), col = "orange", size = 1) +

    # Draw the name of each stomata in the 2mer
    geom_text(data = mer.data$contigable.2mers, aes(x = S1.x, y = S1.y - 15, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$contigable.2mers, aes(x = S2.x, y = S2.y - 15, label = S2), col = "pink1", size = 2) +
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

# plot contig data
plot.contigs <- function(mer.data) {
  create.background.image.plot(mer.data) +

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

# Plot chains after rotation to horizontal. Includes bounding boxes of
# overlapping chains
#
# mer.data - the complete data
plot.rotated.contigs <- function(mer.data) {
  ggplot(mer.data$oriented.contigs) +
    geom_sf(data = mer.data$combined.chain.rectangles, fill = "pink", alpha = 0.6) +
    geom_sf(data = st_polygon(mer.data$oriented.image.bounds), col = "black", fill = NA) +
    geom_sf(data = st_multipolygon(mer.data$rotated.1mers$RotatedPolygon), fill = "darkgreen", alpha = 0.6) +
    geom_segment(
      aes(
        x = RotatedContigStart[, 1],
        y = RotatedContigStart[, 2],
        xend = RotatedContigEnd[, 1],
        yend = RotatedContigEnd[, 2]
      ),
      col = "grey20"
    ) +
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
      dplyr::filter(!is_empty(has.overlap)) |> # may be close but non-overlapping
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
    dplyr::filter(Var1 != Var2) |>
    dplyr::distinct() |>
    merge(mer.data$mer1[, c("stomata", "x.com", "y.com")], by.x = "Var1", by.y = "stomata", all.y = F) |> # add coordinates for S1
    dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>% # rename for clarity
    merge(mer.data$mer1[, c("stomata", "x.com", "y.com")], by.x = "S2", by.y = "stomata", all.y = F) |> # add coordinates for S2
    dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) |> # rename for clarity
    # Apply a length filter to exclude 2mers that are too distant
    # The plurality of the remaining edges will be orientated along the chains
    dplyr::mutate(length = euclidean(S1.x, S1.y, S2.x, S2.y)) |>
    # Calculate the angle of a 2mer based on the locations of the points at either end
    # Wrap at 180 degrees
    dplyr::mutate(angle.of.2mer = angle.to.vertical(S1.x, S1.y, S2.x, S2.y) %% 180) |>
    dplyr::filter(between(length, min.distance, max.distance)) |>
    dplyr::distinct()

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
    dplyr::mutate(kmer = paste0(S1, S2)) |>
    dplyr::group_by(kmer) |>
    dplyr::slice_head(n = 1) |>
    dplyr::mutate(mer2id = row_number())

  if (mer.data$write.debug.images) {
    dup.data <- data.frame("angle.of.2mer" = c(result$angle.of.2mer, result$angle.of.2mer + 180))

    histo.plot <- ggplot(dup.data, aes(x = angle.of.2mer)) +
      geom_histogram(aes(y = after_stat(density)), binwidth = 5) +
      geom_density() +
      geom_vline(xintercept = mer.data$modal.2mer.angle) +
      geom_vline(xintercept = mer.data$modal.2mer.angle + 180) +
      labs(title = sprintf("Modal angle %.2f°", mer.data$modal.2mer.angle)) +
      scale_x_continuous(breaks = seq(0, 360, 45)) +
      theme_bw()
    save.plot(histo.plot, mer.data$image.file, ".diagnostic.mer2.histo.png")
  }

  cat(
    "Created", nrow(result), "2mers from", nrow(mer.data$mer1), "1mers within length bounds",
    min.distance, "-", max.distance, "pixels\n"
  )

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
  # Create line objects and buffer to a rectangle somewhat narrower than a stomata
  median.stomata.radius <- median(mer.data$mer1$perimeter / (2 * pi)) * 0.4

  lines <- lapply(1:nrow(mer.data$mer2), function(i) {
    sf::st_buffer(sf::st_linestring(rbind(
      c(mer.data$mer2$S1.x[i], mer.data$mer2$S1.y[i]),
      c(mer.data$mer2$S2.x[i], mer.data$mer2$S2.y[i])
    )), dist = median.stomata.radius)
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
      geom_histogram(aes(y = after_stat(density)), binwidth = 5) +
      geom_density() +
      scale_x_continuous(breaks = seq(0, 360, 5)) +
      theme_bw()
    save.plot(histo.plot, mer.data$image.file, ".diagnostic.mer3.angle.histo.png")
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
      geom_histogram(aes(y = after_stat(density)), binwidth = 5) +
      geom_density() +
      geom_vline(xintercept = mer.data$modal.mer3.angle) +
      labs(title = sprintf("Modal angle %.2f°", mer.data$modal.mer3.angle)) +
      scale_x_continuous(breaks = seq(0, 360, 5)) +
      theme_bw()
    save.plot(histo.plot, mer.data$image.file, ".diagnostic.mer3.angle.vertical.histo.png")
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


# Join 2mers into a graph based on shared 3mers.
#
# mer.data - the complete data
create.2mer.graph <- function(mer.data) {
  # Find the 2mers that are present in valid 3mers. This assumes the 2mers are
  # consistently ordered (i.e. we cannot find both s1-s2 and s2-s1 in the data).
  unique.2mers.within.3mers <- data.frame(
    kmer = unique(c(mer.data$mer3$merFirst, mer.data$mer3$merLast)),
    source = "In 3mer"
  )

  if (nrow(unique.2mers.within.3mers) == 0) {
    return(make_empty_graph())
  }

  cat(sprintf(
    "Found %i unique 2mers within the remaining %i 3mers\n",
    nrow(unique.2mers.within.3mers), nrow(mer.data$mer3)
  ))

  # Create a table of 2mers that are valid for inclusion in contigs
  mer.data$contigable.2mers <- mer.data$mer2 |>
    dplyr::filter(kmer %in% unique.2mers.within.3mers$kmer)

  # There may still be 2mers remaining that we want to include in the graph.
  # e.g. standalone 2mers that are not part of a 3mer but still part of a chain.
  # (a chain of only 2 stomata, or 2mers connecting other chains)
  standalone.2mers <- mer.data$mer2 |>
    dplyr::filter(
      !(kmer %in% mer.data$contigable.2mers$kmer),
      intersects <= 2,
      angle.is.within.range(
        angle.of.2mer,
        mer.data$modal.mer3.angle,
        ANGLE.DELTA.2MER.STRINGENT
      )
    )

  # Of these potential 2mers, some are invalid - for example, they start at a
  # stomata that is already in a 2mer. Ensure we do not have kmers starting or
  # ending on the same stomata
  standalone.2mers <- standalone.2mers |>
    dplyr::filter(!(S1 %in% mer.data$contigable.2mers$S1)) |>
    dplyr::filter(!(S2 %in% mer.data$contigable.2mers$S2)) |>
    # We may have added in 2mers that start or end on the same stomata. Keep
    # only the shortest 2mers. TODO: replace with angle test
    dplyr::group_by(S1) |>
    dplyr::arrange(length) |>
    dplyr::slice_head(n = 1) |>
    dplyr::group_by(S2) |>
    dplyr::arrange(length) |>
    dplyr::slice_head(n = 1)

  # Add these standalone 2mers to those to be contiged
  mer.data$contigable.2mers <- rbind(mer.data$contigable.2mers, standalone.2mers)

  # Annotate the reason each 2mer was included for debugging
  debug.kmer.selection <- unique.2mers.within.3mers
  if (nrow(standalone.2mers) > 0) {
    debug.kmer.selection <- rbind(unique.2mers.within.3mers, data.frame(
      kmer = standalone.2mers$kmer,
      source = "Standalone"
    ))
  }

  mer.data$contigable.2mers <- mer.data$contigable.2mers |>
    merge(debug.kmer.selection, by = "kmer")

  cat(sprintf(
    "Selected %i 2mers: %i 2mers within valid 3mers, %i 2mers within %.2f° of the modal 3mer angle (%.2f°)\n",
    nrow(mer.data$contigable.2mers), nrow(unique.2mers.within.3mers),
    nrow(standalone.2mers), ANGLE.DELTA.2MER.STRINGENT,
    mer.data$modal.mer3.angle
  ))


  # We can at this point use the 2mers for graph construction with the stomata
  # as nodes
  all.1mers <- unique(c(mer.data$mer2$S1, mer.data$mer2$S2))

  cat(sprintf(
    "Creating graph from %i stomata connected by %i 2mers\n",
    length(all.1mers), nrow(mer.data$contigable.2mers)
  ))

  # Find the vertex id of each 1mer in the graph
  stomata.index <- function(stomata.id) which(all.1mers == stomata.id)
  mer.data$contigable.2mers$S1.id <- sapply(mer.data$contigable.2mers$S1, stomata.index)
  mer.data$contigable.2mers$S2.id <- sapply(mer.data$contigable.2mers$S2, stomata.index)

  # Set stomata as vertices in the network
  mer.data$mer2.graph <- igraph::make_empty_graph()
  mer.data$mer2.graph <- igraph::add_vertices(mer.data$mer2.graph, length(all.1mers))
  igraph::V(mer.data$mer2.graph)$mer1 <- all.1mers

  # Connect by 2mer edges if the 2mer is in the desired subset
  for (i in 1:nrow(mer.data$contigable.2mers)) {
    mer.data$mer2.graph <- igraph::add_edges(mer.data$mer2.graph,
      c(mer.data$contigable.2mers$S1.id[i], mer.data$contigable.2mers$S2.id[i]),
      mer2 = mer.data$contigable.2mers$kmer[i]
    )
  }

  cat(sprintf(
    "Graph contains %i nodes and %i edges\n",
    length(V(mer.data$mer2.graph)), length(E(mer.data$mer2.graph))
  ))

  if (mer.data$write.debug.images) save.plot(plot.graph.2mers(mer.data), mer.data$image.file, ".mer3.6.contigs.png")

  mer.data
}

# Create contigs from stomata connected by valid 2mers
#
# mer.data - the complete data
create.2mer.contigs <- function(mer.data) {
  if (nrow(mer.data$mer2) == 0) {
    return(mer.data$mer2 %>% dplyr::mutate(Contig = list()))
  }

  mer.data <- create.2mer.graph(mer.data)

  # Decompose unlinked contigs
  contig.subgraphs <- igraph::decompose(mer.data$mer2.graph)
  cat(sprintf("Decomposed graph into %s subgraphs\n", length(contig.subgraphs)))

  # Create unique identifier for each contig
  # Get the number of the contig each stomata belongs to
  assign.contig.number <- function(i) {
    subgraph <- contig.subgraphs[[i]]
    stomata.id <- vertex_attr(subgraph, "mer1")
    mer2.id <- edge_attr(subgraph, "mer2")
    data.frame(kmer = mer2.id, Contig = rep(i, length(mer2.id)))
  }

  # Assign each stomata a contig number
  contig.mer2.map <- do.call(rbind, lapply(1:length(contig.subgraphs), assign.contig.number))
  cat(sprintf("Decomposed graph into %i contigs\n", length(unique(contig.mer2.map$Contig))))

  mer.data$contigs <- merge(contig.mer2.map, mer.data$contigable.2mers,
    all.x = TRUE, by = "kmer"
  ) |>
    dplyr::mutate(
      Image = mer.data$image.file,
      Folder = basename(dirname(Image)),
      File = basename(Image)
    )

  cat("Created", length(unique(mer.data$contigs$Contig)), "contigs from 2mers\n")

  if (mer.data$write.chain.image) save.plot(plot.contigs(mer.data), mer.data$image.file, ".result.chains.png")

  mer.data
}

# Create contigs from kmers. Each unlinked group of 3mers is separated to a
# separate contig and the stomata within the contig are ordered. Any unlinked
# 2mers are combined into a contig where possible
#
# mer.data - the complete data
# create.contigs <- function(mer.data) {
#   if (nrow(mer.data$mer3) == 0) {
#     return(mer.data$mer3 %>% dplyr::mutate(Contig = list()))
#   }
#
#   complete.graph <- create.2mer.graph(mer.data)
#
#   # Decompose unlinked contigs
#   contig.subgraphs <- igraph::decompose(complete.graph)
#
#
#   # Get the number of the contig each 2mer belongs to
#   get.contig.number <- function(i) {
#     gph <- contig.subgraphs[[i]]
#     contig <- vertex_attr(gph, "kmer")
#     contig.num <- rep(i, length(contig))
#     names(contig.num) <- contig
#     contig.num
#   }
#
#   # Assign each contig a number
#   contig.numbers <- do.call(c, lapply(1:length(contig.subgraphs), get.contig.number))
#   result <- mer.data$mer3 %>%
#     dplyr::mutate(Contig = map_int(merFirst, function(x) contig.numbers[names(contig.numbers) == x]))
#
#   # Keep only the list of 2mers belonging to each contig
#   result <- result %>%
#     dplyr::select(Contig, merFirst, merLast, Contig) %>%
#     tidyr::pivot_longer(c(merFirst, merLast), names_to = "kmer_type", values_to = "kmer") %>%
#     dplyr::select(Contig, kmer)
#
#   # What about leftover 2mers that should be in chains? Even more stringent
#   # angle filter to get only those that are on the correct orientation
#   remaining.2mers <- mer.data$mer2 %>%
#     dplyr::filter(!(kmer %in% result$kmer))
#
#   # If there are any valid 2mers left, try to bind them in
#   if (nrow(remaining.2mers) > 0) {
#     remaining.2mers <- remaining.2mers |>
#       dplyr::filter(angle.is.within.range(
#         angle.of.2mer,
#         mer.data$modal.2mer.angle,
#         ANGLE.DELTA.2MER.STRINGENT
#       ))
#
#     remaining.2mers <- remaining.2mers |>
#       dplyr::mutate(Contig = max(result$Contig) + row_number()) %>%
#       dplyr::select(Contig, kmer)
#
#     result <- rbind(result, remaining.2mers)
#   }
#
#   mer.data$contigs <- merge(result, mer.data$mer2, all.y = FALSE, by = "kmer") %>%
#     dplyr::mutate(
#       Image = mer.data$image.file,
#       Folder = basename(dirname(Image)),
#       File = basename(Image)
#     )
#
#
#   cat("Created", length(unique(mer.data$contigs$Contig)), "contigs from 3mers\n")
#
#   if (mer.data$write.chain.image) save.plot(plot.contigs(mer.data), mer.data$image.file, ".result.chains.png")
#
#   mer.data
# }

# Orient contigs consistently on the horizontal image axis. Calculates rotations
# needed to visualise stomata with chains horizontal
#
# mer.data - the complete data
orient.contigs <- function(mer.data) {
  # We need the image dimensions and centre as pivot point for rotation
  image.width <- dim(mer.data$img)[2]
  image.height <- dim(mer.data$img)[1]
  image.xcom <- image.width / 2
  image.ycom <- image.height / 2
  image.centre <- c(image.xcom, image.ycom)

  # Create a rectangle defining the image bounds. We will use this later to
  # limit measurement of rotated objects to regions within the original image.
  mer.data$image.bounds <- sf::st_polygon(list(
    matrix(
      c(
        image.width, 0,
        image.width, image.height,
        0, image.height,
        0, 0,
        image.width, 0
      ),
      ncol = 2, byrow = TRUE
    )
  ))

  # Rotate the given xy coordinates around the centre of the current image
  #
  # x - an sf object, a matrix of XY coordinates, or an x coordinate
  # y - a y coordinate. If missing, x is assumed to be an sf or matrix
  # degrees - the angle for rotation
  # Returns: sf_polygon, sf_point or matrix of rotated coordinates
  rotate.about.image.com <- function(x, y = NULL, degrees) {
    # if we have an st_polygon, we only need the first two columns
    if (any(class(x) == "POLYGON")) {
      return(sf::st_polygon(list(autoimage::rotate(x[[1]][, 1:2],
        -deg2rad(degrees),
        pivot = c(image.xcom, image.ycom)
      ))))
    }

    if (any(class(x) == "POINT")) {
      return(sf::st_point(autoimage::rotate(x[1, ],
        -deg2rad(degrees),
        pivot = c(image.xcom, image.ycom)
      )))
    }

    # if we have a point already as a matrix of x and y elements
    if (is.matrix(x)) {
      return(sf::st_point(autoimage::rotate(x,
        -deg2rad(degrees),
        pivot = c(image.xcom, image.ycom)
      )))
    }
    # # Otherwise create a matrix from x and y coordinates and return as st_point
    rotate.about.image.com(x = matrix(c(x, y), nrow = 1), degrees = degrees)
  }

  # Calculate the orientation of the contigs. This proceeds in two steps; first,
  # apply the modal 3mer angle rotation. This makes most of the contigs mostly
  # horizontal. We can then decide which stomata is at the start and end of each
  # contig. The overall contig orientation is the straight line from start-end.
  # This probably does not match the overall 3mer orientation perfectly, because
  # the 3mers wibble around within the contigs. We can calculate the median
  # overall angle of the contigs, and apply this second angle correction. The
  # contigs are now as horizontal as we can make them.
  mer.data$oriented.contigs <- mer.data$contigs |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # Rotate each 2mer coordinate about the contig start point
      RotatedS1 = rotate.about.image.com(
        S1.x, S1.y,
        mer.data$modal.mer3.angle - 90
      ),
      RotatedS2 = rotate.about.image.com(
        S2.x, S2.y,
        mer.data$modal.mer3.angle - 90
      )
    ) |>
    # The contigs are now roughly horizontal.
    # Identify the endpoints of the contig
    dplyr::group_by(Contig) |>
    # What is the contig orientation? Left-right or right-left?
    # Find the endpoints of the contig
    dplyr::mutate(
      Contig.left = min(RotatedS1[, 1]) < min(RotatedS2[, 1]),
      Contig.start.index = ifelse(Contig.left, which.min(RotatedS1[, 1]), which.min(RotatedS2[, 1])),
      Contig.end.index = ifelse(Contig.left, which.max(RotatedS1[, 1]), which.max(RotatedS2[, 1])),
      Contig.start.x = ifelse(Contig.left, RotatedS1[Contig.start.index, 1], RotatedS2[Contig.start.index, 1]),
      Contig.start.y = ifelse(Contig.left, RotatedS1[Contig.start.index, 2], RotatedS2[Contig.start.index, 2]),
      Contig.end.x = ifelse(Contig.left, RotatedS2[Contig.end.index, 1], RotatedS1[Contig.end.index, 1]),
      Contig.end.y = ifelse(Contig.left, RotatedS2[Contig.end.index, 2], RotatedS1[Contig.end.index, 2])
    ) |>
    # Find the endpoints of the contig
    # dplyr::mutate(
    #   Contig.start.index = which.min(RotatedS1[, 1]),
    #   Contig.end.index = which.max(RotatedS2[, 1]),
    #   Contig.start.x = RotatedS1[Contig.start.index, 1],
    #   Contig.start.y = RotatedS1[Contig.start.index, 2],
    #   Contig.end.x = RotatedS2[Contig.end.index, 1],
    #   Contig.end.y = RotatedS2[Contig.end.index, 2]
    # ) |>
    # Now find the remaining offset angle for the contigs. Calculate the median
    # across all contigs.
    dplyr::rowwise() |>
    dplyr::mutate(
      OrientedContigStart = list(sf::st_point(x = c(Contig.start.x, Contig.start.y))),
      OrientedContigEnd = list(sf::st_point(x = c(Contig.end.x, Contig.end.y))),
      Contig.abs.angle = angle.to.horizontal(
        Contig.start.x, Contig.start.y,
        Contig.end.x, Contig.end.y
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(contig.median.angle = median(Contig.abs.angle)) |>
    # Adjust the rotation of 2mers again to make the contigs as close to
    # horizontal as possible.
    dplyr::rowwise() |>
    dplyr::mutate(
      # Rotate each 2mer coordinate again about the image centre
      RotatedS1 = rotate.about.image.com(x = RotatedS1[, 1], y = RotatedS1[, 2], degrees = contig.median.angle),
      RotatedS2 = rotate.about.image.com(x = RotatedS2[, 1], y = RotatedS2[, 2], degrees = contig.median.angle),

      # Rotate the contig again about the image centre
      OrientedContigStart = rotate.about.image.com(x = Contig.start.x, y = Contig.start.y, degrees = contig.median.angle),
      OrientedContigEnd = rotate.about.image.com(x = Contig.end.x, y = Contig.end.y, degrees = contig.median.angle)
    ) |>
    dplyr::select(-(Contig.start.x:Contig.end.y))

  # We now know the total rotation needing to be applied to the image
  total.rotation <- mer.data$modal.mer3.angle - 90 + unique(mer.data$oriented.contigs$contig.median.angle)

  # Create the oriented image bounding box
  mer.data$oriented.image.bounds <- rotate.about.image.com(
    x = mer.data$image.bounds,
    degrees = total.rotation
  )

  # Rotate the individual stomata polygons
  # Make a rotated copy of the polygons for display
  mer.data$rotated.1mers <- mer.data$mer1 |>
    dplyr::rowwise() |>
    dplyr::mutate(
      OrientedPolygon = list(rotate.about.image.com(
        x = polygons,
        degrees = total.rotation
      )),

      # Also note the min and max y position of each rotated polygon. We will use
      # this to establish the shared rectangles covering each chain
      polygon.min.y = min(st_coordinates(OrientedPolygon)[, 2]),
      polygon.max.y = max(st_coordinates(OrientedPolygon)[, 2])
    )

  # mer.data$rotated.1mers$OrientedPolygon <- lapply(
  #   mer.data$rotated.1mers$OrientedCoordinates,
  #   function(x) sf::st_polygon(list(x))
  # )


  # Now create the unioned rectangles along each chain
  mer.data$oriented.chain.rectangles <- mapply(
    function(y.min, y.max) {
      sf::st_polygon(list(
        matrix(
          c(
            min(st_coordinates(mer.data$oriented.image.bounds)[, 1]), y.min,
            max(st_coordinates(mer.data$oriented.image.bounds)[, 1]), y.min,
            max(st_coordinates(mer.data$oriented.image.bounds)[, 1]), y.max,
            min(st_coordinates(mer.data$oriented.image.bounds)[, 1]), y.max,
            min(st_coordinates(mer.data$oriented.image.bounds)[, 1]), y.min
          ),
          ncol = 2, byrow = TRUE
        )
      ))
    },
    y.min = mer.data$rotated.1mers$polygon.min.y,
    y.max = mer.data$rotated.1mers$polygon.max.y,
    SIMPLIFY = FALSE
  )

  # Combine all overlapping chain rectangles
  mer.data$combined.chain.rectangles <- sf::st_simplify(sf::st_combine(do.call("c", mer.data$oriented.chain.rectangles)))

  # Intersect chain rectangles with the image bounds
  mer.data$combined.chain.rectangles <- sf::st_intersection(mer.data$combined.chain.rectangles, mer.data$oriented.image.bounds)


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

  # Create contigs from the 2mers based on overlapping 3mers
  data <- create.2mer.contigs(data)
  data <- orient.contigs(data)

  # The entire dataset is too large to save at scale and mostly not needed.
  # Keep the relevant contig and stomata information
  output.data <- list(
    "oriented.contigs" = data$oriented.contigs,
    "combined.chain.bounds" = data$combined.chain.rectangles,
    "oriented.image.bounds" = data$oriented.image.bounds,
    "mer1" = data$mer1
  )

  saveRDS(output.data, rds.output)

  data
}

#### Run the analysis ####
# fs::dir_create("analysis")

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


# Read the split YOLO files that contain border info for a single image each.
# We assume there that the YOLO outputs and raw images are in ./analysis
# Serialise the detected chains for later use
# input.yolo.files <- list.files(path = ".", pattern = "*.tsv", include.dirs = TRUE, full.names = TRUE, recursive = FALSE)
input.yolo.files <- list.files(path = "analysis", pattern = "*.tsv", include.dirs = TRUE, full.names = TRUE, recursive = TRUE)


# Test on just one
# data <- process.yolo.predictions(input.yolo.files[27], write.chain.image = TRUE, write.debug.images = TRUE)

#### Run parallel ####

# Parallelsugar gives parallel syntax for Windows. Behaves like parallel on
# other platforms
parallelsugar::mclapply(input.yolo.files, process.yolo.predictions,
  mc.cores = 4,
  write.chain.image = TRUE, write.debug.images = FALSE
)
