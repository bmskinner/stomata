# Create chains from point coordinates

#### Imports #####
source("functions.R")

#### Constants #####
# min and maximum distance between stomata
MIN.DISTANCE <- 100
MAX.DISTANCE <- 400

# Max angle difference from 180 degrees
ANGLE.DELTA <- 15

#### Functions #####

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
    geom_text(
      data = mer.data$mer1, aes(x = x.com, y = y.com + 20, label = sprintf("%.0f", abs.angle)),
      col = "red", size = 1.5
    ) +
    # geom_text(data = mer.data$mer1, aes(x = x.com, y = y.com-20, label = sprintf("%.0f", horzAngle)),
    #           col = "red", size = 1.5)+

    # Draw the name of each stomata
    geom_text(data = mer.data$mer1, aes(x = x.com, y = y.com, label = stomata), col = "pink1", size = 2) +
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
        label = sprintf("%.0f", abs.angle)
      ),
      col = "red", size = 3
    ) +

    # Draw the centroid of each stomata in the 2mer
    geom_point(data = mer.data$mer2, aes(x = S1.x, y = S1.y), col = "blue", size = 2) +
    geom_point(data = mer.data$mer2, aes(x = S2.x, y = S2.y), col = "orange", size = 1) +

    # Draw the name of each stomata in the 2mer
    geom_text(data = mer.data$mer2, aes(x = S1.x, y = S1.y, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer2, aes(x = S2.x, y = S2.y, label = S2), col = "pink1", size = 2) +
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
        label = sprintf("%.0f", angle)
      ),
      col = "red", size = 3
    ) +
    geom_text(data = mer.data$mer3, aes(x = S1.x, y = S1.y, label = S1), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer3, aes(x = S2.x, y = S2.y, label = S2), col = "pink1", size = 2) +
    geom_text(data = mer.data$mer3, aes(x = S3.x, y = S3.y, label = S3), col = "pink1", size = 2) +
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


# create a deBruijn graph from kmers
create.debruijn.graph <- function(mer.data) {
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

# create contigs from kmers
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
    dplyr::filter(!(kmer %in% result$kmer)) %>%
    dplyr::filter(between(
      abs.angle, mer.data$modal.mer2.angle - 5,
      mer.data$modal.mer2.angle + 5
    )) %>%
    dplyr::mutate(Contig = max(result$Contig) + row_number()) %>%
    dplyr::select(Contig, kmer)

  result <- rbind(result, remaining.2mers)

  merge(result, mer.data$mer2, all.y = FALSE, by = "kmer")

  # return(result)
}

# From given 1-mers, create 2mers. filter to those within a given distance of
# each other and annotate with absolute angles on image
create.2mers <- function(mer.data, min.distance, max.distance) {
  # Create all pairwise combinations of 1mers
  result <- expand.grid(mer.data$mer1$stomata, mer.data$mer1$stomata, stringsAsFactors = F) %>%
    dplyr::distinct() %>% # remove duplicates
    dplyr::filter(Var1 != Var2) %>%
    merge(., mer.data$mer1[, c("stomata", "x.com", "y.com")], by.x = "Var1", by.y = "stomata", all.y = F) %>% # add coordinates for S1
    dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>% # rename for clarity
    merge(., mer.data$mer1[, c("stomata", "x.com", "y.com")], by.x = "S2", by.y = "stomata", all.y = F) %>% # add coordinates for S2
    dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) # rename for clarity


  # Convert stomata CoMs to list of coordinates
  result$S1.point <- mapply(\(x, y, n) list(X = x, Y = y, name = n), result$S1.x, result$S1.y, result$S1, SIMPLIFY = FALSE)
  result$S2.point <- mapply(\(x, y, n) list(X = x, Y = y, name = n), result$S2.x, result$S2.y, result$S2, SIMPLIFY = FALSE)

  choose.first <- function(p1, p2) {
    if (mer.data$modal.angle.type == "stomataAreHorizontal") {
      return(left(p1, p2))
    } else {
      return(upper(p1, p2))
    }
  }

  choose.last <- function(p1, p2) {
    if (mer.data$modal.angle.type == "stomataAreHorizontal") {
      return(right(p1, p2))
    } else {
      return(lower(p1, p2))
    }
  }

  result$Sf <- mapply(choose.first, result$S1.point, result$S2.point, SIMPLIFY = FALSE)
  result$Sl <- mapply(choose.last, result$S1.point, result$S2.point, SIMPLIFY = FALSE)

  result$abs.angle <- mapply(
    \(pf, pl) mer.data$angle.function(pl$X, pl$Y, pf$X, pf$Y),
    result$Sf, result$Sl
  )

  result$length <- mapply(\(p1, p2) euclidean(p1$X, p1$Y, p2$X, p2$Y), result$Sf, result$Sl)

  result$kmer <- mapply(\(p1, p2) paste0(p1$name, p2$name), result$Sf, result$Sl)

  result <- result %>%
    #
    dplyr::select(Sf, Sl, length, abs.angle, kmer) %>%
    unnest_wider(Sf) %>%
    dplyr::rename(S1 = name, S1.x = X, S1.y = Y) %>%
    unnest_wider(Sl) %>%
    dplyr::rename(S2 = name, S2.x = X, S2.y = Y) %>%
    dplyr::filter(S1 != S2 & between(length, min.distance, max.distance)) %>%
    dplyr::distinct()

  if (nrow(result) > 0) {
    result$mer2id <- 1:nrow(result)
  }
  cat("Created", nrow(result), "2mers from", nrow(mer.data$mer1), "1mers\n")
  return(result)
}

# From given 2mers, create 3mers. Filter to those within an angle delta of 180
# degrees
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
      S1S2.angle = abs.angleA, S2S3.angle = abs.angleB
    ) %>% # rename for clarity

    # Calculate angle of 3mer. Check with stomata are top and bottom
    dplyr::rowwise() %>%
    dplyr::mutate(

      # Internal angle of the 3mer (closeness to straight line)
      angle = LearnGeom::Angle(c(S1.x, S1.y), c(S2.x, S2.y), c(S3.x, S3.y)),

      # From crate.2mers we know that S1-S2-S3 is already ordered L-R or T-B
      abs.angle = ifelse(mer.data$modal.angle.type == "vert",
        angle.to.vertical(S1.x, S1.y, S3.x, S3.y),
        angle.to.horizontal(S3.x, S3.y, S1.x, S1.y)
        # angle.to.horizontal(S.bottom.x, S.bottom.y, S.top.x, S.top.y)
      )
    ) %>%
    dplyr::mutate(
      merFirst = paste0(S1, S2),
      merLast = paste0(S2, S3)
    ) %>%
    # Remove the unused columns
    dplyr::distinct()

  cat("Created", nrow(result), "3mers from", nrow(mer.data$mer2), "2mers\n")
  result
}

#### Main function #####

# Given YOLO 1mers read by read.1mers() for a single image, create contigs
# The 1mers object should contain image, bounding polygon, and stomata orientation
process.yolo.predictions <- function(image.1mers, write.chain.image = FALSE, write.debug.images = FALSE) {
  data <- list()
  data$mer1 <- image.1mers
  data$image.file <- unique(image.1mers$Image)
  cat("Analysing", data$image.file, "\n")

  # Check for overlapping objects and keep only the object with highest
  # circularity
  filter.1mers.for.overlaps <- function() {
    # Look for stomata with centroids closer than 50% of the median stomata diameter
    stomata.avg.diameter <- median(data$mer1$max.feret) / 2
    coms <- dplyr::select(data$mer1, stomata, x.com, y.com, circularity)
    distances <- expand.grid(
      s1 = data$mer1$stomata, s2 = data$mer1$stomata,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(s1 != s2) %>%
      merge(., data$mer1, by.x = "s1", by.y = "stomata") %>%
      merge(., data$mer1, by.x = "s2", by.y = "stomata") %>%
      dplyr::mutate(d1d2 = euclidean(x.com.x, y.com.x, x.com.y, y.com.y)) %>%
      dplyr::filter(d1d2 < stomata.avg.diameter) %>%
      dplyr::select(s1, s2, circularity.x, circularity.y) %>%
      dplyr::mutate(toRemove = ifelse(circularity.x < circularity.y, s1, s2))

    data$mer1 <<- data$mer1[!(data$mer1$stomata %in% distances$toRemove), ]
  }
  filter.1mers.for.overlaps()

  # Decide if absolute angles should be calculated from vertical or horizontal
  find.best.angle.for.filtering <- function() {
    dx <- mean(data$mer1$dx)
    dy <- mean(data$mer1$dy)

    if (dy < dx) {
      data$modal.angle.type <<- "stomataAreHorizontal"
      data$angle.function <<- angle.to.vertical
    } else {
      data$modal.angle.type <<- "stomataAreVertical"
      data$angle.function <<- angle.to.horizontal
    }

    # Now we know the correct orientation for measuring, choose which points
    # in each stomata should be used first and last in the calculations
    data$choose.first <<- function(p1, p2) {
      if (data$modal.angle.type == "stomataAreHorizontal") {
        return(left(p1, p2))
      } else {
        return(upper(p1, p2))
      }
    }

    data$choose.last <<- function(p1, p2) {
      if (data$modal.angle.type == "stomataAreHorizontal") {
        return(right(p1, p2))
      } else {
        return(lower(p1, p2))
      }
    }

    p1 <- mapply(data$choose.first, data$mer1$p1, data$mer1$p2, SIMPLIFY = FALSE)
    p2 <- mapply(data$choose.last, data$mer1$p1, data$mer1$p2, SIMPLIFY = FALSE)

    data$mer1$p1 <<- p1
    data$mer1$p2 <<- p2

    # Calculate the stomata angle given the type of data
    data$mer1$abs.angle <<- mapply(\(p1, p2) data$angle.function(p2$X, p2$Y, p1$X, p1$Y),
      data$mer1$p1, data$mer1$p2,
      SIMPLIFY = TRUE
    )

    # sapply(1:nrow(data$mer1), \(i){
    # data$angle.function(data$mer1$p2[[i]]["X"], data$mer1$p2[[i]]["Y"],
    #                     data$mer1$p1[[i]]["X"], data$mer1$p1[[i]]["Y"])
    # })

    data$modal.stomata.angle <<- find.mode(data$mer1$abs.angle)
  }
  find.best.angle.for.filtering()

  # Read the image for making annotations
  if (!file.exists(data$image.file)) stop(paste("The image file", data$image.file, "was not found"))
  data$img <- OpenImageR::readImage(data$image.file)
  data$img <- OpenImageR::flipImage(data$img, mode = "vertical") # to draw as expected

  if (write.debug.images) save.plot(plot.1mers(data, include.stomata = FALSE), data$image.file, ".mer1.raw.png")

  if (write.debug.images) save.plot(plot.1mers(data), data$image.file, ".mer1.points.png")

  data$mer2 <- create.2mers(data, min.distance = MIN.DISTANCE, max.distance = MAX.DISTANCE)

  if (write.debug.images) save.plot(plot.2mers(data), data$image.file, ".mer2.raw.png")

  # Filter the 2mers to the orientation of stomata in the image
  # Which orientation do we expect? The stomata are oriented with their long
  # diameter aligned with the chain. 2mers should be close to the modal
  # stomata orientation.
  filter.2mers.by.angle <- function() {
    angle.tolerance <- 20

    filt <- data$mer2 %>% dplyr::filter(between(
      abs.angle,
      data$modal.stomata.angle - angle.tolerance,
      data$modal.stomata.angle + angle.tolerance
    ))

    data$modal.mer2.angle <<- find.mode(filt$abs.angle)

    filt <- filt %>% dplyr::filter(between(
      abs.angle,
      data$modal.mer2.angle - angle.tolerance,
      data$modal.mer2.angle + angle.tolerance
    ))

    cat("Filtered from", nrow(data$mer2), "to", nrow(filt), "2mers\n")
    filt
  }

  data$mer2 <- filter.2mers.by.angle()

  # Filter 2mers intersecting a third stomata
  filter.2mers.by.intersections <- function() {
    # Remove edges in the graph that intersect a third stomata
    # can use sf: https://stackoverflow.com/questions/61703791/determine-lines-that-intersect-a-polygon-in-r
    # Create line objects
    lines <- lapply(1:nrow(data$mer2), function(i) {
      sf::st_linestring(rbind(
        c(data$mer2$S1.x[i], data$mer2$S1.y[i]),
        c(data$mer2$S2.x[i], data$mer2$S2.y[i])
      ))
    })

    # Check each 2mer for intersections with a stomata
    # If we intersect, then this is not a valid 2mer
    data$mer2$intersects <- lapply(lines, function(l) sum(unlist(sapply(data$mer1$polygons, sf::st_intersects, y = l))))

    filt <- data$mer2 %>% dplyr::filter(intersects <= 2)

    cat("Filtered from", nrow(data$mer2), "to", nrow(filt), "2mers\n")

    return(filt)
  }

  data$mer2 <- filter.2mers.by.intersections()

  if (write.debug.images) save.plot(plot.2mers(data), data$image.file, ".mer2.filt.png")

  # Join the 2mers to create a 3mer chain
  data$mer3 <- create.3mers(data)
  if (write.debug.images) save.plot(plot.3mers(data), data$image.file, ".mer3.raw.png")

  # Prune the 3mers. If two 3mers share a 2mer, keep the 3mer that is straightest
  filter.3mers.by.straightness <- function() {
    filt <- data$mer3 %>%
      # First keep only the 3mers in straight lines
      dplyr::filter(angle > 180 - ANGLE.DELTA)

    data$modal.mer3.angle <<- find.mode(filt$abs.angle)

    # angle.tolerance <- 20
    # # Now keep only the 3mers that are at the modal angle
    # filt <- filt %>%
    #   dplyr::filter( between(abs.angle,
    #                          data$modal.mer3.angle - angle.tolerance,
    #                          data$modal.mer3.angle + angle.tolerance))
    #
    #
    # # Check the length of the 3mers - any jumping rows will be much longer
    # data$mean.mer3.length <<- mean(filt$S1S2 + filt$S2S3)
    # filt <- filt %>% dplyr::filter(S1S2+S2S3 < (data$mean.mer3.length * 1.5))

    filt <- filt %>%
      # Look for 2mers present more than once
      # Drop the 3mers with the lowest angle.
      dplyr::group_by(merFirst) %>%
      dplyr::arrange(merFirst, desc(angle)) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::group_by(merLast) %>%
      dplyr::arrange(merLast, desc(angle)) %>%
      dplyr::slice_head(n = 1)

    cat("Filtered from", nrow(data$mer3), "to", nrow(filt), "3mers\n")
    filt
  }

  data$mer3 <- filter.3mers.by.straightness()
  if (write.debug.images) save.plot(plot.3mers(data), data$image.file, ".mer3.straight.png")

  filter.3mers.by.branch <- function() {
    # Pruning step - remove any kmers where the endpoint is in the middle of a chain
    # to prevent branches
    mer3.longer <- data$mer3 %>%
      tidyr::pivot_longer(c(S1, S2, S3), names_to = "StomataPosition", values_to = "Stomata") %>%
      dplyr::group_by(StomataPosition, Stomata) %>%
      dplyr::mutate(abs.mer3.angle.diff = abs(data$modal.mer3.angle - abs.angle)) %>%
      dplyr::arrange(StomataPosition, Stomata, abs.mer3.angle.diff) %>%
      dplyr::slice_head(n = 1) %>% # take only the 3mer closest to median angle
      tidyr::pivot_wider(names_from = StomataPosition, values_from = Stomata) %>%
      na.omit() # remove the rows with NAs due to our slice



    cat("Filtered from", nrow(data$mer3), "to", nrow(mer3.longer), "3mers\n")
    return(mer3.longer)
  }

  data$mer3 <- filter.3mers.by.branch()
  if (write.debug.images) save.plot(plot.3mers(data), data$image.file, ".mer3.prune.png")

  filter.3mers.by.terminal <- function() {
    # This still leaves branches when the endpoints meet a terminal kmer.
    # Prune again - cases where the branch meets a terminal 2mer
    # Look for 1mers that are are the S2 of 2 different 2mers
    # TODO: does this need refactor with top/bottom vs left right?
    mer3.terminal <- data$mer3 %>%
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
      dplyr::arrange(AllConnectedStomata, desc(angle)) %>% # order so the one to keep is top
      dplyr::mutate(rownum = row_number()) %>% # scruffy; figure out which row in group
      dplyr::filter(NumberOfConnectedStomata == 1 | rownum == 1) %>% # keep the first from multi2mers, or every row from unique 2mers
      dplyr::ungroup() %>%
      dplyr::select(-c(StomataFirst, StomataLast, AllConnectedStomata, CharactersPerStomata, NumberOfConnectedStomata, rownum)) %>%
      tidyr::pivot_wider(names_from = "mer2Type", values_from = "mer2") %>%
      na.omit()

    # Now do the same for when the branch is at the start, not the end

    mer3.terminal <- mer3.terminal %>%
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
      dplyr::arrange(AllConnectedStomata, desc(angle)) %>% # order so the one to keep is top
      dplyr::mutate(rownum = row_number()) %>% # scruffy; figure out which row in group
      dplyr::filter(NumberOfConnectedStomata == 1 | rownum == 1) %>% # keep the first from multi2mers, or every row from unique 2mers
      dplyr::ungroup() %>%
      dplyr::select(-c(StomataFirst, StomataLast, AllConnectedStomata, CharactersPerStomata, NumberOfConnectedStomata, rownum)) %>%
      tidyr::pivot_wider(names_from = "mer2Type", values_from = "mer2") %>%
      na.omit()

    cat("Filtered from", nrow(data$mer3), "to", nrow(mer3.terminal), "3mers\n")
    return(mer3.terminal)
  }
  data$mer3 <- filter.3mers.by.terminal()
  if (write.debug.images) save.plot(plot.3mers(data), data$image.file, ".mer3.terminal.png")


  # Create contigs from the 3mers
  data$contigs <- create.contigs(data) %>%
    dplyr::mutate(
      Image = data$image.file,
      Folder = basename(dirname(Image)),
      File = basename(Image)
    )

  if (write.chain.image) save.plot(plot.contigs(data), data$image.file, ".chains.png")

  measure.contigs <- function() {
    # Calculate average distances per contig
    # and angle variation within the contig
    dist.contigs <- data$contigs %>%
      dplyr::mutate(File = data$image.file) %>%
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
    data$stomata.in.contigs <<- unique(c(data$contigs$S1, data$contigs$S2))
    data$stomata.unassigned <<- data$mer1$stomata[!(data$mer1$stomata %in% data$stomata.in.contigs)]

    data$nStomata <<- function() length(unique(data$mer1$stomata))
    data$nUnassignedStomata <<- function() length(unique(data$stomata.unassigned))
    data$fUnassignedStomata <<- function() length(unique(data$stomata.unassigned)) / data$nStomata()

    return(dist.contigs)
  }

  data$measured.contigs <- measure.contigs()

  data$metadata <- data.frame(
    Image = data$image.file,
    File = basename(data$image.file),
    Folder = basename(dirname(data$image.file)),
    nStomata = data$nStomata(),
    nUnassignedStomata = data$nUnassignedStomata(),
    fUnassignedStomata = data$fUnassignedStomata(),
    nContigs = length(unique(data$measured.contigs$Contig))
  )

  if (write.debug.images) save.plot(plot.rotated.contigs(data), data$image.file, ".rotated.png")
  data
}


#### Run the analysis ####
fs::dir_create("analysis")

# The complete YOLO output is too large to process in one go - split to a single
# output file per input image. Only needs to be done once for the YOLO output
split.yolo.output.to.single.image.files <- function(file) {
  border.data <- readr::read_tsv(file,
    col_names = TRUE, progress = FALSE,
    col_types = cols()
  ) %>%
    dplyr::group_by(Image) %>%
    dplyr::mutate(
      Folder = basename(dirname(Image)),
      File = basename(Image)
    ) %>%
    dplyr::group_by(Folder) %>%
    dplyr::group_walk(~ fs::dir_create(path = fs::path("analysis", .y$Folder))) %>%
    dplyr::group_by(Folder, File) %>%
    dplyr::group_walk(~ write_tsv(.x, file = fs::path("analysis", .y$Folder, .y$File, ext = "tsv")))
}
# split.yolo.output.to.single.image.files("output.txt")


# Read the split YOLO files that contain border info for a single image each
# Serialise the detected chains for later use
# input.yolo.files <- list.files(path = ".", pattern = "*.tsv", include.dirs = TRUE, full.names = TRUE, recursive = FALSE)
input.yolo.files <- list.files(path = "analysis", pattern = "*.tsv", include.dirs = TRUE, full.names = TRUE, recursive = TRUE)
input.yolo.files <- input.yolo.files[1:2]

process.image.file <- function(image.file, write.chain.image = FALSE, write.debug.images = FALSE) {
  rds.output <- gsub("tsv", "Rds", image.file)
  if (file.exists(rds.output)) {
    cat("Skipping existing file", image.file, "\n")
    return()
  } # skip any already done

  yolo.data <- read.1mers(image.file) %>%
    dplyr::mutate(
      Image = str_replace(Image, "/home/bs19022/projects/stomata/", ""),
      Folder = basename(dirname(Image)),
      File = basename(Image)
    )

  processed.data <- suppressWarnings(process.yolo.predictions(yolo.data,
    write.chain.image = write.chain.image,
    write.debug.images = write.debug.images
  ))

  1 # The entire dataset is too large to save at scale and mostly not needed
  output <- list(
    "contigs" = processed.data$contigs,
    "measurments" = processed.data$measured.contigs,
    "metadata" = processed.data$metadata
  )

  saveRDS(output, rds.output)
}

# Parallelsugar gives parallel syntax for Windows. Behaves like parallel on
# other platforms
parallelsugar::mclapply(input.yolo.files, process.image.file,
  mc.cores = 4,
  write.chain.image = TRUE, write.debug.images = TRUE
)

# Test on just the first
process.image.file(input.yolo.files[1], write.chain.image = TRUE, write.debug.images = TRUE)
