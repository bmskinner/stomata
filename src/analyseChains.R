# Parse Rds output from analysed YOLO stomata coordinates
source("src/functions.R")
fs::dir_create("figure")

# Create a linestring from two st_points
#
# p1 - the first point
# p2 - the second point
# Returns: an st_linestring
st_point2linestring <- function(p1, p2) {
  # cat("Making linestring from points\n")
  sf::st_linestring(matrix(
    c(p1[, 1], p1[, 2], p2[, 1], p2[, 2]),
    byrow = TRUE, ncol = 2
  ))
}

# Project the chain start-end line to the image boundary. Measure the fractional
# length occupied by the chain. This accounts for how much of the chain we could
# image depending on the orientation of the leaf.
#
# p1 - the start of an oriented contig
# p2 - the end of an oriented contig
# image.bounds - the oriented image bounds as an st_polygon
calculate.max.linestring.length <- function(p1, p2, image.bounds) {
  chainLinestring <- st_point2linestring(p1, p2)
  if (sf::st_length(chainLinestring) < 1) {
    cat("Chain zero-length: ", chainLinestring, "\n")
  }
  # cat("Extending linestring\n")
  chainLinestringExt <- extend.linestring(chainLinestring, distance = 6000)
  # cat("Intersecting image bounds\n")
  chainIntersected <- sf::st_intersection(chainLinestringExt, image.bounds)
  sf::st_length(chainIntersected)
}


#### Measurement functions ####

# Calculate straightness, angles and lengths of contigs
#
# mer.data - the complete data
measure.chains <- function(mer.data) {
  mer.data$measured.contigs <- mer.data$oriented.contigs |>
    # Calculate the deviation between the contig line and the individual points
    dplyr::mutate(
      S1.deviance = RotatedS1[, 2] - OrientedContigStart[, 2],
      S2.deviance = RotatedS2[, 2] - OrientedContigStart[, 2]
    ) |>
    dplyr::select(-mer2id) |>
    # Measure brokenness of chains
    dplyr::group_by(Contig) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      chainMaxPossibleLength = calculate.max.linestring.length(
        OrientedContigStart,
        OrientedContigEnd,
        mer.data$oriented.image.bounds
      ),
      chainActualLength = sf::st_length(st_point2linestring(OrientedContigStart, OrientedContigEnd)),
      chainLengthFraction = chainActualLength / chainMaxPossibleLength
    )

  mer.data$chain.measurements <- mer.data$measured.contigs |>
    dplyr::group_by(Folder, File, Contig) |>
    dplyr::summarise(
      nStomataInChain = length(unique(S1, S2)),
      medianInterStomataDistance = median(length),
      medianChain2merAngle = median(angle.of.2mer),
      chainMaxPossibleLength = unique(chainMaxPossibleLength),
      chainActualLength = unique(chainActualLength),
      chainLengthFraction = unique(chainLengthFraction),
      .groups = "drop_last"
    )


  # These measurements are aggregated per-image

  mer.data$image.measurements <- data.frame(
    Folder = unique(mer.data$chain.measurements$Folder),
    File = unique(mer.data$chain.measurements$File),
    nStomataInChains = length(unique(c(mer.data$oriented.contigs$S1, mer.data$oriented.contigs$S2)))
  ) |>
    dplyr::mutate(
      nStomataUnassigned = length(mer.data$mer1$stomata[!(mer.data$mer1$stomata %in% unique(c(mer.data$oriented.contigs$S1, mer.data$oriented.contigs$S2)))]),
      nStomataInImage = length(unique(mer.data$mer1$stomata)),
      fStomataUnassigned = nStomataUnassigned / nStomataInImage,
      fAreaOfChainBounds = sf::st_area(mer.data$combined.chain.bounds) / sf::st_area(mer.data$oriented.image.bounds),
      nChainsInImage = length(unique(mer.data$chain.measurements$Contig))
    )

  # if (mer.data$write.debug.images) save.plot(plot.rotated.contigs(mer.data), mer.data$image.file, ".result.rotated.png")

  mer.data
}


#### Read serialised data from chains ####

input.chain.files <- list.files(
  path = "analysis", pattern = "*.Rds",
  include.dirs = TRUE, full.names = TRUE, recursive = TRUE
)

# Read the given serialised chain data and make chain measurements
read.chain.data <- function(file) {
  input.data <- readRDS(file)
  measure.chains(input.data)
}

input.data <- lapply(input.chain.files, read.chain.data)
stomata.measurements <- do.call(rbind, lapply(input.data, \(x) x$mer1))
image.measurements <- do.call(rbind, lapply(input.data, \(x) x$image.measurements))
chain.measurements <- do.call(rbind, lapply(input.data, \(x) x$chain.measurements))

cat(sprintf(
  "%i stomata analysed from %i images, aggregated into %i chains\n",
  sum(image.measurements$nStomataInImage), nrow(image.measurements), sum(image.measurements$nChainsInImage)
))

#### Summary at image level ####

# Number of stomata per image
save.ggplot(
  ggplot(image.measurements, aes(x = Folder, y = nStomataInImage)) +
    geom_hline(yintercept = median(image.measurements$nStomataInImage)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Stomata per image") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Stomata_per_image.png"
)


# Unassigned stomata per image - those not in chains
save.ggplot(
  ggplot(image.measurements, aes(x = Folder, y = fStomataUnassigned)) +
    geom_hline(yintercept = median(image.measurements$fStomataUnassigned)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Fraction of stomata not in chains") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Stomata_not_in_chains_per_image.png"
)

# Area of image within chain bounds
save.ggplot(
  ggplot(image.measurements, aes(x = Folder, y = fAreaOfChainBounds)) +
    geom_hline(yintercept = median(image.measurements$fAreaOfChainBounds)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Area of image within chain bounds") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Chain_bounds_per_image.png"
)

#### Summary at chain level ####

save.ggplot(
  ggplot(image.measurements, aes(x = Folder, y = nChainsInImage)) +
    geom_hline(yintercept = median(image.measurements$nChainsInImage)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Chains per image") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Chains_per_image.png"
)

# Length of the chains in pixels
save.ggplot(
  ggplot(chain.measurements, aes(x = Folder, y = medianInterStomataDistance)) +
    geom_hline(yintercept = median(chain.measurements$medianInterStomataDistance)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Distance between stomata pairs (pixels)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Distance_between_stomata_pixels.png"
)

# Fractional chain length against the possible in that image
save.ggplot(
  ggplot(chain.measurements, aes(x = Folder, y = chainLengthFraction)) +
    geom_hline(yintercept = median(chain.measurements$chainLengthFraction)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Fraction of maximum measurable chain length") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Fractional_chain_length.png"
)

# Stomata per chain - affected by orientation of the image
save.ggplot(
  ggplot(chain.measurements, aes(x = medianChain2merAngle, y = nStomataInChain)) +
    annotate("rect", xmin = 0, xmax = 45, ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
    annotate("rect", xmin = 135, xmax = 180, ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, 180, 30)) +
    geom_point() +
    labs(x = "Median 2mer angle to vertical within chain (degrees)", y = "Stomata per chain") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Stomata_per_chain.png"
)

#### Summary at stomata level ####

# Area of stomata
save.ggplot(
  ggplot(stomata.measurements, aes(x = Folder, y = area)) +
    geom_hline(yintercept = median(stomata.measurements$area)) +
    geom_violin() +
    geom_boxplot(width = 0.2, alpha = 0) +
    labs(y = "Area of stomata bounding box (pixels)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  "figure/Stomata_area_per_image.png"
)

#### Detailed analysis ####

# Are the distances significantly different?
# kruskal.data <- data.frame("Folder" =  chain.measure$Folder, "length" = chain.measure$length)
# kruskal.test(length ~ Folder, data  = kruskal.data)
# kruskal.data.dunn <- rstatix::dunn_test(kruskal.data, formula = length ~ Folder, detailed = T)


# Look at deviances in each contig
deviance.data <- chain.measure %>%
  dplyr::select(Contig, File, Folder, S1, S2, S1.deviance, S2.deviance) %>%
  tidyr::pivot_longer(cols = c(S1, S2), names_to = "mer1", values_to = "Stomata") %>%
  tidyr::pivot_longer(cols = c(S1.deviance, S2.deviance), names_to = "dev_type", values_to = "Deviance") %>%
  dplyr::filter((mer1 == "S1" & dev_type == "S1.deviance") | (mer1 == "S2" & dev_type == "S2.deviance")) %>%
  dplyr::select(-mer1, -dev_type) %>%
  dplyr::distinct() %>%
  dplyr::group_by(File, Folder, Contig) %>%
  dplyr::summarise(
    nStomata = n(),
    SumDeviance = sum(Deviance),
    SumAbsDeviance = sum(abs(Deviance)),
    MeanDeviance = SumDeviance / nStomata,
    MeanAbsDeviance = SumAbsDeviance / nStomata,
    RootSumSqareDeviance = sqrt(sum(Deviance^2)),
    RootMeanSquareDeviance = sqrt(sum(Deviance^2) / nStomata),
    DevianceRatio = MeanDeviance / MeanAbsDeviance
  )

# Overall levels of deviance from straight line in contig
save.ggplot(
  ggplot(deviance.data, aes(x = Folder, y = MeanDeviance)) +
    geom_hline(yintercept = 0) +
    geom_violin() +
    # geom_beeswarm()+
    geom_boxplot(width = 0.2, alpha = 0) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ),
  "figure/Deviance_per_folder.png"
)

# Mean deviance versus sum of squares - consistency of bends
save.ggplot(
  ggplot(deviance.data, aes(x = MeanAbsDeviance, y = MeanDeviance)) +
    geom_hex(bins = 100) +
    geom_abline(intercept = c(0, 0), slope = 1, col = "grey") + # Consistent bend
    geom_abline(intercept = c(0, 0), slope = -1, col = "grey") + # Consistent bend
    geom_abline(intercept = c(0, 0), slope = 0.5, col = "grey") + # Wibble
    geom_abline(intercept = c(0, 0), slope = -0.5, col = "grey") + # Wibble
    # geom_point()+
    coord_fixed(xlim = c(0, 30), ylim = c(-30, 30), expand = FALSE) +
    scale_fill_viridis_c() +
    labs(fill = "Number of chains", x = "Mean absolute deviance (pixels)", y = "Mean deviance (pixels)") +
    # facet_wrap(~File)+
    theme_bw() +
    theme(legend.position = "top"),
  "figure/Deviance_consistency.png"
)

#### G-function analysis ####

# Run the g-function on each image
# Get the CoMs of each stomata
coms <- chain.contigs %>% dplyr::select(S = S1, x = S1.x, y = S1.y, Image, Folder, File)
coms <- rbind(coms, chain.contigs %>% dplyr::select(S = S2, x = S2.x, y = S2.y, Image, Folder, File))
coms %<>% dplyr::distinct()

g.results <- do.call(rbind, lapply(unique(coms$Image), function(i) g.function(as.matrix(coms[coms$Image == i, c("x", "y")]), i)))

g.results$Folder <- basename(dirname(g.results$Image))

# Plot the per-file plots
save.ggplot(
  ggplot(g.results, aes(x = distance, y = Gd, group = Image)) +
    geom_line(data = g.results[, c("distance", "Gd", "Image")], col = "grey", alpha = 0.1) +
    geom_line(col = "blue", alpha = 0.1) +
    labs(x = "Distance between stomata", y = "Cumulative fraction") +
    facet_wrap(~Folder) +
    theme_bw(),
  "figure/G-function.png",
  height = 170
)

# What is the mean value for each folder?
# Linear interpolation of per-file curve to consistent spacing

windows <- as.data.frame(IRanges(
  start = seq(40, 400, by = 10), # vector of window start positions
  end = seq(50, 410, by = 10)
))

calc.mean <- function(start, end) {
  do.call(rbind, lapply(unique(g.results$Folder), function(folder) {
    subset.data <- g.results[g.results$Folder == folder & g.results$distance > start & g.results$distance <= end, ]
    data.frame(start = start, end = end, Gd = mean(subset.data$Gd), Gd.sd = sd(subset.data$Gd), Folder = folder, distance = (start + end) / 2)
  }))
}

g.summary <- do.call(rbind, mapply(calc.mean, windows$start, windows$end, SIMPLIFY = FALSE))

save.ggplot(
  ggplot(g.summary, aes(x = distance, y = Gd)) +
    geom_line(data = g.results[, c("distance", "Gd", "Image")], aes(group = Image), col = "grey", alpha = 0.1) +
    geom_line(data = g.results, aes(group = Image), col = "lightblue", alpha = 0.3) +
    geom_line(col = "blue") +
    labs(x = "Distance between stomata", y = "Cumulative fraction") +
    facet_wrap(~Folder) +
    theme_bw(),
  "figure/G-function_complete.png",
  height = 170
)

save.ggplot(
  ggplot(g.summary, aes(x = distance, y = Gd, col = Folder)) +
    geom_line() +
    labs(x = "Distance between stomata", y = "Cumulative fraction") +
    theme_bw() +
    theme(legend.position = "top"),
  "figure/G-function_summary.png",
  height = 170
)
