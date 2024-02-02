# Parse the JSON from Anylabelling to extract stomata CoMs

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
library(filesstrings)
library(autoimage)

MAX.DISTANCE <- 300
TRIM.DISTANCE <- 500
SHORT.DISTANCE <- 350
ANGLE.DELTA <- 10

# Calculate centre of mass from JSON
get.border <- function(file){
  data <- jsonlite::read_json(file)
  shape.number = 1
  
  process.shape <- function(shape){
    if(shape$shape_type!="polygon")  return(NULL)
    
    bounds <- split(unlist(shape$points, recursive = T), 1:2) %>% as.data.frame
    colnames(bounds) <- c("X", "Y")
    result <- data.frame("shape" = shape.number,
               "x" = bounds$X,
               "y" = bounds$Y,
               "file" = file)
    shape.number <<- shape.number+1
    return(result)
  }
  
  
  do.call(rbind, lapply(data$shapes, process.shape)) %>% 
    as.data.frame 
}

calc.coms <- function(border.data){
  border.data %>%
    dplyr::group_by(shape, file) %>%
    dplyr::summarise(x.com = mean(x), 
                     y.com = mean(y)) %>%
    dplyr::mutate(stomata  = paste0("s", sprintf("%02d", shape))) %>%
    dplyr::ungroup() %>%
    dplyr::select(stomata, x.com, y.com) %>%
    as.data.frame
  
}

# Calculate distance between two points
euclidean <- function(x1, y1, x2, y2) sqrt( (x1-x2)^2 + (y1-y2)^2)

# convert degrees to radians
deg2rad <- function(deg)(deg * pi) / (180)

# Find the maximum value in the density plot of the given vector
find.mode <-  function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

plot.2mers <- function(mer.data, img, border.data){
  img.grob <- rasterGrob(img, interpolate=TRUE)
  ggplot(mer.data)+
    annotation_custom(img.grob, xmin=0, xmax=dim(img)[2], ymin=0, ymax=dim(img)[1]) +
    coord_fixed(xlim = c(0, dim(img)[2]), ylim = c(0, dim(img)[1]))+
    scale_color_viridis_c()+
    geom_point(data = border.data, aes(x = x, y = y), col = "green", size=0.1)+
    geom_point(aes(x = S1.x, y = S1.y))+
    geom_point(aes(x = S2.x, y = S2.y))+
    geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y, col = S1S2))+
    geom_text(aes(x = (S2.x+S1.x)/2, y = (S2.y+S1.y)/2, label = sprintf("%.0f", S1S2)), col = "blue", size = 2)+
    theme_bw()+
    theme(axis.title = element_blank())
}

plot.3mers <- function(mer.data, img, border.data){
  img.grob <- rasterGrob(img, interpolate=TRUE)
  ggplot(mer.data)+
    annotation_custom(img.grob, xmin=0, xmax=dim(img)[2], ymin=0, ymax=dim(img)[1]) +
    coord_fixed(xlim = c(0, dim(img)[2]), ylim = c(0, dim(img)[1]))+
    scale_color_viridis_c()+
    geom_point(data = border.data, aes(x = x, y = y), col = "green", size=0.1)+
    geom_point(aes(x = S1.x, y = S1.y))+
    geom_point(aes(x = S2.x, y = S2.y))+
    geom_point(aes(x = S3.x, y = S3.y))+
    geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y, col = S1S2))+
    geom_segment(aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y, col = S2S3)) +
    geom_text(aes(x = S2.x, y = S2.y+50, label = sprintf("%.0f", angle)), col = "red", size = 2)+
    geom_text(aes(x = (S2.x+S1.x)/2, y = (S2.y+S1.y)/2, label = sprintf("%.0f", S1S2)), col = "blue", size = 2)+
    geom_text(aes(x = (S2.x+S3.x)/2, y = (S2.y+S3.y)/2, label = sprintf("%.0f", S2S3)), col = "blue", size = 2)+
    geom_text(aes(x = S1.x, y = S1.y, label = S1), col = "purple", size = 2)+
    geom_text(aes(x = S2.x, y = S2.y, label = S2), col = "purple", size = 2)+
    geom_text(aes(x = S3.x, y = S3.y, label = S3), col = "purple", size = 2)+
    theme_bw()+
    theme(axis.title = element_blank())
}

create.debruijn.graph <- function(mer.data){
  uks <- unique(c(mer.data$merL, mer.data$merR)) # unique 2mers
  if(length(uks)==0)
    return(make_empty_graph())
  # 
  # # Assign an id to each kmer
  k.index <- function(kmer) which(uks==kmer)
  mer.data$merLid <- sapply(mer.data$merL, k.index)
  mer.data$merRid <- sapply(mer.data$merR, k.index)
  
  # Make a deBruijn graph from kmers
  g <- make_empty_graph()
  g <- add_vertices(g, length(uks))
  V(g)$kmer <- uks
  
  # Link kmers by edges
  for(i in 1:nrow(mer.data)){
    g <- g + edges( c(mer.data$merLid[i],mer.data$merRid[i]))
  }
  g
}

create.contigs <- function(mer.data){
  if(nrow(mer.data)==0) return(mer.data %>% dplyr::mutate(Contig = list()))
 
   g <- create.debruijn.graph(mer.data)

  # Visualise the graph
  layout <- layout_with_kk(g)
  plot(g, layout = layout)
  
  # Decompose unlinked contigs
  gphs <- decompose.graph(g)
  
  # Get the number of the contig each 2mer belongs to
  get.contig.number <- function(i){
    gph <- gphs[[i]]
    contig <- vertex_attr(gph, "kmer")
    contig.num <- rep(i, length(contig))
    names(contig.num) <- contig
    contig.num
  }
  
  contig.numbers <- do.call(c, sapply(1:length(gphs), get.contig.number))
  
  mer.data %>%
    dplyr::mutate(Contig = map_int(merL, function(x) contig.numbers[names(contig.numbers)==x]))
}

plot.contigs <- function(mer.contigs, img, border.data){
  img.grob <- rasterGrob(img, interpolate=TRUE)
  # visualise the contigs
  ggplot(mer.contigs, aes(col = as.factor(Contig)))+
    annotation_custom(img.grob, xmin=0, xmax=dim(img)[2], ymin=0, ymax=dim(img)[1]) +
    coord_fixed(xlim = c(0, dim(img)[2]), ylim = c(0, dim(img)[1]))+
    geom_point(data = border.data, aes(x = x, y = y), col = "green", size=0.1)+
    geom_point(aes(x = S1.x, y = S1.y), size=2)+
    geom_point(aes(x = S2.x, y = S2.y), size=2)+
    geom_point(aes(x = S3.x, y = S3.y), size=2)+
    geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y), linewidth=1) +
    geom_segment(aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y), linewidth=1) +
    labs(col = "Chain")+
    theme_bw()+
    theme(axis.title = element_blank())
}

# From a table of centres of mass, create 2mers
# filter to those within a given distance of each other
# and annotate with absolute angles on image
create.2mers <- function(coms, min.distance = 100, max.distance = 300){
  data <- expand.grid(coms$stomata, coms$stomata, stringsAsFactors = F) %>%

    dplyr::distinct() %>% # remove duplicates
    dplyr::filter(Var1!=Var2) %>%
    merge(., coms[,1:3], by.x = "Var1", by.y = "stomata", all.y = F) %>% # add coordinates for S1
    dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>% # rename for clarity
    merge(., coms[,1:3], by.x = "S2", by.y = "stomata", all.y = F) %>% # add coordinates for S2
    dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) %>% # rename for clarity
    
    # now reorder S1 and S2 such that S1 is always to the left of S2
    dplyr::rowwise() %>%
    dplyr::mutate(S1t   = ifelse(S1.x < S2.x, S1, S2), 
                  S1t.x = ifelse(S1.x < S2.x, S1.x, S2.x), 
                  S1t.y = ifelse(S1.x < S2.x, S1.y, S2.y), 
                  S2t   = ifelse(S1.x < S2.x, S2, S1),
                  S2t.x = ifelse(S1.x < S2.x, S2.x, S1.x),
                  S2t.y = ifelse(S1.x < S2.x, S2.y, S1.y)) %>% # ensure x asc order
    dplyr::ungroup() %>%
    dplyr::select(S1 = S1t, S1.x = S1t.x, S1.y = S1t.y,
                  S2 = S2t, S2.x = S2t.x, S2.y = S2t.y) %>% # remove unused columns
    dplyr::rowwise() %>%
    dplyr::mutate(S1S2 = euclidean(S1.x, S1.y, S2.x, S2.y),
                  kmer = paste0(S1, S2), 
                  abs.angle = LearnGeom::Angle(c(S1.x, S1.y+10), c(S1.x,  S1.y),  c(S2.x,S2.y))) %>% # angle of 2mer relative to image
    dplyr::filter(S1 != S2 &
                  between(S1S2, min.distance, max.distance))
  data$mer2id <- 1:nrow(data)
  return(data)
}

# From a table of 2mers, create 3mers. Filter to those within
# an angle delta of 180 degrees
create.3mers <- function(mer2){
  merge(mer2, mer2, by.x = c("S2","S2.x", "S2.y" ), by.y = c("S1","S1.x", "S1.y" ), 
        all.y = F, suffixes = c("A", "B")) %>%
    dplyr::select(S1, S1.x, S1.y, S2, S2.x, S2.y, S3 = S2B, S3.x = S2.xB, 
                  S3.y = S2.yB, S1S2 = S1S2A, S2S3 = S1S2B, 
                  S1S2.abs.angle = abs.angleA, S2S3.abs.angle = abs.angleB) %>% # rename for clarity
    dplyr::rowwise() %>%
    dplyr::mutate(angle = LearnGeom::Angle(c(S1.x,  S1.y), c(S2.x, S2.y), c(S3.x,S3.y)),
                  abs.angle = LearnGeom::Angle(c(S1.x, S1.y+10), c(S1.x,  S1.y),  c(S3.x,S3.y))) %>% 
    dplyr::mutate(merL = paste0(S1, S2),
                  merR = paste0(S2, S3))
}

# Given a json file, extract the stomata and 
# create linear chains
process.json.file <- function(file){
  cat("Analysing", file, "\n")
  # Read the image with stomata
  image <- str_replace(file, "json", "jpg")
  img <- readImage(image)
  img <- flipImage (img, mode = "vertical")
  
  border.data <- get.border(file)
  
  mer1 <- calc.coms(border.data) # centres of mass as points

  # Create distance table between pairs of stomata
  # Filter to only those within a given distance
  cat("  Calculating distances\n")
  mer2 <- create.2mers(mer1, min.distance = 100, max.distance = TRIM.DISTANCE)
  mer.2.plot <- plot.2mers(mer2, img, border.data)
  ggsave(str_replace(file, ".json", ".mer2.raw.png"), plot = mer.2.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  # Join the tables to create a 3-mer chain
  mer3 <- create.3mers(mer2)
  mer.3.plot <- plot.3mers(mer3, img, border.data)
  ggsave(str_replace(file, ".json", ".mer3.raw.png"), plot = mer.3.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  # Find the 3mers in straight lines
  mer3.straight <- mer3 %>% dplyr::filter( angle > 180 - ANGLE.DELTA)
  mer.3.straight.plot <- plot.3mers(mer3.straight, img, border.data)
  ggsave(str_replace(file, ".json", ".mer3.straight.png"), plot = mer.3.straight.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  # Filter the straight 3mers to the most common orientation in the image
  modal.angle <- find.mode(mer3.straight$abs.angle)
  mer3.filt <- mer3.straight %>% dplyr::filter(between(abs.angle, modal.angle - ANGLE.DELTA, modal.angle+ANGLE.DELTA))
  
  # Check the distribution of angles in mer3 vs mer2
  density.plot <- ggplot()+
    geom_density(data = mer2, aes(x = abs.angle, col="mer2 raw" ), alpha=0) +
    geom_density(data = mer3, aes(x = abs.angle, col="mer3 raw" ), alpha=0) +
    geom_density(data = mer3.straight, aes(x = abs.angle, col="mer3 straight")) +
    geom_density(data = mer3.filt, aes(x = abs.angle, col="mer3 oriented" )) +
    labs(x = "Angle of 2mer to vertical") +
    scale_color_manual(values = c("black", "red", "green", "blue"))+
    theme_bw()+
    theme(legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          legend.background = element_blank())
  ggsave(str_replace(file, ".json", ".angle_density.png"), plot = density.plot, dpi = 300, units = "mm", width = 100, height = 85)
  
  
  mer.3.plot <- plot.3mers(mer3.filt, img, border.data)
  ggsave(str_replace(file, ".json", ".mer3.filtered.png"), plot = mer.3.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  # Shows we can't make the graph in one step - the distance values between chains and non-chains can overlap
  # What if we link the close kmers, then widen up and elongate chains with more distant kmers?
  
  short.3mers <- mer3.filt %>% dplyr::filter(S1S2 < SHORT.DISTANCE & S2S3 < SHORT.DISTANCE)
  short.mer.3.plot <- plot.3mers(short.3mers, img, border.data)
  ggsave(str_replace(file, ".json", ".mer3.short.png"), plot = short.mer.3.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  short.contigs <- create.contigs(short.3mers) %>% 
    dplyr::group_by(Contig)
  
  
  # Pruning step - remove any kmers where the endpoint is in the middle of a chain
  # to prevent branching
  
  # Look for 1mers present more than 3 times (should be in positions S1, S2, S3 at most)
  # If a 1 mer is present in more than 3 3mers, drop the 3mer with the highest angle.
  # TODO
  
  chain.plot <- plot.contigs(short.contigs, img, border.data)
  
  # Match contigs to the 2mers they contain
  contig.assignment <- short.contigs %>% dplyr::select(Contig, merL, merR) %>%
    tidyr::pivot_longer(c(merL, merR), values_to = "kmer") %>%
    dplyr::select(Contig, kmer) %>%
    dplyr::distinct() %>%
    merge(., mer2, by= "kmer")
  
  # Calculate average distances per contig
  # and angle variation within the contig
  dist.contigs <- contig.assignment %>% 
    dplyr::mutate(File  = file) %>%
    dplyr::group_by(Contig) %>%
    
    # Find the first and last stomata in the contig
    dplyr::mutate(Contig.start.index = which.min(S1.x),
                  Contig.start.x = S1.x[Contig.start.index],
                  Contig.start.y = S1.y[Contig.start.index],
                  Contig.end.index = which.max(S2.x),
                  Contig.end.x = S2.x[Contig.end.index],
                  Contig.end.y = S2.y[Contig.end.index]) %>%
    
    # Calculate the absolute angle of the contig in the image against the vertical
    # Use this to calculate the angle against the horizontal
    dplyr::rowwise() %>%
    dplyr::mutate(Contig.abs.angle = 90 - LearnGeom::Angle(c(Contig.start.x, Contig.start.y+10), 
                                                      c(Contig.start.x,  Contig.start.y),  
                                                      c(Contig.end.x, Contig.end.y)),
                  Contig.abs.length = euclidean(Contig.start.x, Contig.start.y, Contig.end.x, Contig.end.y),
                  Contig.abs.angle.radians = deg2rad(-Contig.abs.angle),
                  
                  # Rotate the end of the contig about the contig start point
                  RotatedContigEnd = autoimage::rotate(matrix(c(Contig.end.x, Contig.end.y), nrow=1), 
                                                       Contig.abs.angle.radians, 
                                                       pivot = c(Contig.start.x, Contig.start.y)),
                  
                  # Rotate the kmer coordinates about the contig start point
                  RotatedS1 = autoimage::rotate(matrix(c(S1.x, S1.y), nrow=1), 
                                              Contig.abs.angle.radians, 
                                              pivot = c(Contig.start.x, Contig.start.y)),
                  RotatedS2 = autoimage::rotate(matrix(c(S2.x, S2.y), nrow=1), 
                                                Contig.abs.angle.radians, 
                                                pivot = c(Contig.start.x, Contig.start.y))
                  ) %>%
    
    # Calculate the deviation between the contig line and the individual points
    dplyr::mutate(S1.deviance = RotatedS1[,2] - Contig.start.y,
                  S2.deviance = RotatedS2[,2] - Contig.start.y) %>%

    # Ensure no duplicate kmers
    dplyr::select(-mer2id) %>%
    dplyr::distinct()
  
  rotated.plot <- ggplot(dist.contigs)+
    geom_segment(aes(x = Contig.start.x, y = Contig.start.y, xend = RotatedContigEnd[,1], yend = RotatedContigEnd[,2]), col = "black")+
    geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y), col = "grey")+
    geom_segment(aes(x = RotatedS1[,1], y = RotatedS1[,2], xend = RotatedS2[,1], yend = RotatedS2[,2]), col = "blue")+
    theme_bw()
  
  ggsave(str_replace(file, ".json", ".chain.rotated.png"), plot = rotated.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  chain.file <- str_replace(file, ".json", ".chains.png")
  cat("  Plotting contigs\n")
  ggsave(chain.file, plot = chain.plot, dpi = 300, units = "mm", width = 170, height = 140)
  return(dist.contigs)
}


# Read all json files
files <- list.files(path = "stomatal_image", pattern = "*.json",  recursive = T, full.names = T)

chain.distances <- do.call(rbind, lapply(files, process.json.file))
chain.distances$Folder <- dirname(chain.distances$File)

# dplyr::group_by(Contig) %>%
#   dplyr::summarise(MeanDistance  = mean(S1S2),
#                    MedianDistance = median(S1S2),
#                    SD = sd(S1S2),
#                    nPairs = n()) %>%


dist.plot <- ggplot(chain.distances, aes(x=Folder, y = S1S2))+
  geom_beeswarm(col="darkgrey")+
  geom_boxplot(alpha=0)+
  labs(y = "Distance between stomata pairs (pixels)")+
  theme_bw()+
  theme(axis.text.x = element_blank())

# Look at deviances in each contig
deviance.data <- chain.distances %>%
  dplyr::select(Contig, File, Folder, S1, S2, S1.deviance, S2.deviance) %>%
  tidyr::pivot_longer(cols = c(S1, S2), names_to = "mer1", values_to = "Stomata")%>%
  tidyr::pivot_longer(cols = c(S1.deviance, S2.deviance), names_to = "dev_type", values_to = "Deviance") %>%
  dplyr::filter( (mer1=="S1" & dev_type=="S1.deviance") | (mer1=="S2" & dev_type=="S2.deviance") ) %>%
  dplyr::select(-mer1, -dev_type) %>%
  dplyr::distinct() %>%
  dplyr::group_by(File, Folder, Contig) %>%
  dplyr::summarise( nStomata = n(),
                    SumDeviance = sum(Deviance),
                    SumAbsDeviance = sum(abs(Deviance)),
                    MeanDeviance = SumDeviance/nStomata,
                    RootSumSqareDeviance = sqrt( sum(Deviance^2)),
                    RootMeanSquareDeviance = sqrt( sum(Deviance^2)/nStomata  )
  )

# Overall levels of deviance from straight line in contig
ggplot(deviance.data, aes(x = File, y = MeanDeviance))+
  geom_hline(yintercept = 0)+
  geom_beeswarm()+
  geom_boxplot(width = 0.2, alpha = 0)+
  theme_bw()

# Mean deviance versus sum of squares - consistency of bends
ggplot(deviance.data, aes(x = RootMeanSquareDeviance, y = MeanDeviance, col = nStomata))+
  geom_point()+
  theme_bw()+
  theme(legend.position = "none")

ggsave("Distance_plot.png", plot = dist.plot, dpi = 300, units = "mm", width = 85, height = 85)

# TODO - find endpoints of each contig, measure deviation from straight line
