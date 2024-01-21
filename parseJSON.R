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

MAX.DISTANCE <- 300
TRIM.DISTANCE <- 500
SHORT.DISTANCE <- 300
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
    geom_text(aes(x = S2.x, y = S2.y+50, label = sprintf("%.0f", angle)))+
    geom_text(aes(x = (S2.x+S1.x)/2, y = (S2.y+S1.y)/2, label = sprintf("%.0f", S1S2)))+
    geom_text(aes(x = (S2.x+S3.x)/2, y = (S2.y+S3.y)/2, label = sprintf("%.0f", S2S3)))
}

create.debruijn.graph <- function(mer.data){
  uks <- unique(c(mer.data$merL, mer.data$merR)) # unique 2mers
  if(length(uks)==0)
    return(make_empty_graph())
  
  # Assign an id to each kmer
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
    coord_fixed()+
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

# Given a json file, extract the stomata and 
# create linear chains
process.json.file <- function(file){
  cat("Analysing", file, "\n")
  # Read the image with stomata
  image <- str_replace(file, "json", "jpg")
  img <- readImage(image)
  img <- flipImage (img, mode = "vertical")
  
  border.data <- get.border(file)
  
  coms <- calc.coms(border.data)

  # Create distance table between pairs of stomata
  # Filter to only those within MAX.DISTANCE
  cat("  Calculating distances\n")
  dist.table <- expand.grid(coms$stomata, coms$stomata, stringsAsFactors = F) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(S1t = min(Var1, Var2), S2t = max(Var1, Var2)) %>% # ensure alphabetical order
    dplyr::ungroup() %>%
    dplyr::select(Var1 = S1t, Var2 = S2t) %>% # remove unused columns
    dplyr::distinct() %>% # remove duplicates
    dplyr::filter(Var1!=Var2) %>%
    merge(., coms[,1:3], by.x = "Var1", by.y = "stomata", all.y = F) %>% # add coordinates for S1
    dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>% # rename for clarity
    merge(., coms[,1:3], by.x = "S2", by.y = "stomata", all.y = F) %>% # add coordinates for S1
    dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) %>% # rename for clarity
    dplyr::rowwise() %>%
    dplyr::mutate(S1S2 = euclidean(S1.x, S1.y, S2.x, S2.y),
                  kmer = paste0(S1, S2),
                  abs.angle = LearnGeom::Angle(c(S1.x, S1.y+10), c(S1.x,  S1.y),  c(S2.x,S2.y))) %>% # angle of 2mer relative to image
    dplyr::filter(S1 != S2, S1S2 < TRIM.DISTANCE)
  
  # Join the tables to create a 3-mer chain
  mer3 <- merge(dist.table, dist.table, by.x = c("S2","S2.x", "S2.y" ), by.y = c("S1","S1.x", "S1.y" ), 
                all.y = F, suffixes = c("A", "B")) %>%
    dplyr::select(S1, S1.x, S1.y, S2, S2.x, S2.y, S3 = S2B, S3.x = S2.xB, 
                  S3.y = S2.yB, S1S2 = S1S2A, S2S3 = S1S2B) %>% # rename for clarity
    dplyr::rowwise() %>%
    dplyr::mutate(angle = LearnGeom::Angle(c(S1.x,  S1.y), c(S2.x, S2.y), c(S3.x,S3.y)),
                  abs.angle = LearnGeom::Angle(c(S1.x, S1.y+10), c(S1.x,  S1.y),  c(S2.x,S2.y))) %>% 
    dplyr::filter( angle > 180 - ANGLE.DELTA)  %>% # rename for clarity
    dplyr::mutate(merL = paste0(S1, S2),
                  merR = paste0(S2, S3))
  
  # Check the distrinbution of 2mer angles in mer3 vs mer2
  ggplot(mer3, aes(x = abs.angle))+
    geom_density()
  
  # Crete density and find max value
  find.mode <-  function(x) {
    d <- density(x)
    d$x[which.max(d$y)]
  }
  
  # Filter the 2mers to those in a reasonable range
  modal.angle <- find.mode(mer3$abs.angle)
  mer3 %<>% dplyr::filter(abs.angle >= modal.angle - 5 & abs.angle <= modal.angle+5)
  
  
  mer.3.plot <- plot.3mers(mer3, img, border.data)
  #ggsave(str_replace(file, ".json", ".mer3.png"), plot = mer.3.plot, dpi = 300, units = "mm", width = 170, height = 140)
  
  # Shows we can't make the graph in one step - the distance values between chains and non-chains can overlap
  # What if we link the close kmers, then widen up and elongate chains with more distant kmers?
  
  short.3mers <- mer3 %>% dplyr::filter(S1S2 < SHORT.DISTANCE & S2S3 < SHORT.DISTANCE)
  
  short.contigs <- create.contigs(short.3mers) %>% 
    dplyr::group_by(Contig) %>%
    dplyr::filter(n()>1)
  
  chain.plot <- plot.contigs(short.contigs, img, border.data)
  
  # Match contigs to the 2mers they contain
  contig.assignment <- short.contigs %>% dplyr::select(Contig, merL, merR) %>%
    tidyr::pivot_longer(c(merL, merR), values_to = "kmer") %>%
    dplyr::select(Contig, kmer) %>%
    dplyr::distinct() %>%
    merge(., dist.table, by= "kmer")
  
  # Calculate average distances per contig
  dist.contigs <- contig.assignment %>% 
    dplyr::mutate(File  = file)
  
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

ggsave("Distance_plot.png", plot = dist.plot, dpi = 300, units = "mm", width = 85, height = 85)

# TODO - find endpoints of each contig, measure deviation from straight line
