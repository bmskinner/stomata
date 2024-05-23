# Create chains from point coordinates

#### Imports #####
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
source("functions.R")

#### Constants #####
# min and maximum distance between stomata
MIN.DISTANCE <- 100
MAX.DISTANCE <- 500

# Max angle difference from 180 degrees
ANGLE.DELTA <- 15

#### Functions #####



# Give a file of stomata locations, create 1mers.
# the input file can be in .json (AnyLabelling) format
# or txt (YOLO bounding box format)
read.1mers <- function(file){
  
  
  read.1mer.from.json <- function(file){

    # Create 1mers from stomata outlines
    create.1mers <- function(border.data){
      data <- border.data %>%
        dplyr::group_by(shape, file) %>%
        dplyr::summarise(x.com = mean(x), 
                         y.com = mean(y), 
                         .groups = 'drop') %>%
        dplyr::mutate(stomata  = paste0("s", sprintf("%02d", shape))) %>%
        dplyr::ungroup() %>%
        dplyr::select(stomata, x.com, y.com, file) %>%
        as.data.frame
      
      # Create polygons from border
      polys <- border.data %>% dplyr::group_by(shape, file) %>%
        dplyr::summarise(poly.matrix = list(matrix(c(x, x[1], y, y[1]), ncol=2, byrow=F)), .groups = 'drop')
      data$polygons <- lapply(polys$poly.matrix, function(x) sf::st_polygon(list(x)))
      data %>% dplyr::select(-file)
    }
    
    border.data <- read.border.from.json(file)
    
    if(nrow(border.data)<=1){
      # return empty dist.contigs dataframe
      cat("  No objects in json file, returning\n")
      return(data.frame())
    }
    
    create.1mers(border.data)
  }
  
  if(str_ends(file, "json")){
    mer1 <- read.1mer.from.json(file)
  } else {
    mer1 <- read.border.from.yolo(file)
  }
  
  mer1
}


# create a deBruijn graph from kmers
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

# create contigs from kmers
create.contigs <- function(mer.data){
  if(nrow(mer.data)==0) return(mer.data %>% dplyr::mutate(Contig = list()))
  
  g <- create.debruijn.graph(mer.data)
  
  # Visualise the graph
  layout <- layout_with_kk(g)
  # plot(g, layout = layout)
  
  # Decompose unlinked contigs
  gphs <- decompose.graph(g)
  
  # Get the number of the contig each 2mer belongs to
  get.contig.number <- function(i){
    # cat("    Getting contig", i, "\n")
    gph <- gphs[[i]]
    contig <- vertex_attr(gph, "kmer")
    contig.num <- rep(i, length(contig))
    names(contig.num) <- contig
    # cat("    Contig", i, "is",contig.num, "\n")
    # cat("    Contig", i, ": ",names(contig.num), "\n")
    # cat("    Contig", i, ": ",str(contig.num), "\n")
    contig.num
  }
  
  # cat("  Graph has", length(gphs), "contigs\n")
  contig.numbers <- do.call(c, lapply(1:length(gphs), get.contig.number))
  
  mer.data %>%
    dplyr::mutate(Contig = map_int(merL, function(x) contig.numbers[names(contig.numbers)==x]))
}

# From a table of centres of mass, create 2mers
# filter to those within a given distance of each other
# and annotate with absolute angles on image
create.2mers <- function(coms, min.distance, max.distance){
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
                    between(S1S2, min.distance, max.distance)) %>%
    dplyr::distinct()
  
  if(nrow(data)>0){
    data$mer2id <- 1:nrow(data)
  }
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
                  merR = paste0(S2, S3)) %>%
    dplyr::distinct()
}

#### Main function #####

# Given a json file, extract the stomata and 
# create linear chains
process.coordinate.file <- function(file){
  cat("Analysing", file, "\n")
  
  # TODO - when this is a single YOLO output file with multiple images
  # Read the image with stomata
  image <- str_replace(file, "json", "jpg")
  img <- OpenImageR::readImage(image)
  img <- OpenImageR::flipImage(img, mode = "vertical") # to draw as expected
  


  
  plot.2mers <- function(mer.data){
    img.grob <- rasterGrob(img, interpolate=TRUE)
    ggplot(mer.data)+
      annotation_custom(img.grob, xmin=0, xmax=dim(img)[2], ymin=0, ymax=dim(img)[1]) +
      coord_fixed(xlim = c(0, dim(img)[2]), ylim = c(0, dim(img)[1]))+
      scale_color_viridis_c()+
      geom_polygon(data = border.data, aes(x = x, y = y, group=shape), fill = "darkgreen", alpha=0.6)+
      geom_point(aes(x = S1.x, y = S1.y))+
      geom_point(aes(x = S2.x, y = S2.y))+
      geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y), col="gray20")+
      geom_text(aes(x = (S2.x+S1.x)/2, y = (S2.y+S1.y)/2, label = sprintf("%.0f", S1S2)), col = "blue", size = 2)+
      geom_text(aes(x = S1.x, y = S1.y, label = S1), col = "pink1", size = 2)+
      geom_text(aes(x = S2.x, y = S2.y, label = S2), col = "pink1", size = 2)+
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(), 
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank()
      )
  }
  
  plot.3mers <- function(mer.data){
    img.grob <- rasterGrob(img, interpolate=TRUE)
    ggplot(mer.data)+
      annotation_custom(img.grob, xmin=0, xmax=dim(img)[2], ymin=0, ymax=dim(img)[1]) +
      coord_fixed(xlim = c(0, dim(img)[2]), ylim = c(0, dim(img)[1]))+
      scale_color_viridis_c()+
      geom_polygon(data = border.data, aes(x = x, y = y, group=shape), fill = "darkgreen", alpha=0.6)+
      geom_point(aes(x = S1.x, y = S1.y))+
      geom_point(aes(x = S2.x, y = S2.y))+
      geom_point(aes(x = S3.x, y = S3.y))+
      geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y), col="gray20")+
      geom_segment(aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y), col="gray20") +
      geom_text(aes(x = S2.x, y = S2.y+50, label = sprintf("%.0f", angle)), col = "red", size = 2)+
      geom_text(aes(x = (S2.x+S1.x)/2, y = (S2.y+S1.y)/2, label = sprintf("%.0f", S1S2)), col = "blue", size = 2)+
      geom_text(aes(x = (S2.x+S3.x)/2, y = (S2.y+S3.y)/2, label = sprintf("%.0f", S2S3)), col = "blue", size = 2)+
      geom_text(aes(x = S1.x, y = S1.y, label = S1), col = "pink1", size = 2)+
      geom_text(aes(x = S2.x, y = S2.y, label = S2), col = "pink1", size = 2)+
      geom_text(aes(x = S3.x, y = S3.y, label = S3), col = "pink1", size = 2)+
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(), 
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank()
      )
  }
  
  # plot contig data
  plot.contigs <- function(mer.contigs){
    img.grob <- rasterGrob(img, interpolate=TRUE)
    # visualise the contigs
    ggplot(mer.contigs)+
      annotation_custom(img.grob, xmin=0, xmax=dim(img)[2], ymin=0, ymax=dim(img)[1]) +
      coord_fixed(xlim = c(0, dim(img)[2]), ylim = c(0, dim(img)[1]))+
      geom_polygon(data = border.data, aes(x = x, y = y, group=shape), fill = "darkgreen", alpha=0.6)+
      geom_point(aes(x = S1.x, y = S1.y, col = as.factor(Contig)), size=2)+
      geom_point(aes(x = S2.x, y = S2.y, col = as.factor(Contig)), size=2)+
      geom_point(aes(x = S3.x, y = S3.y, col = as.factor(Contig)), size=2)+
      geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y, col = as.factor(Contig)), linewidth=1) +
      geom_segment(aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y, col = as.factor(Contig)), linewidth=1) +
      labs(col = "Chain")+
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(), 
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank()
      )
  }
  
  save.plot <- function(plot, filename){
    ggsave(str_replace(file, ".json", filename), plot = plot, dpi = 300, units = "mm", width = 170, height = 140)
    plot
  }
  
  cat("  Creating 1mers\n")
  mer1 <- read.1mers(file)
  
  # Create distance table between pairs of stomata
  # Filter to only those within a given distance
  cat("  Creating 2mers\n")
  mer2 <- create.2mers(mer1, min.distance = MIN.DISTANCE, max.distance = MAX.DISTANCE)
  
  # What happens if no/very few stomata are present in the image that cannot be 
  # contiged? Check there are rows present, and skip if not.
  if(nrow(mer2)<=1){
    # return empty dist.contigs dataframe
    cat("  No 2mers, returning\n")
    return(data.frame())
  }
  
  mer.2.plot <- plot.2mers(mer2)
  save.plot(mer.2.plot, ".mer2.raw.png")
  
  
  ###################
  
  cat("  Removing intersecting 2mers\n")
  # Remove edges in the graph that intersect a third stomata
  # can use sf: https://stackoverflow.com/questions/61703791/determine-lines-that-intersect-a-polygon-in-r
  # Create line objects
  mer2$lines <- lapply(1:nrow(mer2), function(i) sf::st_linestring(rbind(c(mer2$S1.x[i], mer2$S1.y[i]),
                                                                         c(mer2$S2.x[i], mer2$S2.y[i]))) )
  
  # Check each 2mer for intersections with a stomata
  mer2$intersects <- lapply(mer2$lines, function(l) sum(unlist(sapply(mer1$polygons, sf::st_intersects, y=l))))
  
  mer2.filt <- mer2 %>%
    dplyr::filter(intersects <= 2)
  
  if(nrow(mer2.filt)<=1){
    # return empty dist.contigs dataframe
    cat("  No 2mers after filtering, returning\n")
    return(data.frame())
  }
  
  # mer.2.drop.plot <- plot.2mers(mer2[!mer2.keep,], img, border.data)
  # ggsave(str_replace(file, ".json", ".mer2.drop.png"), plot = mer.2.drop.plot, dpi = 300, units = "mm", width = 170, height = 140)
  # 
  mer.2.filt.plot <- plot.2mers(mer2.filt)
  save.plot(mer.2.filt.plot, ".mer2.filt.png")
  
  # TODO: remove the 2mers which are at a different angle to the stomata
  # orientation (use the angle measurements from mer1)
  
  # Join the tables to create a 3-mer chain
  mer3 <- create.3mers(mer2.filt)
  # mer.3.plot <- save.plot(plot.3mers(mer3), ".mer3.raw.png")
  
  if(nrow(mer3)<=1){
    # return empty dist.contigs dataframe
    cat("  No 3mers, returning\n")
    return(data.frame())
  }
  
  # Find the 3mers in straight lines
  mer3.straight <- mer3 %>% dplyr::filter( angle > 180 - ANGLE.DELTA)
  mer.3.straight.plot <- save.plot(plot.3mers(mer3.straight), ".mer3.straight.png")
  
  # cat("Filtering on angle\n")
  
  # Filter the straight 3mers to the most common orientation in the image
  # TODO - this does not work when the 3mers are vertical. 
  # Possible fix - find the mode on mer2 raw instead; seems more consistent
  # Still not perfect - may need to look only at local maxima
  modal.angle <- find.mode(mer2$abs.angle)
  # modal.angle <- find.mode(mer3.straight$abs.angle)
  mer3.filt <- mer3.straight %>% dplyr::filter(between(abs.angle, modal.angle - ANGLE.DELTA, modal.angle+ANGLE.DELTA))
  # mer.3.filt.plot <- save.plot(plot.3mers(mer3.filt), ".mer3.filtered.png")
  
  if(nrow(mer3.filt)<=1){
    # return empty dist.contigs dataframe
    cat("  No 3mers after filtering, returning\n")
    return(data.frame())
  }
  
  
  # cat("Pruning 3mers\n")
  
  # Pruning step
  # Look for 2mers present more than once
  # Drop the 3mers with the lowest angle.
  # TODO
  mer3.pruned <- mer3.filt %>% 
    dplyr::group_by(merL) %>%
    dplyr::arrange(merL, desc(angle)) %>%
    dplyr::slice_head(n=1) %>%
    dplyr::group_by(merR) %>%
    dplyr::arrange(merR, desc(angle)) %>%
    dplyr::slice_head(n=1)
  
  # Pruning step - remove any kmers where the endpoint is in the middle of a chain
  # to prevent branches
  mer3.longer <- mer3.pruned %>%
    pivot_longer(c(S1, S2, S3), names_to = "StomataPosition", values_to = "Stomata") %>%
    dplyr::group_by(StomataPosition, Stomata) %>%
    dplyr::arrange(StomataPosition, Stomata, desc(angle)) %>%
    # dplyr::mutate(Count = n()) %>%
    dplyr::slice_head(n=1) %>% # take only the 3mer with the highest angle
    tidyr::pivot_wider(names_from = StomataPosition, values_from = Stomata) %>%
    na.omit # remove the rows with NAs due to our slice
  
  
  # What happens if no/very few stomata are present in the image that cannot be 
  # contiged? Check there are rows present, and skip if not.
  if(nrow(mer3.longer)<=1){
    # return empty dist.contigs dataframe
    cat("No rows in filtered 3mers, returning\n")
    return(data.frame())
  }
  
  # cat("Removing terminal branches\n")
  # This still leaves branches when the endpoints meet a terminal kmer.
  # Prune again - cases where the branch meets a terminal 2mer
  # Look for 1mers that are are the S2 of 2 different 2mers
  mer3.terminal <- mer3.longer %>%
    tidyr::pivot_longer(c(merL, merR), names_to = "mer2Type", values_to = "mer2") %>%
    dplyr::mutate(SL = ifelse(mer2Type == "merL", S1, S2),
                  SR = ifelse(mer2Type == "merL", S2, S3)) %>%
    dplyr::group_by(SL) %>%
    dplyr::mutate(splits = paste(unique(SR), collapse = ""),
                  nsplits = nchar(splits)) %>% # how many letters in combined mers
    dplyr::group_by(splits) %>%
    dplyr::arrange(splits, desc(angle)) %>% # order so the one to keep is top
    dplyr::mutate(rownum = row_number()) %>% # scruffy; figure out which row in group
    dplyr::filter(nsplits<=3 | rownum ==1) %>% # and keep the first from multi2mers
    dplyr::ungroup() %>%
    dplyr::select(-c(SL, SR, splits, nsplits, rownum)) %>%
    tidyr::pivot_wider(names_from = "mer2Type", values_from = "mer2") %>%
    na.omit
  
  # Now do the same for when the branch is at the start, not the end
  
  mer3.terminal <- mer3.terminal %>%
    tidyr::pivot_longer(c(merL, merR), names_to = "mer2Type", values_to = "mer2") %>%
    dplyr::mutate(SL = ifelse(mer2Type == "merL", S1, S2),
                  SR = ifelse(mer2Type == "merL", S2, S3)) %>%
    dplyr::group_by(SR) %>%
    dplyr::mutate(splits = paste(unique(SL), collapse = ""),
                  nsplits = nchar(splits)) %>% # how many letters in combined mers
    dplyr::group_by(splits) %>%
    dplyr::arrange(splits, desc(angle)) %>% # order so the one to keep is top
    dplyr::mutate(rownum = row_number()) %>% # scruffy; figure out which row in group
    dplyr::filter(nsplits<=3 | rownum ==1) %>% # and keep the first from multi2mers
    dplyr::ungroup() %>%
    dplyr::select(-c(SL, SR, splits, nsplits, rownum)) %>%
    tidyr::pivot_wider(names_from = "mer2Type", values_from = "mer2") %>%
    na.omit
  
  if(nrow(mer3.terminal)<=1){
    # return empty dist.contigs dataframe
    cat("  No 3mers after pruning, returning\n")
    return(data.frame())
  }
  
  
  # mer.3.pruned.plot <- save.plot(plot.3mers(mer3.terminal), ".mer3.pruned.png")
  
  # Check the distribution of angles in mer3 vs mer2
  density.plot <- ggplot()+
    geom_density(data = mer2, aes(x = abs.angle, col="mer2 raw" ), alpha=0) +
    geom_density(data = mer3, aes(x = abs.angle, col="mer3 raw" ), alpha=0) +
    geom_density(data = mer3.straight, aes(x = abs.angle, col="mer3 straight")) +
    geom_density(data = mer3.terminal, aes(x = abs.angle, col="mer3 longer" )) +
    labs(x = "Angle of kmer to vertical") +
    scale_color_manual(values = c("black", "red", "green", "blue"))+
    theme_bw()+
    theme(legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          legend.background = element_blank())
  ggsave(str_replace(file, ".json", "mer2.mer3.angle_density.png"), plot = density.plot, dpi = 300, units = "mm", width = 100, height = 85)
  
  mer3.short.plot <- save.plot(plot.3mers(mer3.terminal), ".mer3.terminal.png")
  
  short.contigs <- create.contigs(mer3.terminal) %>% 
    dplyr::group_by(Contig)
  
  chain.plot <- save.plot(plot.contigs(short.contigs), ".chains.png")
  
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
  
  # How many stomata are not in a chain? TODO - check
  stomata.in.contgs <- unique(c(short.contigs$S1, short.contigs$S2, short.contigs$S3))
  remainder.stomata <- mer1$stomata[ !(mer1$stomata %in% stomata.in.contgs)]
  
  dist.contigs$nStomata <- length(unique(mer1$stomata))
  dist.contigs$unassignedStomata <- length(remainder.stomata)
  dist.contigs$fUnassignedStomata <- length(remainder.stomata)/length(unique(mer1$stomata))
  
  rotated.plot <- ggplot(dist.contigs)+
    geom_segment(aes(x = Contig.start.x, y = Contig.start.y, xend = RotatedContigEnd[,1], yend = RotatedContigEnd[,2]), col = "grey20")+
    # geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y), col = "grey")+
    geom_segment(aes(x = RotatedS1[,1], y = RotatedS1[,2], xend = RotatedS2[,1], yend = RotatedS2[,2]), col = "blue")+
    geom_point(aes(x = RotatedS1[,1], y = RotatedS1[,2]), col="blue", size=2)+
    geom_point(aes(x = RotatedS2[,1], y = RotatedS2[,2]), col="blue", size=2)+
    theme_bw()+
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank())
  
  ggsave(str_replace(file, ".json", ".chain.rotated.png"), plot = rotated.plot, dpi = 300, units = "mm", width = 170, height = 140)
  return(dist.contigs)
}
