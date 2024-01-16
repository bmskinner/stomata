# Parse the JSON from Anylabelling to extract stomata CoMs

library(jsonlite)
library(tidyverse)

MAX.DISTANCE <- 300

# Calculate centre of mass from JSON
calc.coms <- function(file){
  
  data <- jsonlite::read_json(file)
  
  process.shape <- function(shape){
    bounds <- split(unlist(shape$points, recursive = T), 1:2) %>% as.data.frame
    colnames(bounds) <- c("X", "Y")
    
    data.frame("x.com" = mean(bounds$X),
         "y.com" = mean(bounds$Y))
  }
  
  
  do.call(rbind, lapply(data$shapes, process.shape)) %>% 
    as.data.frame %>%
    dplyr::mutate(stomata  = paste0(file, "_s", row_number()),
                  file = file)
  
}

# Read all json files
files <- list.files(path = "stomatal_image", pattern = "*.json",  recursive = T, full.names = T)

# Find CoMs of all objects in one file for now
coms <- do.call(rbind, lapply(files[2], calc.coms)) %>% as.data.frame


# Now add in the chain code

# Create the distance matrix for an image
euclidean <- function(x1, y1, x2, y2) sqrt( (x1-x2)^2 + (y1-y2)^2)

# Create distance table between pairs of stomata
# Filter to only those within MAX.DISTANCE
dist.table <- expand.grid(coms$stomata, coms$stomata) %>%
  merge(., coms[,1:3], by.x = "Var1", by.y = "stomata", all.y = F) %>%
  dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>%
  merge(., coms[,1:3], by.x = "S2", by.y = "stomata", all.y = F) %>%
  dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) %>%
  dplyr::mutate(S1S2 = euclidean(S1.x, S1.y, S2.x, S2.y)) %>%
  dplyr::filter(S1 != S2, S1S2 <= MAX.DISTANCE)

# Join the tables to create a 3-mer chain
kmer <- merge(dist.table, dist.table, by.x = "S2", by.y = "S1", all.y = F)
# TODO
