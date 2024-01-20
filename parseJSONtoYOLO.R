# Parse json file to YOLO format

# File should have object outlines e.g created by SAM

library(jsonlite)
library(tidyverse)
library(magrittr)

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
                         "file" = file,
                         "folder" = dirname(file),
                         "imagename" = str_replace(basename(file), ".json", ".jpg"))
    shape.number <<- shape.number+1
    return(result)
  }
  
  
  do.call(rbind, lapply(data$shapes, process.shape)) %>% 
    as.data.frame 
}

files <- list.files(path = "stomatal_images", pattern = "*.json",  recursive = T, full.names = T)
border.data <- get.border(files[1])

# Annotations for YOLO are in the form of txt files. Each line in a txt file fol YOLO must have the following format:
# image2.jpg 100,94,613,814,0 31,420,220,540,1
# xmin,ymin,xmax,ymax,class_id
class.id <- 0 

bbox.data <- border.data %>%
  dplyr::group_by(shape, folder, imagename) %>%
  dplyr::summarise(xmin = min(x),
                   xmax = max(x),
                   ymin = min(y),
                   ymax = max(y),
                   class_id = 0,
                   outstring = paste(xmin, ymin, xmax, ymax, class_id, sep = ","))


out.data <- bbox.data %>%
  dplyr::group_by(imagename, folder) %>%
  dplyr::summarise(bboxes = paste(outstring, collapse = " "))

write.txt <- function(folder, image, bboxes){
  write.table(paste(image, bboxes), file = paste0(folder, "/", str_replace(image, ".jpg", ".txt")), 
            quote = F, row.names = F, col.names = F, sep = " ")
}

mapply(write.txt, out.data$folder, out.data$imagename, out.data$bboxes)


