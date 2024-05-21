# Parse json file to YOLO format. Make YOLO directory structure

# File should have object outlines e.g created by SAM

library(jsonlite)
library(tidyverse)
library(magrittr)
library(fs)

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
                         "imageName" = str_replace(basename(file), ".json", ".jpg"))
    shape.number <<- shape.number+1
    return(result)
  }
  
  
  do.call(rbind, lapply(data$shapes, process.shape)) %>% 
    as.data.frame %>%
    dplyr::mutate("imageWidth" = data$imageWidth,
                  "imageHeight" = data$imageHeight)
}

files <- list.files(path = "stomatal_image", pattern = "*.json",  recursive = T, full.names = T)
# border.data <- get.border(files[1])
border.data <-  do.call(rbind, lapply(files, get.border))

# Annotations for YOLOv8 are in the form of txt files.
# where each object bound is
# class x_center y_center width height.
# Coordinates should be normalised (to 0-1)

bbox.data <- border.data %>%
  dplyr::group_by(shape, file, folder, imageName, imageWidth, imageHeight) %>%
  dplyr::mutate(xmin = min(x)/imageWidth, # Normalise coordinates
                   xmax = max(x)/imageWidth,
                   ymin = min(y)/imageHeight,
                   ymax = max(y)/imageHeight,
                   width = (xmax-xmin),
                   height = (ymax-ymin),
                   xcentre = xmin+(width/2),
                   ycentre = ymin+(height/2),
                   area = width*height,
                   class_id = 0, # One class only for stomata
                   outstring = paste(class_id, xcentre, ycentre, width, height, sep = " ")) %>%
  dplyr::select(-x, -y) %>%
  dplyr::distinct() %>%
  dplyr::filter(area < 0.01) # remove misclicks that cover the entire image


# Split to training and validation image, by image
set.seed(42)
train.images <- caTools::sample.split(unique(bbox.data$imageName), SplitRatio = 0.7)
names(train.images) <- unique(bbox.data$imageName)



out.data <- bbox.data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(isTraining = train.images[[imageName]]) 

# Write the YOLO bounding box format for a given folder
# dir structure should be 
# images/train/a.jpg
# images/val/b.jpg
# labels/train/a.txt
# labels/val/b.txt

fs::file_delete("data.zip")
fs::dir_delete("data")
fs::dir_create("data/images/train/")
fs::dir_create("data/images/val/")
fs::dir_create("data/labels/train/")
fs::dir_create("data/labels/val/")

write.yolo.v8.bounds <- function(source.image.name){
  yolo.data <- out.data %>% dplyr::filter(imageName==source.image.name)
  
  source.folder <- unique(yolo.data$folder)
  is.Training <- unique(yolo.data$isTraining)
  
  source.image <- paste0(source.folder, "/", source.image.name)
  
  out.label.folder <- ifelse(is.Training, "data/labels/train", "data/labels/val" )
  
  out.label.file <- gsub("jpg", "txt", paste0(out.label.folder, "/", paste0(basename(source.folder),"_", source.image.name)))

  # Write the label file for the image
  write.table(yolo.data$outstring, file = out.label.file,
            quote = F, row.names = F, col.names = F, sep = " ")
  
  # Copy the image to the correct folder
  out.image.folder <- gsub("labels", "images", out.label.folder)
  out.image.file <- paste0(out.image.folder, "/", paste0(basename(source.folder),"_", source.image.name))
  if(!file.exists(out.image.file)) file.copy(source.image, out.image.file)
}

sapply(unique(out.data$imageName), write.yolo.v8.bounds)


YOLO.MODEL <- "yolov8n" # try the tiny version of v9 once paper is published: likely "yolov9t"
YOLO.YAML.FILE  <- paste0("stomata.", YOLO.MODEL, ".yaml")
YOLO.PYTHON.FILE  <- paste0("stomata.", YOLO.MODEL, ".py")
YOLO.PREDICT.FILE  <- paste0("stomata.", YOLO.MODEL, "_predict.py")
YOLO.WRAPPER.FILE  <- paste0("stomata.", YOLO.MODEL, ".sh")


# Create YOLO training YAML
yolo.yaml <- paste0(
  "path: /home/bs19022/projects/stomata/data\n",
  "train: images/train\n",
  "val: images/val\n",
  "\n",
  "# class names\n",
  "names:\n",
  "  0: stomata\n"
)

write_file(yolo.yaml, file = YOLO.YAML.FILE)

# Create python file to run training

yolo.python <- paste0(
"from ultralytics import YOLO

model = YOLO(\"", YOLO.MODEL, ".pt\") # load a pretrained model
  
results = model.train(
  data=\"", YOLO.YAML.FILE, "\",
  imgsz=1280,
  epochs=70, 
  batch=8, 
  name=\"", YOLO.MODEL, "_stomata\"
) # train the model
metrics = model.val() # evaluate model performance on the validation set

# Do we want to tune the model hyperparameters?
# result_grid = model.tune(data=\"stomata.", YOLO.MODEL, ".yaml\", epochs=30, iterations=300, optimizer=\"AdamW\", plots=False, save=False, val=False)

path = model.export(format=\"onnx\")  # export the model to ONNX format
  
"
)

write_file(yolo.python, file = YOLO.PYTHON.FILE)

yolo.predict <- paste0(
"import cv2
from ultralytics import YOLO

# Read the pretrained model
model = YOLO(\"runs/detect/", YOLO.MODEL, "_stomata/weights/best.pt\")

# Read a test image and predict stomata locations
im2 = cv2.imread(\"stomatal_image/August\ 21\ -\ 25\ samples/1012-1-1.jpg\")
results = model.predict(source=im2, save=True, save_txt=True) 
"  
)

write_file(yolo.predict, file = YOLO.PREDICT.FILE)

# Create wrapper script to submit the job to a GPU node
yolo.wrapper <- paste0(
"#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q gpu.q
#$ -l gpu=1
#$ -M b.skinner@essex.ac.uk
#$ -m e
#$ -o /home/bs19022/projects/stomata/runs/", YOLO.MODEL, ".log.txt
#$ -e /home/bs19022/projects/stomata/runs/", YOLO.MODEL, ".log.txt

source /usr/local/gpuallocation.sh
# This is only needed once to create the conda env on a GPU node
#conda create -y -n stomata ultralytics pytorch torchvision

source activate stomata
python ", YOLO.PYTHON.FILE, "
conda deactivate
"
)

write_file(yolo.wrapper, file = YOLO.WRAPPER.FILE)
# Zip the data
zip("data.zip", "data")

