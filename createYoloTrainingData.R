# Parse json file to YOLO format. Make YOLO directory structure

# File should have object outlines e.g created by SAM

library(jsonlite)
library(tidyverse)
library(magrittr)
library(fs)
source("functions.R")

#### Detect CoM in json outline data ####

files <- list.files(path = "stomatal_image", pattern = "*.json",  recursive = T, full.names = T)

border.data <-  do.call(rbind, lapply(files, read.border.from.json))

#### Split images to training and validation ####

# Split to training and validation image, by image
set.seed(42)
train.images <- caTools::sample.split(unique(border.data$imageName), SplitRatio = 0.7)
names(train.images) <- unique(border.data$imageName)

# Write the YOLO bounding box format for a given folder
# dir structure should be 
# images/train/a.jpg
# images/val/b.jpg
# labels/train/a.txt
# labels/val/b.txt

fs::file_delete("data.zip")
fs::dir_delete("data")
fs::dir_create("data/bbox/images/train/")
fs::dir_create("data/bbox/images/val/")
fs::dir_create("data/bbox/labels/train/")
fs::dir_create("data/bbox/labels/val/")
fs::dir_create("data/seg/images/train/")
fs::dir_create("data/seg/images/val/")
fs::dir_create("data/seg/labels/train/")
fs::dir_create("data/seg/labels/val/")

write.yolo.v8.bbox <- function(source.image.name){
  
  # Annotations for YOLOv8 are in the form of txt files.
  # where each object bound is
  # class x_center y_center width height.
  # Coordinates should be normalised (to 0-1)
  
  bbox.data <- border.data %>%
    dplyr::filter(imageName==source.image.name) %>%
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

  out.data <- bbox.data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(isTraining = train.images[[imageName]]) 
  
  yolo.data <- out.data %>% dplyr::filter(imageName==source.image.name)
  
  source.folder <- unique(yolo.data$folder)
  is.Training <- unique(yolo.data$isTraining)
  
  source.image <- paste0(source.folder, "/", source.image.name)
  
  out.label.folder <- ifelse(is.Training, "data/bbox/labels/train", "data/bbox/labels/val" )
  
  out.label.file <- gsub("jpg", "txt", paste0(out.label.folder, "/", paste0(basename(source.folder),"_", source.image.name)))

  # Write the label file for the image
  write.table(yolo.data$outstring, file = out.label.file,
            quote = F, row.names = F, col.names = F, sep = " ")
  
  # Copy the image to the correct folder
  out.image.folder <- gsub("labels", "images", out.label.folder)
  out.image.file <- paste0(out.image.folder, "/", paste0(basename(source.folder),"_", source.image.name))
  if(!file.exists(out.image.file)) file.copy(source.image, out.image.file)
}

sapply(unique(border.data$imageName), write.yolo.v8.bbox)


write.yolo.v8.seg <- function(source.image.name){
  
  # Get the objects for the image, concatenaate normalised xy pairs
  yolo.data <- border.data %>%
    dplyr::filter(imageName==source.image.name) %>%
    dplyr::group_by(shape, file, folder, imageName) %>%
    dplyr::reframe(x = x/imageWidth, # Normalise coordinates
                  y = y/imageHeight,
                  xmin = min(x),
                  xmax = max(x),
                  ymin = min(y),
                  ymax = max(y),
                  width = (xmax-xmin),
                  height = (ymax-ymin),
                  area = width*height,
                  class_id = 0) %>% # One class only for stomata
    dplyr::filter(area < 0.01) %>% # remove misclicks that cover the entire image
    dplyr::select(shape, file, folder, imageName, x, y, class_id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(shape, file, folder, imageName) %>%
    dplyr::reframe(outstring = paste(class_id, paste(x, y, collapse = " "), sep = " ")) %>%
    dplyr::distinct()

  source.folder <- unique(yolo.data$folder)
  is.Training <- train.images[[source.image.name]]
  
  source.image <- paste0(source.folder, "/", source.image.name)
  
  out.label.folder <- ifelse(is.Training, "data/seg/labels/train", "data/seg/labels/val" )
  
  out.label.file <- gsub("jpg", "txt", paste0(out.label.folder, "/", paste0(basename(source.folder),"_", source.image.name)))
  
  # Write the label file for the image
  write.table(yolo.data$outstring, file = out.label.file,
              quote = F, row.names = F, col.names = F, sep = " ")
  
  # Copy the image to the correct folder
  out.image.folder <- gsub("labels", "images", out.label.folder)
  out.image.file <- paste0(out.image.folder, "/", paste0(basename(source.folder),"_", source.image.name))
  if(!file.exists(out.image.file)) file.copy(source.image, out.image.file)
}

sapply(unique(border.data$imageName), write.yolo.v8.seg)


#### Create scripts to train models and predict ####

YOLO.BBOX.MODEL            <- "yolov8n" # try the nano version of v10 once available "yolov10n.pt"
YOLO.BBOX.YAML.FILE        <- paste0("stomata.", YOLO.BBOX.MODEL, ".yaml")
YOLO.TRAIN.BBOX.FILE       <- paste0("stomata.", YOLO.BBOX.MODEL, ".bbox.train.py")
YOLO.PREDICT.BBOX.FILE     <- paste0("stomata.", YOLO.BBOX.MODEL, ".bbox.predict.py")
YOLO.TRAIN.BBOX.WRAPPER    <- paste0("stomata.", YOLO.BBOX.MODEL, ".bbox.train.sh")
YOLO.PREDICT.BBOX.WRAPPER  <- paste0("stomata.", YOLO.BBOX.MODEL, ".bbox.predict.sh")

YOLO.SEG.MODEL            <- "yolov8n-seg"
YOLO.SEG.YAML.FILE        <- paste0("stomata.", YOLO.SEG.MODEL, ".yaml")
YOLO.TRAIN.SEG.FILE       <- paste0("stomata.", YOLO.SEG.MODEL, ".seg.train.py")
YOLO.PREDICT.SEG.FILE     <- paste0("stomata.", YOLO.SEG.MODEL, ".seg.predict.py")
YOLO.TRAIN.SEG.WRAPPER    <- paste0("stomata.", YOLO.SEG.MODEL, ".seg.train.sh")
YOLO.PREDICT.SEG.WRAPPER  <- paste0("stomata.", YOLO.SEG.MODEL, ".seg.predict.sh")

OUTPUT.DIR                <- "runs/segment/single_file"

# Create YOLO bbox training YAML
write_file(paste0(
  "path: /home/bs19022/projects/stomata/data/bbox\n",
  "train: images/train\n",
  "val: images/val\n",
  "\n",
  "# class names\n",
  "names:\n",
  "  0: stomata\n"
), file = YOLO.BBOX.YAML.FILE)

# Create YOLO seg training YAML
write_file(paste0(
  "path: /home/bs19022/projects/stomata/data/seg\n",
  "train: images/train\n",
  "val: images/val\n",
  "\n",
  "# class names\n",
  "names:\n",
  "  0: stomata\n"
), file = YOLO.SEG.YAML.FILE)

# Create python file to run bbox training
write_file( paste0(
  "from ultralytics import YOLO

model = YOLO(\"", YOLO.BBOX.MODEL, ".pt\") # load a pretrained model
  
results = model.train(
  data=\"", YOLO.BBOX.YAML.FILE, "\",
  imgsz=2560, # train at close to full size
  epochs=100, 
  batch=2, 
  name=\"", YOLO.BBOX.MODEL, "_stomata\"
) # train the model
metrics = model.val() # evaluate model performance on the validation set

# Do we want to tune the model hyperparameters?
# result_grid = model.tune(data=\"stomata.", YOLO.BBOX.MODEL, ".yaml\", epochs=30, iterations=300, optimizer=\"AdamW\", plots=False, save=False, val=False)

path = model.export(format=\"onnx\")  # export the model to ONNX format
  
"
), file = YOLO.TRAIN.BBOX.FILE)

# Create python file to run seg training
write_file( paste0(
  "from ultralytics import YOLO

model = YOLO(\"", YOLO.SEG.MODEL, ".pt\") # load a pretrained model
  
results = model.train(
  data=\"", YOLO.SEG.YAML.FILE, "\",
  imgsz=1280, # train at lower resolution
  epochs=100, 
  batch=1, # seg training uses higher memory than bbox 
  name=\"", YOLO.SEG.MODEL, "_stomata\"
) # train the model
metrics = model.val() # evaluate model performance on the validation set

path = model.export(format=\"onnx\")  # export the model to ONNX format
  
"
), file = YOLO.TRAIN.SEG.FILE)


# Script to predict bbox based on val data
write_file(paste0(
  "import cv2
from ultralytics import YOLO

# Read the pretrained model
model = YOLO(\"runs/detect/", YOLO.BBOX.MODEL, "_stomata/weights/best.pt\")

# Test on the original validation group
source = \"data/bbox/images/val/*.jpg\"

# Run inference
results = model(source, stream=True, conf=0.05, imgsz=2560)  # generator of Results objects

# Process results generator
with open(\"runs/detect/predict_bbox/output.txt\", 'a') as f:
  print(\"Image\\tx\\ty\\tw\\th\\tconf\", file=f)
  for result in results:
      boxes = result.boxes  # Boxes object for bounding box outputs

      out_path = result.path.replace(\"data/bbox/images/val\", \"runs/detect/predict_bbox\")
      result.save(filename=out_path, labels=False)  # save annotated image to disk

      for box in boxes:
        xywh = box.xywh.tolist()[0]
        print(result.path, xywh[0], xywh[1], xywh[2], xywh[3], box.conf.item(), sep=\"\\t\", file=f)
        
f.close()
"  
), file = YOLO.PREDICT.BBOX.FILE)

# Script to predict seg based on val data
write_file(paste0(
  "import cv2
import numpy as np
import os
from ultralytics import YOLO

# Read the pretrained model
model = YOLO(\"runs/segment/", YOLO.SEG.MODEL, "_stomata/weights/best.pt\")

# Run on the complete image set
source = \"stomatal_image/*/*.jpg\"

# Run inference
results = model(source, stream=True, conf=0.05, imgsz=1280)  # generator of Results objects

# Process results generator
for result in results:

  # Output file paths
  out_dir = os.path.basename(os.path.dirname(result.path))
  out_file = os.path.basename(result.path)
  
  out_base = os.path.join(\"",OUTPUT.DIR,"\",out_dir)
  os.makedirs(out_base, exist_ok=True)
  
  mask_image_path = os.path.join(out_base, out_file + \".png\")
  mask_coord_path = os.path.join(out_base, out_file + \".txt\")
  
  mask_data = result.cpu().masks.data
  img = (mask_data[0].numpy() * 255).astype(\"uint8\")

  height,width = img.shape
  
  # Combine all masks to one image
  masked = np.zeros((height, width), dtype=\"uint8\")
  num_masks = len(mask_data)
  for i in range(num_masks):
    masked = cv2.add(masked, cv2.min((mask_data[i].numpy() * 255).astype(\"uint8\"), 255))
  cv2.imwrite(mask_image_path, masked)

  # Write the marks coordinates
  with open(mask_coord_path, 'a') as f:
    print(\"Image\\tObject\\tx\\ty\", file=f)
    i=0 # track which object is which in output file
    for mask in result.masks:
      xy = mask.xy[0]
      for c in xy:
        print(result.path, str(i), \"\\t\".join( map(str, c) ), sep=\"\\t\", file=f)
      i+=1
    f.close()
"  
), file = YOLO.PREDICT.SEG.FILE)

# Create wrapper script to submit bbox train to a GPU node
write_file(paste0(
  "#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q gpu.q
#$ -l gpu=1
#$ -M b.skinner@essex.ac.uk
#$ -m e
#$ -o /home/bs19022/projects/stomata/runs/", YOLO.BBOX.MODEL, ".train.txt
#$ -e /home/bs19022/projects/stomata/runs/", YOLO.BBOX.MODEL, ".train.txt

source /usr/local/gpuallocation.sh
# This is only needed once to create the conda env on a GPU node
#conda create -y -n stomata ultralytics pytorch torchvision

source activate stomata
python ", YOLO.TRAIN.BBOX.FILE, "
conda deactivate
"
), file = YOLO.TRAIN.BBOX.WRAPPER)

# Create wrapper script to seg train to a GPU node
write_file(paste0(
  "#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q gpu.q
#$ -l gpu=1
#$ -M b.skinner@essex.ac.uk
#$ -m e
#$ -o /home/bs19022/projects/stomata/runs/", YOLO.SEG.MODEL, ".train.txt
#$ -e /home/bs19022/projects/stomata/runs/", YOLO.SEG.MODEL, ".train.txt

source /usr/local/gpuallocation.sh
source activate stomata
python ", YOLO.TRAIN.SEG.FILE, "
conda deactivate
"
), file = YOLO.TRAIN.SEG.WRAPPER)

# Create wrapper script to submit bbox predict to a GPU node
write_file(paste0(
  "#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q gpu.q
#$ -l gpu=1
#$ -M b.skinner@essex.ac.uk
#$ -m e
#$ -o /home/bs19022/projects/stomata/runs/", YOLO.BBOX.MODEL, ".predict.txt
#$ -e /home/bs19022/projects/stomata/runs/", YOLO.BBOX.MODEL, ".predict.txt

source /usr/local/gpuallocation.sh
source activate stomata
mkdir -p runs/detect/predict_bbox
python ", YOLO.PREDICT.BBOX.FILE, "
conda deactivate
"
), file = YOLO.PREDICT.BBOX.WRAPPER)

write_file(paste0(
  "#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q gpu.q
#$ -l gpu=1
#$ -M b.skinner@essex.ac.uk
#$ -m e
#$ -o /home/bs19022/projects/stomata/runs/", YOLO.SEG.MODEL, ".predict.txt
#$ -e /home/bs19022/projects/stomata/runs/", YOLO.SEG.MODEL, ".predict.txt

source /usr/local/gpuallocation.sh
source activate stomata
mkdir -p ", OUTPUT.DIR, "
python ", YOLO.PREDICT.SEG.FILE, "
conda deactivate
"
), file = YOLO.PREDICT.SEG.WRAPPER)

# Zip the data
zip("data.zip", c("data", 
                  YOLO.BBOX.MODEL, YOLO.BBOX.YAML.FILE, 
                  YOLO.TRAIN.BBOX.FILE, YOLO.PREDICT.BBOX.FILE, 
                  YOLO.TRAIN.BBOX.WRAPPER, YOLO.PREDICT.BBOX.WRAPPER,
                  YOLO.SEG.MODEL, YOLO.SEG.YAML.FILE, 
                  YOLO.TRAIN.SEG.FILE, YOLO.PREDICT.SEG.FILE,
                  YOLO.TRAIN.SEG.WRAPPER, YOLO.PREDICT.SEG.WRAPPER))

