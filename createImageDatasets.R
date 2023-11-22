# Create test and training folders from a mixed set of images.
# Names of images and annotation files will be changed to allow
# multiple folders to be combined into one training/test set of data

# Note that each image should have an accompanying json file containing the 
# keypoint location of the stomata centre:

# {"keypoints":[[[1217, 315, 1]], 
# [[762, 567, 1]], 
# [[524, 132, 1]], 
# [[601, 191, 1]], 
# [[407, 303, 1]], 
# [[1131, 669, 1]], 
# [[123, 534, 1]], 
# [[744, 605, 1]], 
# [[278, 175, 1]]
# ]} 

library(fs)

# Create a json keypoints file from ImageJ manually saved coordinates
make.json <- function(in.file){
  data <- read.csv(in.file)
  data$X <- round(data$X) # correct half pixels
  data$Y <- round(data$Y)

  outstring.start <- "{\"keypoints\":["
  outstring.end <- "]}"
  outstring.mid <- paste(paste0("[[", data$X, ",", data$Y, ",1]]"), collapse=",")

  out.file <- gsub("csv", "json", in.file)
  if(file.exists(out.file)) fs::file_delete(out.file)
  
  file.conn <- file(out.file)
  writeLines( paste0(outstring.start, outstring.mid, outstring.end), file.conn)
  close(file.conn)
}

# Create a new output folder based on the input folder name
# The input folder will contain all the annotated images, possibly in subfolders
in.folder = commandArgs(trailingOnly=T)[1]
if(is.na(in.folder)) stop("Must provide an input directory")

# Delete any previous training folders
out.folder = paste0(in.folder, "_trainable")
if(dir.exists(out.folder)) fs::file_delete(out.folder)
dir.create(out.folder)

# Randomly order the image files and get the corresponding annotation file names
files = sample(list.files(in.folder, full.names = T, pattern = ".jpg", recursive = T))

# Convert csv with manual annotations of keypoints into json format
annot.csv = list.files(in.folder, full.names = T, pattern = ".csv", recursive = T)
sapply(annot.csv, make.json)

annots = list.files(in.folder, full.names = T, pattern = ".json", recursive = T)

# Split training and test images 80-20
pct_index = floor(length(files)*0.8)
test = files[(pct_index+1):length(files)]
train = files[0:pct_index]

# Prepare output folders
test_image_folder = paste0(out.folder, "/test/images")
train_image_folder = paste0(out.folder, "/train/images")
test_annot_folder = paste0(out.folder, "/test/annotations")
train_annot_folder = paste0(out.folder, "/train/annotations")
weights_folder = paste0(out.folder, "/model/weights")

fs::dir_create(c(test_image_folder, test_annot_folder, 
  train_image_folder, train_annot_folder,
  weights_folder), recurse = T)

# Prepare each test image and annotation file
for(file in test){
  # new.name = folder name + file name
  paths       = fs::path_split(fs::path_rel(file, in.folder))
  folder_name =  paths[[1]][1]
  image_name  = paths[[1]][length(paths[[1]])]
  base        = fs::path_ext_remove(image_name)
  
  new.image.path = paste0(test_image_folder, "/", folder_name, "_",base, ".tiff" ) 
  new.annot.path = paste0(test_annot_folder, "/", folder_name, "_",base, ".json" )

  annot.file = fs::path_filter(annots, glob = paste0("*", folder_name, "*", base, ".json"))

  
  # Exported annotations may not contain nuclei in an image - skip these
  if(length(annot.file)==length(new.image.path)){
    fs::file_copy(file, new.image.path, overwrite=T)
    fs::file_copy(annot.file, new.annot.path, overwrite=T)
  }

}

# Prepare each training image and annotation file
for(file in train){
  # new.name = folder name + file name
  paths       = fs::path_split(fs::path_rel(file, in.folder))
  folder_name =  paths[[1]][1]
  image_name  = paths[[1]][length(paths[[1]])]
  base        = fs::path_ext_remove(image_name)
  
  new.image.path = paste0(train_image_folder, "/", folder_name, "_",base, ".tiff" ) 
  new.annot.path = paste0(train_annot_folder, "/", folder_name, "_",base, ".json" )

  annot.file = fs::path_filter(annots, glob = paste0("*", folder_name, "*", base, ".json"))

  
  # Exported annotations may not contain nuclei in an image - skip these
  if(length(annot.file)==length(new.image.path)){
    fs::file_copy(file, new.image.path, overwrite=T)
    fs::file_copy(annot.file, new.annot.path, overwrite=T)
  }
}