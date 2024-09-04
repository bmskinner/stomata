# Stomata detection via YOLO

This repo has analysis for the maize stomata images.

Images were annotated using AnyLabelling with Segment Anything Model (SAM) to generate segmented outlines as JSON.

# Detecting stomata

createYOLOTrainingData.R

- JSON outlines are converted to YOLO format for (a) detection by bounding box and (b) segmentation
- folder structure for training and validation built
- images are split to training and validation sets in `./data`
- python scripts are written for submission onto the cluster GPU nodes to train and predict the detection and segmentation models
- scripts and images are zipped for upload to the cluster

# Analysing stomatal patterns

- createChains.R to contig stomata in straight lines
- analyseChains.R to calculate chain measurements
- process.nma.data.R to analyse NMA outputs from predicted segmented outlines