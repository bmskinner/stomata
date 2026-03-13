# Stomata detection via YOLO

This repo has analysis for maize stomata images. 



# Detecting stomata

- Images were annotated using (X)AnyLabelling with Segment Anything Model (SAM) to generate segmented outlines as JSON.
- JSON outlines were converted to YOLO format for detection by oriented bounding boxes
- OBBs from YOLO inferencing were saved to the `./analysis` directory

# Analysing stomatal patterns

We are interested in how the stomata are arranged into files. To do this, we determine which stomata are in 'chains'; straight lines in a consistent orientation within the image.


- createChains.R to identify stomata in straight lines and serialise the output to a `.Rds` file in the `./analysis` directory
- analyseChains.R to calculate chain measurements from the serialised file