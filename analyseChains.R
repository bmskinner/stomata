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
library(sf)
library(patchwork)
library(fs)
fs::dir_create("figure")

chain.distance.file <- "chain.distances.Rds"
if(!file.exists(chain.distance.file)){
  
  # Remove existing png output files
  unlink(list.files(path="stomatal_image", pattern = "*.mer*.*.png", recursive = T, full.names = T))
  unlink(list.files(path="stomatal_image", pattern = "*chain.*.png", recursive = T, full.names = T))
  
  # Read all json files
  files <- list.files(path = "stomatal_image", pattern = "*.json",  recursive = T, full.names = T)
  
  chain.distances <- do.call(rbind, lapply(files, process.coordinate.file))
  chain.distances$Folder <- dirname(chain.distances$File)
  saveRDS(chain.distances, file="chain.distances.Rds")
}

chain.distances <- readRDS("chain.distances.Rds") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Folder = dirname(File),
                File = basename(File))

dist.plot <- ggplot(chain.distances, aes(x=File, y = length))+
  geom_hline(yintercept = median(chain.distances$length))+
  geom_violin()+
  geom_boxplot(width=0.2, alpha=0)+
  labs(y = "Distance between stomata pairs (pixels)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("figure/Distance_plot.png", plot = dist.plot, dpi = 300, units = "mm", width = 170, height = 170)

# Are the distances significantly different?
kruskal.data <- data.frame("Folder" =  chain.distances$Folder, "S1S2" = chain.distances$length)
kruskal.test(S1S2 ~ Folder, data  = kruskal.data)
kruskal.data.dunn <- rstatix::dunn_test(kruskal.data, formula = length ~ Folder, detailed = T)


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
                    MeanAbsDeviance = SumAbsDeviance / nStomata,
                    RootSumSqareDeviance = sqrt( sum(Deviance^2)),
                    RootMeanSquareDeviance = sqrt( sum(Deviance^2)/nStomata  ),
                    DevianceRatio = MeanDeviance/MeanAbsDeviance
  )

# Overall levels of deviance from straight line in contig
deviance.plot <- ggplot(deviance.data, aes(x = File, y = MeanDeviance))+
  geom_hline(yintercept = 0)+
  geom_violin()+
  # geom_beeswarm()+
  geom_boxplot(width = 0.2, alpha = 0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")
ggsave("figure/deviance_overall.png", deviance.plot, dpi = 300, units = "mm", width = 170, height = 170)

# Mean deviance versus sum of squares - consistency of bends
deviance.consistency.plot <- ggplot(deviance.data, aes(x = MeanAbsDeviance, y = MeanDeviance))+
  # geom_hex(bins = 100)+
  geom_abline(intercept = c(0, 0), slope = 1, col="grey")+ # Consistent bend
  geom_abline(intercept = c(0, 0), slope = -1, col="grey")+ # Consistent bend
  geom_abline(intercept = c(0, 0), slope = 0.5, col="grey")+ # Wibble
  geom_abline(intercept = c(0, 0), slope = -0.5, col="grey")+ # Wibble
  geom_point()+
  coord_fixed(xlim=c(0, 30), ylim=c(-30, 30), expand = FALSE)+
  scale_fill_viridis_c(limits = c(0, 70))+
  labs(fill = "Number of chains", x = "Mean absolute deviance (pixels)", y = "Mean deviance (pixels)")+
  facet_wrap(~File)+
  theme_bw()+
  theme(legend.position = "top")
ggsave("figure/deviance_consistency.png", deviance.consistency.plot, dpi = 300, units = "mm", width = 85, height = 170)


deviance.ratio.plot.data <- deviance.data %>%
  dplyr::mutate(Class = case_when(DevianceRatio > 0.9 ~ "Perfect",
                                  DevianceRatio < -0.9 ~ "Perfect",
                                  .default = "Intermediate"))

deviance.ratio.plot <- ggplot(deviance.ratio.plot.data, aes(x = File, y = DevianceRatio))+
  geom_hline(yintercept = median(deviance.data$DevianceRatio))+
  geom_violin()+
  geom_boxplot(width=0.2, alpha=0)+
  coord_cartesian()+
  facet_wrap(~Class)+
  theme_bw()+
  theme(legend.position = "none")
ggsave("figure/deviance_ratio.png", deviance.ratio.plot, dpi = 300, units = "mm", width = 170, height = 170)

unassigned.stomata <- chain.distances %>%
  dplyr::ungroup() %>%
  dplyr::select(Folder, File, fUnassignedStomata) %>%
  dplyr::distinct()
unass.plot <- ggplot(unassigned.stomata, aes(x=Folder, y=fUnassignedStomata))+
  geom_hline(yintercept = median(unassigned.stomata$fUnassignedStomata))+
  geom_violin()+
  geom_boxplot(width=0.2, alpha=0)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("figure/unass.plot.png", unass.plot, dpi = 300, units = "mm", width = 170, height = 170)

overall.plot <- (dist.plot + deviance.plot )/ (deviance.consistency.plot +unass.plot) + plot_annotation(tag_levels = c("A"))
ggsave("figure/output.png", overall.plot, dpi = 300, units = "mm", width = 170, height = 170)

