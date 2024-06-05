# Parse Rds output from analysed YOLO stomata coordinates
source("functions.R")
fs::dir_create("figure")

input.chain.files <- list.files(path="analysis", pattern = "*.Rds", include.dirs = TRUE, full.names = TRUE, recursive = TRUE)

read.chain.file.metadata <- function(file){
  readRDS(file)$metadata
}

read.chain.file.contigs <- function(file){
  readRDS(file)$contigs
}

read.chain.file.measurements <- function(file){
  readRDS(file)$measurments
}

chain.metadata <- do.call(rbind, lapply(input.chain.files, read.chain.file.metadata))
chain.contigs  <- do.call(rbind, lapply(input.chain.files, read.chain.file.contigs))
chain.measure  <- do.call(rbind, lapply(input.chain.files, read.chain.file.measurements))

cat(sum(chain.metadata$nStomata), "stomata analysed\n")

#### Basic summary info ####

# How many images per folder?
save.ggplot(ggplot(chain.metadata %>% dplyr::group_by(Folder) %>% dplyr::summarise(n = n()), aes(x=Folder, y = n))+
              geom_col()+
              labs(y = "Images per folder")+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1)),
            "figure/Images_per_folder.png")


# How many contigs per image?
save.ggplot(ggplot(chain.metadata, aes(x=Folder, y = nContigs))+
              geom_hline(yintercept = median(chain.metadata$nContigs))+
              geom_violin()+
              geom_boxplot(width=0.2, alpha=0)+
              labs(y = "Chains per image")+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1)),
            "figure/Chains_per_image.png")

# How many stomata per image?
save.ggplot(ggplot(chain.metadata, aes(x=Folder, y = nStomata))+
              geom_hline(yintercept = median(chain.metadata$nStomata))+
              geom_violin()+
              geom_boxplot(width=0.2, alpha=0)+
              labs(y = "Stomata per image")+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1)),
            "figure/Stomata_per_image.png")


# What is the length of the contigs in pixels?
save.ggplot(ggplot(chain.contigs, aes(x=Folder, y = length))+
              geom_hline(yintercept = median(chain.contigs$length))+
              geom_violin()+
              geom_boxplot(width=0.2, alpha=0)+
              labs(y = "Distance between stomata pairs (pixels)")+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1)),
            "figure/Distance_between_stomata_pixels.png")

contig.totals <- chain.contigs %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Contig, Folder, File) %>%
  dplyr::summarise(n = n(),
                   length = sum(length))

# How many stomata per contig?
save.ggplot(ggplot(contig.totals, aes(x=Folder, y = n))+
              geom_hline(yintercept = median(contig.totals$n))+
              # geom_beeswarm()+
              geom_violin()+
              geom_boxplot(width=0.2, alpha=0)+
              labs(y = "Contig length (stomata")+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1)),
            "figure/Stomata_per_contig.png")

# How many unassigned stomata per image?
save.ggplot(ggplot(chain.metadata, aes(x=Folder, y = fUnassignedStomata))+
              geom_hline(yintercept = median(chain.metadata$fUnassignedStomata))+
              # geom_beeswarm()+
              geom_violin()+
              geom_boxplot(width=0.2, alpha=0)+
              labs(y = "Fraction unassigned stomata")+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1)),
            "figure/Unassigned_stomata_per_image.png")

#### Detailed analysis ####

# Are the distances significantly different?
# kruskal.data <- data.frame("Folder" =  chain.measure$Folder, "length" = chain.measure$length)
# kruskal.test(length ~ Folder, data  = kruskal.data)
# kruskal.data.dunn <- rstatix::dunn_test(kruskal.data, formula = length ~ Folder, detailed = T)


# Look at deviances in each contig
deviance.data <- chain.measure %>%
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
save.ggplot(ggplot(deviance.data, aes(x = Folder, y = MeanDeviance))+
              geom_hline(yintercept = 0)+
              geom_violin()+
              # geom_beeswarm()+
              geom_boxplot(width = 0.2, alpha = 0)+
              theme_bw()+
              theme(axis.text.x = element_text(angle = 45, hjust=1),
                    legend.position = "none"),
            "figure/Deviance_per_folder.png")

# Mean deviance versus sum of squares - consistency of bends
save.ggplot(ggplot(deviance.data, aes(x = MeanAbsDeviance, y = MeanDeviance))+
              geom_hex(bins = 100)+
              geom_abline(intercept = c(0, 0), slope = 1, col="grey")+ # Consistent bend
              geom_abline(intercept = c(0, 0), slope = -1, col="grey")+ # Consistent bend
              geom_abline(intercept = c(0, 0), slope = 0.5, col="grey")+ # Wibble
              geom_abline(intercept = c(0, 0), slope = -0.5, col="grey")+ # Wibble
              # geom_point()+
              coord_fixed(xlim=c(0, 30), ylim=c(-30, 30), expand = FALSE)+
              scale_fill_viridis_c()+
              labs(fill = "Number of chains", x = "Mean absolute deviance (pixels)", y = "Mean deviance (pixels)")+
              # facet_wrap(~File)+
              theme_bw()+
              theme(legend.position = "top"),
            "figure/Deviance_consistency.png")

#### G-function analysis ####

# Run the g-function on each image
# Get the CoMs of each stomata
coms <- chain.contigs %>% dplyr::select(S = S1 , x = S1.x, y= S1.y, Image, Folder, File)
coms <- rbind (coms, chain.contigs %>% dplyr::select(S = S2 , x = S2.x, y= S2.y, Image, Folder, File))
coms %<>% dplyr::distinct()

g.results <- do.call(rbind, lapply(unique(coms$Image), function(i) g.function(as.matrix(coms[coms$Image==i, c("x", "y")]), i)))

g.results$Folder <- basename(dirname(g.results$Image))

# Plot the per-file plots
save.ggplot(ggplot(g.results, aes(x=distance, y=Gd, group=Image))+
  geom_line(data =g.results[,c("distance", "Gd", "Image")], col="grey", alpha = 0.1)+
  geom_line(col="blue", alpha = 0.1)+
  labs(x = "Distance between stomata", y = "Cumulative fraction")+
  facet_wrap(~Folder )+
  theme_bw(),
  "figure/G-function.png", height = 170)

# What is the mean value for each folder?
# Linear interpolation of per-file curve to consistent spacing

windows <- as.data.frame(IRanges(start = seq(40, 400, by = 10), # vector of window start positions
                end   = seq(50, 410, by = 10)))

calc.mean <- function(start, end){
  do.call(rbind, lapply(unique(g.results$Folder), function(folder){
    subset.data <- g.results[g.results$Folder==folder & g.results$distance > start & g.results$distance <= end,]
    data.frame(start = start, end = end, Gd = mean(subset.data$Gd), Gd.sd = sd(subset.data$Gd), Folder = folder, distance = (start+end)/2)
  }))
}

g.summary <- do.call(rbind, mapply(calc.mean, windows$start, windows$end, SIMPLIFY = FALSE))

save.ggplot(ggplot(g.summary, aes(x=distance, y=Gd))+
  geom_line(data=g.results[,c("distance", "Gd", "Image")], aes(group=Image), col="grey", alpha = 0.1)+
  geom_line(data=g.results, aes(group=Image), col="lightblue", alpha = 0.3)+
  geom_line(col="blue")+
  labs(x = "Distance between stomata", y = "Cumulative fraction")+
  facet_wrap(~Folder )+
  theme_bw(),
  "figure/G-function_complete.png", height = 170)

save.ggplot(ggplot(g.summary, aes(x=distance, y=Gd, col = Folder))+
  geom_line()+
  labs(x = "Distance between stomata", y = "Cumulative fraction")+
  theme_bw()+
    theme(legend.position = "top"),
  "figure/G-function_summary.png", height = 170)

