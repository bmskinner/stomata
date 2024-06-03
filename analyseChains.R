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
