# Parse the JSON from Anylabelling to extract stomata CoMs

library(jsonlite)
library(tidyverse)
library(magrittr)
library(LearnGeom)
library(text.alignment)

MAX.DISTANCE <- 300
ANGLE.DELTA <- 10

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
    dplyr::mutate(stomata  = paste0("s", sprintf("%02d", row_number())),
                  file = file)
  
}

# Read all json files
files <- list.files(path = "stomatal_image", pattern = "*.json",  recursive = T, full.names = T)

# Find CoMs of all objects in one file for now
coms <- do.call(rbind, lapply(files[2], calc.coms)) %>% as.data.frame


# Now add in the chain code

# Calculate distance between two points
euclidean <- function(x1, y1, x2, y2) sqrt( (x1-x2)^2 + (y1-y2)^2)

# Create distance table between pairs of stomata
# Filter to only those within MAX.DISTANCE
dist.table <- expand.grid(coms$stomata, coms$stomata, stringsAsFactors = F) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(S1t = min(Var1, Var2), S2t = max(Var1, Var2)) %>% # ensure alphabetical order
  dplyr::ungroup() %>%
  dplyr::select(Var1 = S1t, Var2 = S2t) %>% # remove unused columns
  dplyr::distinct() %>% # remove duplicates
  dplyr::filter(Var1!=Var2) %>%
  merge(., coms[,1:3], by.x = "Var1", by.y = "stomata", all.y = F) %>% # add coordinates for S1
  dplyr::select(S1 = Var1, S1.x = x.com, S1.y = y.com, S2 = Var2) %>% # rename for clarity
  merge(., coms[,1:3], by.x = "S2", by.y = "stomata", all.y = F) %>% # add coordinates for S1
  dplyr::select(S1, S1.x, S1.y, S2, S2.x = x.com, S2.y = y.com) %>% # rename for clarity
  dplyr::mutate(S1S2 = euclidean(S1.x, S1.y, S2.x, S2.y)) %>%
  dplyr::filter(S1 != S2, S1S2 <= MAX.DISTANCE)

# Join the tables to create a 3-mer chain
mer3 <- merge(dist.table, dist.table, by.x = c("S2","S2.x", "S2.y" ), by.y = c("S1","S1.x", "S1.y" ), 
              all.y = F, suffixes = c("A", "B")) %>%
  dplyr::select(S1, S1.x, S1.y, S2, S2.x, S2.y, S3 = S2B, S3.x = S2.xB, 
                S3.y = S2.yB, S1S2 = S1S2A, S2S3 = S1S2B) %>% # rename for clarity
  dplyr::rowwise() %>%
  dplyr::mutate(angle = LearnGeom::Angle(c(S1.x,  S1.y), c(S2.x, S2.y), c(S3.x,S3.y))) %>% # calculate angles
  dplyr::filter( angle > 180 - ANGLE.DELTA)


# visualise the 3-mers
ggplot(mer3)+
  geom_point(aes(x = S1.x, y = S1.y))+
  geom_point(aes(x = S2.x, y = S2.y))+
  geom_point(aes(x = S3.x, y = S3.y))+
  geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y))+
  geom_segment(aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y)) + 
  geom_text(aes(x = S2.x, y = S2.y+50, label = sprintf("%.2f", angle)))
  

# Contig the 3mers based on last 2 and first 2 points
mer4 <- merge(mer3, mer3, by.x = c("S2","S3"), by.y = c("S1","S2" ), 
              all.y = F, suffixes = c("A", "B")) %>%
  dplyr::select(S1, S1.x = S1.xA, S1.y = S1.yA, 
                S2, S2.x = S2.xA, S2.y = S2.yA, 
                S3, S3.x = S3.xA, S3.y = S3.yA,
                S4 = S3B,  S4.x = S3.xB, S4.y = S3.yB) # rename for clarity

# visualise the 4-mers
ggplot(mer4)+
  geom_point(aes(x = S1.x, y = S1.y))+
  geom_point(aes(x = S2.x, y = S2.y))+
  geom_point(aes(x = S3.x, y = S3.y))+
  geom_point(aes(x = S4.x, y = S4.y))+
  geom_segment(aes(x = S1.x, y = S1.y, xend = S2.x, yend = S2.y)) +
  geom_segment(aes(x = S2.x, y = S2.y, xend = S3.x, yend = S3.y)) +
  geom_segment(aes(x = S3.x, y = S3.y, xend = S4.x, yend = S4.y)) 


# List all 4mers

mer4 %<>% dplyr::rowwise() %>%
  dplyr::mutate( chain =  paste(S1, S2, S3, S4, sep = "_"))
# 
# # Join 4mers to 6mers
# contig.chains <- function(chains){
#   output.chains <- list()
#   chain.number <- 1
# 
#   for(chain1 in chains){
#     
#     chain1.length = length(chain1)
#     
#     for(chain2 in chains){
#       
#       chain2.length = length(chain2)
#       
#       # Check overlap of end 2 elements
#       if(chain1[[chain2.length-1]] == chain2[[1]] & chain1[[chain2.length]] == chain2[[2]]){
#         print("overlap found")
#         
#         # Join into new chain
#         new.chain =  c(chain1[1:chain1.length-1], chain2[2:chain2.length])
#         
#         output.chains[[chain.number]] <-  new.chain
#         chain.number <- chain.number+1
#       }
#       
#       
#     }
#   }
#   output.chains
# }
# 
# 
# mer6 <- contig.chains(mer4$chain)
# mer8 <- contig.chains(mer6)

# TODO - make SW alignment of contigs
# scores <- needleScores(mer4$chain[1], subjects = mer4$chain, params=defaultNeedleParams)
# scores[scores > 2]
# 
# for(chain in mer4$chain){
#   scores <- needleScores(chain, subjects = mer4$chain, params=defaultNeedleParams)
#   print(paste(chain, " best matches:", paste(names(scores[scores > 6 & names(scores)!=chain]), collapse = "; "),
#               paste(scores[scores > 6 & names(scores)!=chain], collapse = "; ")))
# }

create.contigs <- function(kmers){
  sw.data <- data.frame(doc_id = seq(1:length(kmers)), text = kmers, stringsAsFactors = F)
  sw.result <- do.call(rbind, text.alignment::smith_waterman_pairwise(sw.data, sw.data)) %>%
    as_tibble(.) %>%
    unnest_wider(c(a, b), names_sep = "_") %>% # turn lists into new columns
    dplyr::mutate(sw = unlist(sw), a_doc_id = unlist(a_doc_id), b_doc_id = unlist(b_doc_id)) %>%
    dplyr::filter(a_doc_id != b_doc_id) %>%
    unnest_wider(c(a_alignment, b_alignment), names_sep = "_") %>%
    dplyr::filter(!str_detect(a_alignment_text, "#"),
                  !str_detect(b_alignment_text, "#")) %>% # only allow complete matches
    dplyr::select(sw, a_text, b_text, a_n, b_n, a_alignment_from, a_alignment_to, b_alignment_from, b_alignment_to) %>%
    dplyr::mutate(contig = paste0(str_sub(b_text, start = 1, end = b_alignment_to),
                                  str_sub(a_text, start = a_alignment_to+1, end = a_n)))
  # TODO- ensure we do not increase the contig count
  # only allow contigs if they are not completely embedded within another contig
  unique(sw.result$contig)
  
}

new.contigs <- create.contigs(mer4$chain)
prev.contig.count <- length(new.contigs)
while(TRUE){
  new.contigs <- create.contigs(new.contigs) 
  if(length(new.contigs)==prev.contig.count) break
  prev.contig.count = length(new.contigs)
}

# contigs <- create.contigs(mer4$chain)
# 
# cc.result <- create.contigs(contigs)

# sw.data <- data.frame(doc_id = seq(1:nrow(mer4)), text = mer4$chain, stringsAsFactors = F)
# sw.result <- do.call(rbind, text.alignment::smith_waterman_pairwise(sw.data, sw.data)) %>%
#   as_tibble(.) %>%
#   unnest_wider(c(a, b), names_sep = "_") %>% # turn lists into new columns
#   dplyr::mutate(sw = unlist(sw), a_doc_id = unlist(a_doc_id), b_doc_id = unlist(b_doc_id)) %>%
#   dplyr::filter(a_doc_id != b_doc_id) %>%
#   unnest_wider(c(a_alignment, b_alignment), names_sep = "_") %>%
#   dplyr::filter(a_alignment_n == a_alignment_from + a_alignment_to - 1) %>% # only allow complete matches
#   dplyr::select(sw, a_text, b_text, a_alignment_from, a_alignment_to, b_alignment_from, b_alignment_to) %>%
#   dplyr::mutate(contig = paste0(str_sub(b_text, start = 1, end = b_alignment_to),
#                                 str_sub(a_text, start = a_alignment_to+1, end = length(a_text)+1)))


