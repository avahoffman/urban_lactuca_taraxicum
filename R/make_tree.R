#' We used clonal distance to discern among individuals. According to GenoDive:
#' Clonal – This option calculates the number of mutationsteps that is needed to 
#' transfer the genotype of the first individual into the genotype of the second 
#' individual (Meirmans & Van Tienderen, 2004). Strictly taken, this assumes 
#' clonal reproduction but the distance may also be useful in other cases. 
#' Selecting this option from the distances menu will show a dialog in which you 
#' can choose the mutation model that is used and set some additional parameters.
#' 
#' We choose an infinite allele model here because each mutation creates a new 
#' allele, rather than stepwise changes in copy number.

library(adegenet)      # CRAN v2.1.10 # as.matrix
library(tidyverse)     # CRAN v2.0.0
library(SNPRelate)     # Bioconductor v1.38.1
library(ggtree)        # Bioconductor v3.12.0
# BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse) # [github::YuLab-SMU/treedataverse] v0.0.1


create_tree <- function(spp_, dist_type = "clonal") {
  # Read in distances as calculated by GenoDive
  mat_ <- read_delim(paste0("data/", spp_, "_", dist_type, "_dist.tsv"))
  names_ <- mat_ %>% pull(`Obs.`)
  
  # Clean up genetic sample names so they match Eric's names
  matched_names_ <-
    names_ %>%
    str_replace("\\.BA\\.(?=.*\\.)", ".BAL.") %>% #replace BA with BAL
    str_replace("\\.MN\\.(?=.*\\.)", ".MSP.") %>% #replace BA with BAL
    str_replace("-[^.]*\\.", ".") #look for dash in site name and remove it + spp
  
  # Read in Eric's samples
  selection_ <-
    read_csv("data/LineageList_04.07.2025_EGY.csv") %>%
    mutate(sample = paste0(Species, ".", City, ".", `Site.Mng`, ".", Line)) %>%
    filter(Species == spp_) %>%
    pull(sample)
  
  # Keep only Eric's samples with genetic info
  # Many samples don't have genetic info due to low coverage/quality :(
  selection_with_genetic_info <-
    selection_[selection_ %in% matched_names_]
  
  # Subset accordingly
  mat_ <- as.matrix(mat_ %>% dplyr::select(-1), labels = T)
  dimnames(mat_) <- list(matched_names_, matched_names_)
  mat_selection_ <- mat_[selection_with_genetic_info, selection_with_genetic_info]
  
  # Perform hierarchical cluster analysis on the dissimilarity matrix
  sample.hc <- SNPRelate::snpgdsHCluster(mat_selection_)
  hc <- sample.hc$hclust
  
  # Set up for plotting
  p <- ggtree::ggtree(hc, branch.length = "none")
  d <- data.frame(
    label = hc$labels,
    city = dimnames(mat_selection_)[[1]] %>% str_extract("(?<=\\.)[^.]+(?=\\.)")
  )
  
  # Plotting color params
  colors_ <- viridis::turbo(n = 5)
  my_pal <- setNames(colors_, c("MSP", "BO", "BAL", "LA", "PX"))
  
  # Make plot
  p %<+% d +
    geom_tippoint(
      aes(fill = factor(city), x = x),
      size = 4,
      shape = 21,
      color = 'black'
    ) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 20),
      plot.margin = margin(0, 2, 0, 0, "cm")
    ) +
    geom_tiplab(hjust = -0.1) +
    scale_fill_manual(values = my_pal) +
    guides(fill = guide_legend(title = "")) +
    ggplot2::xlim(0, max(p$data$x) + 6) # margin so labels don't get cut off
  
  # Save plot
  ggsave(paste0(paste0(
    "figures/tree_", spp_, "_", dist_type, ".png"
  )),
  dpi = "print",
  width = 5,
  height = 5)
}

create_tree("LS")
create_tree("TO")

# Test some alternatives:
# create_tree("LS", "SP")
# create_tree("TO", "clonal_stepwise")
