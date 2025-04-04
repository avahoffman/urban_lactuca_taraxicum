library(tidyverse) # CRAN v2.0.0
library(ggh4x)     # CRAN v0.2.8
library(cowplot)   # CRAN v1.1.3

make_structure_plot <- function(spp_,
                                species_name,
                                structure_plot = T,
                                width = 12,
                                height = 3.5) {
  # Read in structure results
  long_df <- read_csv(paste0("data/", spp_, "_structure_results.csv"))
  
  # Clean up genetic sample names so they match Eric's names
  long_df <-
    long_df %>%
    mutate(sample = str_replace(sample, "\\.BA\\.(?=.*\\.)", ".BAL.")) %>% #replace BA with BAL
    mutate(sample = str_replace(sample, "\\.MN\\.(?=.*\\.)", ".MSP.")) %>% #replace BA with BAL
    mutate(sample = str_replace(sample, "\\.[UM]\\.", ".")) %>% #remove managed/unmanaged
    mutate(sample = str_replace(sample, "-[^.]*\\.", ".")) #look for dash in site name and remove it + spp
  
  # Read in Eric's samples
  selection_ <-
    read_csv("data/LineageList.csv") %>%
    mutate(sample = paste0(Species, ".", City, ".", `Site.Line`)) %>%
    filter(Species == spp_) %>%
    pull(sample)
  
  # Filter accordingly
  long_df <-
    long_df %>% 
    filter(sample %in% selection_)
  
  long_df <- long_df %>%
    mutate(city = case_when(
      city == "Minneapolis" ~ "Minneapolis-Saint Paul",
      TRUE ~ city
    ))  %>%
    mutate(sample = fct_reorder(sample, nlcd_urban_pct))
  
  # ----- Reorder cities -----
  
  long_df$city <- factor(long_df$city, levels = c("Minneapolis-Saint Paul", "Boston", "Baltimore", "Los Angeles", "Phoenix"))
  
  # Make plot
  if (structure_plot) {
    gg <-
      ggplot(data = long_df, aes(x = sample, y = value, fill = name)) +
      geom_col(width = 1, color = NA) +
      facet_nested(~ city,
                   scales = "free_x",
                   space = "free") +
      theme_classic() +
      labs(y = species_name) +
      scale_fill_manual(values = c("#30123BFF", "#28BBECFF", "#A2FC3CFF", "#FB8022FF")) + # viridis::turbo(n = 5)
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        ggh4x.facet.nestline = element_line(linetype = 3),
        axis.ticks.length = unit(-1, "inch"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(hjust = 0.1),
        plot.margin = unit(c(0,5,0,5), "pt"),
        text = element_text(size = 8)
      )
    
  } else {
    gg <-
      ggplot(data = long_df, aes(x = sample, y = value, fill = nlcd_urban_pct)) +
      geom_col(width = 1, color = NA) +
      facet_nested(~ city,
                   scales = "free_x",
                   space = "free",
                   switch = "both") +
      theme_void() +
      scale_fill_viridis_c(option = "A") +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        panel.spacing = unit(0, "lines"),
        #ggh4x.facet.nestline = element_line(linetype = 3),
        strip.text.x = element_blank(),
        axis.ticks.length = unit(-1, "inch"),
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x.top = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(15,5,0,5), "pt")
      )
  }
  
  ggsave(
    paste0("figures/structure_", species_name, ".png"),
    dpi = "print",
    width = width,
    height = height
  )
  return(gg)
  
}

make_all <- function(){
  plot_grid(
    make_structure_plot(spp = "LS", structure_plot = F, species_name = "prickly lettuce"),
    make_structure_plot(spp = "LS", species_name = "prickly lettuce"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  ggsave(
    paste0("figures/structure_with_urban_prickly_lettuce.png"),
    dpi = "print",
    width = 5,
    height = 2
  )  
  
  plot_grid(
    make_structure_plot(spp = "TO", structure_plot = F, species_name = "dandelion"),
    make_structure_plot(spp = "TO", species_name = "dandelion"),
    ncol = 1,
    align = "v",
    rel_heights = c(30, 100)
  )
  ggsave(
    paste0("figures/structure_with_urban_dandelion.png"),
    dpi = "print",
    width = 5,
    height = 2
  )  
}
