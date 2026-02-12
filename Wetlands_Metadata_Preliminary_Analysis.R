library(tidyverse)
library(phyloseq)
library(ape)
library(dplyr)
library(ggplot2)
library(GGally)
library(psych)

metafp <- "wetlands_metadata.txt"
meta_wetlands <- read_delim(metafp, delim = "\t")
http://127.0.0.1:23559/graphics/plot_zoom_png?width=1904&height=970
meta_scatter_data <- meta_wetlands %>%
  select(conductivity,
         exchange_capacity,
         potential_exchange_capacity,
         ch4_flux,
         co2_flux,
         soil_no3,
         ph,
         soil_nh4,
         potential_net_n_mineralization,
         potential_net_nitrification) %>%
  mutate(across(everything(),
                ~ suppressWarnings(as.numeric(.))))%>%
  filter(ch4_flux < 0.2)

Scatter_Matrix <- ggpairs(
  meta_scatter_data,
  title = "Scatter Plot Matrix"
)

Scatter_Matrix

ggsave("scatter plot matrix.png", Scatter_Matrix, 
       width = 7, height = 7)
