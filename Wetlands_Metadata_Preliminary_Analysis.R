library(tidyverse)
library(phyloseq)
library(ape)
library(dplyr)
library(ggplot2)
library(GGally)
library(psych)

metafp <- "wetlands_metadata.txt"
meta_wetlands <- read_delim(metafp, delim = "\t")
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
  filter(potential_net_n_mineralization < 10 |
           is.na(potential_net_n_mineralization))
  

Scatter_Matrix <- ggpairs(
  meta_scatter_data,
  title = "Scatter Plot Matrix"
)

Scatter_Matrix

ggsave("scatter plot matrix.png", Scatter_Matrix, 
       width = 7, height = 7)

cond_vs_soilnh4 <- meta_scatter_data %>%
  ggplot(aes(x = conductivity, y = soil_nh4)) +
  geom_point() +
  labs(x = "Conductivity (mS/cm)",
    y = "Soil NH4 (µg N per g dry soil)")
cond_vs_soilnh4

ggsave("Conductivity vs Soil NH4.png", cond_vs_soilnh4, 
       width = 7, height = 7)

cond_vs_n_mineralization <- meta_scatter_data %>%
  ggplot(aes(x = conductivity, y = potential_net_n_mineralization)) +
  geom_point() +
  geom_vline(xintercept = 0.3,
             linetype = "dashed",
             linewidth = 1) +
  labs(x = "Conductivity (mS/cm)",
       y = "N Mineralization (ug N per g dry soil per day)")
cond_vs_n_mineralization

ggsave("Conductivity vs Nitrogen Mineralization.png", cond_vs_n_mineralization, 
       width = 7, height = 7)

