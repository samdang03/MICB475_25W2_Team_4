library(tidyverse)
library(phyloseq)
library(ape)
library(dplyr)
library(ggplot2)
library(GGally)
library(psych)

metafp <- "wetlands_metadata.txt"
meta_wetlands <- read_delim(metafp, delim = "\t")
meta_wetlands_categorical <- meta_wetlands %>%
  filter(conductivity != "Missing: Not provided")%>%
  mutate(conductivity_category = factor(
    if_else(conductivity >= 0.3, "high", "low"),
    levels = c("low", "high")))

write.table(meta_wetlands_categorical, file = "wetlands_metadata_categorical.txt",
    sep = "\t", quote = FALSE, row.names = FALSE)