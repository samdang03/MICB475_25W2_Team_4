library(tidyverse)
library(phyloseq)
library(ape)
library(dplyr)
library(ggplot2)
library(GGally)
library(psych)

meta_wetlands_categorical <- meta_wetlands %>%
  mutate(across(where(~ is.character(.x) | is.factor(.x)),
                ~ na_if(.x, "Missing: Not provided"))) %>%
  mutate(conductivity = as.numeric(conductivity)) %>%
  mutate(conductivity_category = factor(
    if_else(conductivity >= 0.3, "high", "low"),
    levels = c("low", "high")))

write.table(meta_wetlands_categorical, file = "wetlands_metadata_categorical.txt",
    sep = "\t", quote = FALSE, row.names = FALSE)