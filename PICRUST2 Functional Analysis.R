install.packages("devtools")
devtools::install_github("cafferychen777/ggpicrust2")
install.packages('MicrobiomeStat')
BiocManager::install("KEGGREST")
install.packages("GGally")

library(tidyverse)
library(phyloseq)
library(ggpicrust2)

meta <- readRDS('../Datasets/phyloseq_taxonomy.rds') %>% .@sam_data %>% data.frame() %>% 
  rownames_to_column('sample_name')