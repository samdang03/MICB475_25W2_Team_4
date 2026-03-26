library(tidyverse)
library(phyloseq)
library(indicspecies)

#Load Data 
load("phylo_objects/phylo_final.RData")

#### Indicator Species/Taxa Analysis ####

# Filter out rare genus 
phylo_filtered <- filter_taxa(phylo_final, function(x) sum(x > 0) > 5, TRUE)

# glom to Genus
isa_genus <- tax_glom(phylo_filtered, "Genus", NArm = FALSE)
isa_genus_RA <- transform_sample_counts(isa_genus, fun=function(x) x/sum(x))
table(sample_data(isa_genus_RA)$`conductivity levels`)


#ISA
isa_multipatt <- multipatt(t(otu_table(isa_genus_RA)), cluster = sample_data(isa_genus_RA)$`conductivity_category`)
summary(isa_multipatt)
taxtable <- tax_table(phylo_filtered) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Generate ISA Table 
ISA_table <- isa_multipatt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()




