####Load Packages####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

load("Phylo_final.RData")

#### "core" microbiome ####

# Convert to relative abundance
Phylo_RA <- transform_sample_counts(Phylo_final, fun=function(x) x/sum(x))

# Filter dataset by conductivity category
Phylo_high <- subset_samples(Phylo_RA, `conductivity_category`=="high")
Phylo_low <- subset_samples(Phylo_RA, `conductivity_category`=="low")

# What ASVs are found in more than 50% of samples in each conductivity category?
# trying changing the prevalence to see what happens
high_ASVs <- core_members(Phylo_high, detection=0.001, prevalence = 0.5)
low_ASVs <- core_members(Phylo_low, detection=0.001, prevalence = 0.5)

# What are these ASVs (Family Level Interest)? 
tax_table(prune_taxa(high_ASVs,Phylo_final))
tax_table(prune_taxa(low_ASVs,Phylo_final))

# can plot those ASVs' relative abundance
core_plot_0.5p <- prune_taxa(high_ASVs,Phylo_RA) %>% 
  plot_bar(fill="Family") + 
  facet_wrap(.~`conductivity_category`, scales ="free")

# Generate Venn Diagram
high_list <- core_members(Phylo_high, detection=0.001, prevalence = 0.5)
low_list  <- core_members(Phylo_low,  detection=0.001, prevalence = 0.5)

cond_list <- list(High = high_list, Low = low_list)

venn_plot <- ggVennDiagram(cond_list)

ggsave("venn_conductivity.png", venn_plot)
ggsave("core_plot_0.5p.png", core_plot_0.5p)