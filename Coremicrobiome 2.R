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

# What are these ASVs? 
tax_table(prune_taxa(high_ASVs,Phylo_final))

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

####Troubleshooting###
# 1. Search the tax_table directly to see if the family exists
any(tax_table(Phylo_final)[,"Family"] == "f__Desulfobulbaceae")

# 2. If it exists, let's subset specifically for it and plot its abundance
Phylo_Desulfo <- subset_taxa(Phylo_RA, Family == "f__Desulfobulbaceae")

plot_bar(Phylo_Desulfo, fill="Family") + 
  facet_wrap(~conductivity_category, scales="free") +
  theme(legend.position="bottom")

####Try changing parameters to 0.2 prevalence and plot top 10 families based on abundance
## 1. Identify "Core" members at 0.2 prevalence
high_ASVs_02 <- core_members(Phylo_high, detection=0.001, prevalence = 0.2)
core_02_phylo <- prune_taxa(high_ASVs_02, Phylo_RA)

## 2. Get the NAMES of the Top 10 families
# We merge (glom) by Family first to rank them accurately
family_glom <- tax_glom(core_02_phylo, taxrank = "Family")

# Get the IDs of the top 10 most abundant family groups
top10_ids <- names(sort(taxa_sums(family_glom), decreasing = TRUE)[1:10])

# Map those IDs to their actual Family names (e.g., "f__Desulfobulbaceae")
top10_families <- as.character(tax_table(family_glom)[top10_ids, "Family"])

## 3. Plot results
core_plot <- core_02_phylo %>% 
  subset_taxa(Family %in% top10_families) %>%
  plot_bar(fill="Family") + 
  facet_wrap(~conductivity_category, scales ="free") +
  theme_bw() +
  theme(axis.text.x = element_blank(), # Removes messy filename text
        axis.ticks.x = element_blank()) +
  labs(title = "Top 10 Core Families (20% Prevalence)",
       y = "Relative Abundance")
core_plot


