library(phyloseq)
library(ape)
library(tidyverse)
library(ggplot2)
library(picante)

#Initial loading (code adapted from Ruth + Richard)
metafp <- "Phyloseq_Files/wetlands_metadata_NEW 2.txt" 
meta <- read_delim(metafp, delim="\t")

otufp <- "Phyloseq_Files/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "Phyloseq_Files/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "Phyloseq_Files/tree.nwk"
phylotree <- read.tree(phylotreefp)

#Taxonomy formatting
tax_format <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) 

#Matrix/dataframe conversion 

otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`

tax_mat <- as.matrix(tax_format[,-1])
rownames(tax_mat) <- tax$`Feature ID`

met_df <- as.data.frame(meta[,-1])
rownames(met_df) <- meta$`#SampleID`

#Phyloseq component conversion

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
SAMP <- sample_data(met_df)

#Phyloseq creation

phylo <- phyloseq(OTU, SAMP, TAX, phylotree)
phylo

#Removing no conductivity data samples
phylo_filtered <- subset_samples(phylo, !is.na(conductivity_category))

#Removing chloroplasts and eukaryotes
phylo_nmnc <- subset_taxa(phylo_filtered, Domain == "d__Bacteria" & Family != "f__Chloroplast" & Family != "f__Mitochondria")

phylo_notree <- phyloseq(
  otu_table(phylo_nmnc),
  tax_table(phylo_nmnc),
  sample_data(phylo_nmnc))


# Convert to relative abundance (no tree for speed)
phy_rel <- transform_sample_counts(phylo_notree, function(x) x / sum(x))

# Extract OTU + taxonomy
otu_rel <- as.data.frame(otu_table(phy_rel))
tax_rel <- as.data.frame(tax_table(phy_rel))

# Ensure correct orientation
if(!taxa_are_rows(phy_rel)){
  otu_rel <- t(otu_rel)
}

# Identify Desulfobulbaceae
target_row <- rownames(tax_rel)[grepl("Desulfobulbaceae", tax_rel$Family)]

# Sum across all matching taxa
abund_vec <- colSums(otu_rel[target_row, , drop = FALSE])

# Build abundance dataframe
abund_df <- data.frame(
  Sample = names(abund_vec),
  Abundance = as.numeric(abund_vec)
)

# Extract metadata
meta_rel <- data.frame(sample_data(phylo_nmnc))
meta_rel$Sample <- rownames(meta_rel)

# Merge
final_df <- dplyr::left_join(meta_rel, abund_df, by = "Sample")

# Replace missing with 0
final_df$Abundance[is.na(final_df$Abundance)] <- 0

ggplot(final_df, aes(x =conductivity, y = Abundance)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "Conductivity (mS/cm)", y = "Relative Abundance of Desulfobulbaceae",
    title = "Conductivity vs Desulfobulbaceae Relative Abundance") +
  theme_classic()


ggsave("Continuous Conductivity vs Desulfobulbaceae.png",
       width = 8, height = 6, dpi = 300)