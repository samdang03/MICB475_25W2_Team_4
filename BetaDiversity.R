library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggplot2)
library(vegan)

#### Load data ####
metafp <- "phylo_objects/wetlands_metadata_NEW.txt" 
meta <- read_delim(metafp, delim="\t")

otufp <- "phylo_objects/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "phylo_objects/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "phylo_objects/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)
#Convert to phyloseq Object:
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
class(OTU) 

#### Clean Metadata and Align the IDs###
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'#SampleID'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object

phylo <- phyloseq(OTU, SAMP, TAX, phylotree)
phylo

# Filter out samples with "Missing Data" in conductivity_category

phylo_final <- subset_samples(phylo, !is.na(conductivity_category))

#Filter out Chloroplast and Mitochondria 
phylo_final <- subset_taxa(phylo_final, Domain == "d__Bacteria" & Family!= "f__Chloroplast" & Family != "f__Mitochondria")

#Rarefy to the depth of your smallest sample
rarecurve(t(as.data.frame(otu_table(phylo_final))), cex=0.1)
phylo_rare <- rarefy_even_depth(phylo_final, rngseed =123, replace = FALSE)

#### Beta Diversity ####
# Calculate weighted UniFrac distance
wu_dm <- distance(phylo_rare, method = "wunifrac")

# Perform PCoA ordination
pcoa_wu <- ordinate(phylo_rare, method = "PCoA", distance = wu_dm)

# Plot ordination

gg_pcoa <- plot_ordination(phylo_final, pcoa_wu, color = "conductivity_category", shape = "conductivity_category") +
  labs(col = "Conductivity Level",shape = "Conductivity Level")

# Display the plot
gg_pcoa

#Save the Plot 
ggsave("wunifrac_pcoa.png"
       , gg_pcoa
       , height=4, width=5)

#### PERMANOVA analysis on Weighted Unifrac ####
# Extract sampel data from phyloseq 
samp_df <- data.frame(sample_data(phylo_rare))

# Use phyloseq to calculate weighted Unifrac distance matrix
adonis2(wu_dm ~ conductivity_category, data = samp_df)

# re-plot the above PCoA with ellipses to show a significant difference 
#between high and low conductiivty  using ggplot2

gg_pcoa_ellipse <- plot_ordination(
  phylo_rare, pcoa_wu,
  color = "conductivity_category",
  shape = "conductivity_category"
) +
  stat_ellipse(type = "norm") +
  labs(
    color = "Conductivity Level",
    shape = "Conductivity Level",
    title = "Weighted UniFrac PCoA",
    x = paste0("Axis 1 (", round(pcoa_wu$values$Relative_eig[1]*100, 1), "%)"),
    y = paste0("Axis 2 (", round(pcoa_wu$values$Relative_eig[2]*100, 1), "%)")
  ) +
  scale_color_manual(values = c("high" = "#0072B2", "low" = "#D55E00"))

gg_pcoa_ellipse

# Save the Plot 
ggsave("wunifrac_PERMANOVA_pcoa.png"
       , gg_pcoa_ellipse
       , height=4, width=5)

