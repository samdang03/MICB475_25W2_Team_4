library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)

#### Load data ####
metafp <- "Phyloseq_Files/wetlands_metadata_NEW 2.txt" 
meta <- read_delim(metafp, delim="\t")

otufp <- "Phyloseq_Files/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "Phyloseq_Files/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "Phyloseq_Files/tree.nwk"
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

Phylo_object <- phyloseq(OTU, SAMP, TAX, phylotree)
Phylo_object

tax_table(Phylo_final)
#Filter out NAs
Phylo_filtered <- subset_samples(Phylo_object, !is.na(conductivity_category))

#Filter out Chloroplast and Mitochondria
Phylo_final <- subset_taxa(Phylo_filtered, Domain == "d__Bacteria" & Family != "f__Chloroplast" & Family != "f__Mitochondria")

#Rarefaction
#rarecurve(t(as(otu_table(Phylo_final), "matrix")), cex=0.1)
Phylo_rarefied <- rarefy_even_depth(Phylo_final, rngseed = 123, replace = FALSE) 

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(Phylo_rarefied)), phy_tree(Phylo_rarefied),
                 include.root = FALSE)

# add PD to metadata table
sample_data(Phylo_rarefied)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(Phylo_rarefied), aes(conductivity_category, PD)) + 
  geom_boxplot() +
  xlab("Conductivity Category") +
  ylab("Faith's Phylogenetic Diversity")

# view plot
plot.pd

ggsave(filename = "newplot_Faith's.png"
       , plot.pd
       , height=4, width=6)

##stats##

#build dataframe combining alpha measures and metadata
df <- data.frame(sample_data(Phylo_rarefied))
df$sample <- sample_names(Phylo_rarefied)
df$PD <- phylo_dist$PD[match(df$sample, rownames(phylo_dist))]

#run wilcox.test
wilcox.test(PD ~ conductivity_category, data = df)

plot_nostats <- ggplot(df, aes(x = conductivity_category, y = PD)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  xlab("Conductivity Category") +
  ylab("Faith's Phylogenetic Diversity")
plot_nostats

plot_stats <-plot_stats <- ggplot(df, aes(x = conductivity_category, y = PD, fill = conductivity_category)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_signif(
    comparisons = list(c("low", "high")),
    y_position = max(df$PD, na.rm = TRUE) * 1.05,
    annotations = c("**")
  ) +
  scale_fill_manual(values = c(
    "low" = "#D55E00",
    "high" = "#0072B2"
  )) +
  expand_limits(y = max(df$PD, na.rm = TRUE) * 1.15) +
  xlab("Conductivity Category") +
  ylab("Faith's Phylogenetic Diversity")
plot_stats

ggsave(filename = "plot_Faith's_withstats.png"
       , plot_stats
       , height=4, width=6)