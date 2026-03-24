library(phyloseq)
library(ape)
library(tidyverse)
library(ggplot2)
library(picante)

#Initial loading (code adapted from Ruth + Richard)
metafp <- "wetlands_metadata_NEW 2.txt" 
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
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

#Desulfobacteriae phylum only
phylo_desulfo <- subset_taxa(phylo_nmnc, Phylum =="p__Thermodesulfobacteriota")

unique_family <- get_taxa_unique(phylo_desulfo, "Family")

#Thermodesulfobacteriaceae family only

phylo_desulfo_f <- subset_taxa(phylo_nmnc, Family == "f__Thermodesulfobacteriaceae")
get_taxa_unique(phylo_desulfo_f, "Genus")
#There's only one genus (Thermosulfuriphilus), so this will not be continued with
#The final plot will be focused entirely on families within the Thermodesulfobacteriota phylum

#Filtering taxonomy table to top 10 + other
phylo_desulfo_family <- tax_glom(phylo_desulfo, taxrank = "Family")

top_10 <- names(sort(taxa_sums(phylo_desulfo_family), TRUE)[1:10])

tax_table(phylo_desulfo_family)[!taxa_names(phylo_desulfo_family) %in% top_10, "Family"] <- "Other"

#Rarefaction
phylo_desulfo_rarefied <- rarefy_even_depth(phylo_desulfo_family, sample.size = 1000, rngseed = 42, replace = FALSE) 

#Plot creation
phy_plot <- plot_bar(phylo_desulfo_rarefied, fill = "Family") +
  facet_wrap(~conductivity_category, scales = "free_x") +
  theme_classic() +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack", color = NA) +
  labs(title = "Relative Abundance of Microbial Phyla",
    subtitle = "Grouped by Sediment Conductivity Category",
    x = "Individual Samples",
    y = "Relative Abundance (%)",
    caption = "Non-bacterial taxa filtered.")
phy_plot

