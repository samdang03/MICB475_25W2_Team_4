library(phyloseq)
library(ape)
library(tidyverse)
library(ggplot2)
library(picante)
library(ggpubr)
library(dplyr)

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

#Desulfobacteriae phylum only
phylo_desulfo <- subset_taxa(phylo_nmnc, Phylum =="p__Thermodesulfobacteriota")

unique_family <- get_taxa_unique(phylo_desulfo, "Family")
unique_family

#Thermodesulfobacteriaceae family only

phylo_desulfo_f <- subset_taxa(phylo_nmnc, Family == "f__Thermodesulfobacteriaceae")
get_taxa_unique(phylo_desulfo_f, "Genus")
#There's only one genus (Thermosulfuriphilus), so this will not be continued with
#The final plots will be focused entirely on families within the Thermodesulfobacteriota phylum

#Filtering taxonomy table to top 10 + other
phylo_desulfo <- subset_taxa(phylo_nmnc, Phylum =="p__Thermodesulfobacteriota")

ps_fam <- tax_glom(phylo_desulfo, taxrank = "Family")

top10 <- names(sort(taxa_sums(ps_fam), TRUE)[1:10])

top10_names <- as.character(tax_table(ps_fam)[top10, "Family"])
top10_names

tax_table(phylo_desulfo_family)[!taxa_names(phylo_desulfo_family) %in% top_10, "Family"] <- "Other"

#Relative abundance
phylo_desulfo_relative <- transform_sample_counts(phylo_desulfo_family, function(x) x / sum(x))

#Plot creation
phylo_plot_1 <- plot_bar(phylo_desulfo_relative, fill = "Family") +
  facet_wrap(~conductivity_category, scales = "free_x") +
  theme_classic() +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack", color = NA) +
  labs(title = "Relative Abundance of Families in Thermodesulfobacteriaceae Phylum",
       subtitle = "Grouped by Sediment Conductivity Category",
       x = "Individual Samples",
       y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 
phylo_plot_1

ggsave("results/taxonomy_family_plot.png", plot = phylo_plot_1, width = 10, height = 5, dpi = 300)


#Filtering taxonomy table to desulfosbulbaceae + other
phylo_desulfo_family <- tax_glom(phylo_desulfo, taxrank = "Family")

desulfo_only <- names(sort(taxa_sums(phylo_desulfo_family), TRUE)[9])
desulfo_only

tax_table(phylo_desulfo_family)[!taxa_names(phylo_desulfo_family) %in% desulfo_only, "Family"] <- "Other"

#Relative abundance
phylo_desulfo_relative <- transform_sample_counts(phylo_desulfo_family, function(x) x / sum(x))

#Plot creation
phylo_plot_2 <- plot_bar(phylo_desulfo_relative, fill = "Family") +
  facet_wrap(~conductivity_category, scales = "free_x") +
  theme_classic() +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack", color = NA) +
  labs(title = "Relative Abundance of Desulfobulbaceae Family in Thermodesulfobacteriaceae Phylum",
    subtitle = "Grouped by Sediment Conductivity Category",
    x = "Individual Samples",
    y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank()) 

phylo_plot_2

ggsave("results/taxonomy_desulfobulbaceae_relative_plot.png", plot = phylo_plot_2, width = 10, height = 5, dpi = 300)

#Wilcoxon box plot
desulfo_df <- phylo_desulfo_relative %>%
  psmelt()

desulfo_wilcox_plot <- ggplot(desulfo_df, aes(x=conductivity_category, y=Abundance, fill=conductivity_category)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width=0.2, alpha = 0.5) + 
  theme_classic() + 
  scale_fill_manual(values = c("high" = "red", "low" = "blue")) + 
  labs(
    title = "Abundance of Desulfobulbaceae Family by Conductivity Category",
    x = "Conductivity Category",
    y = "Relative Abundance", 
    subtitle = "Non-parametric Wilcoxon Rank-Sum Test") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5)

desulfo_wilcox_plot

ggsave("results/desulfobulbaceae_wilcox_plot.png", plot = desulfo_wilcox_plot, width = 7, height = 5, dpi = 300)

