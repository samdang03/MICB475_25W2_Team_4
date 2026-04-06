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

#The final plots will be focused entirely on families within the Thermodesulfobacteriota phylum

#Filtering taxonomy table to top 10 + other
phylo_desulfo <- subset_taxa(phylo_nmnc, Phylum =="p__Thermodesulfobacteriota")

ps_fam <- tax_glom(phylo_desulfo, taxrank = "Family")

top10 <- names(sort(taxa_sums(ps_fam), TRUE)[1:10])

top10_names <- as.character(tax_table(ps_fam)[top10, "Family"])
top10_names

tax_table(phylo_desulfo_family)[!taxa_names(phylo_desulfo_family) %in% top_10, "Family"] <- "Other"

#There was some interesting correlation with Geobacter, so here are new graphs focusing on those
#Geobacter abundance plot

phylo_desulfo_family <- tax_glom(phylo_desulfo, taxrank = "Family")

geo_only <- names(sort(taxa_sums(phylo_desulfo_family), TRUE)[1])
geo_only

tax_table(phylo_desulfo_family)[!taxa_names(phylo_desulfo_family) %in% geo_only, "Family"] <- "Other"

#Relative abundance
phylo_geo_relative <- transform_sample_counts(phylo_desulfo_family, function(x) x / sum(x))

#Plot creation
phylo_plot_3 <- plot_bar(phylo_geo_relative, fill = "Family") +
  facet_wrap(~conductivity_category, scales = "free_x") +
  theme_classic() +
  geom_bar(aes(fill = Family), stat = "identity", position = "stack", color = NA) +
  labs(title = "Relative Abundance of Desulfobulbaceae Family in Thermodesulfobacteriaceae Phylum",
       subtitle = "Grouped by Sediment Conductivity Category",
       x = "Individual Samples",
       y = "Relative Abundance (%)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

phylo_plot_3

ggsave("results/taxonomy_geobacter_relative_plot.png", plot = phylo_plot_3, width = 10, height = 5, dpi = 300)

#Conversion to relative abundance by sample for Wilcoxon box plot
tax_df <- as.data.frame(as(tax_table(phylo_desulfo), "matrix"))

otu_ids <- rownames(otu_table(phylo_desulfo))

geo_asvs <- rownames(tax_df)[grepl("Geobacteraceae", tax_df$Family, ignore.case = TRUE)]
geo_asvs <- intersect(geo_asvs, otu_ids)

otu_mat <- as.matrix(otu_table(phylo_desulfo))
if (!taxa_are_rows(phylo_desulfo)) {otu_mat <- t(otu_mat) }

target_counts <- otu_mat[geo_asvs, , drop=FALSE]
sample_totals <- colSums(otu_mat)
target_totals <- colSums(target_counts)

geo_stats <- data.frame(
  SampleID = names(target_totals),
  RelativeAbundance = target_totals / sample_totals)

meta_df <- as.data.frame(as(sample_data(phylo_desulfo), "matrix"))
meta_df$SampleID <- rownames(meta_df)
geo_stats <- merge(geo_stats, meta_df, by = "SampleID")

colnames(geo_stats) <- make.unique(names(geo_stats))

#Wilcoxon box plot

geo_wilcox_plot <- ggplot(geo_stats, aes(x = conductivity_category, y = RelativeAbundance, fill = conductivity_category)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_classic() + 
  scale_fill_manual(values = c("high" = "#0072B2", "low" = "#D62728")) + 
  labs(title = "Abundance of Geobacteraceae Family by Conductivity Category",
       subtitle = "Non-parametric Wilcoxon Rank-Sum Test",
       x = "Conductivity Category", 
       y = "Relative Abundance") + 
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) 

geo_wilcox_plot

ggsave("results/geobacter_wilcox_plot.png", plot = geo_wilcox_plot, width = 7, height = 5, dpi = 300)

