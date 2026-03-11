#### 1. Data Transformation ####
# Transform absolute counts to relative abundance (0 to 1) 
# This standardizes samples with different sequencing depths
ps_rel <- transform_sample_counts(ps_rare, function(x) x / sum(x))

#### 2. Visualization: Relative Abundance Bar Plot ####
# Creating a publication-quality bar plot faceted by Conductivity
phy_plot <- plot_bar(ps_rel, fill = "Phylum") +
  facet_wrap(~conductivity_category, scales = "free_x") +
  theme_classic() +
  # Remove the black outlines around every tiny OTU box for a cleaner look
  geom_bar(aes(fill = Phylum), stat = "identity", position = "stack", color = NA) +
  # Professional labels
  labs(
    title = "Relative Abundance of Microbial Phyla",
    subtitle = "Grouped by Sediment Conductivity Category",
    x = "Individual Samples",
    y = "Relative Abundance (%)",
    caption = "Data rarefied to 1000 reads; non-bacterial taxa filtered."
  ) +
  # Clean up X-axis since sample names are often too crowded
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgrey") # Highlights facet headers
  )

# Display the plot
print(phy_plot)

#### 3. Save Output ####
# Good practice to save a high-res version for your repo
ggsave("results/taxonomy_barplot_relative.png", plot = phy_plot, width = 12, height = 8, dpi = 300)

