# Cochliomyia macellaria (Secondary screwworm) adults 16S analysis part 2
# 2025-09-16

# Load custom theme
theme_pub <- function(base_size = 14, base_family = "CMU Sans Serif", legend_pos = "bottom") {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_pos,
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.title = element_text(face = "bold"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
} 

# Load QIIME2 artificats as phyloseq object
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Filter to only Cochliomyia macellaria samples
ps_macellaria <- subset_samples(ps, InsectSpecies == "Cochliomyia macellaria")

# Load required packages for analysis
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(microbiome)
library(patchwork)
library(picante)  # For Faith's PD calculation

# Print summary of the filtered dataset
print("Summary of Cochliomyia macellaria samples:")
ps_macellaria
print(sample_variables(ps_macellaria))

#------------------------------------------------------------------------------
# Comparison of Endogenous vs Exogenous Microbiota in Cochliomyia macellaria
#------------------------------------------------------------------------------

# Check if we have both endogenous and exogenous samples
microbiota_types <- unique(sample_data(ps_macellaria)$Microbiota)
print("Microbiota types in the dataset:")
print(microbiota_types)

# Create separate phyloseq objects for endogenous and exogenous samples
ps_endo <- subset_samples(ps_macellaria, Microbiota == "Endogenous")
ps_exo <- subset_samples(ps_macellaria, Microbiota == "Exogenous")

# Print summary of each dataset
print("Summary of Endogenous samples:")
ps_endo
print("Summary of Exogenous samples:")
ps_exo

#------------------------------------------------------------------------------
# 1. Alpha Diversity Comparison
#------------------------------------------------------------------------------

# Calculate alpha diversity metrics
alpha_div <- microbiome::alpha(ps_macellaria, index = c("shannon", "observed"))
alpha_div$SampleID <- rownames(alpha_div)

# Calculate Faith's Phylogenetic Diversity (Faith's PD)
# Extract the OTU table and convert to presence/absence
otu_table_df <- as.data.frame(t(otu_table(ps_macellaria)))
otu_pa <- ifelse(otu_table_df > 0, 1, 0)

# Extract the phylogenetic tree
tree <- phy_tree(ps_macellaria)

# Calculate Faith's PD
# Make sure sample names are rownames and ASV/OTU names are column names
pd_result <- picante::pd(otu_pa, tree, include.root = FALSE)

# Add SampleID column for merging
pd_result$SampleID <- rownames(pd_result)

# Merge with other alpha diversity metrics
alpha_div <- merge(alpha_div, pd_result, by = "SampleID")

# Print the first few rows of the merged data to verify
print("First few rows of alpha diversity metrics including Faith's PD:")
print(head(alpha_div))

# Merge with sample data
sample_data_df <- data.frame(sample_data(ps_macellaria))
sample_data_df$SampleID <- rownames(sample_data_df)
alpha_div_merged <- merge(alpha_div, sample_data_df, by = "SampleID")

# Statistical comparison of alpha diversity between endogenous and exogenous
# Shannon diversity
shannon_test <- wilcox.test(diversity_shannon ~ Microbiota, data = alpha_div_merged)
print("Wilcoxon test for Shannon diversity between Endogenous and Exogenous:")
print(shannon_test)

# Observed ASVs
observed_test <- wilcox.test(observed ~ Microbiota, data = alpha_div_merged)
print("Wilcoxon test for Observed ASVs between Endogenous and Exogenous:")
print(observed_test)

# Faith's PD
faith_pd_test <- wilcox.test(PD ~ Microbiota, data = alpha_div_merged)
print("Wilcoxon test for Faith's Phylogenetic Diversity between Endogenous and Exogenous:")
print(faith_pd_test)

# Plot alpha diversity metrics
p_shannon <- ggplot(alpha_div_merged, aes(x = Microbiota, y = diversity_shannon, fill = Microbiota)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Shannon Diversity Index") +
  scale_fill_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  theme_pub() +
  theme(legend.position = "none")

p_observed <- ggplot(alpha_div_merged, aes(x = Microbiota, y = observed, fill = Microbiota)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Observed ASVs") +
  scale_fill_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  theme_pub() +
  theme(legend.position = "none")

p_faith_pd <- ggplot(alpha_div_merged, aes(x = Microbiota, y = PD, fill = Microbiota)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Faith's Phylogenetic Diversity") +
  scale_fill_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  theme_pub() +
  theme(legend.position = "none")

# Combine plots
p_alpha_combined <- p_shannon + p_observed + p_faith_pd + plot_layout(ncol = 3)
print(p_alpha_combined)

# Save alpha diversity comparison plot
ggsave("macellaria_endo_exo_alpha_diversity.png", p_alpha_combined, width = 15, height = 9)

#------------------------------------------------------------------------------
# 2. Beta Diversity Comparison
#------------------------------------------------------------------------------

# Calculate Bray-Curtis distance
bc_dist <- phyloseq::distance(ps_macellaria, method = "bray")

# Perform PCoA
pcoa <- ordinate(ps_macellaria, method = "PCoA", distance = bc_dist)

# Extract the variance explained by the first two axes
variance_explained <- pcoa$values$Relative_eig[1:2] * 100

# Plot PCoA
p_pcoa <- plot_ordination(ps_macellaria, pcoa, color = "Microbiota") +
  geom_point(size = 7, alpha = 0.8) +
  stat_ellipse(aes(group = Microbiota), type = "t", linetype = 2) +
  labs(
    x = paste0("PCoA1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PCoA2 (", round(variance_explained[2], 1), "%)"),
    title = "Bray-Curtis PCoA: Endogenous vs Exogenous Microbiota"
  ) +
  scale_color_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  theme_pub()

print(p_pcoa)

# Save beta diversity plot
ggsave("macellaria_endo_exo_beta_diversity_pcoa.png", p_pcoa, width = 15, height = 9)

# Perform PERMANOVA to test for significant differences
metadata_df <- data.frame(sample_data(ps_macellaria))
permanova_result <- adonis2(bc_dist ~ Microbiota, data = metadata_df)
print("PERMANOVA results for Bray-Curtis distances between Endogenous and Exogenous:")
print(permanova_result)

# Test for homogeneity of dispersion (assumption of PERMANOVA)
beta_dispersion <- betadisper(bc_dist, metadata_df$Microbiota)
beta_dispersion_test <- permutest(beta_dispersion)
print("Test for homogeneity of dispersion:")
print(beta_dispersion_test)

#------------------------------------------------------------------------------
# 3. Genus-level Relative Abundance Comparison
#------------------------------------------------------------------------------

# Transform to relative abundance
ps_macellaria_rel <- transform_sample_counts(ps_macellaria, function(x) x / sum(x))

# Melt the phyloseq object for plotting
ps_melt <- psmelt(ps_macellaria_rel)

# Get the top 15 most abundant genera
top_genera <- ps_melt %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  top_n(15) %>%
  pull(Genus)

# Filter for top genera and group others
ps_melt <- ps_melt %>%
  mutate(Genus = if_else(Genus %in% top_genera, as.character(Genus), "Other"))

# Calculate mean abundance by genus and microbiota type
genus_abundance_summary <- ps_melt %>%
  group_by(Microbiota, Genus) %>%
  summarise(
    mean_abundance = mean(Abundance),
    se_abundance = sd(Abundance) / sqrt(n())
  ) %>%
  ungroup() %>%
  arrange(Microbiota, desc(mean_abundance))

# Reorder genera by overall abundance
genus_order <- genus_abundance_summary %>%
  group_by(Genus) %>%
  summarise(overall_abundance = mean(mean_abundance)) %>%
  arrange(desc(overall_abundance)) %>%
  pull(Genus)

# Set the factor levels for genus based on the order
genus_abundance_summary$Genus <- factor(genus_abundance_summary$Genus, levels = genus_order)

# Plot mean relative abundance by genus and microbiota type
p_genus_bar <- ggplot(genus_abundance_summary, aes(x = Genus, y = mean_abundance, fill = Microbiota)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_abundance - se_abundance, ymax = mean_abundance + se_abundance),
    position = position_dodge(width = 0.8), width = 0.25
  ) +
  labs(
    x = "",
    y = "Mean Relative Abundance",
    title = "Genus-level Abundance: Endogenous vs Exogenous"
  ) +
  scale_fill_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_genus_bar)

# Save genus-level abundance comparison plot
ggsave("macellaria_endo_exo_genus_abundance.png", p_genus_bar, width = 12, height = 6)

# Create heatmaps of mean relative abundance by genus and microbiota type
# Use the genus_abundance_summary data frame we created earlier
# Ensure genera are ordered by overall abundance
genus_heatmap_data <- genus_abundance_summary %>%
  select(Microbiota, Genus, mean_abundance) %>%
  # Ensure Genus is ordered by overall abundance
  mutate(Genus = factor(Genus, levels = genus_order))

# Create log-transformed version of the heatmap
genus_heatmap_log <- genus_heatmap_data %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = Genus, y = Microbiota, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Genus",
    y = "Microbiota",
    title = "Mean Relative Abundance by Genus and Microbiota Type - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12)
  )

# Create capped color scale version of the heatmap
# Calculate 95th percentile for capping
cap_value <- quantile(genus_heatmap_data$mean_abundance, 0.95)

genus_heatmap_capped <- genus_heatmap_data %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value)) %>%
  ggplot(aes(x = Genus, y = Microbiota, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Genus",
    y = "Microbiota",
    title = "Mean Relative Abundance by Genus and Microbiota Type - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12)
  )

# Display heatmaps
print(genus_heatmap_log)
print(genus_heatmap_capped)

# Save genus-level heatmaps
ggsave("macellaria_endo_exo_genus_heatmap_log.png", genus_heatmap_log, width = 12, height = 6)
ggsave("macellaria_endo_exo_genus_heatmap_capped.png", genus_heatmap_capped, width = 12, height = 6)

# Create a bubble plot of mean relative abundance by genus and microbiota type
# Similar to the heatmap but using bubbles instead of tiles
p_genus_bubble <- ggplot(genus_heatmap_data, aes(x = Microbiota, y = Genus, size = mean_abundance, color = mean_abundance)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 10), name = "Mean Relative Abundance") +
  scale_color_viridis_c(option = "plasma", name = "Mean Relative Abundance") +
  # Combine size and color legends
  guides(color = guide_legend(override.aes = list(size = 5)), size = guide_legend()) +
  labs(
    x = "Microbiota",
    y = "Genus"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10)
  )

print(p_genus_bubble)

# Save genus-level bubble plot
ggsave("macellaria_endo_exo_genus_bubble.png", p_genus_bubble, width = 12, height = 6)

# Create a bubble plot with sample-level data
# Get the original melted data with individual samples
sample_genus_data <- ps_melt %>%
  filter(Genus %in% top_genera) %>%
  # Ensure Genus is ordered by overall abundance
  mutate(Genus = factor(Genus, levels = genus_order))

# Create bubble plot with samples on x-axis, genus on y-axis, faceted by microbiota
p_sample_genus_bubble <- ggplot(sample_genus_data, aes(x = Sample, y = Genus, size = Abundance, color = Abundance)) +
  geom_point(alpha = 0.7) +
  facet_grid(. ~ Microbiota, scales = "free_x", space = "free_x") +
  scale_size_continuous(range = c(0.5, 8), name = "Relative Abundance") +
  scale_color_viridis_c(option = "plasma", name = "Relative Abundance") +
  # Combine size and color legends
  guides(color = guide_legend(override.aes = list(size = 5)), size = guide_legend()) +
  labs(
    x = "Sample",
    y = "Genus",
    title = "Sample-level Genus Abundance by Microbiota Type"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "#f0f0f0", size = 0.2)
  )

print(p_sample_genus_bubble)

# Save sample-level genus bubble plot
ggsave("macellaria_endo_exo_sample_genus_bubble.png", p_sample_genus_bubble, width = 14, height = 8)

#------------------------------------------------------------------------------
# 4. Differential Abundance Analysis
#------------------------------------------------------------------------------

# Identify differentially abundant genera between endogenous and exogenous samples
# Using ANCOM-BC2 for differential abundance analysis

# Load ANCOMBC library
library(ANCOMBC)

# Agglomerate taxa at genus level
ps_genus <- tax_glom(ps_macellaria, taxrank = "Genus")

# Run ANCOM-BC2 analysis
cat("\nRunning ANCOM-BC2 for Microbiota (Endogenous vs Exogenous) at Genus level...\n")

# Set seed for reproducibility
set.seed(123)

# Run ANCOM-BC2 following the tutorial approach
ancombc_result <- ancombc2(
  data = ps_genus,
  tax_level = "Genus",
  fix_formula = "Microbiota",
  p_adj_method = "holm",
  pseudo_sens = TRUE,  # Enable sensitivity analysis for pseudo-count addition
  prv_cut = 0.10,      # Filter features with low prevalence
  lib_cut = 1000,      # Minimum library size
  s0_perc = 0.05,      # Percentile for variance regularization
  group = "Microbiota",
  struc_zero = TRUE,   # Detect and handle structural zeros
  neg_lb = TRUE,       # Use negative lower bound for log-fold change
  alpha = 0.05,        # Significance level
  global = FALSE,      # No global test needed for two-group comparison
  pairwise = FALSE     # No pairwise test needed for two-group comparison
)

# Extract results
res_prim <- ancombc_result$res

# Create a data frame for visualization
df_microbiota <- res_prim %>%
  dplyr::select(taxon, contains("Microbiota")) 

# Filter for significant results
df_sig <- df_microbiota %>%
  dplyr::filter(diff_MicrobiotaExogenous == 1) %>%
  dplyr::arrange(desc(lfc_MicrobiotaExogenous))

# Print summary
if (nrow(df_sig) > 0) {
  cat(paste0("\nFound ", nrow(df_sig), " differentially abundant genera between Endogenous and Exogenous\n"))
  
  # Save significant results to CSV
  write.csv(df_sig, "ancombc2_Genus_macellaria_MicrobiotaExogenous.csv", row.names = FALSE)
  
  # Get taxonomy information for plotting
  tax_info <- as.data.frame(tax_table(ps_genus))
  tax_info$taxon <- rownames(tax_info)  # Add row names as a column
  
  # Create a data frame for plotting
  plot_data <- df_sig %>%
    dplyr::left_join(tax_info, by = "taxon") %>%
    dplyr::mutate(
      # Use taxon ID as the plot name (guaranteed to be unique)
      PlotName = taxon,
      # Add genus name for display
      DisplayName = ifelse(is.na(Genus), "Unknown", Genus),
      # Replace underscores with spaces
      DisplayName = gsub("_", " ", DisplayName),
      # Create a combined name for display
      PlotName = paste0(DisplayName, " (", taxon, ")")
    )
  
  # Order the data by log fold change
  plot_data <- plot_data %>%
    dplyr::arrange(lfc_MicrobiotaExogenous)
  
  # Create factor with ordered levels
  plot_data$PlotName <- factor(plot_data$PlotName, levels = unique(plot_data$PlotName))
  
  # Create bar plot of log fold changes
  p_sig_genera <- ggplot(plot_data, aes(x = PlotName, y = lfc_MicrobiotaExogenous, fill = lfc_MicrobiotaExogenous > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(
      values = c("TRUE" = "#1B9E77", "FALSE" = "#D95F02"),
      labels = c("TRUE" = "Higher in Exogenous", "FALSE" = "Higher in Endogenous"),
      name = "Abundance"
    ) +
    labs(
      x = "Genus",
      y = "Log Fold Change",
      title = "ANCOM-BC2: Differentially Abundant Genera",
      subtitle = "Exogenous vs Endogenous Microbiota"
    ) +
    theme_pub() +
    theme(axis.text.y = element_text(size = 10))
  
  # Display and save bar plot
  print(p_sig_genera)
  ggsave("macellaria_endo_exo_ancombc_Genus.png", p_sig_genera, width = 10, height = 8)
  
  # Create volcano plot
  # Check if sensitivity analysis results are available
  has_sensitivity <- "diff_robust_MicrobiotaExogenous" %in% colnames(df_microbiota)
  
  volcano_data <- df_microbiota %>%
    dplyr::mutate(
      Significant = diff_MicrobiotaExogenous == 1,
      # Only use sensitivity results if available
      SensitivityPass = if(has_sensitivity) {
        diff_robust_MicrobiotaExogenous == 1
      } else {
        # If sensitivity results not available, set all to TRUE
        TRUE
      },
      direction = case_when(
        Significant & lfc_MicrobiotaExogenous > 0 ~ "Higher in Exogenous",
        Significant & lfc_MicrobiotaExogenous < 0 ~ "Higher in Endogenous",
        TRUE ~ "Not Significant"
      )
    ) %>%
    dplyr::left_join(tax_info, by = "taxon") %>%
    dplyr::mutate(
      label = ifelse(Significant, Genus, NA)
    )
  
  # Create the volcano plot
  p_volcano <- ggplot(volcano_data, aes(x = lfc_MicrobiotaExogenous, y = -log10(q_MicrobiotaExogenous), color = direction)) +
    geom_point(aes(shape = SensitivityPass), alpha = 0.7, size = 3) +
    # Only include shape legend if sensitivity analysis was performed
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1), 
      name = "Passed Sensitivity",
      drop = FALSE,
      guide = if(has_sensitivity) guide_legend() else "none"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray") +
    scale_color_manual(values = c(
      "Higher in Exogenous" = "#1B9E77",
      "Higher in Endogenous" = "#D95F02",
      "Not Significant" = "gray80"
    )) +
    labs(
      x = "Log Fold Change",
      y = "-Log10 Adjusted P-value",
      title = "ANCOM-BC2: Differential Abundance",
      subtitle = "Exogenous vs Endogenous Microbiota"
    ) +
    ggrepel::geom_text_repel(
      aes(label = label),
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "gray50"
    ) +
    theme_pub()
  
  # Display and save volcano plot
  print(p_volcano)
  ggsave("macellaria_endo_exo_ancombc_volcano_Genus.png", p_volcano, width = 10, height = 8)
  
  # Create a heatmap of significant genera
  if (nrow(plot_data) > 1) {
    # Get bias-corrected abundances
    bias_correct_log_table <- ancombc_result$bias_correct_log_table
    
    # Replace NAs with zeros (equivalent to adding pseudo-count of 1 to zero counts)
    bias_correct_log_table[is.na(bias_correct_log_table)] <- 0
    
    # Extract only significant taxa
    sig_taxa_names <- plot_data$taxon
    bias_correct_sig <- bias_correct_log_table[rownames(bias_correct_log_table) %in% sig_taxa_names, ]
    
    # Prepare data for heatmap
    heatmap_data <- as.data.frame(t(bias_correct_sig))
    heatmap_data$Sample <- rownames(heatmap_data)
    
    # Add metadata
    sample_data_df <- data.frame(sample_data(ps_genus))
    sample_data_df$Sample <- rownames(sample_data_df)
    heatmap_data <- merge(heatmap_data, sample_data_df[, c("Sample", "Microbiota")], by = "Sample")
    
    # Reshape for ggplot
    heatmap_long <- heatmap_data %>%
      tidyr::pivot_longer(cols = -c(Sample, Microbiota), names_to = "Taxon", values_to = "Abundance")
    
    # Map taxa names to genus names
    taxon_to_genus <- plot_data %>% 
      dplyr::select(taxon, PlotName) %>% 
      dplyr::rename(Taxon = taxon, Genus = PlotName)
    
    heatmap_long <- heatmap_long %>%
      dplyr::left_join(taxon_to_genus, by = "Taxon")
    
    # Create heatmap
    p_heatmap <- ggplot(heatmap_long, aes(x = Sample, y = Genus, fill = Abundance)) +
      geom_tile() +
      facet_grid(. ~ Microbiota, scales = "free_x", space = "free_x") +
      scale_fill_viridis_c(option = "plasma", name = "Bias-corrected\nLog Abundance") +
      labs(
        x = "",
        y = "",
        title = "Bias-corrected Abundances of Differentially Abundant Genera"
      ) +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 10)
      )
    
    # Display and save heatmap
    print(p_heatmap)
    ggsave("macellaria_endo_exo_ancombc_heatmap_Genus.png", p_heatmap, width = 12, height = 8)
  }
  
  # Create a very basic heatmap of top 15 most differentially abundant taxa
  # Get the top 15 taxa by absolute log fold change
  top_diff_taxa <- df_microbiota %>%
    dplyr::mutate(abs_lfc = abs(lfc_MicrobiotaExogenous)) %>%
    dplyr::arrange(desc(abs_lfc)) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::pull(taxon)
  
  # Get taxonomy information for these taxa
  top_tax_info <- tax_info %>%
    dplyr::filter(taxon %in% top_diff_taxa) %>%
    dplyr::mutate(
      # Create a clean display name for each taxon using the most specific taxonomic level available
      DisplayName = case_when(
        !is.na(Genus) ~ paste0(gsub("_", " ", Genus), " (", taxon, ")"),
        !is.na(Family) ~ paste0("Family ", gsub("_", " ", Family), " (", taxon, ")"),
        !is.na(Order) ~ paste0("Order ", gsub("_", " ", Order), " (", taxon, ")"),
        !is.na(Class) ~ paste0("Class ", gsub("_", " ", Class), " (", taxon, ")"),
        !is.na(Phylum) ~ paste0("Phylum ", gsub("_", " ", Phylum), " (", taxon, ")"),
        TRUE ~ paste0("Unknown (", taxon, ")")
      ),
      # Add significance information
      Significant = taxon %in% df_sig$taxon
    )
  
  # Get bias-corrected abundances for these taxa
  bias_correct_log_table <- ancombc_result$bias_correct_log_table
  bias_correct_log_table[is.na(bias_correct_log_table)] <- 0
  
  # Extract the subset of the bias-corrected abundance matrix for the top taxa
  bias_correct_top <- bias_correct_log_table[rownames(bias_correct_log_table) %in% top_diff_taxa, ]
  
  # Create a data frame with sample metadata
  sample_data_df <- data.frame(sample_data(ps_genus))
  sample_data_df$Sample <- rownames(sample_data_df)
  
  # Create a simple plot showing the top taxa as a bar plot
  # This is a fallback in case the heatmap doesn't work
  top_taxa_lfc <- df_microbiota %>%
    dplyr::filter(taxon %in% top_diff_taxa) %>%
    dplyr::arrange(desc(abs(lfc_MicrobiotaExogenous))) %>%
    dplyr::left_join(top_tax_info, by = "taxon") %>%
    dplyr::mutate(
      Direction = ifelse(lfc_MicrobiotaExogenous > 0, "Higher in Exogenous", "Higher in Endogenous"),
      DisplayName = factor(DisplayName, levels = DisplayName[order(lfc_MicrobiotaExogenous)])
    )
  
  # Create a bar plot of log fold changes for the top 15 taxa
  p_top_lfc <- ggplot(top_taxa_lfc, aes(x = DisplayName, y = lfc_MicrobiotaExogenous, fill = Direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(
      values = c("Higher in Exogenous" = "#1B9E77", "Higher in Endogenous" = "#D95F02")
    ) +
    labs(
      x = "",
      y = "Log Fold Change",
      title = "Top 15 Most Differentially Abundant Genera"
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 10)
    )
  
  # Display and save the bar plot
  print(p_top_lfc)
  ggsave("macellaria_endo_exo_ancombc_top15_lfc_Genus.png", p_top_lfc, width = 10, height = 8)
  
  # Try to create a simple heatmap with base R
  # This is a fallback in case ggplot2 heatmap doesn't work
  png("macellaria_endo_exo_ancombc_top15_heatmap_Genus.png", width = 12, height = 8, units = "in", res = 300)
  
  # Get display names for taxa
  row_names <- sapply(rownames(bias_correct_top), function(taxon) {
    info <- top_tax_info[top_tax_info$taxon == taxon, ]
    if (nrow(info) > 0) {
      return(info$DisplayName[1])
    } else {
      return(paste0("Unknown (", taxon, ")"))
    }
  })
  
  # Split samples by microbiota type
  endo_samples <- rownames(sample_data_df)[sample_data_df$Microbiota == "Endogenous"]
  exo_samples <- rownames(sample_data_df)[sample_data_df$Microbiota == "Exogenous"]
  
  # Create separate matrices
  endo_matrix <- bias_correct_top[, colnames(bias_correct_top) %in% endo_samples]
  exo_matrix <- bias_correct_top[, colnames(bias_correct_top) %in% exo_samples]
  
  # Set up the plot layout
  layout(matrix(c(1, 2), nrow = 1), widths = c(1, 1))
  
  # Plot endogenous samples
  heatmap(
    as.matrix(endo_matrix),
    Rowv = NA,  # Don't cluster rows
    Colv = NA,  # Don't cluster columns
    scale = "none",  # Don't scale the data
    col = viridis::viridis(100, option = "plasma"),  # Use viridis color palette
    main = "Endogenous Microbiota",
    xlab = "",
    ylab = "",
    labRow = row_names,
    margins = c(8, 10)  # Adjust margins for labels
  )
  
  # Plot exogenous samples
  heatmap(
    as.matrix(exo_matrix),
    Rowv = NA,  # Don't cluster rows
    Colv = NA,  # Don't cluster columns
    scale = "none",  # Don't scale the data
    col = viridis::viridis(100, option = "plasma"),  # Use viridis color palette
    main = "Exogenous Microbiota",
    xlab = "",
    ylab = "",
    labRow = row_names,
    margins = c(8, 10)  # Adjust margins for labels
  )
  
  dev.off()
  
  # Create a simple bar plot of log fold changes
  # This approach doesn't rely on the abundance data and should be more robust
  
  # Create a data frame with the top 15 taxa and their log fold changes
  top_taxa_df <- df_microbiota %>%
    dplyr::filter(taxon %in% top_diff_taxa) %>%
    dplyr::mutate(
      Significant = taxon %in% df_sig$taxon,
      Direction = ifelse(lfc_MicrobiotaExogenous > 0, "Higher in Exogenous", "Higher in Endogenous")
    ) %>%
    dplyr::left_join(tax_info, by = "taxon")
  
  # Just use the taxon ID as the display name
  top_taxa_df <- top_taxa_df %>%
    dplyr::mutate(
      DisplayName = taxon
    )
  
  # Order by log fold change
  top_taxa_df <- top_taxa_df %>%
    dplyr::arrange(lfc_MicrobiotaExogenous) %>%
    dplyr::mutate(
      DisplayName = factor(DisplayName, levels = DisplayName)
    )
  
  # Create a bar plot of log fold changes
  p_lfc_bar <- ggplot(top_taxa_df, aes(x = DisplayName, y = lfc_MicrobiotaExogenous, fill = Direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(
      values = c("Higher in Exogenous" = "#1B9E77", "Higher in Endogenous" = "#D95F02")
    ) +
    # Add a border to significant taxa
    geom_bar(
      data = top_taxa_df %>% filter(Significant),
      stat = "identity",
      fill = NA,
      color = "black",
      size = 1
    ) +
    # Set y-axis limits from -2 to 2
    scale_y_continuous(limits = c(-2, 2)) +
    labs(
      x = "",
      y = "Log Fold Change",
      title = "Top 15 Most Differentially Abundant Taxa"
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 10)
    )
  
  # Save the bar plot
  ggsave("macellaria_endo_exo_ancombc_top15_lfc_bar_Genus.png", p_lfc_bar, width = 10, height = 8)
  
  # Create a horizontal dot plot of log fold changes
  p_lfc_dot <- ggplot(top_taxa_df, aes(x = lfc_MicrobiotaExogenous, y = DisplayName, color = Direction)) +
    geom_point(size = 4) +
    # Add a border to significant taxa
    geom_point(
      data = top_taxa_df %>% filter(Significant),
      shape = 21,
      fill = NA,
      color = "black",
      size = 8,
      stroke = 1
    ) +
    scale_color_manual(
      values = c("Higher in Exogenous" = "#1B9E77", "Higher in Endogenous" = "#D95F02")
    ) +
    # Set x-axis limits from -2 to 2 (same as the bar plot)
    scale_x_continuous(limits = c(-2, 2)) +
    labs(
      x = "Log Fold Change",
      y = "",
      title = "Top 15 Most Differentially Abundant Taxa"
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 10)
    )
  
  # Save the dot plot
  ggsave("macellaria_endo_exo_ancombc_top15_lfc_dot_Genus.png", p_lfc_dot, width = 10, height = 8)
  
  # Store results for summary
  genus_results <- list(
    sig_taxa = df_sig,
    plot_data = plot_data,
    p_sig_genera = p_sig_genera,
    p_volcano = p_volcano
  )
} else {
  cat("\nNo differentially abundant genera found between Endogenous and Exogenous using ANCOM-BC2.\n")
  genus_results <- NULL
}

# Print summary of the comparison
cat("\nComparison of Endogenous vs Exogenous Microbiota in Cochliomyia macellaria completed.\n")
cat("Output files:\n")
cat("- macellaria_endo_exo_alpha_diversity.png\n")
cat("- macellaria_endo_exo_beta_diversity_pcoa.png\n")
cat("- macellaria_endo_exo_genus_abundance.png\n")
cat("- macellaria_endo_exo_genus_heatmap_log.png\n")
cat("- macellaria_endo_exo_genus_heatmap_capped.png\n")
cat("- macellaria_endo_exo_genus_bubble.png\n")
cat("- macellaria_endo_exo_sample_genus_bubble.png\n")
if (!is.null(genus_results)) {
  cat("- ancombc2_Genus_macellaria_MicrobiotaExogenous.csv\n")
  cat("- macellaria_endo_exo_ancombc_Genus.png\n")
  cat("- macellaria_endo_exo_ancombc_volcano_Genus.png\n")
  cat("- macellaria_endo_exo_ancombc_heatmap_Genus.png\n")
  cat("- macellaria_endo_exo_ancombc_top15_heatmap_Genus.png\n")
}
