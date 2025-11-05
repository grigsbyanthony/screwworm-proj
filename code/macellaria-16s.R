# Cochliomyia macellaria (Secondary screwworm) adults 16S analysis
# 2025-09-05

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

# Print summary of the filtered dataset
print("Summary of Cochliomyia macellaria samples:")
ps_macellaria
print(sample_variables(ps_macellaria))

#------------------------------------------------------------------------------
# 1. Alpha Diversity Analysis
#------------------------------------------------------------------------------

# Calculate alpha diversity metrics
alpha_div <- microbiome::alpha(ps_macellaria, index = c("shannon", "simpson", "observed"))
alpha_div$SampleID <- rownames(alpha_div)

# Print column names to check the actual names
print("Alpha diversity column names:")
print(colnames(alpha_div))

# Merge with sample data
sample_data_df <- data.frame(sample_data(ps_macellaria))
sample_data_df$SampleID <- rownames(sample_data_df)
alpha_div_merged <- merge(alpha_div, sample_data_df, by = "SampleID")

# Plot alpha diversity (using the correct column names from microbiome::alpha)
p_shannon <- ggplot(alpha_div_merged, aes(x = 1, y = diversity_shannon)) +
  geom_boxplot(width = 0.5, fill = "#A4D3EE") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Shannon Index") +
  facet_wrap(~ Microbiota) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub()

p_observed <- ggplot(alpha_div_merged, aes(x = 1, y = observed)) +
  geom_boxplot(width = 0.5, fill = "#A4D3EE") +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Number of ASVs") +
  facet_wrap(~ Microbiota) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub()

# Combine plots
p_alpha <- p_shannon + p_observed
print(p_alpha)

# Save alpha diversity plot
ggsave("macellaria_alpha_diversity.pdf", p_alpha, width = 10, height = 5)

#------------------------------------------------------------------------------
# 2. Beta Diversity Analysis
#------------------------------------------------------------------------------

# Calculate Bray-Curtis distance
bc_dist <- phyloseq::distance(ps_macellaria, method = "bray")

# Perform PCoA
pcoa <- ordinate(ps_macellaria, method = "PCoA", distance = bc_dist)

# Plot PCoA
p_pcoa <- plot_ordination(ps_macellaria, pcoa, color = "SampleID") +
  geom_point(size = 3) +
  labs(x = "Axis 1", y = "Axis 2") +
  facet_wrap(~ Microbiota) +
  theme_pub()

print(p_pcoa)

# Save beta diversity plot
ggsave("macellaria_beta_diversity_pcoa.pdf", p_pcoa, width = 8, height = 6)

#------------------------------------------------------------------------------
# 3. Taxonomic Composition Analysis
#------------------------------------------------------------------------------

# Transform to relative abundance
ps_macellaria_rel <- transform_sample_counts(ps_macellaria, function(x) x / sum(x))

# Melt the phyloseq object for plotting
ps_melt <- psmelt(ps_macellaria_rel)

# Print column names to check the actual names
print("Melted phyloseq column names:")
print(colnames(ps_melt))

# Get the top 10 most abundant phyla
top_phyla <- ps_melt %>%
  group_by(Phylum) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  top_n(10) %>%
  pull(Phylum)

# Filter for top phyla and group others
ps_melt <- ps_melt %>%
  mutate(Phylum = if_else(Phylum %in% top_phyla, as.character(Phylum), "Other"))

# Plot phylum-level composition
p_phylum <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Relative Abundance") +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_phylum)

# Get the top 20 most abundant genera
top_genera <- ps_melt %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  top_n(20) %>%
  pull(Genus)

# Filter for top genera and group others
ps_melt <- ps_melt %>%
  mutate(Genus = if_else(Genus %in% top_genera, as.character(Genus), "Other"))

# Plot genus-level composition
p_genus <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Relative Abundance") +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_genus)

# Save taxonomic composition plots
ggsave("macellaria_phylum_composition.pdf", p_phylum, width = 10, height = 6)
ggsave("macellaria_genus_composition.pdf", p_genus, width = 10, height = 6)

#------------------------------------------------------------------------------
# 4. Core Microbiome Analysis
#------------------------------------------------------------------------------

# Identify core microbiome (present in at least 50% of samples)
core_microbiome <- core_members(ps_macellaria_rel, detection = 0.001, prevalence = 0.5)
print("Core microbiome ASVs (present in at least 50% of samples):")
print(length(core_microbiome))

# Extract taxonomy of core ASVs
core_taxa <- tax_table(ps_macellaria)[core_microbiome, ]
print("Taxonomy of core ASVs:")
print(core_taxa)

# Create a phyloseq object with only core ASVs
ps_core <- prune_taxa(core_microbiome, ps_macellaria_rel)

# Plot core microbiome composition
ps_core_melt <- psmelt(ps_core)

# Plot genus-level core composition
p_core_genus <- ggplot(ps_core_melt, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Relative Abundance") +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_core_genus)

# Save core microbiome plot
ggsave("macellaria_core_microbiome.pdf", p_core_genus, width = 10, height = 6)

#------------------------------------------------------------------------------
# 5. Assign Providencia to Core Taxa with NA Taxonomic Assignments
#------------------------------------------------------------------------------

# Identify core ASVs with NA in genus-level taxonomy
core_taxa_df <- as.data.frame(core_taxa)
na_genus_asv <- rownames(core_taxa_df)[is.na(core_taxa_df$Genus)]

# Print the number of core ASVs with NA genus
cat("\nNumber of core ASVs with NA genus: ", length(na_genus_asv), "\n")

# If there are ASVs with NA genus, assign them as Providencia
if (length(na_genus_asv) > 0) {
  cat("\nCore ASVs with NA genus (before assignment):\n")
  
  # Extract taxonomy information for ASVs with NA genus
  na_genus_taxa <- core_taxa_df[na_genus_asv, ]
  
  # Print the taxonomy information before assignment
  print(na_genus_taxa)
  
  # Save the original taxonomy information to a file
  write.csv(na_genus_taxa, "macellaria_core_na_genus_taxa_original.csv")
  cat("Original taxonomy information for core ASVs with NA genus saved to: macellaria_core_na_genus_taxa_original.csv\n")
  
  # Assign Providencia to the ASVs with NA genus
  cat("\nAssigning genus 'Providencia' to ASVs with NA genus...\n")
  
  # Get the tax table from the phyloseq object
  tax_table_df <- as.data.frame(tax_table(ps_macellaria))
  
  # Assign Providencia to the ASVs with NA genus
  tax_table_df[na_genus_asv, "Genus"] <- "Providencia"
  
  # Update the tax table in the phyloseq object
  tax_table(ps_macellaria) <- tax_table(as.matrix(tax_table_df))
  
  # Update the core taxa object
  core_taxa <- tax_table(ps_macellaria)[core_microbiome, ]
  
  # Print the updated taxonomy information
  cat("\nUpdated taxonomy for previously NA genus ASVs (now assigned as Providencia):\n")
  print(core_taxa[na_genus_asv, ])
  
  # Save the updated taxonomy information to a file
  write.csv(as.data.frame(core_taxa[na_genus_asv, ]), "macellaria_core_providencia_taxa.csv")
  cat("Updated taxonomy information for core ASVs now assigned as Providencia saved to: macellaria_core_providencia_taxa.csv\n")
  
  # Update the relative abundance object with the new taxonomy
  ps_macellaria_rel <- transform_sample_counts(ps_macellaria, function(x) x / sum(x))
  
  # Update the core microbiome object
  ps_core <- prune_taxa(core_microbiome, ps_macellaria_rel)
  
  # Re-melt the phyloseq object for plotting with updated taxonomy
  ps_core_melt <- psmelt(ps_core)
  
  # Re-plot the core microbiome composition with updated taxonomy
  p_core_genus_updated <- ggplot(ps_core_melt, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity") +
    labs(x = "Sample", y = "Relative Abundance") +
    facet_wrap(~ Microbiota, scales = "free_x") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_core_genus_updated)
  
  # Save the updated core microbiome plot
  ggsave("macellaria_core_microbiome_with_providencia.pdf", p_core_genus_updated, width = 10, height = 6)
  cat("Updated core microbiome plot with Providencia assignment saved to: macellaria_core_microbiome_with_providencia.pdf\n")
} else {
  cat("No core ASVs with NA genus found.\n")
}

# Print summary of analysis
cat("\nAnalysis of Cochliomyia macellaria 16S data completed.\n")
cat("Output files:\n")
cat("- macellaria_alpha_diversity.pdf\n")
cat("- macellaria_beta_diversity_pcoa.pdf\n")
cat("- macellaria_phylum_composition.pdf\n")
cat("- macellaria_genus_composition.pdf\n")
cat("- macellaria_core_microbiome.pdf\n")
if (length(na_genus_asv) > 0) {
  cat("- macellaria_core_na_genus_taxa_original.csv\n")
  cat("- macellaria_core_providencia_taxa.csv\n")
  cat("- macellaria_core_microbiome_with_providencia.pdf\n")
}
