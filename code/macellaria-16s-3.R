# Comparison of Cochliomyia macellaria to other fly species 16S analysis
# 2025-09-17

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

# Load required packages for analysis
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(microbiome)
library(patchwork)
library(ANCOMBC)
library(viridis)
library(ggsignif)  # For adding significance bars to plots

# Print summary of the dataset
print("Summary of all samples:")
ps
print(sample_variables(ps))

# Check what fly species we have in the dataset
fly_species <- unique(sample_data(ps)$InsectSpecies)
print("Fly species in the dataset:")
print(fly_species)

# Check if we have both endogenous and exogenous samples
microbiota_types <- unique(sample_data(ps)$Microbiota)
print("Microbiota types in the dataset:")
print(microbiota_types)

# Print sample counts by species and microbiota type
sample_counts <- table(sample_data(ps)$InsectSpecies, sample_data(ps)$Microbiota)
print("Sample counts by species and microbiota type:")
print(sample_counts)

# Filter out species with only one sample
species_sample_counts <- table(sample_data(ps)$InsectSpecies)
species_to_keep <- names(species_sample_counts[species_sample_counts > 1])
ps <- subset_samples(ps, InsectSpecies %in% species_to_keep)

# Filter out species with no Exogenous samples
exo_sample_counts <- table(sample_data(ps)[sample_data(ps)$Microbiota == "Exogenous", ]$InsectSpecies)
species_with_exo <- names(exo_sample_counts[exo_sample_counts > 0])
ps <- subset_samples(ps, InsectSpecies %in% species_with_exo)

# Filter out Hydrotaea aenescens and Necrodes surinamensis
ps <- subset_samples(ps, InsectSpecies != "Hydrotaea aenescens")

ps <- subset_samples(ps, InsectSpecies != "Necrodes surinamensis")


# Print updated sample counts after filtering
updated_sample_counts <- table(sample_data(ps)$InsectSpecies, sample_data(ps)$Microbiota)
print("Updated sample counts after filtering:")
print(updated_sample_counts)

# Print summary of filtered dataset
print("Summary of filtered dataset:")
print(paste("Number of species:", length(unique(sample_data(ps)$InsectSpecies))))
print(paste("Total number of samples:", nsamples(ps)))
print(paste("Number of ASVs:", ntaxa(ps)))



#------------------------------------------------------------------------------
# 1. Alpha Diversity Comparison Across Fly Species, Faceted by Microbiota
#------------------------------------------------------------------------------

# Calculate alpha diversity metrics
alpha_div <- microbiome::alpha(ps, index = c("shannon", "observed"))
alpha_div$SampleID <- rownames(alpha_div)

# Calculate Faith's Phylogenetic Diversity (Faith's PD)
# Extract the OTU table and convert to presence/absence
otu_table_df <- as.data.frame(t(otu_table(ps)))
otu_pa <- ifelse(otu_table_df > 0, 1, 0)

# Extract the phylogenetic tree
tree <- phy_tree(ps)

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
sample_data_df <- data.frame(sample_data(ps))
sample_data_df$SampleID <- rownames(sample_data_df)
alpha_div_merged <- merge(alpha_div, sample_data_df, by = "SampleID")

# Statistical comparison of alpha diversity between species
# Shannon diversity
shannon_test_endo <- kruskal.test(diversity_shannon ~ InsectSpecies, 
                                  data = alpha_div_merged[alpha_div_merged$Microbiota == "Endogenous", ])
print("Kruskal-Wallis test for Shannon diversity between species (Endogenous):")
print(shannon_test_endo)

shannon_test_exo <- kruskal.test(diversity_shannon ~ InsectSpecies, 
                                 data = alpha_div_merged[alpha_div_merged$Microbiota == "Exogenous", ])
print("Kruskal-Wallis test for Shannon diversity between species (Exogenous):")
print(shannon_test_exo)

# Observed ASVs
observed_test_endo <- kruskal.test(observed ~ InsectSpecies, 
                                   data = alpha_div_merged[alpha_div_merged$Microbiota == "Endogenous", ])
print("Kruskal-Wallis test for Observed ASVs between species (Endogenous):")
print(observed_test_endo)

observed_test_exo <- kruskal.test(observed ~ InsectSpecies, 
                                  data = alpha_div_merged[alpha_div_merged$Microbiota == "Exogenous", ])
print("Kruskal-Wallis test for Observed ASVs between species (Exogenous):")
print(observed_test_exo)

# Faith's PD
faith_pd_test_endo <- kruskal.test(PD ~ InsectSpecies, 
                                   data = alpha_div_merged[alpha_div_merged$Microbiota == "Endogenous", ])
print("Kruskal-Wallis test for Faith's PD between species (Endogenous):")
print(faith_pd_test_endo)

faith_pd_test_exo <- kruskal.test(PD ~ InsectSpecies, 
                                  data = alpha_div_merged[alpha_div_merged$Microbiota == "Exogenous", ])
print("Kruskal-Wallis test for Faith's PD between species (Exogenous):")
print(faith_pd_test_exo)

# Perform post-hoc pairwise tests (regardless of Kruskal-Wallis results)
# Shannon diversity - Endogenous
shannon_posthoc_endo <- pairwise.wilcox.test(
  alpha_div_merged$diversity_shannon[alpha_div_merged$Microbiota == "Endogenous"],
  alpha_div_merged$InsectSpecies[alpha_div_merged$Microbiota == "Endogenous"],
  p.adjust.method = "holm"
)
print("Post-hoc pairwise tests for Shannon diversity (Endogenous):")
print(shannon_posthoc_endo)

# Shannon diversity - Exogenous
shannon_posthoc_exo <- pairwise.wilcox.test(
  alpha_div_merged$diversity_shannon[alpha_div_merged$Microbiota == "Exogenous"],
  alpha_div_merged$InsectSpecies[alpha_div_merged$Microbiota == "Exogenous"],
  p.adjust.method = "holm"
)
print("Post-hoc pairwise tests for Shannon diversity (Exogenous):")
print(shannon_posthoc_exo)

# Observed ASVs - Endogenous
observed_posthoc_endo <- pairwise.wilcox.test(
  alpha_div_merged$observed[alpha_div_merged$Microbiota == "Endogenous"],
  alpha_div_merged$InsectSpecies[alpha_div_merged$Microbiota == "Endogenous"],
  p.adjust.method = "holm"
)
print("Post-hoc pairwise tests for Observed ASVs (Endogenous):")
print(observed_posthoc_endo)

# Observed ASVs - Exogenous
observed_posthoc_exo <- pairwise.wilcox.test(
  alpha_div_merged$observed[alpha_div_merged$Microbiota == "Exogenous"],
  alpha_div_merged$InsectSpecies[alpha_div_merged$Microbiota == "Exogenous"],
  p.adjust.method = "holm"
)
print("Post-hoc pairwise tests for Observed ASVs (Exogenous):")
print(observed_posthoc_exo)

# Faith's PD - Endogenous
faith_pd_posthoc_endo <- pairwise.wilcox.test(
  alpha_div_merged$PD[alpha_div_merged$Microbiota == "Endogenous"],
  alpha_div_merged$InsectSpecies[alpha_div_merged$Microbiota == "Endogenous"],
  p.adjust.method = "holm"
)
print("Post-hoc pairwise tests for Faith's PD (Endogenous):")
print(faith_pd_posthoc_endo)

# Faith's PD - Exogenous
faith_pd_posthoc_exo <- pairwise.wilcox.test(
  alpha_div_merged$PD[alpha_div_merged$Microbiota == "Exogenous"],
  alpha_div_merged$InsectSpecies[alpha_div_merged$Microbiota == "Exogenous"],
  p.adjust.method = "holm"
)
print("Post-hoc pairwise tests for Faith's PD (Exogenous):")
print(faith_pd_posthoc_exo)

# Function to extract significant pairs from pairwise test results
extract_sig_pairs <- function(pwtest, alpha = 0.05) {
  # Get the p-value matrix
  pvals <- pwtest$p.value
  
  # Initialize list to store significant pairs
  sig_pairs <- list()
  
  # Loop through the matrix to find significant pairs
  for (i in 1:nrow(pvals)) {
    for (j in 1:ncol(pvals)) {
      if (!is.na(pvals[i, j]) && pvals[i, j] < alpha) {
        # Add the pair to the list
        pair <- c(rownames(pvals)[i], colnames(pvals)[j])
        sig_pairs[[length(sig_pairs) + 1]] <- pair
      }
    }
  }
  
  return(sig_pairs)
}

# Extract significant pairs
shannon_sig_pairs_endo <- extract_sig_pairs(shannon_posthoc_endo)
shannon_sig_pairs_exo <- extract_sig_pairs(shannon_posthoc_exo)
observed_sig_pairs_endo <- extract_sig_pairs(observed_posthoc_endo)
observed_sig_pairs_exo <- extract_sig_pairs(observed_posthoc_exo)
faith_pd_sig_pairs_endo <- extract_sig_pairs(faith_pd_posthoc_endo)
faith_pd_sig_pairs_exo <- extract_sig_pairs(faith_pd_posthoc_exo)

# Print significant pairs
cat("Significant pairs for Shannon diversity (Endogenous):\n")
print(shannon_sig_pairs_endo)
cat("Significant pairs for Shannon diversity (Exogenous):\n")
print(shannon_sig_pairs_exo)
cat("Significant pairs for Observed ASVs (Endogenous):\n")
print(observed_sig_pairs_endo)
cat("Significant pairs for Observed ASVs (Exogenous):\n")
print(observed_sig_pairs_exo)
cat("Significant pairs for Faith's PD (Endogenous):\n")
print(faith_pd_sig_pairs_endo)
cat("Significant pairs for Faith's PD (Exogenous):\n")
print(faith_pd_sig_pairs_exo)

# Function to format significant pairs for ggsignif
format_sig_pairs <- function(sig_pairs) {
  if (length(sig_pairs) == 0) {
    return(NULL)
  }
  
  # Convert list of pairs to data frame
  result <- data.frame(
    group1 = sapply(sig_pairs, function(x) x[1]),
    group2 = sapply(sig_pairs, function(x) x[2]),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

# Format significant pairs for ggsignif
shannon_sig_pairs_endo_df <- format_sig_pairs(shannon_sig_pairs_endo)
shannon_sig_pairs_exo_df <- format_sig_pairs(shannon_sig_pairs_exo)
observed_sig_pairs_endo_df <- format_sig_pairs(observed_sig_pairs_endo)
observed_sig_pairs_exo_df <- format_sig_pairs(observed_sig_pairs_exo)
faith_pd_sig_pairs_endo_df <- format_sig_pairs(faith_pd_sig_pairs_endo)
faith_pd_sig_pairs_exo_df <- format_sig_pairs(faith_pd_sig_pairs_exo)

# Plot alpha diversity metrics faceted by Microbiota
# Faith's PD - Endogenous
p_faith_pd_endo <- ggplot(alpha_div_merged[alpha_div_merged$Microbiota == "Endogenous", ], 
                         aes(x = InsectSpecies, y = PD, fill = InsectSpecies)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Faith's Phylogenetic Diversity", title = "Faith's PD (Endogenous)") +
  scale_fill_viridis_d(option = "plasma") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Add significance bars if there are significant pairs
if (!is.null(faith_pd_sig_pairs_endo_df)) {
  p_faith_pd_endo <- p_faith_pd_endo +
    geom_signif(
      comparisons = lapply(1:nrow(faith_pd_sig_pairs_endo_df), function(i) {
        c(faith_pd_sig_pairs_endo_df$group1[i], faith_pd_sig_pairs_endo_df$group2[i])
      }),
      map_signif_level = TRUE,
      tip_length = 0.01,
      y_position = max(alpha_div_merged$PD[alpha_div_merged$Microbiota == "Endogenous"]) * 
        (1.05 + 0.1 * (1:nrow(faith_pd_sig_pairs_endo_df) - 1))
    )
}

# Faith's PD - Exogenous
p_faith_pd_exo <- ggplot(alpha_div_merged[alpha_div_merged$Microbiota == "Exogenous", ], 
                        aes(x = InsectSpecies, y = PD, fill = InsectSpecies)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Faith's Phylogenetic Diversity", title = "Faith's PD (Exogenous)") +
  scale_fill_viridis_d(option = "plasma") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Add significance bars if there are significant pairs
if (!is.null(faith_pd_sig_pairs_exo_df)) {
  p_faith_pd_exo <- p_faith_pd_exo +
    geom_signif(
      comparisons = lapply(1:nrow(faith_pd_sig_pairs_exo_df), function(i) {
        c(faith_pd_sig_pairs_exo_df$group1[i], faith_pd_sig_pairs_exo_df$group2[i])
      }),
      map_signif_level = TRUE,
      tip_length = 0.01,
      y_position = max(alpha_div_merged$PD[alpha_div_merged$Microbiota == "Exogenous"]) * 
        (1.05 + 0.1 * (1:nrow(faith_pd_sig_pairs_exo_df) - 1))
    )
}

# Combine Faith's PD plots
p_faith_pd <- p_faith_pd_endo + p_faith_pd_exo + plot_layout(ncol = 2)

# Shannon diversity - Endogenous
p_shannon_endo <- ggplot(alpha_div_merged[alpha_div_merged$Microbiota == "Endogenous", ], 
                         aes(x = InsectSpecies, y = diversity_shannon, fill = InsectSpecies)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Shannon Diversity Index", title = "Shannon Diversity (Endogenous)") +
  scale_fill_viridis_d(option = "plasma") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Add significance bars if there are significant pairs
if (!is.null(shannon_sig_pairs_endo_df)) {
  p_shannon_endo <- p_shannon_endo +
    geom_signif(
      comparisons = lapply(1:nrow(shannon_sig_pairs_endo_df), function(i) {
        c(shannon_sig_pairs_endo_df$group1[i], shannon_sig_pairs_endo_df$group2[i])
      }),
      map_signif_level = TRUE,
      tip_length = 0.01,
      y_position = max(alpha_div_merged$diversity_shannon[alpha_div_merged$Microbiota == "Endogenous"]) * 
        (1.05 + 0.1 * (1:nrow(shannon_sig_pairs_endo_df) - 1))
    )
}

# Shannon diversity - Exogenous
p_shannon_exo <- ggplot(alpha_div_merged[alpha_div_merged$Microbiota == "Exogenous", ], 
                        aes(x = InsectSpecies, y = diversity_shannon, fill = InsectSpecies)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Shannon Diversity Index", title = "Shannon Diversity (Exogenous)") +
  scale_fill_viridis_d(option = "plasma") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Add significance bars if there are significant pairs
if (!is.null(shannon_sig_pairs_exo_df)) {
  p_shannon_exo <- p_shannon_exo +
    geom_signif(
      comparisons = lapply(1:nrow(shannon_sig_pairs_exo_df), function(i) {
        c(shannon_sig_pairs_exo_df$group1[i], shannon_sig_pairs_exo_df$group2[i])
      }),
      map_signif_level = TRUE,
      tip_length = 0.01,
      y_position = max(alpha_div_merged$diversity_shannon[alpha_div_merged$Microbiota == "Exogenous"]) * 
        (1.05 + 0.1 * (1:nrow(shannon_sig_pairs_exo_df) - 1))
    )
}

# Combine Shannon plots
p_shannon <- p_shannon_endo + p_shannon_exo + plot_layout(ncol = 2)

# Observed ASVs - Endogenous
p_observed_endo <- ggplot(alpha_div_merged[alpha_div_merged$Microbiota == "Endogenous", ], 
                          aes(x = InsectSpecies, y = observed, fill = InsectSpecies)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Observed ASVs", title = "Observed ASVs (Endogenous)") +
  scale_fill_viridis_d(option = "plasma") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Add significance bars if there are significant pairs
if (!is.null(observed_sig_pairs_endo_df)) {
  p_observed_endo <- p_observed_endo +
    geom_signif(
      comparisons = lapply(1:nrow(observed_sig_pairs_endo_df), function(i) {
        c(observed_sig_pairs_endo_df$group1[i], observed_sig_pairs_endo_df$group2[i])
      }),
      map_signif_level = TRUE,
      tip_length = 0.01,
      y_position = max(alpha_div_merged$observed[alpha_div_merged$Microbiota == "Endogenous"]) * 
        (1.05 + 0.1 * (1:nrow(observed_sig_pairs_endo_df) - 1))
    )
}

# Observed ASVs - Exogenous
p_observed_exo <- ggplot(alpha_div_merged[alpha_div_merged$Microbiota == "Exogenous", ], 
                         aes(x = InsectSpecies, y = observed, fill = InsectSpecies)) +
  geom_boxplot(width = 0.5, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(x = "", y = "Observed ASVs", title = "Observed ASVs (Exogenous)") +
  scale_fill_viridis_d(option = "plasma") +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Add significance bars if there are significant pairs
if (!is.null(observed_sig_pairs_exo_df)) {
  p_observed_exo <- p_observed_exo +
    geom_signif(
      comparisons = lapply(1:nrow(observed_sig_pairs_exo_df), function(i) {
        c(observed_sig_pairs_exo_df$group1[i], observed_sig_pairs_exo_df$group2[i])
      }),
      map_signif_level = TRUE,
      tip_length = 0.01,
      y_position = max(alpha_div_merged$observed[alpha_div_merged$Microbiota == "Exogenous"]) * 
        (1.05 + 0.1 * (1:nrow(observed_sig_pairs_exo_df) - 1))
    )
}

# Combine Observed ASVs plots
p_observed <- p_observed_endo + p_observed_exo + plot_layout(ncol = 2)

# Display plots
print(p_shannon)
print(p_observed)
print(p_faith_pd)

# Save alpha diversity comparison plots
ggsave("fly_species_alpha_diversity_shannon.png", p_shannon, width = 12, height = 6)
ggsave("fly_species_alpha_diversity_observed.png", p_observed, width = 12, height = 6)
ggsave("fly_species_alpha_diversity_faith_pd.png", p_faith_pd, width = 12, height = 6)

# Combine plots
p_alpha_combined <- p_shannon / p_faith_pd
print(p_alpha_combined)
ggsave("fly_species_alpha_diversity_combined.png", p_alpha_combined, width = 15, height = 10)

#------------------------------------------------------------------------------
# 2. Beta Diversity Comparison Across Fly Species, Faceted by Microbiota
#------------------------------------------------------------------------------

# Calculate Bray-Curtis distance
bc_dist <- phyloseq::distance(ps, method = "bray")

# Perform PCoA
pcoa <- ordinate(ps, method = "PCoA", distance = bc_dist)

# Extract the variance explained by the first two axes
variance_explained <- pcoa$values$Relative_eig[1:2] * 100

# Plot PCoA faceted by Microbiota
p_pcoa <- plot_ordination(ps, pcoa, color = "InsectSpecies") +
  geom_point(size = 6, alpha = 0.8) +
  stat_ellipse(aes(group = InsectSpecies), type = "t", linetype = 2, alpha = 0.5) +
  facet_wrap(~ Microbiota) +
  labs(
    x = paste0("PCoA1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PCoA2 (", round(variance_explained[2], 1), "%)"),
    title = "Bray-Curtis PCoA by Fly Species"
  ) +
  scale_color_viridis_d(option = "plasma") +
  theme_pub()

print(p_pcoa)

# Save beta diversity plot
ggsave("fly_species_beta_diversity_pcoa.png", p_pcoa, width = 15, height = 9)

# Perform PERMANOVA to test for significant differences
# For Endogenous samples
metadata_df_endo <- data.frame(sample_data(ps))[sample_data(ps)$Microbiota == "Endogenous", ]
bc_dist_endo <- phyloseq::distance(subset_samples(ps, Microbiota == "Endogenous"), method = "bray")
permanova_result_endo <- adonis2(bc_dist_endo ~ InsectSpecies, data = metadata_df_endo)
print("PERMANOVA results for Bray-Curtis distances between species (Endogenous):")
print(permanova_result_endo)

# For Exogenous samples
metadata_df_exo <- data.frame(sample_data(ps))[sample_data(ps)$Microbiota == "Exogenous", ]
bc_dist_exo <- phyloseq::distance(subset_samples(ps, Microbiota == "Exogenous"), method = "bray")
permanova_result_exo <- adonis2(bc_dist_exo ~ InsectSpecies, data = metadata_df_exo)
print("PERMANOVA results for Bray-Curtis distances between species (Exogenous):")
print(permanova_result_exo)

# Test for homogeneity of dispersion (assumption of PERMANOVA)
# For Endogenous samples
beta_dispersion_endo <- betadisper(bc_dist_endo, metadata_df_endo$InsectSpecies)
beta_dispersion_test_endo <- permutest(beta_dispersion_endo)
print("Test for homogeneity of dispersion (Endogenous):")
print(beta_dispersion_test_endo)

# For Exogenous samples
beta_dispersion_exo <- betadisper(bc_dist_exo, metadata_df_exo$InsectSpecies)
beta_dispersion_test_exo <- permutest(beta_dispersion_exo)
print("Test for homogeneity of dispersion (Exogenous):")
print(beta_dispersion_test_exo)

library(pairwiseAdonis)

# If PERMANOVA is significant, perform pairwise comparisons
if (permanova_result_endo$`Pr(>F)`[1] < 0.05) {
  # Pairwise PERMANOVA for Endogenous samples
  pairwise_permanova_endo <- pairwise.adonis2(bc_dist_endo ~ InsectSpecies, data = metadata_df_endo)
  print("Pairwise PERMANOVA results for Endogenous samples:")
  print(pairwise_permanova_endo)
}

if (permanova_result_exo$`Pr(>F)`[1] < 0.05) {
  # Pairwise PERMANOVA for Exogenous samples
  pairwise_permanova_exo <- pairwise.adonis2(bc_dist_exo ~ InsectSpecies, data = metadata_df_exo)
  print("Pairwise PERMANOVA results for Exogenous samples:")
  print(pairwise_permanova_exo)
}

# Define the pairwise.adonis2 function if it's not available
pairwise.adonis2 <- function(dist, formula, data, p.adjust.m = "holm") {
  co <- combn(unique(as.character(data[,as.character(formula[[2]])])), 2)
  pairs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  for(i in 1:ncol(co)) {
    sub_data <- subset(data, get(as.character(formula[[2]])) %in% c(co[1,i], co[2,i]))
    sub_dist <- as.dist(as.matrix(dist)[rownames(sub_data), rownames(sub_data)])
    
    ad <- adonis2(sub_dist ~ get(as.character(formula[[2]])), data = sub_data)
    
    pairs <- c(pairs, paste(co[1,i], "vs", co[2,i]))
    F.Model <- c(F.Model, ad$F[1])
    R2 <- c(R2, ad$R2[1])
    p.value <- c(p.value, ad$`Pr(>F)`[1])
  }
  
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  
  results <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
  return(results)
}

#------------------------------------------------------------------------------
# 3. Heatmap of Top 20 Most Abundant Genera Across All Fly Species
#------------------------------------------------------------------------------

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Melt the phyloseq object for plotting
ps_melt <- psmelt(ps_rel)

# Get the top 20 most abundant genera across all samples
top_genera_all <- ps_melt %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  top_n(20) %>%
  pull(Genus)

# Filter for top genera
ps_melt_top <- ps_melt %>%
  filter(Genus %in% top_genera_all)

# Calculate mean abundance by genus, species, and microbiota type
genus_abundance_summary <- ps_melt_top %>%
  group_by(Microbiota, InsectSpecies, Genus) %>%
  summarise(
    mean_abundance = mean(Abundance),
    se_abundance = sd(Abundance) / sqrt(n())
  ) %>%
  ungroup()

# Create a heatmap of mean relative abundance by genus, species, and microbiota type
# Ensure genera are ordered by overall abundance
genus_order <- ps_melt_top %>%
  group_by(Genus) %>%
  summarise(overall_abundance = mean(Abundance)) %>%
  arrange(desc(overall_abundance)) %>%
  pull(Genus)

# Set the factor levels for genus based on the order
genus_abundance_summary$Genus <- factor(genus_abundance_summary$Genus, levels = genus_order)

# Create separate heatmaps for Endogenous and Exogenous
# For Endogenous - Log-transformed version
genus_heatmap_endo <- genus_abundance_summary %>%
  filter(Microbiota == "Endogenous") %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 20 Genera Abundance (Endogenous) - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# For Endogenous - Capped color scale version
# Calculate 95th percentile for capping
cap_value_endo <- quantile(genus_abundance_summary$mean_abundance[genus_abundance_summary$Microbiota == "Endogenous"], 0.95)

genus_heatmap_endo_capped <- genus_abundance_summary %>%
  filter(Microbiota == "Endogenous") %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value_endo)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 20 Genera Abundance (Endogenous) - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# For Exogenous - Log-transformed version
genus_heatmap_exo <- genus_abundance_summary %>%
  filter(Microbiota == "Exogenous") %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 20 Genera Abundance (Exogenous) - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# For Exogenous - Capped color scale version
# Calculate 95th percentile for capping
cap_value_exo <- quantile(genus_abundance_summary$mean_abundance[genus_abundance_summary$Microbiota == "Exogenous"], 0.95)

genus_heatmap_exo_capped <- genus_abundance_summary %>%
  filter(Microbiota == "Exogenous") %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value_exo)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 20 Genera Abundance (Exogenous) - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# Display heatmaps
print(genus_heatmap_endo)
print(genus_heatmap_exo)

# Save heatmaps
ggsave("fly_species_top20_genera_heatmap_endogenous.png", genus_heatmap_endo, width = 10, height = 8)
ggsave("fly_species_top20_genera_heatmap_exogenous.png", genus_heatmap_exo, width = 10, height = 8)

# Combined heatmap with faceting - Log-transformed version
genus_heatmap_combined <- genus_abundance_summary %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 20 Genera Abundance by Fly Species and Microbiota Type - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# Combined heatmap with faceting - Capped color scale version
# Calculate 95th percentile for capping across all data
cap_value_combined <- quantile(genus_abundance_summary$mean_abundance, 0.95)

genus_heatmap_combined_capped <- genus_abundance_summary %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value_combined)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 20 Genera Abundance by Fly Species and Microbiota Type - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

print(genus_heatmap_combined)
ggsave("fly_species_top20_genera_heatmap_combined.png", genus_heatmap_combined, width = 14, height = 8)

#------------------------------------------------------------------------------
# 4. Heatmap of Top 5 Most Abundant Genera from Each Fly Species
#------------------------------------------------------------------------------

# Get the top 5 genera for each fly species
top_genera_by_species <- ps_melt %>%
  group_by(InsectSpecies, Genus) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(InsectSpecies, desc(mean_abundance)) %>%
  group_by(InsectSpecies) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  pull(Genus) %>%
  unique()

# Filter for these genera
ps_melt_top_by_species <- ps_melt %>%
  filter(Genus %in% top_genera_by_species)

# Calculate mean abundance by genus, species, and microbiota type
genus_abundance_by_species <- ps_melt_top_by_species %>%
  group_by(Microbiota, InsectSpecies, Genus) %>%
  summarise(
    mean_abundance = mean(Abundance),
    se_abundance = sd(Abundance) / sqrt(n())
  ) %>%
  ungroup()

# Order genera by overall abundance
genus_order_by_species <- ps_melt_top_by_species %>%
  group_by(Genus) %>%
  summarise(overall_abundance = mean(Abundance)) %>%
  arrange(desc(overall_abundance)) %>%
  pull(Genus)

# Set the factor levels for genus based on the order
genus_abundance_by_species$Genus <- factor(genus_abundance_by_species$Genus, levels = genus_order_by_species)

# Create separate heatmaps for Endogenous and Exogenous
# For Endogenous - Log-transformed version
genus_heatmap_by_species_endo <- genus_abundance_by_species %>%
  filter(Microbiota == "Endogenous") %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 5 Genera from Each Species (Endogenous) - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# For Endogenous - Capped color scale version
# Calculate 95th percentile for capping
cap_value_by_species_endo <- quantile(genus_abundance_by_species$mean_abundance[genus_abundance_by_species$Microbiota == "Endogenous"], 0.95)

genus_heatmap_by_species_endo_capped <- genus_abundance_by_species %>%
  filter(Microbiota == "Endogenous") %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value_by_species_endo)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 5 Genera from Each Species (Endogenous) - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# For Exogenous - Log-transformed version
genus_heatmap_by_species_exo <- genus_abundance_by_species %>%
  filter(Microbiota == "Exogenous") %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 5 Genera from Each Species (Exogenous) - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# For Exogenous - Capped color scale version
# Calculate 95th percentile for capping
cap_value_by_species_exo <- quantile(genus_abundance_by_species$mean_abundance[genus_abundance_by_species$Microbiota == "Exogenous"], 0.95)

genus_heatmap_by_species_exo_capped <- genus_abundance_by_species %>%
  filter(Microbiota == "Exogenous") %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value_by_species_exo)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 5 Genera from Each Species (Exogenous) - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# Display heatmaps
print(genus_heatmap_by_species_endo)
print(genus_heatmap_by_species_exo)

# Save heatmaps
ggsave("fly_species_top5_by_species_heatmap_endogenous.png", genus_heatmap_by_species_endo, width = 10, height = 8)
ggsave("fly_species_top5_by_species_heatmap_exogenous.png", genus_heatmap_by_species_exo, width = 10, height = 8)

# Combined heatmap with faceting - Log-transformed version
genus_heatmap_by_species_combined <- genus_abundance_by_species %>%
  # Add a small pseudocount and log-transform
  mutate(log_abundance = log10(mean_abundance + 0.0001)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = log_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_fill_viridis_c(option = "plasma", name = "Log10 Mean\nRelative Abundance") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 5 Genera from Each Fly Species by Microbiota Type - Log-transformed"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

# Combined heatmap with faceting - Capped color scale version
# Calculate 95th percentile for capping across all data
cap_value_by_species_combined <- quantile(genus_abundance_by_species$mean_abundance, 0.95)

genus_heatmap_by_species_combined_capped <- genus_abundance_by_species %>%
  # Cap values at 95th percentile
  mutate(capped_abundance = pmin(mean_abundance, cap_value_by_species_combined)) %>%
  ggplot(aes(x = InsectSpecies, y = Genus, fill = capped_abundance)) +
  geom_tile(color = "white", size = 0.5) +
  facet_wrap(~ Microbiota, scales = "free_x") +
  scale_fill_viridis_c(option = "plasma", name = "Mean Relative\nAbundance\n(capped at 95th percentile)") +
  labs(
    x = "Fly Species",
    y = "Genus",
    title = "Top 5 Genera from Each Fly Species by Microbiota Type - Capped Scale"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )

print(genus_heatmap_by_species_combined)
ggsave("fly_species_top5_by_species_heatmap_combined.png", genus_heatmap_by_species_combined, width = 14, height = 8)

#------------------------------------------------------------------------------
# 5. ANCOM-BC Comparing Species, Split by Microbiota Type
#------------------------------------------------------------------------------

# Agglomerate taxa at genus level
ps_genus <- tax_glom(ps, taxrank = "Genus")

# Split the dataset by microbiota type
ps_endo <- subset_samples(ps_genus, Microbiota == "Endogenous")
ps_exo <- subset_samples(ps_genus, Microbiota == "Exogenous")

# Run ANCOM-BC2 analysis for Endogenous samples
cat("\nRunning ANCOM-BC2 for InsectSpecies (Endogenous) at Genus level...\n")

# Set seed for reproducibility
set.seed(123)

# Run ANCOM-BC2 for Endogenous samples
ancombc_result_endo <- ancombc2(
  data = ps_endo,
  tax_level = "Genus",
  fix_formula = "InsectSpecies",
  p_adj_method = "holm",
  pseudo_sens = TRUE,  # Enable sensitivity analysis for pseudo-count addition
  prv_cut = 0.10,      # Filter features with low prevalence
  lib_cut = 1000,      # Minimum library size
  s0_perc = 0.05,      # Percentile for variance regularization
  group = "InsectSpecies",
  struc_zero = TRUE,   # Detect and handle structural zeros
  neg_lb = TRUE,       # Use negative lower bound for log-fold change
  alpha = 0.05,        # Significance level
  global = TRUE,       # Global test for multi-group comparison
  pairwise = TRUE      # Pairwise test for multi-group comparison
)

# Extract results
res_prim_endo <- ancombc_result_endo$res_pair

# Create a data frame for visualization
df_species_endo <- res_prim_endo %>%
  dplyr::select(taxon, contains("lfc_"), contains("diff_"), contains("q_"))

# Filter for significant results
df_sig_endo <- df_species_endo %>%
  dplyr::filter(if_any(starts_with("diff_"), ~ . == 1)) %>%
  dplyr::arrange(taxon)

# Print summary
if (nrow(df_sig_endo) > 0) {
  cat(paste0("\nFound ", nrow(df_sig_endo), " differentially abundant genera between fly species (Endogenous)\n"))
  
  # Save significant results to CSV
  write.csv(df_sig_endo, "ancombc2_Genus_fly_species_Endogenous.csv", row.names = FALSE)
  
  # Get taxonomy information for plotting
  tax_info <- as.data.frame(tax_table(ps_genus))
  tax_info$taxon <- rownames(tax_info)  # Add row names as a column
  
  # Create a data frame for plotting
  plot_data_endo <- df_sig_endo %>%
    dplyr::left_join(tax_info, by = "taxon") %>%
    dplyr::mutate(
      # Use taxon ID as the display name
      DisplayName = taxon
    )
  
  # Create a heatmap of log fold changes for significant taxa
  # Reshape the data for heatmap
  lfc_cols <- grep("^lfc_", names(plot_data_endo), value = TRUE)
  diff_cols <- grep("^diff_", names(plot_data_endo), value = TRUE)
  
  # Create a long format data frame for log fold changes
  lfc_data <- plot_data_endo %>%
    dplyr::select(DisplayName, all_of(lfc_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(lfc_cols),
      names_to = "Comparison",
      values_to = "LogFoldChange"
    ) %>%
    dplyr::mutate(
      Comparison = gsub("lfc_", "", Comparison)
    )
  
  # Create a long format data frame for significance
  diff_data <- plot_data_endo %>%
    dplyr::select(DisplayName, all_of(diff_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(diff_cols),
      names_to = "Comparison",
      values_to = "Significant"
    ) %>%
    dplyr::mutate(
      Comparison = gsub("diff_", "", Comparison)
    )
  
  # Merge the two data frames
  heatmap_data_endo <- lfc_data %>%
    dplyr::left_join(diff_data, by = c("DisplayName", "Comparison"))
  
  # Create a heatmap of log fold changes
  p_lfc_heatmap_endo <- ggplot(heatmap_data_endo, aes(x = Comparison, y = DisplayName, fill = LogFoldChange)) +
    geom_tile() +
    # Add a border to significant taxa
    geom_tile(
      data = heatmap_data_endo %>% filter(Significant == 1),
      fill = NA,
      color = "black",
      size = 0.5
    ) +
    scale_fill_gradient2(
      low = "#D95F02",
      mid = "white",
      high = "#1B9E77",
      midpoint = 0,
      name = "Log Fold Change"
    ) +
    labs(
      x = "",
      y = "",
      title = "ANCOM-BC2: Differentially Abundant Genera (Endogenous)",
      subtitle = "Black borders indicate statistically significant differences"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # Display and save heatmap
  print(p_lfc_heatmap_endo)
  ggsave("fly_species_ancombc_heatmap_Genus_Endogenous.png", p_lfc_heatmap_endo, width = 14, height = 10)
  
  # Create a dot plot of log fold changes for significant taxa
  # This is an alternative visualization that may be easier to interpret
  p_lfc_dot_endo <- ggplot(heatmap_data_endo, aes(x = Comparison, y = DisplayName)) +
    # Add background tiles for better visibility
    geom_tile(aes(fill = LogFoldChange), alpha = 0.7) +
    # Add points to indicate significance
    geom_point(aes(size = as.factor(Significant)), color = "black", alpha = 0.8) +
    scale_fill_gradient2(
      low = "#D95F02",
      mid = "#FFFFFF",
      high = "#1B9E77",
      midpoint = 0,
      name = "Log Fold Change"
    ) +
    scale_size_manual(
      values = c("0" = 0, "1" = 3),
      name = "Significant",
      labels = c("0" = "No", "1" = "Yes")
    ) +
    labs(
      x = "",
      y = "",
      title = "ANCOM-BC2: Differentially Abundant Genera (Endogenous)",
      subtitle = "Larger points indicate statistically significant differences"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # Display and save dot plot
  print(p_lfc_dot_endo)
  ggsave("fly_species_ancombc_dotplot_Genus_Endogenous.png", p_lfc_dot_endo, width = 14, height = 10)
  
  # Store results for summary
  genus_results_endo <- list(
    sig_taxa = df_sig_endo,
    plot_data = plot_data_endo,
    p_lfc_heatmap = p_lfc_heatmap_endo,
    p_lfc_dot = p_lfc_dot_endo
  )
} else {
  cat("\nNo differentially abundant genera found between fly species (Endogenous) using ANCOM-BC2.\n")
  genus_results_endo <- NULL
}

# Run ANCOM-BC2 analysis for Exogenous samples
cat("\nRunning ANCOM-BC2 for InsectSpecies (Exogenous) at Genus level...\n")

# Set seed for reproducibility
set.seed(123)

# Run ANCOM-BC2 for Exogenous samples
ancombc_result_exo <- ancombc2(
  data = ps_exo,
  tax_level = "Genus",
  fix_formula = "InsectSpecies",
  p_adj_method = "holm",
  pseudo_sens = TRUE,  # Enable sensitivity analysis for pseudo-count addition
  prv_cut = 0.10,      # Filter features with low prevalence
  lib_cut = 1000,      # Minimum library size
  s0_perc = 0.05,      # Percentile for variance regularization
  group = "InsectSpecies",
  struc_zero = TRUE,   # Detect and handle structural zeros
  neg_lb = TRUE,       # Use negative lower bound for log-fold change
  alpha = 0.05,        # Significance level
  global = TRUE,       # Global test for multi-group comparison
  pairwise = TRUE      # Pairwise test for multi-group comparison
)

# Extract results
res_prim_exo <- ancombc_result_exo$res_pair

# Create a data frame for visualization
df_species_exo <- res_prim_exo %>%
  dplyr::select(taxon, contains("lfc_"), contains("diff_"), contains("q_"))

# Filter for significant results
df_sig_exo <- df_species_exo %>%
  dplyr::filter(if_any(starts_with("diff_"), ~ . == 1)) %>%
  dplyr::arrange(taxon)

# Print summary
if (nrow(df_sig_exo) > 0) {
  cat(paste0("\nFound ", nrow(df_sig_exo), " differentially abundant genera between fly species (Exogenous)\n"))
  
  # Save significant results to CSV
  write.csv(df_sig_exo, "ancombc2_Genus_fly_species_Exogenous.csv", row.names = FALSE)
  
  # Get taxonomy information for plotting
  tax_info <- as.data.frame(tax_table(ps_genus))
  tax_info$taxon <- rownames(tax_info)  # Add row names as a column
  
  # Create a data frame for plotting
  plot_data_exo <- df_sig_exo %>%
    dplyr::left_join(tax_info, by = "taxon") %>%
    dplyr::mutate(
      # Use taxon ID as the display name
      DisplayName = taxon
    )
  
  # Create a heatmap of log fold changes for significant taxa
  # Reshape the data for heatmap
  lfc_cols <- grep("^lfc_", names(plot_data_exo), value = TRUE)
  diff_cols <- grep("^diff_", names(plot_data_exo), value = TRUE)
  
  # Create a long format data frame for log fold changes
  lfc_data <- plot_data_exo %>%
    dplyr::select(DisplayName, all_of(lfc_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(lfc_cols),
      names_to = "Comparison",
      values_to = "LogFoldChange"
    ) %>%
    dplyr::mutate(
      Comparison = gsub("lfc_", "", Comparison)
    )
  
  # Create a long format data frame for significance
  diff_data <- plot_data_exo %>%
    dplyr::select(DisplayName, all_of(diff_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(diff_cols),
      names_to = "Comparison",
      values_to = "Significant"
    ) %>%
    dplyr::mutate(
      Comparison = gsub("diff_", "", Comparison)
    )
  
  # Merge the two data frames
  heatmap_data_exo <- lfc_data %>%
    dplyr::left_join(diff_data, by = c("DisplayName", "Comparison"))
  
  # Create a heatmap of log fold changes
  p_lfc_heatmap_exo <- ggplot(heatmap_data_exo, aes(x = Comparison, y = DisplayName, fill = LogFoldChange)) +
    geom_tile() +
    # Add a border to significant taxa
    geom_tile(
      data = heatmap_data_exo %>% filter(Significant == 1),
      fill = NA,
      color = "black",
      size = 0.5
    ) +
    scale_fill_gradient2(
      low = "#D95F02",
      mid = "white",
      high = "#1B9E77",
      midpoint = 0,
      name = "Log Fold Change"
    ) +
    labs(
      x = "",
      y = "",
      title = "ANCOM-BC2: Differentially Abundant Genera (Exogenous)",
      subtitle = "Black borders indicate statistically significant differences"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # Display and save heatmap
  print(p_lfc_heatmap_exo)
  ggsave("fly_species_ancombc_heatmap_Genus_Exogenous.png", p_lfc_heatmap_exo, width = 14, height = 10)
  
  # Create a dot plot of log fold changes for significant taxa
  # This is an alternative visualization that may be easier to interpret
  p_lfc_dot_exo <- ggplot(heatmap_data_exo, aes(x = Comparison, y = DisplayName)) +
    # Add background tiles for better visibility
    geom_tile(aes(fill = LogFoldChange), alpha = 0.7) +
    # Add points to indicate significance
    geom_point(aes(size = as.factor(Significant)), color = "black", alpha = 0.8) +
    scale_fill_gradient2(
      low = "#D95F02",
      mid = "#FFFFFF",
      high = "#1B9E77",
      midpoint = 0,
      name = "Log Fold Change"
    ) +
    scale_size_manual(
      values = c("0" = 0, "1" = 3),
      name = "Significant",
      labels = c("0" = "No", "1" = "Yes")
    ) +
    labs(
      x = "",
      y = "",
      title = "ANCOM-BC2: Differentially Abundant Genera (Exogenous)",
      subtitle = "Larger points indicate statistically significant differences"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  # Display and save dot plot
  print(p_lfc_dot_exo)
  ggsave("fly_species_ancombc_dotplot_Genus_Exogenous.png", p_lfc_dot_exo, width = 14, height = 10)
  
  # Store results for summary
  genus_results_exo <- list(
    sig_taxa = df_sig_exo,
    plot_data = plot_data_exo,
    p_lfc_heatmap = p_lfc_heatmap_exo,
    p_lfc_dot = p_lfc_dot_exo
  )
} else {
  cat("\nNo differentially abundant genera found between fly species (Exogenous) using ANCOM-BC2.\n")
  genus_results_exo <- NULL
}

# Print summary of the comparison
cat("\nComparison of fly species completed.\n")
cat("Output files:\n")
cat("- fly_species_alpha_diversity_shannon.png\n")
cat("- fly_species_alpha_diversity_observed.png\n")
cat("- fly_species_alpha_diversity_faith_pd.png\n")
cat("- fly_species_alpha_diversity_combined.png\n")
cat("- fly_species_beta_diversity_pcoa.png\n")
cat("- fly_species_top20_genera_heatmap_endogenous.png\n")
cat("- fly_species_top20_genera_heatmap_exogenous.png\n")
cat("- fly_species_top20_genera_heatmap_combined.png\n")
cat("- fly_species_top5_by_species_heatmap_endogenous.png\n")
cat("- fly_species_top5_by_species_heatmap_exogenous.png\n")
cat("- fly_species_top5_by_species_heatmap_combined.png\n")
if (!is.null(genus_results_endo)) {
  cat("- ancombc2_Genus_fly_species_Endogenous.csv\n")
  cat("- fly_species_ancombc_heatmap_Genus_Endogenous.png\n")
  cat("- fly_species_ancombc_dotplot_Genus_Endogenous.png\n")
}
if (!is.null(genus_results_exo)) {
  cat("- ancombc2_Genus_fly_species_Exogenous.csv\n")
  cat("- fly_species_ancombc_heatmap_Genus_Exogenous.png\n")
  cat("- fly_species_ancombc_dotplot_Genus_Exogenous.png\n")
}
