# Extraction of filtered OTU tables for Cochliomyia macellaria analysis
# 2025-09-18

# Load required packages
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(microbiome)
library(biomformat)  # For BIOM format export

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

# Load QIIME2 artifacts as phyloseq object
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

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

#------------------------------------------------------------------------------
# 1. Extract OTU table with JUST Cochliomyia macellaria samples
#------------------------------------------------------------------------------

# Filter to keep only Cochliomyia macellaria samples
ps_macellaria <- subset_samples(ps, InsectSpecies == "Cochliomyia macellaria")

# Remove ASVs with zero counts after filtering
ps_macellaria <- prune_taxa(taxa_sums(ps_macellaria) > 0, ps_macellaria)

# Print summary of the filtered dataset
print("Summary of Cochliomyia macellaria dataset:")
ps_macellaria
print(paste("Number of samples:", nsamples(ps_macellaria)))
print(paste("Number of ASVs:", ntaxa(ps_macellaria)))

# Print sample counts by microbiota type
macellaria_counts <- table(sample_data(ps_macellaria)$Microbiota)
print("Sample counts by microbiota type for Cochliomyia macellaria:")
print(macellaria_counts)

# Extract OTU table
otu_table_macellaria <- as.data.frame(otu_table(ps_macellaria))
print("Dimensions of Cochliomyia macellaria OTU table:")
print(dim(otu_table_macellaria))

# Save OTU table to file
write.csv(otu_table_macellaria, "otu_table_macellaria_only.csv")

# Save taxonomy table
tax_table_macellaria <- as.data.frame(tax_table(ps_macellaria))
write.csv(tax_table_macellaria, "tax_table_macellaria_only.csv")

# Save sample data
sample_data_macellaria <- as.data.frame(sample_data(ps_macellaria))
write.csv(sample_data_macellaria, "sample_data_macellaria_only.csv")

# Save phyloseq object
saveRDS(ps_macellaria, "ps_macellaria_only.rds")

# Convert to BIOM format and save
# A simpler approach to avoid metadata dimension issues
# First, ensure the OTU table is in the correct orientation (samples as columns, taxa as rows)
otu_mat <- t(as(otu_table(ps_macellaria), "matrix"))

# Create a simple BIOM object without metadata
biom_data <- make_biom(data = otu_mat)

# Write BIOM file
write_biom(biom_data, "otu_table_macellaria_only.biom")

# Print confirmation
cat("BIOM file for Cochliomyia macellaria created successfully.\n")

#------------------------------------------------------------------------------
# 2. Extract OTU table with selected fly species
#------------------------------------------------------------------------------

# Define the species to keep
selected_species <- c("Cochliomyia macellaria", "Chrysomya rufifacies", 
                      "Lucilia coeruleiviridis", "Chrysomya megacephala")

# Filter to keep only the selected species
ps_selected <- subset_samples(ps, InsectSpecies %in% selected_species)

# Remove ASVs with zero counts after filtering
ps_selected <- prune_taxa(taxa_sums(ps_selected) > 0, ps_selected)

# Print summary of the filtered dataset
print("Summary of selected species dataset:")
ps_selected
print(paste("Number of samples:", nsamples(ps_selected)))
print(paste("Number of ASVs:", ntaxa(ps_selected)))

# Print sample counts by species and microbiota type
selected_counts <- table(sample_data(ps_selected)$InsectSpecies, sample_data(ps_selected)$Microbiota)
print("Sample counts by species and microbiota type for selected species:")
print(selected_counts)

# Extract OTU table
otu_table_selected <- as.data.frame(otu_table(ps_selected))
print("Dimensions of selected species OTU table:")
print(dim(otu_table_selected))

# Save OTU table to file
write.csv(otu_table_selected, "otu_table_selected_species.csv")

# Save taxonomy table
tax_table_selected <- as.data.frame(tax_table(ps_selected))
write.csv(tax_table_selected, "tax_table_selected_species.csv")

# Save sample data
sample_data_selected <- as.data.frame(sample_data(ps_selected))
write.csv(sample_data_selected, "sample_data_selected_species.csv")

# Save phyloseq object
saveRDS(ps_selected, "ps_selected_species.rds")

# Convert to BIOM format and save
# A simpler approach to avoid metadata dimension issues
# First, ensure the OTU table is in the correct orientation (samples as columns, taxa as rows)
otu_mat <- t(as(otu_table(ps_selected), "matrix"))

# Create a simple BIOM object without metadata
biom_data <- make_biom(data = otu_mat)

# Write BIOM file
write_biom(biom_data, "otu_table_selected_species.biom")

# Print confirmation
cat("BIOM file for selected species created successfully.\n")

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

cat("\nExtraction of filtered OTU tables completed.\n")
cat("Output files:\n")
cat("- otu_table_macellaria_only.csv\n")
cat("- otu_table_macellaria_only.biom\n")
cat("- tax_table_macellaria_only.csv\n")
cat("- sample_data_macellaria_only.csv\n")
cat("- ps_macellaria_only.rds\n")
cat("- otu_table_selected_species.csv\n")
cat("- otu_table_selected_species.biom\n")
cat("- tax_table_selected_species.csv\n")
cat("- sample_data_selected_species.csv\n")
cat("- ps_selected_species.rds\n")
