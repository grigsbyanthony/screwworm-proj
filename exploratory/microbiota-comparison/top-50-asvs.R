# ==================================================================================
# Load required libraries
# ==================================================================================
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(ggsignif)
library(rstatix)

# ==================================================================================
# Load custom ggplot2 theme
# ==================================================================================
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
      legend.key.size = grid::unit(0.2, "cm"),
      legend.title = element_text(face = "bold"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
}

# ==================================================================================
# Import data from QIIME2 artifacts
# ==================================================================================
ps <- qiime2R::qza_to_phyloseq(
  features = "merged-table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Print original phyloseq object summary
cat("Original phyloseq object:\n")
print(ps)

# ==================================================================================
# Create species-specific phyloseq objects and extract top 50 ASVs
# ==================================================================================

# Function to categorize microbiota as Endogenous or Exogenous
categorize_microbiota <- function(genus, family) {
  # Define endogenous (obligate symbionts) genera
  endogenous_genera <- c("Wolbachia", "Spiroplasma", "Rickettsia", "Buchnera", 
                        "Blochmannia", "Carsonella", "Portiera", "Hamiltonella",
                        "Regiella", "Serratia", "Arsenophonus")
  
  # Define endogenous families (for cases where genus is NA)
  endogenous_families <- c("Anaplasmataceae", "Rickettsiaceae", "Mycoplasmataceae")
  
  # Check if genus or family indicates endogenous microbiota
  if (!is.na(genus) && genus %in% endogenous_genera) {
    return("Endogenous")
  } else if (!is.na(family) && family %in% endogenous_families) {
    return("Endogenous") 
  } else {
    return("Exogenous")
  }
}

# Get unique species from metadata
metadata <- as.data.frame(sample_data(ps))
species_list <- unique(metadata$InsectSpecies)
cat("Processing", length(species_list), "insect species:\n")
print(species_list)

# Create directory for output files
if (!dir.exists("top_asvs_by_species")) {
  dir.create("top_asvs_by_species")
}

# Initialize list to store all results for combined table
all_species_results <- list()

# Process each species individually
for (i in seq_along(species_list)) {
  species <- species_list[i]
  cat("\n", paste0(strrep("=", 80)), "\n")
  cat("PROCESSING SPECIES", i, "of", length(species_list), ":", species, "\n")
  cat(paste0(strrep("=", 80)), "\n")
  
  # Create species-specific phyloseq object
  species_samples <- metadata[metadata$InsectSpecies == species, ]
  species_sample_names <- rownames(species_samples)
  
  # Subset phyloseq object to species-specific samples
  ps_species <- prune_samples(species_sample_names, ps)
  
  # Remove ASVs with zero abundance in this species
  ps_species <- prune_taxa(taxa_sums(ps_species) > 0, ps_species)
  
  cat("Species-specific phyloseq object:\n")
  print(ps_species)
  
  # Convert to relative abundances
  ps_species_rel <- transform_sample_counts(ps_species, function(x) x / sum(x))
  
  # Calculate mean relative abundance for each ASV
  otu_table_species <- as.data.frame(otu_table(ps_species_rel))
  tax_table_species <- as.data.frame(tax_table(ps_species_rel))
  
  # Calculate mean abundance across samples for this species
  mean_abundances <- rowMeans(otu_table_species)
  
  # Get top 50 ASVs
  top_50_indices <- order(mean_abundances, decreasing = TRUE)[1:min(50, length(mean_abundances))]
  top_50_asvs <- names(mean_abundances)[top_50_indices]
  top_50_abundances <- mean_abundances[top_50_indices]
  
  # Get taxonomy for top 50 ASVs
  top_50_taxonomy <- tax_table_species[top_50_asvs, , drop = FALSE]
  
  # Create species table with microbiota categorization
  species_data <- data.frame(
    Rank = 1:length(top_50_asvs),
    ASV_ID = top_50_asvs,
    Kingdom = top_50_taxonomy$Kingdom,
    Phylum = top_50_taxonomy$Phylum,
    Class = top_50_taxonomy$Class,
    Order = top_50_taxonomy$Order,
    Family = top_50_taxonomy$Family,
    Genus = top_50_taxonomy$Genus,
    Species = top_50_taxonomy$Species,
    Mean_Abundance_Percent = round(top_50_abundances * 100, 4),
    stringsAsFactors = FALSE
  )
  
  # Add microbiota category
  species_data$Microbiota_Type <- mapply(categorize_microbiota, 
                                       species_data$Genus, 
                                       species_data$Family)
  
  # Separate into Endogenous and Exogenous
  endogenous_data <- species_data[species_data$Microbiota_Type == "Endogenous", ]
  exogenous_data <- species_data[species_data$Microbiota_Type == "Exogenous", ]
  
  # Print summary
  cat("\nSUMMARY FOR", species, ":\n")
  cat("Total ASVs found:", nrow(species_data), "\n")
  cat("Endogenous microbiota ASVs:", nrow(endogenous_data), "\n")
  cat("Exogenous microbiota ASVs:", nrow(exogenous_data), "\n")
  
  # Print top endogenous ASVs
  if (nrow(endogenous_data) > 0) {
    cat("\nTOP ENDOGENOUS MICROBIOTA:\n")
    print(endogenous_data[1:min(10, nrow(endogenous_data)), c("Rank", "Genus", "Family", "Mean_Abundance_Percent")])
  }
  
  # Print top exogenous ASVs  
  if (nrow(exogenous_data) > 0) {
    cat("\nTOP EXOGENOUS MICROBIOTA:\n")
    print(exogenous_data[1:min(10, nrow(exogenous_data)), c("Rank", "Genus", "Family", "Mean_Abundance_Percent")])
  }
  
  # Save individual CSV files
  filename_all <- paste0("top_asvs_by_species/", gsub("[^A-Za-z0-9]", "_", species), "_top50_all_asvs.csv")
  filename_endo <- paste0("top_asvs_by_species/", gsub("[^A-Za-z0-9]", "_", species), "_endogenous_asvs.csv")
  filename_exo <- paste0("top_asvs_by_species/", gsub("[^A-Za-z0-9]", "_", species), "_exogenous_asvs.csv")
  
  write.csv(species_data, filename_all, row.names = FALSE)
  write.csv(endogenous_data, filename_endo, row.names = FALSE)
  write.csv(exogenous_data, filename_exo, row.names = FALSE)
  
  cat("Files saved:\n")
  cat("- All ASVs:", filename_all, "\n")
  cat("- Endogenous:", filename_endo, "\n")
  cat("- Exogenous:", filename_exo, "\n")
  
  # Add species name and store for combined table
  species_data$Insect_Species <- species
  all_species_results[[i]] <- species_data
}

# Create combined table with all species
combined_table <- do.call(rbind, all_species_results)
combined_table <- combined_table[, c("Insect_Species", "Rank", "ASV_ID", "Kingdom", "Phylum", 
                                   "Class", "Order", "Family", "Genus", "Species", 
                                   "Mean_Abundance_Percent", "Microbiota_Type")]

# Save combined table
write.csv(combined_table, "top_asvs_by_species/combined_all_species_top50_asvs.csv", row.names = FALSE)

# Create summary statistics
summary_stats <- combined_table %>%
  group_by(Insect_Species, Microbiota_Type) %>%
  summarise(
    n_ASVs = n(),
    total_abundance = sum(Mean_Abundance_Percent),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Microbiota_Type, 
              values_from = c(n_ASVs, total_abundance),
              values_fill = 0)

write.csv(summary_stats, "top_asvs_by_species/microbiota_summary_by_species.csv", row.names = FALSE)

cat("\n", paste0(strrep("=", 80)), "\n")
cat("FINAL SUMMARY:\n")
cat("Processed", length(species_list), "insect species\n")
cat("Files created in 'top_asvs_by_species/' directory:\n")
cat("- Individual species tables (all, endogenous, exogenous)\n")
cat("- Combined table: combined_all_species_top50_asvs.csv\n")
cat("- Summary statistics: microbiota_summary_by_species.csv\n")
cat(paste0(strrep("=", 80)), "\n")