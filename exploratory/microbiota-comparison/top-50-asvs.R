# ==================================================================================
# Clean approach: Subset phyloseq by species, then extract top taxa
# ==================================================================================

library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)

# ==================================================================================
# Load custom ggplot2 theme (keeping your original theme)
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
# Import data from QIIME2 artifacts (uncomment if needed)
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
# Configuration - Choose what to extract
# ==================================================================================

# Choose what taxonomic level to analyze (you can change this)
ANALYSIS_LEVEL <- "ASV"  # Options: "ASV", "Genus", "Family", "Order", "Class", "Phylum"
N_TOP_TAXA <- 50         # Number of top taxa to extract
USE_RELATIVE <- TRUE     # TRUE for relative abundance, FALSE for raw counts

cat("\nAnalysis Configuration:\n")
cat("- Taxonomic level:", ANALYSIS_LEVEL, "\n")
cat("- Number of top taxa:", N_TOP_TAXA, "\n")
cat("- Use relative abundance:", USE_RELATIVE, "\n")

# ==================================================================================
# Get species list and subset phyloseq objects
# ==================================================================================

# Get metadata and clean species names
metadata <- as.data.frame(sample_data(ps))
metadata$InsectSpecies <- trimws(as.character(metadata$InsectSpecies))

# Get unique species
species_list <- unique(metadata$InsectSpecies)
species_list <- species_list[!is.na(species_list) & species_list != ""]

cat("\nFound species:\n")
for(i in 1:length(species_list)) {
  cat(i, ".", species_list[i], "\n")
}
cat("Total species:", length(species_list), "\n")

# ==================================================================================
# Function to extract top taxa from a phyloseq object
# ==================================================================================

extract_top_taxa <- function(ps_subset, n_taxa = 50, level = "ASV", use_relative = TRUE) {
  
  # Transform to relative abundance if requested
  if(use_relative) {
    ps_subset <- transform_sample_counts(ps_subset, function(x) x / sum(x))
  }
  
  # Get appropriate data based on analysis level
  if(level == "ASV") {
    # Use ASV-level data
    otu_table <- as.data.frame(otu_table(ps_subset))
    tax_table <- as.data.frame(tax_table(ps_subset))
    
    # Calculate mean abundance across samples for each ASV
    mean_abundances <- rowMeans(otu_table)
    
    # Get top N ASVs
    top_indices <- order(mean_abundances, decreasing = TRUE)[1:min(n_taxa, length(mean_abundances))]
    top_asvs <- names(mean_abundances)[top_indices]
    
    # Create results table
    results <- data.frame(
      Rank = 1:length(top_asvs),
      ASV_ID = top_asvs,
      Mean_Abundance = mean_abundances[top_asvs],
      stringsAsFactors = FALSE
    )
    
    # Add taxonomy information
    if(nrow(tax_table) > 0) {
      tax_subset <- tax_table[top_asvs, , drop = FALSE]
      results <- cbind(results, tax_subset)
    }
    
  } else {
    # Aggregate to specified taxonomic level
    ps_glom <- tax_glom(ps_subset, taxrank = level, NArm = TRUE)
    
    # Get aggregated abundance table
    otu_table_glom <- as.data.frame(otu_table(ps_glom))
    tax_table_glom <- as.data.frame(tax_table(ps_glom))
    
    # Calculate mean abundance
    mean_abundances <- rowMeans(otu_table_glom)
    
    # Get top N taxa
    top_indices <- order(mean_abundances, decreasing = TRUE)[1:min(n_taxa, length(mean_abundances))]
    top_taxa <- rownames(otu_table_glom)[top_indices]
    
    # Create results table
    results <- data.frame(
      Rank = 1:length(top_taxa),
      Taxa_ID = top_taxa,
      Mean_Abundance = mean_abundances[top_taxa],
      stringsAsFactors = FALSE
    )
    
    # Add taxonomy information
    if(nrow(tax_table_glom) > 0) {
      tax_subset <- tax_table_glom[top_taxa, , drop = FALSE]
      results <- cbind(results, tax_subset)
    }
  }
  
  # Convert abundance to percentage if using relative abundance
  if(use_relative) {
    results$Mean_Abundance_Percent <- round(results$Mean_Abundance * 100, 4)
    results$Mean_Abundance <- NULL  # Remove the decimal version
  } else {
    results$Mean_Abundance <- round(results$Mean_Abundance, 2)
  }
  
  return(results)
}

# ==================================================================================
# Process each species individually
# ==================================================================================

# Create output directory
if (!dir.exists("top_taxa_by_species")) {
  dir.create("top_taxa_by_species")
}

cat("\n", paste0(strrep("=", 80)), "\n")
cat("PROCESSING SPECIES INDIVIDUALLY\n")
cat(paste0(strrep("=", 80)), "\n")

all_results <- list()
successful_species <- 0

for(species in species_list) {
  
  cat("\n--- Processing:", species, "---\n")
  
  # Subset phyloseq object to current species
  species_samples <- rownames(metadata)[metadata$InsectSpecies == species]
  
  if(length(species_samples) == 0) {
    cat("WARNING: No samples found for species:", species, "\n")
    next
  }
  
  cat("Found", length(species_samples), "samples for", species, "\n")
  
  # Create subset phyloseq object
  ps_species <- prune_samples(species_samples, ps)
  
  # Remove ASVs that are not present in any sample of this species
  ps_species <- prune_taxa(taxa_sums(ps_species) > 0, ps_species)
  
  cat("Phyloseq subset: ", ntaxa(ps_species), "ASVs,", nsamples(ps_species), "samples\n")
  
  # Skip if no taxa remain
  if(ntaxa(ps_species) == 0) {
    cat("WARNING: No taxa found for species:", species, "\n")
    next
  }
  
  # Extract top taxa
  top_taxa_table <- extract_top_taxa(
    ps_species, 
    n_taxa = N_TOP_TAXA, 
    level = ANALYSIS_LEVEL,
    use_relative = USE_RELATIVE
  )
  
  cat("Extracted", nrow(top_taxa_table), "top", ANALYSIS_LEVEL, "for", species, "\n")
  
  # Save individual CSV file
  clean_species_name <- gsub("[^A-Za-z0-9_]", "_", species)
  clean_species_name <- gsub("_{2,}", "_", clean_species_name)
  
  abundance_type <- if(USE_RELATIVE) "relative" else "raw"
  filename <- paste0("top_taxa_by_species/", clean_species_name, "_top_", N_TOP_TAXA, "_", ANALYSIS_LEVEL, "_", abundance_type, ".csv")
  
  write.csv(top_taxa_table, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
  
  # Store results for combined table
  top_taxa_table$InsectSpecies <- species
  all_results[[species]] <- top_taxa_table
  
  # Print preview
  cat("Top 3 entries:\n")
  print(head(top_taxa_table[, 1:min(5, ncol(top_taxa_table))], 3))
  
  successful_species <- successful_species + 1
}

# ==================================================================================
# Create combined table
# ==================================================================================

if(length(all_results) > 0) {
  
  cat("\n--- Creating Combined Table ---\n")
  
  # Combine all results
  combined_results <- do.call(rbind, all_results)
  
  # Reorder columns to put species first
  col_order <- c("InsectSpecies", setdiff(names(combined_results), "InsectSpecies"))
  combined_results <- combined_results[, col_order]
  
  # Rename for clarity
  names(combined_results)[1] <- "Insect_Species"
  
  # Save combined table
  abundance_type <- if(USE_RELATIVE) "relative" else "raw"
  combined_filename <- paste0("top_taxa_by_species/COMBINED_all_species_top_", N_TOP_TAXA, "_", ANALYSIS_LEVEL, "_", abundance_type, ".csv")
  
  write.csv(combined_results, combined_filename, row.names = FALSE)
  cat("Combined table saved:", combined_filename, "\n")
  
  # Summary by species
  species_summary <- combined_results %>%
    group_by(Insect_Species) %>%
    summarise(
      n_taxa = n(),
      .groups = "drop"
    )
  
  cat("\nSummary by species:\n")
  print(species_summary)
}

# ==================================================================================
# Final Summary
# ==================================================================================

cat("\n", paste0(strrep("=", 80)), "\n")
cat("FINAL SUMMARY\n")
cat(paste0(strrep("=", 80)), "\n")
cat("Analysis level:", ANALYSIS_LEVEL, "\n")
cat("Abundance type:", if(USE_RELATIVE) "Relative abundance (%)" else "Raw counts", "\n")
cat("Number of top taxa extracted:", N_TOP_TAXA, "\n")
cat("Total species in dataset:", length(species_list), "\n")
cat("Successfully processed species:", successful_species, "\n")
cat("Output directory: top_taxa_by_species/\n")

if(successful_species < length(species_list)) {
  failed_species <- setdiff(species_list, names(all_results))
  cat("\nSpecies that couldn't be processed:\n")
  for(sp in failed_species) {
    cat("  -", sp, "\n")
  }
}

cat("\nFiles created:\n")
cat("- Individual CSV files for each species\n")
cat("- Combined CSV file with all species\n")
cat(paste0(strrep("=", 80)), "\n")

# ==================================================================================
# Quick verification
# ==================================================================================

cat("\nQuick verification - files created:\n")
files_created <- list.files("top_taxa_by_species", pattern = "\\.csv$")
cat("Number of files:", length(files_created), "\n")
for(file in files_created) {
  cat(" -", file, "\n")
}