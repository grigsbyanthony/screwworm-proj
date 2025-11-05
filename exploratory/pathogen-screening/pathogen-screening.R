# ==================================================================================
# Load required libraries
# ==================================================================================
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(vegan)
library(ape)

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
      legend.key.size = grid::unit(0.4, "cm"),
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
# Filter phyloseq object to Cochliomyia macellaria samples only
# ==================================================================================
ps_cmac <- subset_samples(ps, InsectSpecies == "Cochliomyia macellaria")

# Print filtered phyloseq object summary
cat("\nFiltered phyloseq object (Cochliomyia macellaria only):\n")
print(ps_cmac)

# ==================================================================================
# Split phyloseq object by microbiota type
# ==================================================================================
ps_endogenous <- subset_samples(ps_cmac, Microbiota == "Endogenous")
ps_exogenous <- subset_samples(ps_cmac, Microbiota == "Exogenous")

# Print summary of split phyloseq objects
cat("\nEndogenous microbiota phyloseq object:\n")
print(ps_endogenous)

cat("\nExogenous microbiota phyloseq object:\n")
print(ps_exogenous)

# ==================================================================================
# Screen for potentially pathogenic bacteria
# ==================================================================================

# Define list of known pathogenic genera commonly found in flies
pathogenic_genera <- c(
  "Salmonella",         # Foodborne pathogen
  "Escherichia",        # E. coli and relatives
  "Shigella",           # Dysentery
  "Campylobacter",      # Gastroenteritis
  "Vibrio",             # Cholera and related diseases
  "Yersinia",           # Plague and gastroenteritis
  "Clostridium",        # Various toxin-mediated diseases
  "Staphylococcus",     # Skin infections, food poisoning
  "Streptococcus",      # Various infections
  "Bacillus",           # B. anthracis and others
  "Listeria",           # Listeriosis
  "Mycobacterium",      # Tuberculosis and related
  "Corynebacterium",    # Diphtheria and opportunistic infections
  "Enterococcus",       # Healthcare-associated infections
  "Pseudomonas",        # Opportunistic pathogen
  "Acinetobacter",      # Hospital-acquired infections
  "Klebsiella",         # Pneumonia and UTIs
  "Proteus",            # UTIs and wound infections
  "Enterobacter",       # Healthcare-associated infections
  "Citrobacter",        # Opportunistic pathogen
  "Serratia",           # Opportunistic pathogen
  "Morganella",         # UTIs and opportunistic infections
  "Providencia"         # Opportunistic pathogen
)

# Function to screen for pathogenic genera
screen_pathogens <- function(ps_obj, object_name) {
  cat("\n========================================\n")
  cat("Pathogen screening for:", object_name, "\n")
  cat("========================================\n")
  
  # Get taxonomy table
  tax_table_df <- as.data.frame(tax_table(ps_obj))
  
  # Check if Genus column exists
  if (!"Genus" %in% colnames(tax_table_df)) {
    cat("Warning: No 'Genus' column found in taxonomy table\n")
    return(NULL)
  }
  
  # Find pathogenic genera present
  pathogenic_otus <- tax_table_df[tax_table_df$Genus %in% pathogenic_genera & 
                                  !is.na(tax_table_df$Genus), ]
  
  if (nrow(pathogenic_otus) > 0) {
    cat("\n*** POTENTIAL PATHOGENS DETECTED ***\n")
    cat("Number of pathogenic OTUs:", nrow(pathogenic_otus), "\n\n")
    
    # Display pathogenic genera found
    pathogenic_genera_found <- unique(pathogenic_otus$Genus[!is.na(pathogenic_otus$Genus)])
    cat("Pathogenic genera detected:\n")
    for (genus in pathogenic_genera_found) {
      genus_count <- sum(pathogenic_otus$Genus == genus, na.rm = TRUE)
      cat("-", genus, "(", genus_count, "OTUs )\n")
    }
    
    # Get abundance data for pathogenic OTUs
    pathogenic_otu_names <- rownames(pathogenic_otus)
    pathogenic_subset <- prune_taxa(pathogenic_otu_names, ps_obj)
    
    # Calculate relative abundances
    pathogenic_rel <- transform_sample_counts(pathogenic_subset, function(x) x/sum(x))
    
    # Summary statistics
    cat("\nAbundance summary for pathogenic taxa:\n")
    otu_sums <- taxa_sums(pathogenic_subset)
    cat("Total reads:", sum(otu_sums), "\n")
    cat("Mean reads per pathogenic OTU:", round(mean(otu_sums), 2), "\n")
    cat("Max reads for single pathogenic OTU:", max(otu_sums), "\n")
    
    # Return the pathogenic subset for further analysis
    return(list(
      pathogenic_otus = pathogenic_otus,
      pathogenic_subset = pathogenic_subset,
      pathogenic_rel = pathogenic_rel,
      genera_found = pathogenic_genera_found
    ))
    
  } else {
    cat("\n*** NO KNOWN PATHOGENIC GENERA DETECTED ***\n")
    cat("This is expected for healthy fly microbiomes\n")
    return(NULL)
  }
}

# Screen endogenous microbiota
pathogen_results_endo <- screen_pathogens(ps_endogenous, "Endogenous microbiota")

# Screen exogenous microbiota  
pathogen_results_exo <- screen_pathogens(ps_exogenous, "Exogenous microbiota")

# ==================================================================================
# Visualizations for pathogen screening results
# ==================================================================================

# Function to create pathogen abundance plots
plot_pathogen_abundance <- function(pathogen_results, microbiota_type) {
  if (is.null(pathogen_results)) {
    cat("\nNo pathogens detected in", microbiota_type, "- skipping visualization\n")
    return(NULL)
  }
  
  # Get relative abundance data
  ps_rel <- pathogen_results$pathogenic_rel
  
  # Convert to data frame for plotting
  rel_df <- psmelt(ps_rel)
  
  # Create bar plot of pathogenic genera by sample
  p1 <- ggplot(rel_df, aes(x = Sample, y = Abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Pathogenic Genera Abundance -", microbiota_type),
         x = "Sample", y = "Relative Abundance",
         fill = "Pathogenic Genus") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_viridis_d()
  
  # Create box plot of pathogenic genera abundance
  p2 <- ggplot(rel_df, aes(x = Genus, y = Abundance, fill = Genus)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(title = paste("Pathogenic Genera Distribution -", microbiota_type),
         x = "Pathogenic Genus", y = "Relative Abundance") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_fill_viridis_d()
  
  return(list(abundance_plot = p1, distribution_plot = p2))
}

# Function to create summary comparison plot
plot_pathogen_summary <- function(pathogen_results_endo, pathogen_results_exo) {
  
  # Prepare summary data
  summary_data <- data.frame(
    Microbiota = character(),
    Genus = character(),
    OTU_Count = numeric(),
    Total_Abundance = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(pathogen_results_endo)) {
    endo_summary <- pathogen_results_endo$pathogenic_otus %>%
      group_by(Genus) %>%
      summarise(OTU_Count = n(), .groups = 'drop') %>%
      mutate(Microbiota = "Endogenous")
    
    # Get abundance data
    endo_abundance <- taxa_sums(pathogen_results_endo$pathogenic_subset)
    endo_genus_abundance <- aggregate(endo_abundance, 
                                     by = list(pathogen_results_endo$pathogenic_otus$Genus), 
                                     FUN = sum)
    names(endo_genus_abundance) <- c("Genus", "Total_Abundance")
    
    endo_summary <- merge(endo_summary, endo_genus_abundance, by = "Genus")
    summary_data <- rbind(summary_data, endo_summary)
  }
  
  if (!is.null(pathogen_results_exo)) {
    exo_summary <- pathogen_results_exo$pathogenic_otus %>%
      group_by(Genus) %>%
      summarise(OTU_Count = n(), .groups = 'drop') %>%
      mutate(Microbiota = "Exogenous")
    
    # Get abundance data
    exo_abundance <- taxa_sums(pathogen_results_exo$pathogenic_subset)
    exo_genus_abundance <- aggregate(exo_abundance, 
                                    by = list(pathogen_results_exo$pathogenic_otus$Genus), 
                                    FUN = sum)
    names(exo_genus_abundance) <- c("Genus", "Total_Abundance")
    
    exo_summary <- merge(exo_summary, exo_genus_abundance, by = "Genus")
    summary_data <- rbind(summary_data, exo_summary)
  }
  
  if (nrow(summary_data) == 0) {
    cat("\nNo pathogens detected in either microbiota type - creating empty summary plot\n")
    
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No pathogenic genera detected\nin either endogenous or exogenous microbiota", 
               size = 16, hjust = 0.5, vjust = 0.5) +
      labs(title = "Pathogen Screening Summary - Cochliomyia macellaria") +
      theme_pub() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
    
    return(empty_plot)
  }
  
  # Create comparison plots
  p1 <- ggplot(summary_data, aes(x = Genus, y = OTU_Count, fill = Microbiota)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Pathogenic OTU Counts by Microbiota Type",
         x = "Pathogenic Genus", y = "Number of OTUs",
         fill = "Microbiota Type") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  p2 <- ggplot(summary_data, aes(x = Genus, y = log10(Total_Abundance + 1), fill = Microbiota)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Pathogenic Taxa Abundance by Microbiota Type",
         x = "Pathogenic Genus", y = "Log10(Total Abundance + 1)",
         fill = "Microbiota Type") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  return(list(otu_comparison = p1, abundance_comparison = p2))
}

# Generate visualizations
cat("\n\n========================================\n")
cat("Generating pathogen screening visualizations...\n")
cat("========================================\n")

# Individual microbiota plots
endo_plots <- plot_pathogen_abundance(pathogen_results_endo, "Endogenous")
exo_plots <- plot_pathogen_abundance(pathogen_results_exo, "Exogenous")

# Summary comparison plots
summary_plots <- plot_pathogen_summary(pathogen_results_endo, pathogen_results_exo)

# Display plots
if (!is.null(endo_plots)) {
  print(endo_plots$abundance_plot)
  print(endo_plots$distribution_plot)
}

if (!is.null(exo_plots)) {
  print(exo_plots$abundance_plot)
  print(exo_plots$distribution_plot)
}

if (is.list(summary_plots)) {
  print(summary_plots$otu_comparison)
  print(summary_plots$abundance_comparison)
} else {
  print(summary_plots)  # This would be the empty plot
}

# ==================================================================================
# Specific Bacillus screening and classification analysis
# ==================================================================================

# Function to screen for Bacillus genus and species
screen_bacillus <- function(ps_obj, object_name) {
  cat("\n========================================\n")
  cat("Bacillus screening for:", object_name, "\n")
  cat("========================================\n")
  
  # Get taxonomy table
  tax_table_df <- as.data.frame(tax_table(ps_obj))
  
  # Check if Genus column exists
  if (!"Genus" %in% colnames(tax_table_df)) {
    cat("Warning: No 'Genus' column found in taxonomy table\n")
    return(NULL)
  }
  
  # Find Bacillus OTUs
  bacillus_otus <- tax_table_df[tax_table_df$Genus == "Bacillus" & 
                                !is.na(tax_table_df$Genus), ]
  
  if (nrow(bacillus_otus) > 0) {
    cat("\n*** BACILLUS DETECTED ***\n")
    cat("Number of Bacillus OTUs:", nrow(bacillus_otus), "\n\n")
    
    # Display full taxonomy for each Bacillus OTU
    cat("Detailed taxonomy for Bacillus OTUs:\n")
    cat("=====================================\n")
    
    for (i in 1:nrow(bacillus_otus)) {
      otu_id <- rownames(bacillus_otus)[i]
      cat("\nOTU", i, "(", otu_id, "):\n")
      
      # Print each taxonomic level
      tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      for (level in tax_levels) {
        if (level %in% colnames(bacillus_otus)) {
          value <- bacillus_otus[i, level]
          if (!is.na(value) && value != "") {
            cat("  ", level, ":", value, "\n")
          } else {
            cat("  ", level, ": [Unclassified]\n")
          }
        }
      }
    }
    
    # Get abundance data for Bacillus OTUs
    bacillus_otu_names <- rownames(bacillus_otus)
    bacillus_subset <- prune_taxa(bacillus_otu_names, ps_obj)
    
    # Calculate relative abundances
    bacillus_rel <- transform_sample_counts(bacillus_subset, function(x) x/sum(x))
    
    # Summary statistics
    cat("\n\nAbundance summary for Bacillus OTUs:\n")
    cat("====================================\n")
    otu_sums <- taxa_sums(bacillus_subset)
    cat("Total reads:", sum(otu_sums), "\n")
    cat("Mean reads per Bacillus OTU:", round(mean(otu_sums), 2), "\n")
    cat("Max reads for single Bacillus OTU:", max(otu_sums), "\n")
    cat("Min reads for single Bacillus OTU:", min(otu_sums), "\n")
    
    # Species-level summary
    species_info <- bacillus_otus$Species
    species_classified <- sum(!is.na(species_info) & species_info != "")
    cat("\nSpecies classification summary:\n")
    cat("OTUs with species-level classification:", species_classified, "out of", nrow(bacillus_otus), "\n")
    
    if (species_classified > 0) {
      cat("\nSpecies identified:\n")
      unique_species <- unique(species_info[!is.na(species_info) & species_info != ""])
      for (species in unique_species) {
        species_count <- sum(species_info == species, na.rm = TRUE)
        cat("-", species, "(", species_count, "OTUs )\n")
      }
    }
    
    # Return the Bacillus subset for visualization
    return(list(
      bacillus_otus = bacillus_otus,
      bacillus_subset = bacillus_subset,
      bacillus_rel = bacillus_rel,
      species_classified = species_classified,
      total_otus = nrow(bacillus_otus)
    ))
    
  } else {
    cat("\n*** NO BACILLUS DETECTED ***\n")
    cat("No OTUs classified as Bacillus genus found\n")
    return(NULL)
  }
}

# Function to visualize Bacillus results
plot_bacillus_results <- function(bacillus_results_endo, bacillus_results_exo) {
  
  # Create summary data
  summary_data <- data.frame(
    Microbiota = character(),
    Total_OTUs = numeric(),
    Species_Classified = numeric(),
    Total_Abundance = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(bacillus_results_endo)) {
    endo_abundance <- sum(taxa_sums(bacillus_results_endo$bacillus_subset))
    summary_data <- rbind(summary_data, data.frame(
      Microbiota = "Endogenous",
      Total_OTUs = bacillus_results_endo$total_otus,
      Species_Classified = bacillus_results_endo$species_classified,
      Total_Abundance = endo_abundance
    ))
  }
  
  if (!is.null(bacillus_results_exo)) {
    exo_abundance <- sum(taxa_sums(bacillus_results_exo$bacillus_subset))
    summary_data <- rbind(summary_data, data.frame(
      Microbiota = "Exogenous",
      Total_OTUs = bacillus_results_exo$total_otus,
      Species_Classified = bacillus_results_exo$species_classified,
      Total_Abundance = exo_abundance
    ))
  }
  
  if (nrow(summary_data) == 0) {
    cat("\nNo Bacillus detected in either microbiota type - creating summary message\n")
    
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No Bacillus genus detected\nin either endogenous or exogenous microbiota", 
               size = 16, hjust = 0.5, vjust = 0.5) +
      labs(title = "Bacillus Screening Summary - Cochliomyia macellaria") +
      theme_pub() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
    
    return(empty_plot)
  }
  
  # Create plots
  p1 <- ggplot(summary_data, aes(x = Microbiota, y = Total_OTUs, fill = Microbiota)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Total_OTUs), vjust = -0.3, size = 5) +
    labs(title = "Bacillus OTU Counts by Microbiota Type",
         x = "Microbiota Type", y = "Number of Bacillus OTUs") +
    theme_pub() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  p2 <- ggplot(summary_data, aes(x = Microbiota, y = Species_Classified, fill = Microbiota)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste(Species_Classified, "/", Total_OTUs)), vjust = -0.3, size = 4) +
    labs(title = "Bacillus Species-Level Classification",
         x = "Microbiota Type", y = "OTUs with Species Classification") +
    theme_pub() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  p3 <- ggplot(summary_data, aes(x = Microbiota, y = log10(Total_Abundance + 1), fill = Microbiota)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Total_Abundance), vjust = -0.3, size = 4) +
    labs(title = "Bacillus Total Abundance",
         x = "Microbiota Type", y = "Log10(Total Abundance + 1)") +
    theme_pub() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  return(list(otu_counts = p1, species_classification = p2, abundance = p3))
}

# Run Bacillus screening
cat("\n\n========================================\n")
cat("BACILLUS-SPECIFIC SCREENING\n")
cat("========================================\n")

bacillus_results_endo <- screen_bacillus(ps_endogenous, "Endogenous microbiota")
bacillus_results_exo <- screen_bacillus(ps_exogenous, "Exogenous microbiota")

# Generate Bacillus visualizations
cat("\n\nGenerating Bacillus-specific visualizations...\n")
bacillus_plots <- plot_bacillus_results(bacillus_results_endo, bacillus_results_exo)

# Display Bacillus plots
if (is.list(bacillus_plots)) {
  print(bacillus_plots$otu_counts)
  print(bacillus_plots$species_classification)
  print(bacillus_plots$abundance)
} else {
  print(bacillus_plots)  # This would be the empty plot
}

# ==================================================================================
# Screen for known pathogens that specifically affect blow flies
# ==================================================================================

# Define pathogens known to affect blow flies and other Diptera
fly_specific_pathogens <- list(
  # Bacterial pathogens
  bacteria = c(
    "Serratia",           # S. marcescens - causes septicemia in flies
    "Pseudomonas",        # P. aeruginosa - opportunistic pathogen in insects
    "Bacillus",           # B. thuringiensis - insecticide, B. sphaericus
    "Paenibacillus",      # P. larvae and related - affect various insects
    "Xenorhabdus",        # Entomopathogenic bacteria (symbiotic with nematodes)
    "Photorhabdus"        # Entomopathogenic bacteria (symbiotic with nematodes)
  ),
  
  # Potential fungal-related bacteria (some fungi have bacterial symbionts)
  fungal_associated = c(
    "Burkholderia",       # Often associated with fungal pathogens
    "Ralstonia"           # Some species associated with plant/insect pathogens
  ),
  
  # Opportunistic pathogens that can cause issues in stressed flies
  opportunistic = c(
    "Erwinia",            # Plant pathogen that can affect insects
    "Pantoea",            # Opportunistic pathogen in insects
    "Enterobacter",       # Can cause septicemia in insects
    "Acinetobacter"       # Opportunistic pathogen in various hosts
  )
)

# Function to screen for fly-specific pathogens
screen_fly_pathogens <- function(ps_obj, object_name) {
  cat("\n========================================\n")
  cat("Fly-specific pathogen screening for:", object_name, "\n")
  cat("========================================\n")
  
  # Get taxonomy table
  tax_table_df <- as.data.frame(tax_table(ps_obj))
  
  # Check if Genus column exists
  if (!"Genus" %in% colnames(tax_table_df)) {
    cat("Warning: No 'Genus' column found in taxonomy table\n")
    return(NULL)
  }
  
  # Combine all fly pathogen genera
  all_fly_pathogens <- unique(c(fly_specific_pathogens$bacteria, 
                               fly_specific_pathogens$fungal_associated,
                               fly_specific_pathogens$opportunistic))
  
  # Find fly pathogenic genera present
  fly_pathogenic_otus <- tax_table_df[tax_table_df$Genus %in% all_fly_pathogens & 
                                      !is.na(tax_table_df$Genus), ]
  
  if (nrow(fly_pathogenic_otus) > 0) {
    cat("\n*** FLY-SPECIFIC PATHOGENS DETECTED ***\n")
    cat("Number of fly pathogenic OTUs:", nrow(fly_pathogenic_otus), "\n\n")
    
    # Categorize found pathogens
    results_by_category <- list()
    
    for (category in names(fly_specific_pathogens)) {
      category_genera <- fly_specific_pathogens[[category]]
      category_otus <- fly_pathogenic_otus[fly_pathogenic_otus$Genus %in% category_genera, ]
      
      if (nrow(category_otus) > 0) {
        cat("=== ", toupper(category), " PATHOGENS ===\n")
        unique_genera <- unique(category_otus$Genus[!is.na(category_otus$Genus)])
        
        for (genus in unique_genera) {
          genus_count <- sum(category_otus$Genus == genus, na.rm = TRUE)
          cat("-", genus, "(", genus_count, "OTUs )\n")
          
          # Show species if available
          genus_otus <- category_otus[category_otus$Genus == genus, ]
          if ("Species" %in% colnames(genus_otus)) {
            species_info <- genus_otus$Species[!is.na(genus_otus$Species) & genus_otus$Species != ""]
            if (length(species_info) > 0) {
              unique_species <- unique(species_info)
              cat("  Species detected:", paste(unique_species, collapse = ", "), "\n")
            }
          }
        }
        cat("\n")
        results_by_category[[category]] <- category_otus
      }
    }
    
    # Get abundance data for all fly pathogenic OTUs
    fly_pathogenic_otu_names <- rownames(fly_pathogenic_otus)
    fly_pathogenic_subset <- prune_taxa(fly_pathogenic_otu_names, ps_obj)
    
    # Calculate relative abundances
    fly_pathogenic_rel <- transform_sample_counts(fly_pathogenic_subset, function(x) x/sum(x))
    
    # Summary statistics
    cat("Abundance summary for fly-specific pathogenic taxa:\n")
    cat("=================================================\n")
    otu_sums <- taxa_sums(fly_pathogenic_subset)
    cat("Total reads:", sum(otu_sums), "\n")
    cat("Mean reads per pathogenic OTU:", round(mean(otu_sums), 2), "\n")
    cat("Max reads for single pathogenic OTU:", max(otu_sums), "\n")
    
    # Calculate abundance by category
    cat("\nAbundance by pathogen category:\n")
    for (category in names(results_by_category)) {
      category_otu_names <- rownames(results_by_category[[category]])
      if (length(category_otu_names) > 0) {
        category_abundance <- sum(otu_sums[names(otu_sums) %in% category_otu_names])
        cat("-", category, ":", category_abundance, "reads\n")
      }
    }
    
    # Return results for further analysis
    return(list(
      fly_pathogenic_otus = fly_pathogenic_otus,
      fly_pathogenic_subset = fly_pathogenic_subset,
      fly_pathogenic_rel = fly_pathogenic_rel,
      results_by_category = results_by_category,
      total_otus = nrow(fly_pathogenic_otus)
    ))
    
  } else {
    cat("\n*** NO FLY-SPECIFIC PATHOGENS DETECTED ***\n")
    cat("No known fly pathogenic genera found - this is generally good for fly health!\n")
    return(NULL)
  }
}

# Function to create fly pathogen visualizations
plot_fly_pathogen_results <- function(fly_results_endo, fly_results_exo) {
  
  # Prepare summary data
  summary_data <- data.frame(
    Microbiota = character(),
    Category = character(),
    OTU_Count = numeric(),
    Total_Abundance = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process endogenous results
  if (!is.null(fly_results_endo)) {
    endo_abundance <- taxa_sums(fly_results_endo$fly_pathogenic_subset)
    
    for (category in names(fly_results_endo$results_by_category)) {
      category_otus <- fly_results_endo$results_by_category[[category]]
      category_otu_names <- rownames(category_otus)
      category_abundance <- sum(endo_abundance[names(endo_abundance) %in% category_otu_names])
      
      summary_data <- rbind(summary_data, data.frame(
        Microbiota = "Endogenous",
        Category = category,
        OTU_Count = nrow(category_otus),
        Total_Abundance = category_abundance
      ))
    }
  }
  
  # Process exogenous results
  if (!is.null(fly_results_exo)) {
    exo_abundance <- taxa_sums(fly_results_exo$fly_pathogenic_subset)
    
    for (category in names(fly_results_exo$results_by_category)) {
      category_otus <- fly_results_exo$results_by_category[[category]]
      category_otu_names <- rownames(category_otus)
      category_abundance <- sum(exo_abundance[names(exo_abundance) %in% category_otu_names])
      
      summary_data <- rbind(summary_data, data.frame(
        Microbiota = "Exogenous",
        Category = category,
        OTU_Count = nrow(category_otus),
        Total_Abundance = category_abundance
      ))
    }
  }
  
  if (nrow(summary_data) == 0) {
    cat("\nNo fly-specific pathogens detected - creating summary message\n")
    
    empty_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No known fly-specific pathogens detected\nThis suggests healthy fly microbiomes!", 
               size = 14, hjust = 0.5, vjust = 0.5) +
      labs(title = "Fly-Specific Pathogen Screening - Cochliomyia macellaria") +
      theme_pub() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank())
    
    return(empty_plot)
  }
  
  # Create visualizations
  p1 <- ggplot(summary_data, aes(x = Category, y = OTU_Count, fill = Microbiota)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = OTU_Count), position = position_dodge(width = 0.9), vjust = -0.3) +
    labs(title = "Fly-Specific Pathogenic OTUs by Category",
         x = "Pathogen Category", y = "Number of OTUs",
         fill = "Microbiota Type") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  p2 <- ggplot(summary_data, aes(x = Category, y = log10(Total_Abundance + 1), fill = Microbiota)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Fly-Specific Pathogen Abundance by Category",
         x = "Pathogen Category", y = "Log10(Total Abundance + 1)",
         fill = "Microbiota Type") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("Endogenous" = "#2E8B57", "Exogenous" = "#CD853F"))
  
  # Create a detailed breakdown if data exists
  if (!is.null(fly_results_endo) || !is.null(fly_results_exo)) {
    
    # Prepare genus-level data
    genus_data <- data.frame(
      Microbiota = character(),
      Genus = character(),
      Category = character(),
      OTU_Count = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Process both result sets for genus-level breakdown
    for (result_set in list(list(data = fly_results_endo, name = "Endogenous"),
                           list(data = fly_results_exo, name = "Exogenous"))) {
      if (!is.null(result_set$data)) {
        for (category in names(result_set$data$results_by_category)) {
          category_otus <- result_set$data$results_by_category[[category]]
          if (nrow(category_otus) > 0) {
            genus_summary <- table(category_otus$Genus)
            for (genus in names(genus_summary)) {
              genus_data <- rbind(genus_data, data.frame(
                Microbiota = result_set$name,
                Genus = genus,
                Category = category,
                OTU_Count = as.numeric(genus_summary[genus])
              ))
            }
          }
        }
      }
    }
    
    if (nrow(genus_data) > 0) {
      p3 <- ggplot(genus_data, aes(x = Genus, y = OTU_Count, fill = Category)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_wrap(~ Microbiota, scales = "free") +
        labs(title = "Fly-Specific Pathogenic Genera by Category",
             x = "Genus", y = "Number of OTUs",
             fill = "Pathogen Category") +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_viridis_d()
      
      return(list(category_otus = p1, category_abundance = p2, genus_breakdown = p3))
    }
  }
  
  return(list(category_otus = p1, category_abundance = p2))
}

# Run fly-specific pathogen screening
cat("\n\n========================================\n")
cat("FLY-SPECIFIC PATHOGEN SCREENING\n")
cat("========================================\n")

fly_pathogen_results_endo <- screen_fly_pathogens(ps_endogenous, "Endogenous microbiota")
fly_pathogen_results_exo <- screen_fly_pathogens(ps_exogenous, "Exogenous microbiota")

# Generate fly pathogen visualizations
cat("\n\nGenerating fly-specific pathogen visualizations...\n")
fly_pathogen_plots <- plot_fly_pathogen_results(fly_pathogen_results_endo, fly_pathogen_results_exo)

# Display fly pathogen plots
if (is.list(fly_pathogen_plots)) {
  print(fly_pathogen_plots$category_otus)
  print(fly_pathogen_plots$category_abundance)
  if ("genus_breakdown" %in% names(fly_pathogen_plots)) {
    print(fly_pathogen_plots$genus_breakdown)
  }
} else {
  print(fly_pathogen_plots)  # This would be the empty plot
}

