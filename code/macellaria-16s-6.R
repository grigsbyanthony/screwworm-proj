# Cochliomyia macellaria (Secondary screwworm) Microbiota Network Analysis
# Network generation and visualization for endogenous and exogenous microbiota

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

# Filter to only Cochliomyia macellaria samples
ps_macellaria <- subset_samples(ps, InsectSpecies == "Cochliomyia macellaria")

# Load required packages for network analysis
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(SpiecEasi)  # For SparCC analysis
library(igraph)     # For network analysis
library(Matrix)     # For sparse matrix operations
library(RColorBrewer)  # For color palettes
library(ggnetwork)  # For network visualization
library(intergraph) # For converting between igraph and network objects
library(VennDiagram) # For comparing networks
library(patchwork)  # For combining plots

# Print summary of the filtered dataset
print("Summary of Cochliomyia macellaria samples:")
ps_macellaria
print(sample_variables(ps_macellaria))

# Create subsets for endogenous and exogenous microbiota
ps_macellaria_endo <- subset_samples(ps_macellaria, Microbiota == "Endogenous")
ps_macellaria_endo <- prune_taxa(taxa_sums(ps_macellaria_endo) > 0, ps_macellaria_endo)  # Remove taxa with zero counts

ps_macellaria_exo <- subset_samples(ps_macellaria, Microbiota == "Exogenous")
ps_macellaria_exo <- prune_taxa(taxa_sums(ps_macellaria_exo) > 0, ps_macellaria_exo)  # Remove taxa with zero counts

# Print summary of the subsets
cat("\nEndogenous microbiota subset:", nsamples(ps_macellaria_endo), "samples,", ntaxa(ps_macellaria_endo), "taxa\n")
cat("Exogenous microbiota subset:", nsamples(ps_macellaria_exo), "samples,", ntaxa(ps_macellaria_exo), "taxa\n")

#------------------------------------------------------------------------------
# 1. Data preparation for network analysis
#------------------------------------------------------------------------------

# Function to filter low-abundance taxa
filter_taxa_by_prevalence <- function(ps, prevalence_threshold = 0.05, abundance_threshold = 0) {
  # Calculate prevalence of each taxon
  prevalence <- apply(otu_table(ps) > abundance_threshold, 1, sum) / nsamples(ps)
  
  # Filter taxa based on prevalence threshold
  keep_taxa <- names(prevalence[prevalence >= prevalence_threshold])
  ps_filtered <- prune_taxa(keep_taxa, ps)
  
  cat("Filtered from", ntaxa(ps), "to", ntaxa(ps_filtered), "taxa based on prevalence threshold of", 
      prevalence_threshold, "\n")
  
  return(ps_filtered)
}

# Filter taxa with low prevalence for network analysis
ps_macellaria_endo_network <- filter_taxa_by_prevalence(ps_macellaria_endo, prevalence_threshold = 0.1)
ps_macellaria_exo_network <- filter_taxa_by_prevalence(ps_macellaria_exo, prevalence_threshold = 0.1)

#------------------------------------------------------------------------------
# 2. SparCC correlation analysis
#------------------------------------------------------------------------------

# Function to run SparCC analysis
run_sparcc_analysis <- function(ps, iterations = 20, threshold = 0.3) {
  cat("\nRunning SparCC analysis...\n")
  
  # Extract OTU table
  otu_table_sparcc <- t(otu_table(ps))
  
  # Run SparCC
  cat("Calculating correlations with", iterations, "iterations...\n")
  sparcc_results <- sparcc(otu_table_sparcc, iter = iterations)
  
  # Extract correlation matrix
  cor_matrix <- sparcc_results$Cor
  
  # Apply threshold
  cat("Applying correlation threshold of", threshold, "...\n")
  cor_matrix[abs(cor_matrix) < threshold] <- 0
  
  # Create adjacency matrix
  adj_matrix <- cor_matrix
  diag(adj_matrix) <- 0  # Remove self-correlations
  
  # Create igraph object - use mode="max" to handle asymmetry in the matrix
  g <- graph_from_adjacency_matrix(
    adjmatrix = adj_matrix,
    mode = "max",  # Use max mode to handle potential asymmetry
    weighted = TRUE,
    diag = FALSE
  )
  
  # Add taxa names as vertex attributes
  taxa_names <- tax_table(ps)[, "Genus"]
  taxa_names <- make.unique(as.character(taxa_names))
  V(g)$name <- taxa_names
  
  # Add correlation method as graph attribute
  g$method <- "SparCC"
  
  cat("Created SparCC network with", vcount(g), "nodes and", 
      ecount(g), "edges\n")
  
  return(list(network = g, correlation = cor_matrix))
}

# Run SparCC analysis for endogenous and exogenous microbiota
cat("\n=== Running SparCC analysis for endogenous microbiota ===\n")
sparcc_results_endo <- run_sparcc_analysis(ps_macellaria_endo_network, iterations = 20, threshold = 0.3)
sparcc_network_endo <- sparcc_results_endo$network

cat("\n=== Running SparCC analysis for exogenous microbiota ===\n")
sparcc_results_exo <- run_sparcc_analysis(ps_macellaria_exo_network, iterations = 20, threshold = 0.3)
sparcc_network_exo <- sparcc_results_exo$network

#------------------------------------------------------------------------------
# 3. Network visualization
#------------------------------------------------------------------------------

# Function to clean and prepare network for visualization
prepare_network_for_viz <- function(g, min_degree = 1, remove_na = TRUE) {
  # Remove nodes with "NA" names if specified
  if (remove_na) {
    na_nodes <- which(grepl("^NA", V(g)$name))
    if (length(na_nodes) > 0) {
      g <- delete_vertices(g, na_nodes)
    }
  }
  
  # Remove isolated nodes and nodes with degree < min_degree
  # Use igraph:: prefix to explicitly call the igraph package's degree function
  degrees <- igraph::degree(g)
  low_degree_nodes <- which(degrees < min_degree)
  if (length(low_degree_nodes) > 0) {
    g <- delete_vertices(g, low_degree_nodes)
  }
  
  # Detect communities using the Louvain method
  communities <- cluster_louvain(g, weights = abs(E(g)$weight))
  
  # Add community membership to vertices
  V(g)$community <- communities$membership
  
  # Create color palette for communities
  n_communities <- max(communities$membership)
  community_colors <- brewer.pal(min(n_communities, 9), "Set1")
  if (n_communities > 9) {
    community_colors <- colorRampPalette(community_colors)(n_communities)
  }
  V(g)$color <- community_colors[V(g)$community]
  
  # Calculate node properties
  V(g)$degree <- igraph::degree(g)
  
  # Calculate betweenness using absolute weights to handle negative correlations
  # Create a copy of the graph with absolute weights for betweenness calculation
  g_abs <- g
  E(g_abs)$weight <- abs(E(g)$weight)
  V(g)$betweenness <- betweenness(g_abs, weights = E(g_abs)$weight, normalized = TRUE)
  
  # Set vertex sizes based on degree
  V(g)$size <- pmin(10, 2 + 3 * log1p(V(g)$degree))
  
  # Set edge properties
  E(g)$width <- abs(E(g)$weight) * 2
  E(g)$color <- ifelse(E(g)$weight > 0, 
                      rgb(0, 0, 1, alpha = 0.6),  # Blue for positive
                      rgb(1, 0, 0, alpha = 0.6))  # Red for negative
  
  # Store communities object in graph
  g$communities <- communities
  
  return(g)
}

# Prepare networks for visualization
sparcc_network_endo_viz <- prepare_network_for_viz(sparcc_network_endo, min_degree = 2)
sparcc_network_exo_viz <- prepare_network_for_viz(sparcc_network_exo, min_degree = 2)

# Function to visualize network using ggplot2
visualize_network_ggplot <- function(g, title = "") {
  # Check if network is empty
  if (vcount(g) == 0 || ecount(g) == 0) {
    warning("Network is empty, cannot visualize")
    # Return an empty plot with a message
    return(ggplot() + 
             annotate("text", x = 0, y = 0, label = "Empty network - no visualization possible") +
             theme_void())
  }
  
  # Ensure all required vertex attributes exist
  if (!"community" %in% vertex_attr_names(g)) {
    V(g)$community <- 1  # Default all to same community
  }
  if (!"size" %in% vertex_attr_names(g)) {
    V(g)$size <- 5  # Default size
  }
  if (!"name" %in% vertex_attr_names(g)) {
    V(g)$name <- paste0("Node", 1:vcount(g))  # Default names
  }
  if (!"color" %in% vertex_attr_names(g)) {
    V(g)$color <- "lightblue"  # Default color
  }
  
  # Ensure all required edge attributes exist
  if (!"width" %in% edge_attr_names(g)) {
    E(g)$width <- 1  # Default width
  }
  
  # Create a layout manually to ensure consistency
  set.seed(123)  # For reproducibility
  
  # Create a copy of the graph with absolute weights for layout calculation
  g_layout <- g
  E(g_layout)$weight <- abs(E(g)$weight)
  
  # Use the absolute weight graph for layout
  layout <- layout_with_fr(g_layout)
  
  # Create a data frame for nodes
  nodes_df <- data.frame(
    id = 1:vcount(g),
    name = V(g)$name,
    community = V(g)$community,
    size = V(g)$size,
    color = V(g)$color,
    x = layout[, 1],
    y = layout[, 2],
    stringsAsFactors = FALSE
  )
  
  # Create a data frame for edges - using a safer approach
  edge_list <- get.edgelist(g, names = FALSE)
  
  edges_df <- data.frame(
    from = edge_list[, 1],
    to = edge_list[, 2],
    weight = E(g)$weight,
    width = E(g)$width,
    color = E(g)$color,
    stringsAsFactors = FALSE
  )
  
  # Add coordinates to edges
  edges_df$x <- layout[edges_df$from, 1]
  edges_df$y <- layout[edges_df$from, 2]
  edges_df$xend <- layout[edges_df$to, 1]
  edges_df$yend <- layout[edges_df$to, 2]
  
  # Create plot
  p <- ggplot() +
    # Add edges
    geom_segment(data = edges_df, 
                aes(x = x, y = y, xend = xend, yend = yend, 
                    color = weight, linewidth = width),  # Using linewidth instead of size
                alpha = 0.7) +
    # Add nodes
    geom_point(data = nodes_df,
              aes(x = x, y = y, size = size, fill = as.factor(community)),
              shape = 21, color = "black", alpha = 0.8) +
    # Add node labels
    ggrepel::geom_text_repel(data = nodes_df,
                           aes(x = x, y = y, label = name),
                           size = 3, max.overlaps = 20,
                           family = "CMU Sans Serif") +
    # Set scales
    scale_fill_brewer(palette = "Set1", name = "Community") +
    scale_color_gradient2(
      low = "red", high = "blue", mid = "gray", midpoint = 0,
      name = "Correlation"
    ) +
    scale_size_continuous(name = "Node Size") +
    scale_linewidth_continuous(name = "Edge Width") +
    # Set labels
    labs(
      title = title,
      subtitle = paste0(vcount(g), " taxa, ", ecount(g), " correlations")
    ) +
    # Set theme
    theme_void() +
    theme(
      text = element_text(family = "CMU Sans Serif"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, family = "CMU Sans Serif"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, family = "CMU Sans Serif"),
      legend.title = element_text(family = "CMU Sans Serif"),
      legend.text = element_text(family = "CMU Sans Serif"),
      legend.position = "right"
    )
  
  return(p)
}

# Create network visualizations
sparcc_network_endo_plot <- visualize_network_ggplot(
  sparcc_network_endo_viz, 
  title = "Cochliomyia macellaria Endogenous Microbiota Network"
)

sparcc_network_exo_plot <- visualize_network_ggplot(
  sparcc_network_exo_viz, 
  title = "Cochliomyia macellaria Exogenous Microbiota Network"
)

# Display the plots
print(sparcc_network_endo_plot)
print(sparcc_network_exo_plot)

# Save the plots
ggsave("macellaria_endogenous_network.png", 
       plot = sparcc_network_endo_plot,
       width = 12, height = 10, 
       dpi = 300, 
       bg = "white")

ggsave("macellaria_exogenous_network.png", 
       plot = sparcc_network_exo_plot,
       width = 12, height = 10, 
       dpi = 300, 
       bg = "white")

#------------------------------------------------------------------------------
# 3.1 Create reduced networks with only top 75 most important nodes
#------------------------------------------------------------------------------

# Function to create a reduced network with only the top N most important nodes
create_reduced_network <- function(g, top_n = 75) {
  cat("\nCreating reduced network with top", top_n, "most important nodes...\n")
  
  # Calculate centrality measures for all nodes
  # Create a copy of the graph with absolute weights for centrality calculations
  g_abs <- g
  E(g_abs)$weight <- abs(E(g)$weight)
  
  # Calculate centrality measures
  degree_centrality <- igraph::degree(g)
  betweenness_centrality <- igraph::betweenness(g_abs, weights = E(g_abs)$weight, normalized = TRUE)
  closeness_centrality <- igraph::closeness(g_abs, weights = E(g_abs)$weight, normalized = TRUE)
  eigenvector_centrality <- igraph::eigen_centrality(g_abs, weights = E(g_abs)$weight)$vector
  
  # Get vertex names and handle missing values
  vertex_names <- V(g)$name
  # Replace NA values with placeholder names
  if(any(is.na(vertex_names))) {
    cat("Warning: Some vertex names are NA. Replacing with placeholder names.\n")
    na_indices <- which(is.na(vertex_names))
    vertex_names[na_indices] <- paste0("Taxon_", na_indices)
  }
  
  # Create a data frame with node IDs and centrality measures
  # Explicitly set row.names = NULL to avoid issues with missing values
  node_centrality <- data.frame(
    node_id = 1:vcount(g),
    name = vertex_names,
    degree = degree_centrality,
    betweenness = betweenness_centrality,
    closeness = closeness_centrality,
    eigenvector = eigenvector_centrality,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # Calculate overall importance score (mean of normalized centrality measures)
  # Handle potential errors in scaling
  tryCatch({
    node_centrality$importance <- apply(node_centrality[, c("degree", "betweenness", "closeness", "eigenvector")], 1, 
                                       function(x) mean(scale(x)))
  }, error = function(e) {
    cat("Warning: Error in scaling centrality measures. Using unscaled average instead.\n")
    node_centrality$importance <- rowMeans(node_centrality[, c("degree", "betweenness", "closeness", "eigenvector")])
  })
  
  # Sort by importance score and get top N nodes
  top_nodes <- node_centrality %>%
    arrange(desc(importance)) %>%
    head(top_n) %>%
    pull(node_id)
  
  # Create a subgraph with only the top nodes
  g_reduced <- induced_subgraph(g, top_nodes)
  
  # Return the reduced network
  return(g_reduced)
}

# Create reduced networks
sparcc_network_endo_reduced <- create_reduced_network(sparcc_network_endo, top_n = 75)
sparcc_network_exo_reduced <- create_reduced_network(sparcc_network_exo, top_n = 75)

# Prepare reduced networks for visualization
sparcc_network_endo_reduced_viz <- prepare_network_for_viz(sparcc_network_endo_reduced, min_degree = 1)
sparcc_network_exo_reduced_viz <- prepare_network_for_viz(sparcc_network_exo_reduced, min_degree = 1)

# Create visualizations for reduced networks
sparcc_network_endo_reduced_plot <- visualize_network_ggplot(
  sparcc_network_endo_reduced_viz, 
  title = "Cochliomyia macellaria Endogenous Microbiota Network (Top 75 Taxa)"
)

sparcc_network_exo_reduced_plot <- visualize_network_ggplot(
  sparcc_network_exo_reduced_viz, 
  title = "Cochliomyia macellaria Exogenous Microbiota Network (Top 75 Taxa)"
)

# Display the reduced network plots
print(sparcc_network_endo_reduced_plot)
print(sparcc_network_exo_reduced_plot)

# Save the reduced network plots
ggsave("macellaria_endogenous_network_top75.png", 
       plot = sparcc_network_endo_reduced_plot,
       width = 12, height = 10, 
       dpi = 300, 
       bg = "white")

ggsave("macellaria_exogenous_network_top75.png", 
       plot = sparcc_network_exo_reduced_plot,
       width = 12, height = 10, 
       dpi = 300, 
       bg = "white")

#------------------------------------------------------------------------------
# 4. Key taxa identification for each network
#------------------------------------------------------------------------------

# Function to identify keystone species based on centrality measures
identify_key_taxa <- function(g, top_n = 10) {
  # Calculate centrality measures
  degree_centrality <- igraph::degree(g)
  
  # Create a copy of the graph with absolute weights for centrality calculations
  g_abs <- g
  E(g_abs)$weight <- abs(E(g)$weight)
  
  # Calculate centrality measures using absolute weights
  betweenness_centrality <- igraph::betweenness(g_abs, weights = E(g_abs)$weight, normalized = TRUE)
  closeness_centrality <- igraph::closeness(g_abs, weights = E(g_abs)$weight, normalized = TRUE)
  eigenvector_centrality <- igraph::eigen_centrality(g_abs, weights = E(g_abs)$weight)$vector
  
  # Create a sequence of IDs to use as row names
  ids <- seq_len(vcount(g))
  
  # Get vertex names and handle missing values
  vertex_names <- V(g)$name
  # Replace NA or NULL values with placeholder names
  if(any(is.na(vertex_names))) {
    cat("Warning: Some vertex names are NA. Replacing with placeholder names.\n")
    na_indices <- which(is.na(vertex_names))
    vertex_names[na_indices] <- paste0("Taxon_", na_indices)
  }
  
  # Create a data frame with explicit row names
  centrality_df <- data.frame(
    row.names = ids,  # Use explicit row names
    Taxon = vertex_names,
    Degree = degree_centrality,
    Betweenness = betweenness_centrality,
    Closeness = closeness_centrality,
    Eigenvector = eigenvector_centrality,
    stringsAsFactors = FALSE,
    check.names = FALSE  # Don't modify column names
  )
  
  # Calculate overall importance score (mean of normalized centrality measures)
  # Handle potential errors in scaling
  tryCatch({
    centrality_df$Importance <- apply(centrality_df[, c("Degree", "Betweenness", "Closeness", "Eigenvector")], 1, 
                                     function(x) mean(scale(x)))
  }, error = function(e) {
    cat("Warning: Error in scaling centrality measures. Using unscaled average instead.\n")
    centrality_df$Importance <<- rowMeans(centrality_df[, c("Degree", "Betweenness", "Closeness", "Eigenvector")])
  })
  
  # Sort by importance score
  centrality_df <- centrality_df[order(-centrality_df$Importance), ]
  
  # Return top N keystone species
  return(centrality_df[1:min(top_n, nrow(centrality_df)), ])
}

# Identify key taxa for each network
key_taxa_endo <- identify_key_taxa(sparcc_network_endo, top_n = 10)
key_taxa_exo <- identify_key_taxa(sparcc_network_exo, top_n = 10)

# Print key taxa
cat("\nTop 10 key taxa in the endogenous microbiota network:\n")
print(key_taxa_endo[, c("Taxon", "Degree", "Betweenness", "Importance")])

cat("\nTop 10 key taxa in the exogenous microbiota network:\n")
print(key_taxa_exo[, c("Taxon", "Degree", "Betweenness", "Importance")])

# Create bar plots of key taxa importance
key_taxa_endo_plot <- ggplot(key_taxa_endo, aes(x = reorder(Taxon, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#1B9E77") +
  coord_flip() +
  labs(
    title = "Key Taxa in Endogenous Microbiota Network",
    x = "Taxon",
    y = "Importance Score"
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure bars start at y-axis
  theme_pub()

key_taxa_exo_plot <- ggplot(key_taxa_exo, aes(x = reorder(Taxon, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "#D95F02") +
  coord_flip() +
  labs(
    title = "Key Taxa in Exogenous Microbiota Network",
    x = "Taxon",
    y = "Importance Score"
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Ensure bars start at y-axis
  theme_pub()

# Display the plots
print(key_taxa_endo_plot)
print(key_taxa_exo_plot)

# Save the plots
ggsave("macellaria_endogenous_key_taxa.png", 
       plot = key_taxa_endo_plot,
       width = 10, height = 6, 
       dpi = 300, 
       bg = "white")

ggsave("macellaria_exogenous_key_taxa.png", 
       plot = key_taxa_exo_plot,
       width = 10, height = 6, 
       dpi = 300, 
       bg = "white")

# Save key taxa to CSV
write.csv(key_taxa_endo, "macellaria_endogenous_key_taxa.csv", row.names = FALSE)
write.csv(key_taxa_exo, "macellaria_exogenous_key_taxa.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# 5. Community structure analysis for each network
#------------------------------------------------------------------------------

# Function to analyze community structure
analyze_communities <- function(g) {
  # Get communities from graph
  communities <- g$communities
  
  # If communities not already calculated, detect them
  if (is.null(communities)) {
    communities <- igraph::cluster_louvain(g, weights = abs(igraph::E(g)$weight))
  }
  
  # Calculate modularity
  modularity_score <- igraph::modularity(communities)
  
  # Get community membership
  membership <- communities$membership
  
  # Create a copy of the graph with absolute weights for betweenness calculation
  g_abs <- g
  igraph::E(g_abs)$weight <- abs(igraph::E(g)$weight)
  
  # Get vertex names and handle missing values
  vertex_names <- igraph::V(g)$name
  # Replace NA or NULL values with placeholder names
  if(any(is.na(vertex_names))) {
    cat("Warning: Some vertex names are NA in analyze_communities. Replacing with placeholder names.\n")
    na_indices <- which(is.na(vertex_names))
    vertex_names[na_indices] <- paste0("Taxon_", na_indices)
  }
  
  # Create a sequence of IDs to use as row names
  ids <- seq_len(igraph::vcount(g))
  
  # Create community dataframe with explicit row names
  community_df <- data.frame(
    row.names = ids,  # Use explicit row names
    Taxon = vertex_names,
    Community = membership,
    Degree = igraph::degree(g),
    Betweenness = igraph::betweenness(g_abs, weights = igraph::E(g_abs)$weight, normalized = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE  # Don't modify column names
  )
  
  # Calculate community sizes
  community_sizes <- table(membership)
  
  # Calculate statistics for each community
  community_stats <- data.frame(
    Community = as.integer(names(community_sizes)),
    Size = as.integer(community_sizes),
    stringsAsFactors = FALSE
  )
  
  # For each community, find the hub (node with highest degree)
  community_stats$Hub <- sapply(community_stats$Community, function(comm) {
    comm_nodes <- community_df$Taxon[community_df$Community == comm]
    comm_degrees <- community_df$Degree[community_df$Community == comm]
    if(length(comm_nodes) == 0 || length(comm_degrees) == 0) {
      return(NA)  # Handle empty communities
    }
    return(comm_nodes[which.max(comm_degrees)])
  })
  
  # Sort by community size
  community_stats <- community_stats[order(-community_stats$Size), ]
  
  return(list(
    communities = communities,
    modularity = modularity_score,
    community_df = community_df,
    community_stats = community_stats
  ))
}

# Analyze community structure for each network
community_analysis_endo <- analyze_communities(sparcc_network_endo_viz)
community_analysis_exo <- analyze_communities(sparcc_network_exo_viz)

# Print community statistics
cat("\nCommunity structure analysis for endogenous microbiota:\n")
cat("Modularity:", round(community_analysis_endo$modularity, 3), "\n")
cat("Number of communities:", length(unique(community_analysis_endo$community_df$Community)), "\n\n")
cat("Community statistics:\n")
print(community_analysis_endo$community_stats)

cat("\nCommunity structure analysis for exogenous microbiota:\n")
cat("Modularity:", round(community_analysis_exo$modularity, 3), "\n")
cat("Number of communities:", length(unique(community_analysis_exo$community_df$Community)), "\n\n")
cat("Community statistics:\n")
print(community_analysis_exo$community_stats)

# Create bar plots of community sizes
community_size_endo_plot <- ggplot(community_analysis_endo$community_stats, 
                             aes(x = reorder(as.factor(Community), -Size), y = Size)) +
  geom_bar(stat = "identity", fill = "#1B9E77") +
  labs(
    title = "Community Sizes - Endogenous Microbiota",
    x = "Community",
    y = "Number of Taxa"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

community_size_exo_plot <- ggplot(community_analysis_exo$community_stats, 
                             aes(x = reorder(as.factor(Community), -Size), y = Size)) +
  geom_bar(stat = "identity", fill = "#D95F02") +
  labs(
    title = "Community Sizes - Exogenous Microbiota",
    x = "Community",
    y = "Number of Taxa"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plots
print(community_size_endo_plot)
print(community_size_exo_plot)

# Save the plots
ggsave("macellaria_endogenous_community_sizes.png", 
       plot = community_size_endo_plot,
       width = 8, height = 6, 
       dpi = 300, 
       bg = "white")

ggsave("macellaria_exogenous_community_sizes.png", 
       plot = community_size_exo_plot,
       width = 8, height = 6, 
       dpi = 300, 
       bg = "white")

# Save community data
write.csv(community_analysis_endo$community_df, "macellaria_endogenous_community_membership.csv", row.names = FALSE)
write.csv(community_analysis_endo$community_stats, "macellaria_endogenous_community_statistics.csv", row.names = FALSE)
write.csv(community_analysis_exo$community_df, "macellaria_exogenous_community_membership.csv", row.names = FALSE)
write.csv(community_analysis_exo$community_stats, "macellaria_exogenous_community_statistics.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# 6. Network statistics for each network
#------------------------------------------------------------------------------

# Function to calculate network statistics
calculate_network_stats <- function(g, community_analysis) {
  # Basic network properties
  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  density <- edge_density(g)
  
  # Create a copy of the graph with absolute weights for calculations
  g_abs <- g
  E(g_abs)$weight <- abs(E(g)$weight)
  
  # Check if network is connected
  is_connected <- is.connected(g)
  
  # Calculate diameter and average path length if connected
  if (is_connected) {
    diameter <- diameter(g, weights = NA)
    avg_path_length <- mean_distance(g, weights = NA)
  } else {
    # For disconnected networks, calculate for largest component
    components <- decompose(g)
    largest_component <- components[[which.max(sapply(components, vcount))]]
    diameter <- diameter(largest_component, weights = NA)
    avg_path_length <- mean_distance(largest_component, weights = NA)
  }
  
  # Clustering coefficient
  transitivity <- transitivity(g, type = "global")
  
  # Degree distribution
  degree_dist <- degree_distribution(g)
  
  # Positive vs negative correlations
  positive_edges <- sum(E(g)$weight > 0)
  negative_edges <- sum(E(g)$weight < 0)
  
  # Create a data frame of statistics
  stats_df <- data.frame(
    Metric = c("Nodes", "Edges", "Density", "Diameter", "Average Path Length",
               "Clustering Coefficient", "Positive Correlations", "Negative Correlations",
               "Modularity"),
    Value = c(n_nodes, n_edges, density, diameter, avg_path_length,
              transitivity, positive_edges, negative_edges,
              community_analysis$modularity)
  )
  
  return(stats_df)
}

# Calculate network statistics for each network
network_stats_endo <- calculate_network_stats(sparcc_network_endo_viz, community_analysis_endo)
network_stats_exo <- calculate_network_stats(sparcc_network_exo_viz, community_analysis_exo)

# Print network statistics
cat("\nNetwork statistics for endogenous microbiota:\n")
print(network_stats_endo)

cat("\nNetwork statistics for exogenous microbiota:\n")
print(network_stats_exo)

# Combine statistics for comparison
network_stats_combined <- data.frame(
  Metric = network_stats_endo$Metric,
  Endogenous = network_stats_endo$Value,
  Exogenous = network_stats_exo$Value
)

# Print combined statistics
cat("\nComparison of network statistics:\n")
print(network_stats_combined)

# Save network statistics
write.csv(network_stats_endo, "macellaria_endogenous_network_statistics.csv", row.names = FALSE)
write.csv(network_stats_exo, "macellaria_exogenous_network_statistics.csv", row.names = FALSE)
write.csv(network_stats_combined, "macellaria_network_statistics_comparison.csv", row.names = FALSE)

#------------------------------------------------------------------------------
# 7. Network comparison between endogenous and exogenous microbiota
#------------------------------------------------------------------------------

# Function to compare taxa between networks
compare_network_taxa <- function(g1, g2, g1_name = "Network 1", g2_name = "Network 2") {
  # Get taxa names from each network
  taxa1 <- V(g1)$name
  taxa2 <- V(g2)$name
  
  # Handle NA values in taxa names
  if(any(is.na(taxa1))) {
    cat("Warning: NA values found in", g1_name, "taxa names. Replacing with placeholder names.\n")
    na_indices1 <- which(is.na(taxa1))
    taxa1[na_indices1] <- paste0(g1_name, "_Taxon_", na_indices1)
  }
  
  if(any(is.na(taxa2))) {
    cat("Warning: NA values found in", g2_name, "taxa names. Replacing with placeholder names.\n")
    na_indices2 <- which(is.na(taxa2))
    taxa2[na_indices2] <- paste0(g2_name, "_Taxon_", na_indices2)
  }
  
  # Find common and unique taxa
  common_taxa <- intersect(taxa1, taxa2)
  unique_taxa1 <- setdiff(taxa1, taxa2)
  unique_taxa2 <- setdiff(taxa2, taxa1)
  
  # Calculate percentages
  total_unique <- length(unique(c(taxa1, taxa2)))
  percent_common <- round(length(common_taxa) / total_unique * 100, 1)
  percent_unique1 <- round(length(unique_taxa1) / total_unique * 100, 1)
  percent_unique2 <- round(length(unique_taxa2) / total_unique * 100, 1)
  
  # Create summary
  cat("\nComparison of taxa between", g1_name, "and", g2_name, "networks:\n")
  cat("Total unique taxa across both networks:", total_unique, "\n")
  cat("Common taxa:", length(common_taxa), "(", percent_common, "%)\n")
  cat("Taxa unique to", g1_name, ":", length(unique_taxa1), "(", percent_unique1, "%)\n")
  cat("Taxa unique to", g2_name, ":", length(unique_taxa2), "(", percent_unique2, "%)\n")
  
  # Create Venn diagram data
  venn_data <- list(
    Network1 = taxa1,
    Network2 = taxa2
  )
  names(venn_data) <- c(g1_name, g2_name)
  
  # Create Venn diagram with error handling
  tryCatch({
    venn_plot <- venn.diagram(
      x = venn_data,
      filename = NULL,
      fill = c("#1B9E77", "#D95F02"),  # Endogenous and Exogenous colors
      alpha = 0.5,
      label.col = "black",
      cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontfamily = "sans",
      cat.col = c("#1B9E77", "#D95F02"),
      cat.fontface = "bold",
      margin = 0.1
    )
    
    # Return results with venn plot
    return(list(
      common_taxa = common_taxa,
      unique_taxa1 = unique_taxa1,
      unique_taxa2 = unique_taxa2,
      percent_common = percent_common,
      venn_plot = venn_plot
    ))
  }, error = function(e) {
    cat("Error creating Venn diagram:", conditionMessage(e), "\n")
    cat("Returning results without Venn diagram.\n")
    
    # Return results without venn plot
    return(list(
      common_taxa = common_taxa,
      unique_taxa1 = unique_taxa1,
      unique_taxa2 = unique_taxa2,
      percent_common = percent_common,
      venn_plot = NULL
    ))
  })
}

# Compare endogenous and exogenous networks
cat("\n=== Comparing Endogenous and Exogenous Microbiota Networks ===\n")
network_comparison <- compare_network_taxa(
  sparcc_network_endo, 
  sparcc_network_exo, 
  g1_name = "Endogenous", 
  g2_name = "Exogenous"
)

# Display Venn diagram if available
if (!is.null(network_comparison$venn_plot)) {
  grid::grid.newpage()
  grid::grid.draw(network_comparison$venn_plot)
  
  # Save Venn diagram
  png("macellaria_endo_vs_exo_network_taxa_venn.png", width = 8, height = 6, units = "in", res = 300)
  grid::grid.draw(network_comparison$venn_plot)
  dev.off()
}

# Compare network properties
cat("\n=== Network Property Comparison: Endogenous vs Exogenous ===\n")
network_property_comparison <- data.frame(
  Property = c("Nodes", "Edges", "Density", "Modularity", "Clustering Coefficient",
               "Positive Correlations", "Negative Correlations"),
  Endogenous = c(vcount(sparcc_network_endo_viz), ecount(sparcc_network_endo_viz),
          edge_density(sparcc_network_endo_viz), community_analysis_endo$modularity,
          transitivity(sparcc_network_endo_viz),
          sum(E(sparcc_network_endo_viz)$weight > 0), sum(E(sparcc_network_endo_viz)$weight < 0)),
  Exogenous = c(vcount(sparcc_network_exo_viz), ecount(sparcc_network_exo_viz),
             edge_density(sparcc_network_exo_viz), community_analysis_exo$modularity,
             transitivity(sparcc_network_exo_viz),
             sum(E(sparcc_network_exo_viz)$weight > 0), sum(E(sparcc_network_exo_viz)$weight < 0))
)

# Print comparison
print(network_property_comparison)

# Save comparison
write.csv(network_property_comparison, "macellaria_network_property_comparison.csv", row.names = FALSE)

# Create a bar plot comparing network properties
network_property_long <- network_property_comparison %>%
  tidyr::pivot_longer(cols = c("Endogenous", "Exogenous"), names_to = "Network", values_to = "Value")

# Plot only numeric properties that make sense to compare
plot_properties <- c("Nodes", "Edges", "Density", "Modularity", "Clustering Coefficient")
network_property_plot <- ggplot(
  network_property_long %>% filter(Property %in% plot_properties),
  aes(x = Property, y = Value, fill = Network)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  labs(
    title = "Network Property Comparison: Endogenous vs Exogenous",
    x = "Property",
    y = "Value"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(network_property_plot)

# Save the plot
ggsave("macellaria_network_property_comparison.png", 
       plot = network_property_plot,
       width = 10, height = 6, 
       dpi = 300, 
       bg = "white")

# Compare key taxa between networks
cat("\n=== Key Taxa Comparison: Endogenous vs Exogenous ===\n")

# Find common key taxa
common_key_taxa <- intersect(key_taxa_endo$Taxon, key_taxa_exo$Taxon)
cat("Common key taxa:", length(common_key_taxa), "\n")
if(length(common_key_taxa) > 0) {
  cat("Common key taxa list:", paste(common_key_taxa, collapse = ", "), "\n")
}

# Create a combined key taxa data frame for visualization
key_taxa_combined <- bind_rows(
  key_taxa_endo %>% mutate(Network = "Endogenous"),
  key_taxa_exo %>% mutate(Network = "Exogenous")
)

# Create a plot comparing key taxa importance
key_taxa_comparison_plot <- ggplot(key_taxa_combined, 
                                 aes(x = reorder(Taxon, Importance), y = Importance, fill = Network)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Endogenous" = "#1B9E77", "Exogenous" = "#D95F02")) +
  coord_flip() +
  labs(
    title = "Key Taxa Importance: Endogenous vs Exogenous",
    x = "Taxon",
    y = "Importance Score"
  ) +
  scale_y_continuous(expand = c(0,0)) +
  theme_pub()

# Display the plot
print(key_taxa_comparison_plot)

# Save the plot
ggsave("macellaria_key_taxa_comparison.png", 
       plot = key_taxa_comparison_plot,
       width = 12, height = 8, 
       dpi = 300, 
       bg = "white")

#------------------------------------------------------------------------------
# 8. Export network data for Gephi visualization
#------------------------------------------------------------------------------

# Function to export network data for Gephi
export_network_for_gephi <- function(g, output_prefix = "network_for_gephi") {
  cat("\nExporting network data for Gephi visualization...\n")
  
  # Check which attributes exist in the graph
  vertex_attrs <- igraph::vertex_attr_names(g)
  
  # Get vertex names and handle missing values
  vertex_names <- igraph::V(g)$name
  # Replace NA or NULL values with placeholder names
  if(any(is.na(vertex_names))) {
    cat("Warning: Some vertex names are NA in export_network_for_gephi. Replacing with placeholder names.\n")
    na_indices <- which(is.na(vertex_names))
    vertex_names[na_indices] <- paste0("Taxon_", na_indices)
    # Update the graph with the new names
    igraph::V(g)$name <- vertex_names
  }
  
  # Create node IDs vector
  node_ids <- 1:igraph::vcount(g)
  
  # Calculate betweenness if needed
  g_abs <- g
  igraph::E(g_abs)$weight <- abs(igraph::E(g)$weight)
  betweenness_values <- igraph::betweenness(g_abs, weights = igraph::E(g_abs)$weight, normalized = TRUE)
  
  # Create community values (default to 1 if not present)
  community_values <- if("community" %in% vertex_attrs) igraph::V(g)$community else rep(1, igraph::vcount(g))
  
  # Create size values (default to 5 if not present)
  size_values <- if("size" %in% vertex_attrs) igraph::V(g)$size else rep(5, igraph::vcount(g))
  
  # Create degree values
  degree_values <- igraph::degree(g)
  
  # Create nodes data frame column by column to avoid row name issues
  nodes_df <- data.frame(
    Id = node_ids,
    Label = as.character(vertex_names),  # Ensure character type
    Degree = as.numeric(degree_values),  # Ensure numeric type
    Community = as.numeric(community_values),  # Ensure numeric type
    Betweenness = as.numeric(betweenness_values),  # Ensure numeric type
    Size = as.numeric(size_values),  # Ensure numeric type
    stringsAsFactors = FALSE
  )
  
  # Create edge data
  edge_list <- igraph::get.edgelist(g, names = FALSE)
  edge_weights <- igraph::E(g)$weight
  edge_types <- ifelse(edge_weights > 0, "Positive", "Negative")
  
  # Create edges data frame column by column
  edges_df <- data.frame(
    Source = as.integer(edge_list[, 1]),
    Target = as.integer(edge_list[, 2]),
    Weight = as.numeric(edge_weights),
    Type = as.character(edge_types),
    stringsAsFactors = FALSE
  )
  
  # Save nodes and edges to CSV files
  nodes_file <- paste0(output_prefix, "_nodes.csv")
  edges_file <- paste0(output_prefix, "_edges.csv")
  
  tryCatch({
    write.csv(nodes_df, nodes_file, row.names = FALSE)
    write.csv(edges_df, edges_file, row.names = FALSE)
    
    cat("Nodes file saved:", nodes_file, "\n")
    cat("Edges file saved:", edges_file, "\n")
  }, error = function(e) {
    cat("Error saving CSV files:", conditionMessage(e), "\n")
  })
}

# Export networks for Gephi
export_network_for_gephi(sparcc_network_endo_viz, output_prefix = "macellaria_endogenous_network")
export_network_for_gephi(sparcc_network_exo_viz, output_prefix = "macellaria_exogenous_network")

# Export reduced networks for Gephi
export_network_for_gephi(sparcc_network_endo_reduced_viz, output_prefix = "macellaria_endogenous_network_top75")
export_network_for_gephi(sparcc_network_exo_reduced_viz, output_prefix = "macellaria_exogenous_network_top75")

# Print summary of analysis
cat("\nNetwork analysis of Cochliomyia macellaria microbiota completed.\n")
cat("Output files:\n")
cat("- macellaria_endogenous_network.png\n")
cat("- macellaria_exogenous_network.png\n")
cat("- macellaria_endogenous_network_top75.png\n")
cat("- macellaria_exogenous_network_top75.png\n")
cat("- macellaria_endogenous_key_taxa.png\n")
cat("- macellaria_exogenous_key_taxa.png\n")
cat("- macellaria_endogenous_community_sizes.png\n")
cat("- macellaria_exogenous_community_sizes.png\n")
cat("- macellaria_endo_vs_exo_network_taxa_venn.png\n")
cat("- macellaria_network_property_comparison.png\n")
cat("- macellaria_key_taxa_comparison.png\n")
cat("- macellaria_endogenous_key_taxa.csv\n")
cat("- macellaria_exogenous_key_taxa.csv\n")
cat("- macellaria_endogenous_community_membership.csv\n")
cat("- macellaria_endogenous_community_statistics.csv\n")
cat("- macellaria_exogenous_community_membership.csv\n")
cat("- macellaria_exogenous_community_statistics.csv\n")
cat("- macellaria_endogenous_network_statistics.csv\n")
cat("- macellaria_exogenous_network_statistics.csv\n")
cat("- macellaria_network_statistics_comparison.csv\n")
cat("- macellaria_network_property_comparison.csv\n")
cat("- macellaria_endogenous_network_nodes.csv\n")
cat("- macellaria_endogenous_network_edges.csv\n")
cat("- macellaria_exogenous_network_nodes.csv\n")
cat("- macellaria_exogenous_network_edges.csv\n")
cat("- macellaria_endogenous_network_top75_nodes.csv\n")
cat("- macellaria_endogenous_network_top75_edges.csv\n")
cat("- macellaria_exogenous_network_top75_nodes.csv\n")
cat("- macellaria_exogenous_network_top75_edges.csv\n")
