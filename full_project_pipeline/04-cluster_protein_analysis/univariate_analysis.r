prot_df <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")
output_dir <- "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output"

### Function to perform cluster-specific analysis
analyze_cluster <- function(cluster_num) {
    # Get the appropriate binary outcome
    binary_col <- paste0("cluster", cluster_num, "_binary")
    
    # Initialize results dataframe with numeric columns having higher precision
    results <- data.frame(
        protein = character(),
        effect_size = numeric(),
        ci_lower = numeric(),
        ci_upper = numeric(),
        p_value = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Get protein column names (exclude index, Cluster, and binary columns)
    protein_cols <- names(prot_df)[!names(prot_df) %in% c("index", "Cluster", 
                                                         "cluster1_binary", "cluster2_binary", 
                                                         "cluster3_binary", "cluster4_binary")]
    
    # Add tracking variables
    n_proteins <- length(protein_cols)
    failed_proteins <- character()
    
    # Perform linear regression for each protein
    for(protein in protein_cols) {
        tryCatch({
            # Fit model with protein as outcome
            model <- lm(prot_df[[protein]] ~ prot_df[[binary_col]])
            
            # Extract results
            coef_summary <- summary(model)$coefficients[2,]
            
            # Calculate confidence intervals
            ci <- confint(model)[2,]
            
            # Store results with full precision
            results <- rbind(results, data.frame(
                protein = protein,
                effect_size = coef_summary[1],    # This is now the mean difference between groups
                ci_lower = ci[1],
                ci_upper = ci[2],
                p_value = coef_summary[4]
            ))
        }, error = function(e) {
            warning(sprintf("Error analyzing protein %s: %s", protein, e$message))
            failed_proteins <- c(failed_proteins, protein)
        })
    }
    
    # Report on analysis completion
    cat(sprintf("\nCluster %d Analysis Summary:\n", cluster_num))
    cat(sprintf("Total proteins analyzed: %d\n", n_proteins))
    cat(sprintf("Successful analyses: %d\n", nrow(results)))
    if(length(failed_proteins) > 0) {
        cat("Failed proteins:\n")
        print(failed_proteins)
    }
    
    # Modify p-value adjustments to maintain precision
    results$p_adj_fdr <- p.adjust(results$p_value, method = "fdr")
    results$p_adj_bonf <- p.adjust(results$p_value, method = "bonferroni")
    
    # Sort by p-value maintaining precision
    results <- results[order(results$p_value), ]
    
    return(results)
}

# Update Volcano plot function to use effect size instead of odds ratio
Volcano <- function(results, cluster_num, thr = 0.05, colors = c("darkviolet", "navy", "darkgrey")) {
    par(mar = c(4.5, 4.5, 1, 1))
    
    # Handle p-values of 0 by replacing with smallest possible double
    safe_p <- ifelse(results$p_value == 0, .Machine$double.xmin, results$p_value)
    
    # Calculate -log10 of p-values and use effect size directly
    log_pval <- -log10(safe_p)
    effect_sizes <- results$effect_size
    
    # Create the plot with specified colors
    plot(effect_sizes, log_pval,
         pch = 19,
         las = 1,
         cex = 0.5,
         xlab = "Effect Size (Mean Difference)",
         ylab = expression(-log[10](p[value])),
         col = ifelse(results$p_adj_fdr < thr, 
                     yes = colors[1],
                     no = colors[3]))
    
    # Add threshold line for significance after Bonferroni correction
    abline(h = -log10(thr / nrow(results)), lty = 2, col = "gray")
    
    # Add vertical line at x = 0
    abline(v = 0, lty = 2, col = "gray")
    
    # Label ALL significant proteins with color based on direction
    if(sum(results$p_adj_fdr < thr) > 0) {
        # Get all significant hits
        sig_hits <- results[results$p_adj_fdr < thr, ]
        
        # Split significant proteins into higher and lower
        sig_higher <- sig_hits$effect_size > 0
        sig_lower <- sig_hits$effect_size < 0
        
        # Add labels for higher proteins (red)
        if(any(sig_higher)) {
            text(sig_hits$effect_size[sig_higher], 
                 -log10(ifelse(sig_hits$p_value[sig_higher] == 0, 
                              .Machine$double.xmin, 
                              sig_hits$p_value[sig_higher])), 
                 labels = sig_hits$protein[sig_higher],
                 cex = 0.5,
                 pos = 3,
                 col = "red")
        }
        
        # Add labels for lower proteins (blue)
        if(any(sig_lower)) {
            text(sig_hits$effect_size[sig_lower], 
                 -log10(ifelse(sig_hits$p_value[sig_lower] == 0, 
                              .Machine$double.xmin, 
                              sig_hits$p_value[sig_lower])), 
                 labels = sig_hits$protein[sig_lower],
                 cex = 0.5,
                 pos = 3,
                 col = "blue")
        }
    }
    
    # Updated legend
    legend("topleft", 
           legend = c(paste0("Bonferroni threshold at ", thr), 
                     paste0("FDR significant hits (endotype ", cluster_num, ")"), 
                     "Not significant",
                     paste0("Higher in endotype ", cluster_num),
                     paste0("Lower in endotype ", cluster_num)),
           pch = c(NA, 19, 19, 19, 19),
           col = c("gray", colors[1], colors[3], "red", "blue"),
           lty = c(2, NA, NA, NA, NA))
}

# Define color schemes for each cluster
cluster_colors <- list(
    cluster1 = c("darkviolet", "navy", "darkgrey"),
    cluster2 = c("firebrick", "darkred", "darkgrey"),
    cluster3 = c("forestgreen", "darkgreen", "darkgrey"),
    cluster4 = c("dodgerblue", "navy", "darkgrey")
)

# Function to create summary statistics
create_summary_stats <- function(results) {
    summary_stats <- list(
        total_proteins = nrow(results),
        sig_p05 = sum(results$p_value < 0.05),
        sig_fdr = sum(results$p_adj_fdr < 0.05),
        sig_bonf = sum(results$p_adj_bonf < 0.05),
        top_10_proteins = head(results, 10)
    )
    return(summary_stats)
}

# Function to save analysis results and plot
save_cluster_analysis <- function(results, cluster_num) {
    # Save results dataframe
    saveRDS(results, 
            file.path(output_dir, sprintf("univariate_results_cluster%d.rds", cluster_num)))
    
    # Create and save summary statistics
    summary_stats <- create_summary_stats(results)
    saveRDS(summary_stats, 
            file.path(output_dir, sprintf("univariate_summary_stats_cluster%d.rds", cluster_num)))
    
    # Create volcano plot with cluster-specific colors and pass cluster number
    png(file.path(output_dir, sprintf("volcano_cluster%d.png", cluster_num)),
        width = 1000, height = 800, res = 100)
    Volcano(results, cluster_num, colors = cluster_colors[[paste0("cluster", cluster_num)]])
    dev.off()
    
    # Print analysis summary with scientific notation
    cat(sprintf("\n=== Cluster %d Analysis ===\n", cluster_num))
    cat("Number of proteins significant at p < 0.05:", summary_stats$sig_p05)
    cat("\nNumber of proteins significant after FDR correction:", summary_stats$sig_fdr)
    cat("\nNumber of proteins significant after Bonferroni correction:", summary_stats$sig_bonf)
    
    cat("\n\nTop 10 most significant proteins:\n")
    # Format p-values in scientific notation with higher precision
    top_10 <- summary_stats$top_10_proteins
    top_10$p_value <- format(top_10$p_value, scientific = TRUE, digits = 16)
    top_10$p_adj_fdr <- format(top_10$p_adj_fdr, scientific = TRUE, digits = 16)
    top_10$p_adj_bonf <- format(top_10$p_adj_bonf, scientific = TRUE, digits = 16)
    print(top_10)
}

### Perform analysis
# Analyze each cluster
results_cluster1 <- analyze_cluster(1)
results_cluster2 <- analyze_cluster(2)
results_cluster3 <- analyze_cluster(3)
results_cluster4 <- analyze_cluster(4)

# Save and plot results for each cluster
save_cluster_analysis(results_cluster1, 1)
save_cluster_analysis(results_cluster2, 2)
save_cluster_analysis(results_cluster3, 3)
save_cluster_analysis(results_cluster4, 4)

# Save all results in a single list object
all_cluster_results <- list(
    cluster1 = list(results = results_cluster1, summary = create_summary_stats(results_cluster1)),
    cluster2 = list(results = results_cluster2, summary = create_summary_stats(results_cluster2)),
    cluster3 = list(results = results_cluster3, summary = create_summary_stats(results_cluster3)),
    cluster4 = list(results = results_cluster4, summary = create_summary_stats(results_cluster4))
)

# Save combined results
saveRDS(all_cluster_results, 
        "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/univariate_all_cluster_results.rds")