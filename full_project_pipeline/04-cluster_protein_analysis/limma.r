# Read in the preprocessed data
prot_df <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")

# Create output directory if it doesn't exist
output_dir <- "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/limma/limma_output"
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# Load required libraries with error handling
if (!require("limma")) {
    stop("The limma package is required but not installed")
}

### Function to perform cluster-specific analysis using limma
analyze_cluster_limma <- function(cluster_num) {
    # Get the appropriate binary outcome
    binary_col <- paste0("cluster", cluster_num, "_binary")
    
    # Verify binary column exists
    if (!binary_col %in% names(prot_df)) {
        stop(sprintf("Binary column %s not found in data", binary_col))
    }
    
    # Verify binary column coding
    binary_values <- unique(prot_df[[binary_col]])
    if (!all(binary_values %in% c(0, 1))) {
        stop(sprintf("Binary column %s contains values other than 0/1: %s", 
                    binary_col, paste(binary_values, collapse=", ")))
    }
    
    # Get protein column names (exclude index, Cluster, and binary columns)
    protein_cols <- names(prot_df)[!names(prot_df) %in% c("index", "Cluster", 
                                                         "cluster1_binary", "cluster2_binary", 
                                                         "cluster3_binary", "cluster4_binary")]
    
    # Verify we have protein columns
    if (length(protein_cols) == 0) {
        stop("No protein columns found in data")
    }
    
    # Create expression matrix and transpose it for limma (proteins should be rows)
    expr_matrix <- t(as.matrix(prot_df[, protein_cols]))
    
    # Create design matrix more explicitly
    design_df <- data.frame(binary = prot_df[[binary_col]])
    design <- model.matrix(~ binary, data = design_df)
    colnames(design) <- c("Intercept", "ClusterEffect")
    
    # Verify sample alignment
    if (ncol(expr_matrix) != nrow(design)) {
        stop("Number of samples in expression matrix does not match design matrix")
    }
    
    # Fit linear model
    fit <- lmFit(expr_matrix, design)
    
    # Apply empirical Bayes moderation
    fit <- eBayes(fit)
    
    # Get results table
    results_table <- topTable(fit, coef="ClusterEffect", number=Inf)
    
    # Calculate proper standard errors using moderated t-statistics
    se <- abs(results_table$logFC) / abs(results_table$t)
    
    # Format results with corrected confidence intervals
    results <- data.frame(
        protein = rownames(results_table),
        effect_size = results_table$logFC,
        t_statistic = results_table$t,
        std_error = se,
        ci_lower = results_table$logFC - 1.96 * se,
        ci_upper = results_table$logFC + 1.96 * se,
        p_value = results_table$P.Value,
        p_adj_fdr = results_table$adj.P.Val,
        stringsAsFactors = FALSE
    )
    
    # Verify results integrity
    if (nrow(results) == 0) {
        stop("No results generated from the analysis")
    }
    if (any(is.na(results$effect_size))) {
        warning("Some effect sizes are NA")
    }
    if (any(is.na(results$p_value))) {
        warning("Some p-values are NA")
    }
    
    # Add Bonferroni correction
    results$p_adj_bonf <- p.adjust(results$p_value, method = "bonferroni")
    
    # Report analysis completion
    cat(sprintf("\nCluster %d Analysis Summary:\n", cluster_num))
    cat(sprintf("Total proteins analyzed: %d\n", nrow(results)))
    cat(sprintf("Proteins significant at FDR < 0.05: %d\n", sum(results$p_adj_fdr < 0.05)))
    cat(sprintf("Proteins significant at Bonferroni < 0.05: %d\n", sum(results$p_adj_bonf < 0.05)))
    
    return(results)
}

# Volcano plot function (reusing from univariate analysis with minor modifications)
Volcano <- function(results, cluster_num, thr = 0.05, colors = c("darkviolet", "navy", "darkgrey")) {
    par(mar = c(4.5, 4.5, 1, 1))
    
    # Handle p-values of 0 by replacing with smallest possible double
    safe_p <- ifelse(results$p_value == 0, .Machine$double.xmin, results$p_value)
    
    # Calculate -log10 of p-values and use effect size directly
    log_pval <- -log10(safe_p)
    effect_sizes <- results$effect_size
    
    # Filter points to plot (exclude those above 290)
    plot_threshold <- 290
    plot_idx <- log_pval <= plot_threshold
    
    # Create the plot with specified colors and ylim
    plot(effect_sizes[plot_idx], log_pval[plot_idx],
         pch = 19,
         las = 1,
         cex = 0.5,
         xlab = "Log2 Fold Change",
         ylab = expression(-log[10](p[value])),
         ylim = c(0, 300),
         col = ifelse(results$p_adj_fdr[plot_idx] < thr, 
                     yes = colors[1],
                     no = colors[3]))
    
    # Add threshold line for significance after Bonferroni correction
    abline(h = -log10(thr / nrow(results)), lty = 2, col = "gray")
    
    # Add vertical line at x = 0
    abline(v = 0, lty = 2, col = "gray")
    
    # Label significant proteins (only those below the plot threshold)
    if(sum(results$p_adj_fdr < thr) > 0) {
        sig_hits <- results[results$p_adj_fdr < thr, ]
        sig_log_pval <- -log10(ifelse(sig_hits$p_value == 0, 
                                     .Machine$double.xmin, 
                                     sig_hits$p_value))
        sig_below_threshold <- sig_log_pval <= plot_threshold
        
        sig_higher <- sig_hits$effect_size > 0 & sig_below_threshold
        sig_lower <- sig_hits$effect_size < 0 & sig_below_threshold
        
        if(any(sig_higher)) {
            text(sig_hits$effect_size[sig_higher], 
                 sig_log_pval[sig_higher], 
                 labels = sig_hits$protein[sig_higher],
                 cex = 0.5,
                 pos = 3,
                 col = "red")
        }
        
        if(any(sig_lower)) {
            text(sig_hits$effect_size[sig_lower], 
                 sig_log_pval[sig_lower], 
                 labels = sig_hits$protein[sig_lower],
                 cex = 0.5,
                 pos = 3,
                 col = "blue")
        }
    }
    
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
        top_10_proteins = head(results[order(results$p_value),], 10)
    )
    return(summary_stats)
}

# Function to save analysis results and plot
save_cluster_analysis <- function(results, cluster_num) {
    # Save results dataframe
    saveRDS(results, 
            file.path(output_dir, sprintf("results_cluster%d_limma.rds", cluster_num)))
    
    # Create and save summary statistics
    summary_stats <- create_summary_stats(results)
    saveRDS(summary_stats, 
            file.path(output_dir, sprintf("summary_stats_cluster%d_limma.rds", cluster_num)))
    
    # Create volcano plot
    png(file.path(output_dir, sprintf("volcano_cluster%d_limma.png", cluster_num)),
        width = 1000, height = 800, res = 100)
    Volcano(results, cluster_num, colors = cluster_colors[[paste0("cluster", cluster_num)]])
    dev.off()
    
    # Print analysis summary
    cat(sprintf("\n=== Cluster %d Analysis ===\n", cluster_num))
    cat("Number of proteins significant at p < 0.05:", summary_stats$sig_p05)
    cat("\nNumber of proteins significant after FDR correction:", summary_stats$sig_fdr)
    cat("\nNumber of proteins significant after Bonferroni correction:", summary_stats$sig_bonf)
    
    cat("\n\nTop 10 most significant proteins:\n")
    top_10 <- summary_stats$top_10_proteins
    print(top_10)
}

### Perform analysis
# Analyze each cluster
results_cluster1 <- analyze_cluster_limma(1)
results_cluster2 <- analyze_cluster_limma(2)
results_cluster3 <- analyze_cluster_limma(3)
results_cluster4 <- analyze_cluster_limma(4)

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
        file.path(output_dir, "all_cluster_results_limma.rds"))
