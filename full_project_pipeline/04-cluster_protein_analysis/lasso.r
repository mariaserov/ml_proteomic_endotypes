# Read in the preprocessed data
prot_df <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")

# Load required libraries
library(glmnet)

# Create output directory if it doesn't exist
output_dir <- "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/lasso/lasso_output"
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# Function to perform LASSO regression for a specific cluster
lasso_cluster <- function(cluster_num) {
    cat(sprintf("\n==== LASSO Regression for Cluster %d ====\n", cluster_num))
    
    # Get binary outcome for this cluster
    binary_col <- paste0("cluster", cluster_num, "_binary")
    y <- prot_df[[binary_col]]
    
    # Get protein columns (excluding index, Cluster, and binary columns)
    exclude_cols <- c("index", "Cluster", "cluster1_binary", "cluster2_binary", 
                      "cluster3_binary", "cluster4_binary")
    protein_cols <- names(prot_df)[!names(prot_df) %in% exclude_cols]
    
    # Create matrix of predictors
    X <- as.matrix(prot_df[, protein_cols])
    
    # Perform LASSO regression with cross-validation
    set.seed(123) # For reproducibility
    cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial", nfolds = 10)
    
    # Extract optimal lambda values
    lambda_min <- cv_lasso$lambda.min
    lambda_1se <- cv_lasso$lambda.1se
    
    cat(sprintf("Optimal lambda (min): %f\n", lambda_min))
    cat(sprintf("Optimal lambda (1se): %f\n", lambda_1se))
    
    # Fit LASSO model with optimal lambda (min)
    lasso_model <- glmnet(X, y, alpha = 1, family = "binomial", lambda = lambda_min)
    
    # Extract coefficients
    coefs <- as.matrix(coef(lasso_model))
    coefs_df <- data.frame(
        feature = rownames(coefs),
        coefficient = coefs[,1]
    )
    
    # Sort by absolute coefficient value
    coefs_df <- coefs_df[order(abs(coefs_df$coefficient), decreasing = TRUE),]
    
    # Create coefficient path plot
    png(file.path(output_dir, sprintf("cluster%d_coefficient_path.png", cluster_num)),
        width = 1000, height = 800, res = 100)
    plot(glmnet(X, y, alpha = 1, family = "binomial"), xvar = "lambda", 
         main = sprintf("Coefficient Path for Cluster %d", cluster_num),
         sub = "LASSO Regression")
    dev.off()
    
    # Create cross-validation curve
    png(file.path(output_dir, sprintf("cluster%d_cv_curve.png", cluster_num)),
        width = 1000, height = 800, res = 100)
    plot(cv_lasso, 
         main = sprintf("Cross-Validation Curve for Cluster %d", cluster_num))
    abline(v = log(lambda_min), col = "red", lty = 2)
    abline(v = log(lambda_1se), col = "blue", lty = 2)
    legend("topright", 
           legend = c("lambda.min", "lambda.1se"), 
           col = c("red", "blue"), 
           lty = 2)
    dev.off()
    
    # Print non-zero coefficients (LASSO's feature selection)
    cat("\nNon-zero coefficients (selected features):\n")
    nonzero_coefs <- coefs_df[coefs_df$coefficient != 0,]
    print(nonzero_coefs)
    
    # Calculate model performance metrics
    probs <- predict(lasso_model, X, type = "response")
    preds <- ifelse(probs > 0.5, 1, 0)
    accuracy <- mean(preds == y)
    
    # Create summary statistics
    model_summary <- list(
        lambda_min = lambda_min,
        lambda_1se = lambda_1se,
        accuracy = accuracy,
        nonzero_features = nonzero_coefs,
        full_coefficients = coefs_df,
        cv_object = cv_lasso,
        model_object = lasso_model
    )
    
    # Save model summary
    saveRDS(model_summary, 
            file.path(output_dir, sprintf("lasso_results_cluster%d.rds", cluster_num)))
    
    # Save full coefficient table as CSV
    write.csv(coefs_df, 
              file.path(output_dir, sprintf("lasso_coefficients_cluster%d.csv", cluster_num)),
              row.names = FALSE)
    
    return(model_summary)
}

# Perform LASSO regression for each cluster
lasso_results_1 <- lasso_cluster(1)
lasso_results_2 <- lasso_cluster(2)
lasso_results_3 <- lasso_cluster(3)
lasso_results_4 <- lasso_cluster(4)

# Combine all results
all_results <- list(
    cluster1 = lasso_results_1,
    cluster2 = lasso_results_2,
    cluster3 = lasso_results_3,
    cluster4 = lasso_results_4
)

# Save combined results
saveRDS(all_results, file.path(output_dir, "all_lasso_results.rds"))

# Create comparison plot of selected features across clusters
# Extract non-zero features from each cluster
get_selected_features <- function(results) {
    features <- results$nonzero_features$feature
    features[features != "(Intercept)"]
}

selected_features_1 <- get_selected_features(lasso_results_1)
selected_features_2 <- get_selected_features(lasso_results_2)
selected_features_3 <- get_selected_features(lasso_results_3)
selected_features_4 <- get_selected_features(lasso_results_4)

# Find unique features across all clusters
all_selected_features <- unique(c(selected_features_1, selected_features_2, 
                                selected_features_3, selected_features_4))

# Create a summary table
feature_comparison <- data.frame(
    feature = all_selected_features,
    stringsAsFactors = FALSE
)

# Extract coefficients for each cluster
for (i in 1:4) {
    results <- get(paste0("lasso_results_", i))
    coefs <- results$full_coefficients
    
    # Add coefficients to comparison table
    feature_comparison[, paste0("cluster", i)] <- NA
    for (f in all_selected_features) {
        idx <- which(coefs$feature == f)
        if (length(idx) > 0) {
            feature_comparison[feature_comparison$feature == f, paste0("cluster", i)] <- 
                coefs$coefficient[idx]
        }
    }
}

# Save feature comparison
write.csv(feature_comparison, 
          file.path(output_dir, "selected_features_comparison.csv"),
          row.names = FALSE)

cat("\nLASSO regression analysis complete. Results saved to:", output_dir, "\n")
