# Read in the preprocessed data
prot_df <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")

# Load required libraries
library(glmnet)

# Create output directory if it doesn't exist
output_dir <- "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/multiple_ridge_reg/ridge_reg_output"
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
}

# Function to perform ridge regression for a specific cluster
ridge_cluster <- function(cluster_num) {
    cat(sprintf("\n==== Ridge Regression for Cluster %d ====\n", cluster_num))
    
    # Get binary outcome for this cluster
    binary_col <- paste0("cluster", cluster_num, "_binary")
    y <- prot_df[[binary_col]]
    
    # Get protein columns (excluding index, Cluster, and binary columns)
    exclude_cols <- c("index", "Cluster", "cluster1_binary", "cluster2_binary", 
                      "cluster3_binary", "cluster4_binary")
    protein_cols <- names(prot_df)[!names(prot_df) %in% exclude_cols]
    
    # Create matrix of predictors
    X <- as.matrix(prot_df[, protein_cols])
    
    # Perform ridge regression with cross-validation
    set.seed(123) # For reproducibility
    cv_ridge <- cv.glmnet(X, y, alpha = 0, family = "binomial", nfolds = 10)
    
    # Extract optimal lambda values
    lambda_min <- cv_ridge$lambda.min
    lambda_1se <- cv_ridge$lambda.1se
    
    cat(sprintf("Optimal lambda (min): %f\n", lambda_min))
    cat(sprintf("Optimal lambda (1se): %f\n", lambda_1se))
    
    # Fit ridge model with optimal lambda (min)
    ridge_model <- glmnet(X, y, alpha = 0, family = "binomial", lambda = lambda_min)
    
    # Extract coefficients
    coefs <- as.matrix(coef(ridge_model))
    coefs_df <- data.frame(
        feature = rownames(coefs),
        coefficient = coefs[,1]
    )
    
    # Sort by absolute coefficient value
    coefs_df <- coefs_df[order(abs(coefs_df$coefficient), decreasing = TRUE),]
    
    # Create coefficient path plot
    png(file.path(output_dir, sprintf("cluster%d_coefficient_path.png", cluster_num)),
        width = 1000, height = 800, res = 100)
    plot(glmnet(X, y, alpha = 0, family = "binomial"), xvar = "lambda", 
         main = "",
         sub = "")
    dev.off()
    
    # Create cross-validation curve
    png(file.path(output_dir, sprintf("cluster%d_cv_curve.png", cluster_num)),
        width = 1000, height = 800, res = 100)
    plot(cv_ridge, 
         main = "")
    abline(v = log(lambda_min), col = "red", lty = 2)
    abline(v = log(lambda_1se), col = "blue", lty = 2)
    legend("topright", 
           legend = c("lambda.min", "lambda.1se"), 
           col = c("red", "blue"), 
           lty = 2)
    dev.off()
    
    # Print top 20 features by absolute coefficient value
    cat("\nTop 20 features by absolute coefficient value:\n")
    print(head(coefs_df, 20))
    
    # Calculate model performance metrics
    probs <- predict(ridge_model, X, type = "response")
    preds <- ifelse(probs > 0.5, 1, 0)
    accuracy <- mean(preds == y)
    
    # Create summary statistics
    model_summary <- list(
        lambda_min = lambda_min,
        lambda_1se = lambda_1se,
        accuracy = accuracy,
        top_features = head(coefs_df, 50),
        full_coefficients = coefs_df,
        cv_object = cv_ridge,
        model_object = ridge_model
    )
    
    # Save model summary
    saveRDS(model_summary, 
            file.path(output_dir, sprintf("ridge_results_cluster%d.rds", cluster_num)))
    
    # Save full coefficient table as CSV
    write.csv(coefs_df, 
              file.path(output_dir, sprintf("ridge_coefficients_cluster%d.csv", cluster_num)),
              row.names = FALSE)
    
    return(model_summary)
}

# Perform ridge regression for each cluster
ridge_results_1 <- ridge_cluster(1)
ridge_results_2 <- ridge_cluster(2)
ridge_results_3 <- ridge_cluster(3)
ridge_results_4 <- ridge_cluster(4)

# Combine all results
all_results <- list(
    cluster1 = ridge_results_1,
    cluster2 = ridge_results_2,
    cluster3 = ridge_results_3,
    cluster4 = ridge_results_4
)

# Save combined results
saveRDS(all_results, file.path(output_dir, "all_ridge_results.rds"))

# Create comparison plot of top features across clusters
# Extract top 10 features from each cluster
get_top_features <- function(results, n=10) {
    top <- head(results$top_features, n)
    top$feature[top$feature != "(Intercept)"]
}

top_features_1 <- get_top_features(ridge_results_1)
top_features_2 <- get_top_features(ridge_results_2)
top_features_3 <- get_top_features(ridge_results_3)
top_features_4 <- get_top_features(ridge_results_4)

# Find unique features across all clusters
all_top_features <- unique(c(top_features_1, top_features_2, top_features_3, top_features_4))

# Create a summary table
feature_comparison <- data.frame(
    feature = all_top_features,
    stringsAsFactors = FALSE
)

# Extract coefficients for each cluster
for (i in 1:4) {
    results <- get(paste0("ridge_results_", i))
    coefs <- results$full_coefficients
    
    # Add coefficients to comparison table
    feature_comparison[, paste0("cluster", i)] <- NA
    for (f in all_top_features) {
        idx <- which(coefs$feature == f)
        if (length(idx) > 0) {
            feature_comparison[feature_comparison$feature == f, paste0("cluster", i)] <- 
                coefs$coefficient[idx]
        }
    }
}

# Save feature comparison
write.csv(feature_comparison, 
          file.path(output_dir, "top_features_comparison.csv"),
          row.names = FALSE)

cat("\nRidge regression analysis complete. Results saved to:", output_dir, "\n")