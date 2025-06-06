#FINAL - with cleaned -3 responses and n_cat = 2 and with PFER_thr as 1.

library(fake)
library(sharp)
library(glmnet)
library(tidyverse)
library(pheatmap)

cluster_df <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/CGGM_Recoded_Cluster_Labels/cluster_assignments_k4.csv")
matched_df <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/Gianna/matched.csv")
met_score_df <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/Gianna/physical_activity.csv")
diet_score_df <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/Gianna/diet.csv")

merged_df <- cluster_df %>%
  inner_join(matched_df, by='index') %>%
  left_join(met_score_df, by='index') %>%
  left_join(diet_score_df, by='index')

#directory for results
dir.create("lasso_results_renamed_reversed_cleaned_ncat2", showWarnings = FALSE)

#combining England/Wales/Scotland scores as to avoid NA's
preprocess_regional_scores <- function(df) {
  #combine education scores
  df$education_score <- NA
  df$education_score[!is.na(df$education_score_england.0.0)] <- df$education_score_england.0.0[!is.na(df$education_score_england.0.0)]
  df$education_score[!is.na(df$education_score_wales.0.0)] <- df$education_score_wales.0.0[!is.na(df$education_score_wales.0.0)]
  df$education_score[!is.na(df$education_score_scotland.0.0)] <- df$education_score_scotland.0.0[!is.na(df$education_score_scotland.0.0)]
  
  #combine deprivation indices
  df$deprivation_index <- NA
  df$deprivation_index[!is.na(df$index_of_multiple_deprivation_england.0.0)] <- df$index_of_multiple_deprivation_england.0.0[!is.na(df$index_of_multiple_deprivation_england.0.0)]
  df$deprivation_index[!is.na(df$index_of_multiple_deprivation_wales.0.0)] <- df$index_of_multiple_deprivation_wales.0.0[!is.na(df$index_of_multiple_deprivation_wales.0.0)]
  df$deprivation_index[!is.na(df$index_of_multiple_deprivation_scotland.0.0)] <- df$index_of_multiple_deprivation_scotland.0.0[!is.na(df$index_of_multiple_deprivation_scotland.0.0)]
  
  return(df)
}

merged_df <- preprocess_regional_scores(merged_df)

# Remove all "prefer not to answer" (-3) responses
merged_df$alcohol_intake_frequency.0.0[merged_df$alcohol_intake_frequency.0.0 == -3] <- NA
merged_df$smoking_status.0.0[merged_df$smoking_status.0.0 == -3] <- NA
merged_df$sleep_duration.0.0[merged_df$sleep_duration.0.0 == -3] <- NA
merged_df$sleeplessness_insomnia.0.0[merged_df$sleeplessness_insomnia.0.0 == -3] <- NA
merged_df$time_spent_watching_television_tv.0.0[merged_df$time_spent_watching_television_tv.0.0 == -3] <- NA

# Reverse alcohol coding for interpretability (so 6=Never becomes 1, 1=Daily becomes 6)
# This makes higher numbers mean more frequent alcohol consumption
merged_df$alcohol_intake_frequency.0.0 <- 7 - merged_df$alcohol_intake_frequency.0.0

# No need to modify smoking_status.0.0 as it's already in increasing order
# (0=Never, 1=Previous, 2=Current)

#updated exposure variables 
exposure_vars <- c(
  'met_score',
  'diet_score',
  'smoking_status.0.0',           # Kept as 0,1,2
  'alcohol_intake_frequency.0.0',  # Now 6=highest frequency
  'education_score',
  'deprivation_index',
  'sleep_duration.0.0',
  'sleeplessness_insomnia.0.0',
  'time_spent_watching_television_tv.0.0'
)

#remove region-specific scores and add combined scores
exposure_vars <- exposure_vars[!grepl("education_score_|index_of_multiple_deprivation_", exposure_vars)]
exposure_vars <- c(exposure_vars, 'education_score', 'deprivation_index')

#named vector for renaming #new name = old name
new_names <- c(
  'Met Score' = 'met_score',  
  'Diet Score' = 'diet_score',
  'Smoking Status' = 'smoking_status.0.0',
  'Alcohol Intake Frequency' = 'alcohol_intake_frequency.0.0',
  'Education Score' = 'education_score',
  'Deprivation Index' = 'deprivation_index',
  'Sleep Duration' = 'sleep_duration.0.0',
  'Insomnia' = 'sleeplessness_insomnia.0.0',
  'Sedentary Behaviour' = 'time_spent_watching_television_tv.0.0'
)

# Select and rename columns in X
X <- merged_df %>% 
  select(all_of(exposure_vars)) %>%
  rename(!!!new_names)

#check for NA
print(colSums(is.na(X)))

#impute NAs with median
for(col in names(X)) {
  X[[col]][is.na(X[[col]])] <- median(X[[col]], na.rm = TRUE)
}

#check no missing values
print(colSums(is.na(X)))

#run stability selection for each cluster (y is clusters/outcome)
run_stability_selection <- function(X, cluster_labels, target_cluster) {
  #binary outcome (1 for target cluster, 0 for others)
  y <- as.numeric(cluster_labels == target_cluster)
  
  #matrix and scale
  X <- as.matrix(X)
  X <- scale(X)
  
  #printing so I can see
  cat(sprintf("\nRunning stability selection for Cluster %d vs Others\n", target_cluster))
  cat(sprintf("Number of samples in cluster: %d\n", sum(y)))
  
  t0 <- Sys.time()
  out <- VariableSelection(
    xdata = X,
    ydata = y,
    family = "binomial",  
    PFER_thr = 1,         #number of false positives allowed
    n_cat = 2,            # Added this parameter for binary outcomes
    verbose = TRUE        
  )
  t1 <- Sys.time()
  cat(sprintf("Time taken: %s\n", format(t1 - t0)))
  
  #selection proportions
  selprop <- SelectionProportions(out) #freq of each var
  hat_params <- Argmax(out) #for threshold
  
  #visualization
  pdf(sprintf("lasso_results_renamed_reversed_cleaned_ncat2/stability_selection_plot_cluster_%d.pdf", target_cluster))
  par(mar = c(15, 5, 1, 1))
  plot(selprop,
       type = "h", lwd = 3, las = 1, 
       xlab = "", ylab = "Selection Proportion", 
       xaxt = "n",  # Don't plot x axis
       col = ifelse(selprop >= hat_params[2], "red", "grey"),
       cex.lab = 1.5,
       main = sprintf("Stability Selection Results for Cluster %d vs Others", target_cluster)
  )
  
  #threshold line
  abline(h = hat_params[2], lty = 2, col = "darkred")
  
  #x labels
  for (i in seq_along(selprop)) {
    axis(1, at = i, labels = names(selprop)[i], las = 2,
         col.axis = ifelse(selprop[i] >= hat_params[2], "red", "grey"))
  }
  
  dev.off()
  
  #calibration plot
  pdf(sprintf("lasso_results_renamed_reversed_cleaned_ncat2/calibration_plot_cluster_%d.pdf", target_cluster))
  CalibrationPlot(out)
  dev.off()
  
  #save results to CSV
  results <- data.frame(
    Feature = names(selprop),
    Selection_Frequency = selprop
  ) %>% arrange(desc(Selection_Frequency))
  
  write.csv(results, 
            file = sprintf("lasso_results_renamed_reversed_cleaned_ncat2/stability_selection_results_cluster_%d.csv", target_cluster),
            row.names = FALSE)
  
  #print selected features
  selected_features <- names(selprop)[selprop >= hat_params[2]]
  cat("\nSelected features:\n")
  print(selected_features)
  cat(sprintf("\nSelection threshold: %.3f\n", hat_params[2]))
  
  return(list(
    results = results,
    selected_features = selected_features,
    threshold = hat_params[2],
    selection_proportions = selprop
  ))
}

#stability selection for each cluster
cat("Running stability selection (one cluster vs all others)...\n")
all_results <- list()

for(cluster in 1:4) {
  cat(sprintf("\n===== ANALYZING CLUSTER %d vs ALL OTHERS =====\n", cluster))
  all_results[[cluster]] <- run_stability_selection(X, merged_df$Cluster, cluster)
}

#summary table
results_table <- data.frame(
  Model = paste("Cluster", 1:4, "vs Rest"),
  'Penalty parameter (λ)' = sapply(all_results, function(x) x$threshold),
  'Selection proportion threshold (π)' = sapply(all_results, function(x) x$threshold),
  'Number of stably selected features' = sapply(all_results, function(x) length(x$selected_features)),
  'Most stably selected features' = sapply(all_results, function(x) {
    if(length(x$selected_features) > 0) {
      paste(paste0(x$selected_features, " (", round(x$results$Selection_Frequency[match(x$selected_features, x$results$Feature)], 2), ")"),
            collapse = "\n")
    } else {
      "None"
    }
  })
)

# save the summary table
write.csv(results_table, "lasso_results_renamed_reversed_cleaned_ncat2/summary_table.csv", row.names = FALSE)

#combined selection proportion plots
pdf("lasso_results_renamed_reversed_cleaned_ncat2/combined_selection_proportions.pdf", width = 12, height = 8)
par(mfrow = c(2,2), mar = c(10, 4, 4, 2))
for(cluster in 1:4) {
  selprop <- all_results[[cluster]]$selection_proportions
  threshold <- all_results[[cluster]]$threshold
  
  if(length(selprop) > 0) {
    plot(selprop,
         type = "h", lwd = 3, las = 2,
         xlab = "", ylab = "Selection Proportion",
         xaxt = "n",
         col = ifelse(selprop >= threshold, "red", "grey"),
         main = sprintf("Cluster %d vs Others", cluster))
    
    abline(h = threshold, lty = 2, col = "darkred")
    
    for (i in 1:length(selprop)) {
      axis(1, at = i, labels = names(selprop)[i], las = 2,
           col.axis = ifelse(selprop[i] >= threshold, "red", "grey"))
    }
  } else {
    plot(0, 0, type = "n", xlab = "", ylab = "",
         main = sprintf("Cluster %d vs Others - No features selected", cluster))
  }
}
dev.off()

#final summary labelled
cat("\nFinal Summary:\n")
for(cluster in 1:4) {
  cat(sprintf("\nCluster %d:\n", cluster))
  cat(sprintf("Number of samples: %d\n", sum(merged_df$Cluster == cluster)))
  cat(sprintf("Number of selected features: %d\n", 
              length(all_results[[cluster]]$selected_features)))
  
  if(length(all_results[[cluster]]$selected_features) > 0) {
    cat("Selected features (with selection frequencies):\n")
    selected_idx <- match(all_results[[cluster]]$selected_features,
                          all_results[[cluster]]$results$Feature)
    print(data.frame(
      Feature = all_results[[cluster]]$selected_features,
      Selection_Frequency = all_results[[cluster]]$results$Selection_Frequency[selected_idx]
    ))
  } else {
    cat("None.\n")
  }
}
