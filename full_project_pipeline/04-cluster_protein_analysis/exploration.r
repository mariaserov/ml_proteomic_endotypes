# Load the results
results_cluster1 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/results_cluster1.rds")
results_cluster2 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/results_cluster2.rds")
results_cluster3 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/results_cluster3.rds")
results_cluster4 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/results_cluster4.rds")

# Load the summary stats
summary_stats_cluster1 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/summary_stats_cluster1.rds")
summary_stats_cluster2 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/summary_stats_cluster2.rds")
summary_stats_cluster3 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/summary_stats_cluster3.rds")
summary_stats_cluster4 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/univariate_log_reg/univariate_output/summary_stats_cluster4.rds")


# Load the RDS file
results <- readRDS('/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/multiple_ridge_reg/ridge_reg_output/all_ridge_results.rds')

# Extract top features as a table
top_features_table <- data.frame(
  feature = results$cluster1$top_features$feature,
  coefficient = results$cluster1$top_features$coefficient
)

# View the table
head(top_features_table)  # Shows first 6 rows
View(top_features_table)  # Opens in RStudio viewer