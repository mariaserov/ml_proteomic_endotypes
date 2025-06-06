########### Preprocessing
# Read in cluster assignments and main dataset
cluster_assignments <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/CGGM_Recoded_Cluster_Labels/cluster_assignments_k4.csv")

main_dataset <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/Data/main_hypertension_dataset.csv")

# Check how many cluster assignment rows are in main dataset
cluster_in_main <- sum(cluster_assignments$index %in% main_dataset$Row.names)
cat("Number of cluster assignments found in main dataset:", cluster_in_main, "\n")
cat("Total number of cluster assignments:", nrow(cluster_assignments), "\n")

# Check how many main dataset rows are in cluster assignments
main_in_cluster <- sum(main_dataset$Row.names %in% cluster_assignments$index)
cat("Number of main dataset rows found in cluster assignments:", main_in_cluster, "\n")
cat("Total number of main dataset rows:", nrow(main_dataset), "\n")

# Check for any mismatches
missing_from_main <- cluster_assignments$index[!cluster_assignments$index %in% main_dataset$Row.names]
missing_from_clusters <- main_dataset$Row.names[!main_dataset$Row.names %in% cluster_assignments$index]

if(length(missing_from_main) > 0) {
    cat("\nIndices in cluster assignments but not in main dataset:\n")
    print(missing_from_main)
}

if(length(missing_from_clusters) > 0) {
    cat("\nIndices in main dataset but not in cluster assignments:\n")
    print(missing_from_clusters)
}

# Merge the datasets
# First, ensure Row.names is treated as character to match with index
main_dataset$Row.names <- as.character(main_dataset$Row.names)
cluster_assignments$index <- as.character(cluster_assignments$index)

# Merge datasets using index and Row.names as keys
merged_dataset <- merge(main_dataset, 
                       cluster_assignments, 
                       by.x = "Row.names", 
                       by.y = "index", 
                       all.x = FALSE, 
                       all.y = FALSE)

# Print dimensions of merged dataset to verify
cat("\nDimensions of merged dataset:", dim(merged_dataset)[1], "rows,", dim(merged_dataset)[2], "columns\n")

# Drop Row.names column and rename X to index
merged_dataset$Row.names <- NULL  # Remove Row.names column
names(merged_dataset)[names(merged_dataset) == "X"] <- "index"  # Rename X to index

# Print first few column names to verify changes
cat("\nFirst few column names after modifications:\n")
print(head(names(merged_dataset)))

# Find the position of 'AARSD1' column
aarsd1_pos <- which(names(merged_dataset) == "AARSD1")

# Create protein dataset including index and all columns from AARSD1 onwards
prot_df <- merged_dataset[, c("index", names(merged_dataset)[aarsd1_pos:ncol(merged_dataset)])]

# Print dimensions to verify
cat("\nDimensions of protein dataset:", dim(prot_df)[1], "rows,", dim(prot_df)[2], "columns\n")
# Print first few column names to verify
cat("\nFirst few protein columns:\n")
print(head(names(prot_df)))

### Add binary columns for all clusters
prot_df$cluster1_binary <- as.numeric(prot_df$Cluster == 1)
prot_df$cluster2_binary <- as.numeric(prot_df$Cluster == 2)
prot_df$cluster3_binary <- as.numeric(prot_df$Cluster == 3)
prot_df$cluster4_binary <- as.numeric(prot_df$Cluster == 4)

# Save the protein dataframe with binary columns for use in univariate analysis
saveRDS(prot_df, "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")

# Print completion message
cat("\nPreprocessing complete. Data saved for univariate analysis.\n")
cat("To run the analysis, use: Rscript univariate_analysis.r\n")


