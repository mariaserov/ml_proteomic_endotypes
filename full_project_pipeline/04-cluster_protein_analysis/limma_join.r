full_dataset_1<-read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data/full_results_cluster1_with_uniprot.csv")
full_dataset_2<-read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data/full_results_cluster2_with_uniprot.csv")
full_dataset_3<-read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data/full_results_cluster3_with_uniprot.csv")
full_dataset_4<-read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data/full_results_cluster4_with_uniprot.csv")

limma_results_1 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/limma/limma_output/results_cluster1_limma.rds")
limma_results_2 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/limma/limma_output/results_cluster2_limma.rds")
limma_results_3 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/limma/limma_output/results_cluster3_limma.rds")
limma_results_4 <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/limma/limma_output/results_cluster4_limma.rds")

# Function to add prefix to limma columns and join datasets
process_and_join_cluster <- function(limma_results, full_dataset, cluster_num) {
    # Verify both datasets have the 'protein' column
    if(!"protein" %in% names(limma_results) || !"protein" %in% names(full_dataset)) {
        stop(sprintf("Missing protein column in datasets for cluster %d", cluster_num))
    }
    
    # Get columns to rename (all except 'protein')
    cols_to_rename <- setdiff(names(limma_results), "protein")
    
    # Create new names with 'limma_' prefix
    new_names <- paste0("limma_", cols_to_rename)
    
    # Rename columns in limma results
    names(limma_results)[match(cols_to_rename, names(limma_results))] <- new_names
    
    # Verify no duplicate column names before joining
    duplicate_cols <- intersect(names(limma_results), names(full_dataset))
    if(length(duplicate_cols) > 1) { # > 1 because 'protein' should be shared
        stop(sprintf("Duplicate columns found in cluster %d: %s", 
                    cluster_num, 
                    paste(setdiff(duplicate_cols, "protein"), collapse=", ")))
    }
    
    # Perform the join
    merged_data <- merge(full_dataset, limma_results, by="protein", all.x=TRUE)
    
    # Verify the merge
    if(nrow(merged_data) != nrow(full_dataset)) {
        warning(sprintf("Row count changed after merge for cluster %d. Before: %d, After: %d", 
                       cluster_num, nrow(full_dataset), nrow(merged_data)))
    }
    
    # Check for any NA values introduced by the merge
    na_counts <- colSums(is.na(merged_data[, new_names]))
    if(any(na_counts > 0)) {
        warning(sprintf("NAs introduced in limma columns for cluster %d: %s", 
                       cluster_num, 
                       paste(names(na_counts[na_counts > 0]), collapse=", ")))
    }
    
    return(merged_data)
}

# Process each cluster
output_dir <- "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data"

# Cluster 1
result1 <- process_and_join_cluster(limma_results_1, full_dataset_1, 1)
write.csv(result1, 
          file.path(output_dir, "full_results_limma_cluster_1.csv"), 
          row.names=FALSE)

# Cluster 2
result2 <- process_and_join_cluster(limma_results_2, full_dataset_2, 2)
write.csv(result2, 
          file.path(output_dir, "full_results_limma_cluster_2.csv"), 
          row.names=FALSE)

# Cluster 3
result3 <- process_and_join_cluster(limma_results_3, full_dataset_3, 3)
write.csv(result3, 
          file.path(output_dir, "full_results_limma_cluster_3.csv"), 
          row.names=FALSE)

# Cluster 4
result4 <- process_and_join_cluster(limma_results_4, full_dataset_4, 4)
write.csv(result4, 
          file.path(output_dir, "full_results_limma_cluster_4.csv"), 
          row.names=FALSE)

# Print summary of results
for(i in 1:4) {
    result <- get(paste0("result", i))
    cat(sprintf("\nCluster %d Summary:\n", i))
    cat("Total rows:", nrow(result), "\n")
    cat("Total columns:", ncol(result), "\n")
    cat("Limma columns added:", sum(grepl("^limma_", names(result))), "\n")
    cat("Sample of limma columns:", paste(head(grep("^limma_", names(result), value=TRUE)), collapse=", "), "\n\n")
}