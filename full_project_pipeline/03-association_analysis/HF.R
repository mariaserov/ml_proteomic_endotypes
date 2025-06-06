#PACKAGES
#install.packages("tidyverse")


#1) PREPARE DATA
# set up 
rm(list = ls()) 
setwd("/rds/general/project/hda_24-25/live/TDS/Group09")
library(tidyverse)

# import data - ACTION REQUIRED - must select folder name
outcome <- readRDS("./outcome_definition/Outputs_saved/HF/output_final.rds") 
ukb_htn <- readRDS("./Data/main_hypertension_dataset.rds")
k4 <- read.csv("./CGGM_Recoded_Cluster_Labels/cluster_assignments_k4.csv")

# format proteomic data 
selected_cols <- c("sex.0.0", "age_at_recruitment.0.0", "diabetes_diagnosed_by_doctor.0.0", "smoking_status.0.0") # define relevant columns 
df <- ukb_htn[ , selected_cols] # create new data frame with selected variables
colSums(is.na(df)) #check NAs
df <- df[!(is.na(df$diabetes_diagnosed_by_doctor.0.0) & is.na(df$smoking_status.0.0)), ] # drop 30 people with missing data 
colSums(is.na(df)) #now no missing data 

df$smoking_status_clean <- df$smoking_status.0.0
df$smoking_status_clean[df$smoking_status_clean < 0] <- NA
df$smoking_status_clean <- factor(df$smoking_status_clean, levels = c(0, 1, 2))
table(df$smoking_status_clean)

df$diabetes_clean <- df$diabetes_diagnosed_by_doctor.0.0
df$diabetes_clean[df$diabetes_clean < 0] <- NA
df$diabetes_clean <- factor(df$diabetes_clean, levels = c(0, 1))
table(df$diabetes_clean)

df$sex.0.0 <- as.factor(df$sex.0.0)

df <- df %>% select(-smoking_status.0.0, -diabetes_diagnosed_by_doctor.0.0) #drop old columns

# make sure eid is correctly formatted 
df$eid <- rownames(df)  # Make row names an explicit column
outcome$eid <- rownames(outcome)  # Do the same for AAA
k4 <- k4 %>% rename(eid = index)  # Rename 'index' to 'eid' in k4

# ensure eid is numeric 
df$eid <- as.numeric(df$eid)
outcome$eid <- as.numeric(outcome$eid)
k4$eid <- as.numeric(k4$eid)

# merge
df <- merge(df, outcome, by = "eid", all.x = TRUE)
df <- merge(df, k4, by = "eid", all.x = TRUE)
rownames(df) <- NULL
df <- df[!is.na(df$incident_case), ] #drop values for which we do not have outcome data

# ensure cluster is a factor
df$Cluster <- as.factor(df$Cluster)

# 2) ANALYSIS 
#see what outcome data looks like
table(df$incident_case)

# Initialize a list to store each cluster's row for the summary dataframe
results_list <- list()

# --- Loop over each cluster ---
clusters <- sort(unique(as.character(df$Cluster)))  # e.g. "A", "B", "C", "D"

# Load dplyr for data frame operations
library(dplyr)

# Initialize list for full model results and list for summary table rows
model_results <- list()
results_list <- list()

# --- Loop over each cluster ---
clusters <- sort(unique(as.character(df$Cluster)))  # e.g. "A", "B", "C", "D"

for (i in seq_along(clusters)) {
  clust <- clusters[i]
  
  # Create binary variable for cluster membership
  var_name <- paste0("cluster_", clust)
  df[[var_name]] <- factor(ifelse(as.character(df$Cluster) == clust, 1, 0),
                           levels = c(0, 1),
                           labels = c(paste("Not", clust), clust))
  
  # Logistic regression model
  model_binary <- glm(
    incident_case ~ get(var_name) + age_at_recruitment.0.0 + sex.0.0 +
      diabetes_clean + smoking_status_clean,
    data = df,
    family = binomial
  )
  
  # Save the full model in the list for later use
  model_results[[paste0("Cluster ", i)]] <- list(
    model = model_binary,
    summary = summary(model_binary)
  )
  
  # Extract coefficients for the cluster variable (2nd row)
  coef_table <- summary(model_binary)$coefficients
  coef_row <- coef_table[2, ]  # 2nd row corresponds to the cluster variable
  
  # Calculate OR and 95% CI
  estimate <- coef_row["Estimate"]
  std_error <- coef_row["Std. Error"]
  z_value <- coef_row["z value"]
  p_value <- coef_row["Pr(>|z|)"]
  OR <- exp(estimate)
  CI_lower <- exp(estimate - 1.96 * std_error)
  CI_upper <- exp(estimate + 1.96 * std_error)
  
  # Create a data frame row for this cluster
  result_row <- data.frame(
    Cluster     = paste0("Cluster ", i),
    Estimate    = estimate,
    Std_Error   = std_error,
    Z_value     = z_value,
    P_value     = p_value,
    OR          = OR,
    CI_lower    = CI_lower,
    CI_upper    = CI_upper
  )
  
  # Append to the summary list
  results_list[[i]] <- result_row
}

# Summary results

# Combine rows into a single data frame
results_df <- bind_rows(results_list)
results_df <- results_df %>% mutate(Significant = ifelse(P_value < 0.05, "Yes", "No")) # Add a column to flag significant clusters (p < 0.05)
print(results_df)

#Print full results
for (clust in names(model_results)) {
  
  cat("\n=============================\n")
  cat("Results for", clust, "\n")
  cat("=============================\n\n")
  
  # Print the model summary
  print(model_results[[clust]]$summary)
  
  # Extract coefficients
  coef_table <- model_results[[clust]]$summary$coefficients
  
  # Calculate Odds Ratios (ORs) and 95% CIs
  ORs <- exp(coef_table[, "Estimate"])
  CI_lower <- exp(coef_table[, "Estimate"] - 1.96 * coef_table[, "Std. Error"])
  CI_upper <- exp(coef_table[, "Estimate"] + 1.96 * coef_table[, "Std. Error"])
  
  # Combine into a table
  OR_table <- cbind(OR = ORs, CI_lower = CI_lower, CI_upper = CI_upper)
  
  cat("\nOdds Ratios and 95% Confidence Intervals:\n")
  print(round(OR_table, 3))
}

#export results
write.csv(results_df, "/rds/general/user/jdc124/projects/hda_24-25/live/TDS/Group09/full_pipeline/4_outcome_association//Results/HF_df.csv", row.names = FALSE)




