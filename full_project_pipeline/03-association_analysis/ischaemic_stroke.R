# PACKAGES
# install.packages("tidyverse")

# 1) Set up
rm(list = ls()) 
setwd("/rds/general/project/hda_24-25/live/TDS/Group09")

library(tidyverse)
library(broom)

# 2) Prepare data 
outcome <- readRDS("./outcome_definition/Outputs_saved/ischaemic_stroke/output_final.rds") 
ukb_htn <- readRDS("./Data/main_hypertension_dataset.rds")
k4 <- read.csv("./CGGM_Recoded_Cluster_Labels/cluster_assignments_k4.csv")

# Select relevant columns (cholesterol removed)
selected_cols <- c("sex.0.0", "age_at_recruitment.0.0", 
                   "diabetes_diagnosed_by_doctor.0.0", 
                   "smoking_status.0.0") 

df <- ukb_htn[ , selected_cols] 

# Inspect missingness
colSums(is.na(df))

# Remove rows with both diabetes and smoking missing
df <- df[!(is.na(df$diabetes_diagnosed_by_doctor.0.0) & is.na(df$smoking_status.0.0)), ] 
colSums(is.na(df))

# Clean variables
df$smoking_status_clean <- df$smoking_status.0.0
df$smoking_status_clean[df$smoking_status_clean < 0] <- NA 
df$diabetes_clean <- df$diabetes_diagnosed_by_doctor.0.0
df$diabetes_clean[df$diabetes_clean < 0] <- NA 

df <- df %>% filter(!is.na(smoking_status_clean), !is.na(diabetes_clean)) 

# Convert to factors
df$sex.0.0 <- as.factor(df$sex.0.0)
df$smoking_status_clean <- factor(df$smoking_status_clean, levels = c(0, 1, 2)) 
df$diabetes_clean <- factor(df$diabetes_clean, levels = c(0, 1))

# Drop old columns
df <- df %>% select(-smoking_status.0.0, -diabetes_diagnosed_by_doctor.0.0)

# Add eid and merge datasets
df$eid <- rownames(df)
outcome$eid <- rownames(outcome)
k4 <- k4 %>% rename(eid = index)

# Ensure eid is numeric
df$eid <- as.numeric(df$eid)
outcome$eid <- as.numeric(outcome$eid)
k4$eid <- as.numeric(k4$eid)

# Merge datasets
df <- merge(df, outcome, by = "eid", all.x = TRUE)
df <- merge(df, k4, by = "eid", all.x = TRUE)
rownames(df) <- NULL

# Keep only individuals with incident outcome data
df <- df[!is.na(df$incident_case), ]

# Convert Cluster to factor
df$Cluster <- as.factor(df$Cluster)

# Outcome distribution
table(df$incident_case)

# Drop unused columns
df <- df %>% select(-date_diagnosis, -date_death, -time_to_diagnosis, -prevalent_case)

# 3) Logistic regression for each cluster (no imputation)

# Initialize lists to store pooled results
pooled_models <- list()
pooled_summaries <- list()

# Get unique clusters
clusters <- sort(unique(as.character(df$Cluster)))

# Loop through each cluster and run logistic regression
for (clust in clusters) {
  
  # Create binary cluster indicator
  var_name <- paste0("cluster_", clust)
  df[[var_name]] <- factor(ifelse(df$Cluster == clust, 1, 0),
                           levels = c(0, 1),
                           labels = c(paste("Not", clust), clust))
  
  # Run logistic regression (cholesterol removed)
  model <- glm(
    incident_case ~ get(var_name) + age_at_recruitment.0.0 + sex.0.0 +
      diabetes_clean + smoking_status_clean,
    data = df,
    family = binomial
  )
  
  # Store results
  pooled_models[[clust]] <- model
  pooled_summaries[[clust]] <- summary(model)
  
  # Print model summary
  cat("\n=============================\n")
  cat("Results for Cluster", clust, "\n")
  cat("=============================\n\n")
  print(summary(model))
  
  # Odds ratios and 95% CIs
  coef_table <- summary(model)$coefficients
  ORs <- exp(coef_table[, "Estimate"])
  CI_lower <- exp(coef_table[, "Estimate"] - 1.96 * coef_table[, "Std. Error"])
  CI_upper <- exp(coef_table[, "Estimate"] + 1.96 * coef_table[, "Std. Error"])
  
  OR_table <- cbind(OR = ORs, CI_lower = CI_lower, CI_upper = CI_upper)
  
  cat("\nOdds Ratios and 95% Confidence Intervals:\n")
  print(round(OR_table, 3))
}

# 4) Combine cluster results into summary table

results_list <- lapply(clusters, function(cl) {
  
  # Get model summary
  model_sum <- pooled_summaries[[cl]]
  
  # Extract coefficient for cluster indicator
  coef_row <- model_sum$coefficients[2, ]  # 2nd row
  
  estimate <- coef_row["Estimate"]
  std_error <- coef_row["Std. Error"]
  z_value <- coef_row["z value"]
  p_value <- coef_row["Pr(>|z|)"]
  
  OR <- exp(estimate)
  CI_lower <- exp(estimate - 1.96 * std_error)
  CI_upper <- exp(estimate + 1.96 * std_error)
  
  significant <- ifelse(p_value < 0.05, "Yes", "No")
  
  data.frame(
    Cluster     = paste("Cluster", cl),
    Estimate    = estimate,
    Std_Error   = std_error,
    Z_value     = z_value,
    P_value     = p_value,
    OR          = OR,
    CI_lower    = CI_lower,
    CI_upper    = CI_upper,
    Significant = significant
  )
})

results_df <- bind_rows(results_list)

# Print results table
print(results_df)

# 5) Export results
write.csv(results_df, "/rds/general/user/jdc124/projects/hda_24-25/live/TDS/Group09/full_pipeline/4_outcome_association//Results/ischaemic_stroke_df.csv", row.names = FALSE)
