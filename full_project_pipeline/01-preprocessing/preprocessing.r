rm(list=ls())
# Read in the datasets
hypertension_data <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/Data/main_hypertension_dataset.csv")
AAA_data <- readRDS('/rds/general/project/hda_24-25/live/TDS/Group09/outcome_definition/Outputs_saved/AAA/output_final.rds')
protein_data <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")
protein_labels <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data/cluster_3_dea.csv")

# Core variable selection
selected_cols_ukb <- c(
  "Row.names", "sex.0.0", "age_at_recruitment.0.0", "ethnic_background.0.0",
  "waist_circumference.0.0", "hip_circumference.0.0", "body_mass_index_bmi.1.0",
  "diabetes_diagnosed_by_doctor.0.0", "systolic_blood_pressure_automated_reading.0.0",
  "diastolic_blood_pressure_automated_reading.0.0", "age_high_blood_pressure_diagnosed.0.0",
  "medication_for_cholesterol_blood_pressure_or_diabetes.0.0",
  "current_tobacco_smoking.0.0", "past_tobacco_smoking.0.0",
  "alcohol_intake_frequency.0.0", "number_of_daysweek_of_moderate_physical_activity_10_minutes.0.0",
  "number_of_daysweek_of_vigorous_physical_activity_10_minutes.0.0",
  "total_cholesterol.0.0", "hdl_cholesterol_2.0.0", "ldl_cholesterol.0.0", "triglycerides.0.0"
)

# Filter and prepare hypertension data
hypertension_data <- hypertension_data[, selected_cols_ukb]
hypertension_data$Row.names <- as.numeric(hypertension_data$Row.names)
AAA_data$eid <- as.numeric(AAA_data$eid)

# Merge hypertension and AAA data
merged_data <- merge(hypertension_data, AAA_data, 
                    by.x = "Row.names", 
                    by.y = "eid", 
                    all = FALSE)

# Process protein data
not_found_proteins <- protein_labels[grep("Not found: ", protein_labels$uniprot_id), ]
proteins_to_drop <- not_found_proteins$protein
protein_data <- protein_data[, !(colnames(protein_data) %in% proteins_to_drop)]
protein_data$index <- as.numeric(protein_data$index)

# Create final dataset
full_data <- merge(merged_data, protein_data,
                       by.x = "Row.names",
                       by.y = "index",
                       all = FALSE)

# Drop Row.names and unwanted cluster columns
columns_to_drop <- c("Row.names", "date_diagnosis", "date_death", "time_to_diagnosis", "body_mass_index_bmi.1.0", "cluster1_binary", "cluster2_binary", "cluster4_binary", "Cluster")
full_data <- full_data[, !(names(full_data) %in% columns_to_drop)]

# Print all column names
cat("All columns in dataset:\n")
print(colnames(full_data))

# Create summary of NAs in each column
na_counts <- colSums(is.na(full_data))
na_columns <- na_counts[na_counts > 0]

# Print columns with NAs and their counts
cat("\nColumns with NA values and their counts:\n")
print(na_columns)

# Print total number of complete cases
cat("\nNumber of complete cases (rows with no NAs):", sum(complete.cases(full_data)), "\n")
cat("Total number of rows:", nrow(full_data), "\n")
cat("Percentage of complete cases:", round(sum(complete.cases(full_data))/nrow(full_data) * 100, 2), "%\n")

# Drop variables with high missingness
high_missing_vars <- c(
    "age_high_blood_pressure_diagnosed.0.0",                    # 15002 NAs (~51%)
    "medication_for_cholesterol_blood_pressure_or_diabetes.0.0", # 14379 NAs (~48%)
    "total_cholesterol.0.0",                                    # 12539 NAs (~42%)
    "hdl_cholesterol_2.0.0",                                    # 12539 NAs (~42%)
    "ldl_cholesterol.0.0"                                       # 12539 NAs (~42%)
)

# Create reduced dataset
reduced_data <- full_data[, !(names(full_data) %in% high_missing_vars)]

# Check NA situation in reduced dataset
na_counts_reduced <- colSums(is.na(reduced_data))
na_columns_reduced <- na_counts_reduced[na_counts_reduced > 0]

# Print summary of NAs in reduced dataset
cat("Columns with NA values in reduced dataset:\n")
print(na_columns_reduced)
cat("\nNumber of complete cases in reduced dataset:", sum(complete.cases(reduced_data)), "\n")
cat("Total number of rows:", nrow(reduced_data), "\n")
cat("Percentage of complete cases:", round(sum(complete.cases(reduced_data))/nrow(reduced_data) * 100, 2), "%\n")

# Remove prevalent cases
reduced_data <- reduced_data[reduced_data$prevalent_case == 0, ]

# Fit logistic regression
model <- glm(incident_case ~ sex.0.0 + age_at_recruitment.0.0 + IL17A + TXNDC15 + IL17A*TXNDC15, 
             family = binomial(link = "logit"), 
             data = reduced_data)
             
             
model2 <- glm(incident_case ~ sex.0.0 + age_at_recruitment.0.0 + CD34 + ITM2A + CD34*ITM2A, 
              family = binomial(link = "logit"), 
              data = reduced_data)

# View summary of the model
summary(model2)

# Check number of complete cases used in the model
cat("\nNumber of observations used in the model:", length(model$fitted.values))

# Perform multiple imputation on reduced dataset
library(mice)
imp <- mice(reduced_data, m=5, maxit=50, seed=123)

# Save imputed dataset
complete_data <- complete(imp)

# Save complete imputed dataset
write.csv(complete_data, 
          "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/joint_effect_analysis/complete_data.csv",
          row.names = FALSE)
