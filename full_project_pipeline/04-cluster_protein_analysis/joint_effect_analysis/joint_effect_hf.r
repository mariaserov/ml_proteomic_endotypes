setwd("/rds/general/project/hda_24-25/live/TDS/Group09")
library(tidyverse)

# import data - ACTION REQUIRED - must select folder name
outcome <- readRDS("./outcome_definition/Outputs_saved/HF/output_final.rds") 
ukb_htn <- readRDS("./Data/main_hypertension_dataset.rds")

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

# ensure eid is numeric 
df$eid <- as.numeric(df$eid)
outcome$eid <- as.numeric(outcome$eid)

# merge
df <- merge(df, outcome, by = "eid", all.x = TRUE)
rownames(df) <- NULL
df <- df[!is.na(df$incident_case), ] #drop values for which we do not have outcome data

# Add these lines:
# Read protein data
protein_data <- readRDS("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/data/protein_data_with_clusters.rds")
protein_labels <- read.csv("/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/pathway_analysis_data/cluster_3_dea.csv")

# Process protein data
not_found_proteins <- protein_labels[grep("Not found: ", protein_labels$uniprot_id), ]
proteins_to_drop <- not_found_proteins$protein
protein_data <- protein_data[, !(colnames(protein_data) %in% proteins_to_drop)]
protein_data$index <- as.numeric(protein_data$index)

# Merge with protein data
df_hf <- merge(df, protein_data,
            by.x = "eid",
            by.y = "index",
            all = FALSE)

df_hf <- df_hf %>% 
  select(-Cluster, -cluster1_binary, -cluster2_binary, -cluster4_binary)

library(glmnet)

# Get protein columns
start_idx <- which(colnames(df_hf) == "AARSD1")
end_idx <- ncol(df_hf) - 1  # Second last column
protein_cols <- colnames(df_hf)[start_idx:end_idx]

# Check for NAs in columns we'll use in the models
model_cols <- c("incident_case", "sex.0.0", "age_at_recruitment.0.0", 
                "smoking_status_clean", "diabetes_clean", "cluster3_binary", protein_cols)
na_counts <- colSums(is.na(df_hf[, model_cols]))

# Print NA counts if any exist
if(any(na_counts > 0)) {
    cat("\nMissing values in model columns:\n")
    print(na_counts[na_counts > 0])
    
    # Remove rows with NAs in these columns only
    df_hf <- df_hf[complete.cases(df_hf[, model_cols]), ]
    cat("\nRows remaining after removing NAs:", nrow(df_hf), "\n")
}

model1 <- glm(incident_case ~ cluster3_binary + sex.0.0 + age_at_recruitment.0.0 + 
              smoking_status_clean + diabetes_clean,
              family = binomial(link = "logit"),
              data = df_hf)

# Create model matrix
X2 <- model.matrix(~ . - 1, data = df_hf[, c("sex.0.0", "age_at_recruitment.0.0", 
                                            "smoking_status_clean", "diabetes_clean", protein_cols)])
y <- df_hf$incident_case

# Fit ridge regression
set.seed(123) # for reproducibility
cv_fit2 <- cv.glmnet(X2, y, family = "binomial", alpha = 0, nfolds = 10)
model2 <- glmnet(X2, y, family = "binomial", alpha = 0, lambda = cv_fit2$lambda.min)

# Calculate model 2 performance
probs2 <- predict(model2, X2, type = "response")
preds2 <- ifelse(probs2 > 0.5, 1, 0)
accuracy2 <- mean(preds2 == y)

#######################
# Model 3: Combined Ridge Regression
#######################
# Create model matrix including cluster3_binary
X3 <- model.matrix(~ . - 1, data = df_hf[, c("sex.0.0", "age_at_recruitment.0.0", 
                                            "smoking_status_clean", "diabetes_clean",
                                            protein_cols, "cluster3_binary")])

# Fit combined ridge regression
cv_fit3 <- cv.glmnet(X3, y, family = "binomial", alpha = 0, nfolds = 10)
model3 <- glmnet(X3, y, family = "binomial", alpha = 0, lambda = cv_fit3$lambda.min)

# Calculate model 3 performance
probs3 <- predict(model3, X3, type = "response")
preds3 <- ifelse(probs3 > 0.5, 1, 0)
accuracy3 <- mean(preds3 == y)

# Extract and print cluster3_binary coefficient from model 3
cluster3_coef <- coef(model3)[which(rownames(coef(model3)) == "cluster3_binary")]
cat("\nCluster3_binary coefficient in combined model:", cluster3_coef, "\n")

# Add AUC calculations and handle imbalanced data
library(pROC)

# For model 1
probs1 <- predict(model1, type = "response")
roc1 <- roc(df_hf$incident_case, probs1)
auc1 <- auc(roc1)

# For model 2
roc2 <- roc(y, probs2)
auc2 <- auc(roc2)

# For model 3
roc3 <- roc(y, probs3)
auc3 <- auc(roc3)

# Create output directory if it doesn't exist
output_dir <- "/rds/general/project/hda_24-25/live/TDS/Group09/cluster_protein_analysis/joint_effect_analysis/Results"
if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Create log file
log_file <- file.path(output_dir, "hf_results.log")
sink(log_file)

# Print all results to log file
cat("=== Analysis Results for HF ===\n\n")

cat("Sample size after removing NAs:", nrow(df_hf), "\n\n")

cat("Class distribution in outcome:\n")
print(table(df_hf$incident_case))
cat("\n")

cat("Model 1: Basic Logistic Regression Results\n")
print(summary(model1))
cat("\n")

cat("Model 2: Ridge Regression Results\n")
cat("Optimal lambda (min):", cv_fit2$lambda.min, "\n")
cat("Optimal lambda (1se):", cv_fit2$lambda.1se, "\n")
cat("Number of non-zero coefficients:", sum(coef(model2) != 0), "\n")
cat("Accuracy:", accuracy2, "\n")

cat("Model 3: Combined Ridge Regression Results\n")
cat("Optimal lambda (min):", cv_fit3$lambda.min, "\n")
cat("Optimal lambda (1se):", cv_fit3$lambda.1se, "\n")
cat("Number of non-zero coefficients:", sum(coef(model3) != 0), "\n")
cat("Accuracy:", accuracy3, "\n")
cat("Cluster3_binary coefficient:", cluster3_coef, "\n")

cat("AUC Scores:\n")
cat("Model 1 (Logistic with cluster3_binary only):", auc1, "\n")
cat("Model 2 (Ridge with proteins only):", auc2, "\n")
cat("Model 3 (Ridge with proteins + cluster3_binary):", auc3, "\n")
cat("AUC improvement (Model 3 - Model 2):", auc3 - auc2, "\n")

# Close the log file
sink()

# Save models and results
results <- list(
    model1 = model1,
    model2 = model2,
    model3 = model3,
    cv_fit2 = cv_fit2,
    cv_fit3 = cv_fit3,
    performance = list(
        accuracy_model2 = accuracy2,
        accuracy_model3 = accuracy3,
        auc_model1 = auc1,
        auc_model2 = auc2,
        auc_model3 = auc3,
        class_distribution = table(df_hf$incident_case)
    )
)
saveRDS(results, file.path(output_dir, "hf_models.rds"))

# Create and save plots
pdf(file.path(output_dir, "hf_ridge_plots.pdf"))
par(mfrow = c(2,2))
plot(cv_fit2, main = "Model 2: Cross-validation")
plot(cv_fit3, main = "Model 3: Cross-validation")
plot(model2, xvar = "lambda", main = "Model 2: Coefficient Paths")
plot(model3, xvar = "lambda", main = "Model 3: Coefficient Paths")
dev.off()