rm(list = ls())
setwd("/rds/general/project/hda_24-25/live/TDS/Group09")
library(tidyverse)
library(sharp)
library(pROC)
library(caret)

# Function to read and prepare data for each outcome
prepare_data <- function(outcome) {
  df <- read.csv(paste0("prediction/data/prediction_data_", outcome, ".csv"))
  
  # Create cluster variable (assuming we need to recreate it)
  df$Cluster <- factor(df$Cluster)
  
  return(df)
}

# Function to run models and get predictions
run_models <- function(train_data, test_data) {
  # Model 1: Covariates only
  model1 <- glm(incident_case ~ sex.0.0 + age_at_recruitment.0.0 + 
                smoking_status_clean + diabetes_clean,
                data = train_data, family = "binomial")
  
  # Model 2: Cluster only
  model2 <- glm(incident_case ~ Cluster,
                data = train_data, family = "binomial")
  
  # Model 3: Combined
  model3 <- glm(incident_case ~ sex.0.0 + age_at_recruitment.0.0 + 
                smoking_status_clean + diabetes_clean + Cluster,
                data = train_data, family = "binomial")
  
  # Get predictions for test set
  pred1 <- predict(model1, newdata = test_data, type = "response")
  pred2 <- predict(model2, newdata = test_data, type = "response")
  pred3 <- predict(model3, newdata = test_data, type = "response")
  
  return(list(pred1 = pred1, pred2 = pred2, pred3 = pred3))
}

# Function to plot ROC curves
plot_roc_curves <- function(test_labels, predictions, outcome_name) {
  roc1 <- roc(test_labels, predictions$pred1)
  roc2 <- roc(test_labels, predictions$pred2)
  roc3 <- roc(test_labels, predictions$pred3)
  
  plot(roc1, col = "red", main = paste("ROC Curve -", outcome_name),
       xlab = "False Positive Rate", ylab = "True Positive Rate")
  lines(roc2, col = "blue")
  lines(roc3, col = "darkgreen")
  
  legend("bottomright", 
         legend = c(paste("Covariates (AUC =", round(auc(roc1), 3), ")"),
                   paste("Cluster (AUC =", round(auc(roc2), 3), ")"),
                   paste("Combined (AUC =", round(auc(roc3), 3), ")")),
         col = c("red", "blue", "darkgreen"),
         lwd = 2)
}

# List of outcomes
outcomes <- c("aa", "ad", "cad", "hf", "is", "pad")

# Run analysis for each outcome
for(outcome in outcomes) {
  # Read data
  data <- prepare_data(outcome)
  
  # Set up resampling
  set.seed(123)
  train_control <- trainControl(method = "cv", 
                              number = 5,
                              classProbs = TRUE,
                              savePredictions = TRUE)
  
  # Initialize lists to store results
  cv_predictions <- list()
  cv_test_labels <- list()
  
  # Perform cross-validation
  folds <- createFolds(data$incident_case, k = 5)
  
  for(i in 1:5) {
    test_indices <- folds[[i]]
    train_data <- data[-test_indices, ]
    test_data <- data[test_indices, ]
    
    # Run models and get predictions
    predictions <- run_models(train_data, test_data)
    cv_predictions[[i]] <- predictions
    cv_test_labels[[i]] <- test_data$incident_case
  }
  
  # Combine predictions from all folds
  combined_preds <- list(
    pred1 = unlist(lapply(cv_predictions, function(x) x$pred1)),
    pred2 = unlist(lapply(cv_predictions, function(x) x$pred2)),
    pred3 = unlist(lapply(cv_predictions, function(x) x$pred3))
  )
  combined_labels <- unlist(cv_test_labels)
  
  # Plot ROC curves
  png(paste0("prediction/outputs/roc_curves_", outcome, ".png"), 
      width = 800, height = 600)
  plot_roc_curves(combined_labels, combined_preds, 
                 toupper(outcome))
  dev.off()
  
  # Calculate and save performance metrics
  metrics <- data.frame(
    Outcome = outcome,
    Model = c("Covariates", "Cluster", "Combined"),
    AUC = c(
      auc(roc(combined_labels, combined_preds$pred1)),
      auc(roc(combined_labels, combined_preds$pred2)),
      auc(roc(combined_labels, combined_preds$pred3))
    )
  )
  
  write.csv(metrics, 
            paste0("prediction/outputs/metrics_", outcome, ".csv"),
            row.names = FALSE)
}
