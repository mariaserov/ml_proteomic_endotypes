rm(list = ls())
# setwd("/rds/general/project/hda_24-25/live/TDS/Group09")
library(tidyverse)
library(fake)
library(sharp)
library(pROC)

# Function to read and prepare data for each outcome
prepare_data <- function(outcome) {
  df <- read.csv(paste0("prediction/data/prediction_data_", outcome, ".csv"))
  df$Cluster <- factor(df$Cluster)
  return(df)
}

# Function to plot ROC curves - simplified version
plot_roc_curves <- function(test_labels, predictions, outcome_name) {
  # Ensure all inputs are proper numeric vectors
  labels <- as.numeric(test_labels)
  pred1 <- as.numeric(predictions$pred1)
  pred2 <- as.numeric(predictions$pred2)
  pred3 <- as.numeric(predictions$pred3)
  
  # Calculate ROC curves manually
  calc_roc <- function(labels, preds) {
    # Sort by predictions
    ord <- order(preds, decreasing = TRUE)
    labels <- labels[ord]
    
    # Calculate TPR and FPR
    tp <- cumsum(labels)
    fp <- cumsum(1 - labels)
    tpr <- tp / sum(labels)
    fpr <- fp / sum(1 - labels)
    
    # Calculate AUC
    auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1))) / 2
    
    list(fpr = fpr, tpr = tpr, auc = auc)
  }
  
  # Calculate ROC curves
  roc1 <- calc_roc(labels, pred1)
  roc2 <- calc_roc(labels, pred2)
  roc3 <- calc_roc(labels, pred3)
  
  # Create plot
  plot(0:1, 0:1, type = "n", 
       main = paste("ROC Curves -", outcome_name),
       xlab = "False Positive Rate", 
       ylab = "True Positive Rate")
  
  # Add curves
  lines(roc1$fpr, roc1$tpr, col = "red")
  lines(roc2$fpr, roc2$tpr, col = "blue")
  lines(roc3$fpr, roc3$tpr, col = "darkgreen")
  abline(0, 1, lty = 2, col = "gray")
  
  # Add legend
  legend("bottomright", 
         legend = c(paste("Covariates (AUC =", round(roc1$auc, 3), ")"),
                    paste("Cluster (AUC =", round(roc2$auc, 3), ")"),
                    paste("Combined (AUC =", round(roc3$auc, 3), ")")),
         col = c("red", "blue", "darkgreen"),
         lwd = 2)
  
  # Return AUCs for metrics
  list(auc1 = roc1$auc, auc2 = roc2$auc, auc3 = roc3$auc)
}

# List of outcomes
outcomes <- c("aa", "ad", "cad", "hf", "is", "pad")

# Run analysis for each outcome
for(outcome in outcomes) {
  cat("Starting analysis for", toupper(outcome), "\n")
  
  # Prepare data
  data <- prepare_data(outcome)
  
  # Split data into training and testing sets for proper evaluation
  set.seed(42)
  train_indices <- Split(data = data$incident_case, family = "binomial", tau = c(0.8, 0.2))
  train_data <- data[train_indices[[1]], ]
  test_data <- data[train_indices[[2]], ]
  
  # Define models and run them using stability selection
  # Model 1: Covariates only
  model1 <- VariableSelection(
    xdata = as.matrix(model.matrix(~ sex.0.0 + age_at_recruitment.0.0 + 
                                     smoking_status_clean + diabetes_clean - 1, train_data)),
    ydata = train_data$incident_case,
    family = "binomial",
    K = 100  # Number of resampling iterations
  )
  
  # Model 2: Cluster only
  model2 <- VariableSelection(
    xdata = as.matrix(model.matrix(~ Cluster - 1, train_data)),
    ydata = train_data$incident_case,
    family = "binomial",
    K = 100
  )
  
  # Model 3: Combined
  model3 <- VariableSelection(
    xdata = as.matrix(model.matrix(~ sex.0.0 + age_at_recruitment.0.0 + 
                                     smoking_status_clean + diabetes_clean + 
                                     Cluster - 1, train_data)),
    ydata = train_data$incident_case,
    family = "binomial",
    K = 100
  )
  
  # Generate predictions for test data
  pred1 <- as.numeric(predict(
    object = model1, 
    xdata = as.matrix(model.matrix(~ sex.0.0 + age_at_recruitment.0.0 + 
                                     smoking_status_clean + diabetes_clean - 1, train_data)),
    ydata = train_data$incident_case,
    newdata = as.matrix(model.matrix(~ sex.0.0 + age_at_recruitment.0.0 + 
                                       smoking_status_clean + diabetes_clean - 1, test_data)),
    method = "refit",
    type = "response"
  ))
  
  pred2 <- as.numeric(predict(
    object = model2, 
    xdata = as.matrix(model.matrix(~ Cluster - 1, train_data)),
    ydata = train_data$incident_case,
    newdata = as.matrix(model.matrix(~ Cluster - 1, test_data)),
    method = "refit",
    type = "response"
  ))
  
  pred3 <- as.numeric(predict(
    object = model3, 
    xdata = as.matrix(model.matrix(~ sex.0.0 + age_at_recruitment.0.0 + 
                                     smoking_status_clean + diabetes_clean + 
                                     Cluster - 1, train_data)),
    ydata = train_data$incident_case,
    newdata = as.matrix(model.matrix(~ sex.0.0 + age_at_recruitment.0.0 + 
                                       smoking_status_clean + diabetes_clean + 
                                       Cluster - 1, test_data)),
    method = "refit",
    type = "response"
  ))
  
  # Calculate and plot ROC curves
  png(paste0("/rds/general/project/hda_24-25/live/TDS/Group09/prediction/sharp_output/roc_curves_", outcome, ".png"), 
      width = 800, height = 600)
  
  aucs <- plot_roc_curves(test_data$incident_case, 
                          list(
                            pred1 = pred1,
                            pred2 = pred2,
                            pred3 = pred3
                          ),
                          toupper(outcome))
  
  dev.off()
  
  # Create metrics without using pROC
  metrics <- data.frame(
    Outcome = outcome,
    Model = c("Covariates", "Cluster", "Combined"),
    AUC = c(aucs$auc1, aucs$auc2, aucs$auc3)
  )
  
  # Save metrics
  write.csv(metrics, 
            paste0("/rds/general/project/hda_24-25/live/TDS/Group09/prediction/sharp_output/metrics_", outcome, ".csv"),
            row.names = FALSE)
  
  # Print progress
  cat("Completed analysis for", toupper(outcome), "\n")
}