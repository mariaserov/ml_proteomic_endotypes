library(readr)


matched_data_extracted <- read_csv("/rds/general/project/hda_24-25/live/TDS/Group09/extraction_and_recoding/outputs/matched_data_extracted.csv")
df <- matched_data_extracted

# Data Cleaning ----
## missing value ----
missing_values_per_column <- colSums(is.na(df))
missing_values_per_column[missing_values_per_column > 0]
missing_values_per_column[missing_values_per_column == 0]

# checking missing value percentage
missing_percentage <- colMeans(is.na(df)) * 100
missing_percentage[missing_percentage > 0] %>% sort(decreasing = TRUE)

# delete the column which the percentage of missing value is bigger than 50%
# I'm not sure wheather this is correct, because almost all the column from 'ukb_extracted.rds' have a high percentage of missing value
df_less_missing <- df[, colMeans(is.na(df)) < 0.5]

# fit the NA using mean
df_less_missing$`Systolic blood pressure, automated reading.0.1`[is.na(df_less_missing$`Systolic blood pressure, automated reading.0.1`)] <- 
  mean(df_less_missing$`Systolic blood pressure, automated reading.0.1`, na.rm = TRUE)

df_less_missing$`Systolic blood pressure, automated reading.0.0`[is.na(df_less_missing$`Systolic blood pressure, automated reading.0.0`)] <- 
  mean(df_less_missing$`Systolic blood pressure, automated reading.0.0`, na.rm = TRUE)


## checking duplicate ----
# this will running for a long time, 0 duplicate.
duplicates <- sum(duplicated(df_less_missing))
print(paste("Duplicate Rows: ", duplicates))
df_less_missing <- df_less_missing %>% distinct()

df_no_missing <- df_less_missing

## outlier
# using box plot to check the outlier
library(ggplot2)
ggplot(df_no_missing, aes(x = "", y = `Systolic blood pressure, automated reading.0.0`)) +
  geom_boxplot() +
  theme_minimal()

# using IQR to delete the outlier
Q1 <- quantile(df_no_missing$`Systolic blood pressure, automated reading.0.0`, 0.25, na.rm = TRUE)
Q3 <- quantile(df_no_missing$`Systolic blood pressure, automated reading.0.0`, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

df_no_missing_outlier <- df_no_missing %>%
  filter(`Systolic blood pressure, automated reading.0.0` >= (Q1 - 1.5 * IQR) & `Systolic blood pressure, automated reading.0.0` <= (Q3 + 1.5 * IQR))

Q1_01 <- quantile(df_no_missing$`Systolic blood pressure, automated reading.0.1`, 0.25, na.rm = TRUE)
Q3_01 <- quantile(df_no_missing$`Systolic blood pressure, automated reading.0.1`, 0.75, na.rm = TRUE)
IQR_01 <- Q3_01 - Q1_01

df_no_missing_outlier <- df_no_missing_outlier %>%
  filter(`Systolic blood pressure, automated reading.0.1` >= (Q1_01 - 1.5 * IQR_01) & `Systolic blood pressure, automated reading.0.1` <= (Q3_01 + 1.5 * IQR_01))


# EDA ----
summary(df_no_missing_outlier)

library(ggcorrplot)

## heat map ----
cor_matrix <- cor(df_no_missing_outlier[, 4:10], use = "pairwise.complete.obs")  # chose the protein data
ggcorrplot(cor_matrix, lab = FALSE, colors = c("blue", "white", "red"))


df_no_missing_outlier %>%
  select(`Systolic blood pressure, automated reading.0.0`, `Systolic blood pressure, automated reading.0.1`) %>%  # chose numerical variable
  gather(key = "variable", value = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()



# PCA ----
library(ggfortify)

protein_data <- df_no_missing_outlier[, 4:1466]  # chose protein data
protein_scaled <- scale(protein_data)

pca_result <- prcomp(protein_scaled, center = TRUE, scale. = TRUE)

autoplot(pca_result, data = df_no_missing_outlier, colour = "gender")  # 按性别分色






