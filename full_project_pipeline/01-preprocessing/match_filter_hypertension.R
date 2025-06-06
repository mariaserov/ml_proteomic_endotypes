library(dplyr)
library(tibble)

ukb <- readRDS("../Kieran/extraction_and_recoding/outputs/ukb_extracted.rds")
prots <- readRDS("../Proteins/Data/Proteins_imputed.rds")

ukb_filtered <- ukb %>% rownames_to_column("index") %>% filter(row.names(ukb) %in% rownames(prots)) %>% 
  column_to_rownames("index")

dim(ukb)
dim(ukb_filtered)
dim(prots)


df <- merge(ukb_filtered, prots, by = "row.names", all = FALSE)
rownames(df) <- df$Row.names

print("merged dim:")

df_subset <- subset(df, systolic_blood_pressure_automated_reading.0.0 >= 140 | 
                      diastolic_blood_pressure_automated_reading.0.0 >= 90 | 
                      !is.na(age_high_blood_pressure_diagnosed.0.0) |
                      medication_for_cholesterol_blood_pressure_or_diabetes.0.0 ==2)

write.csv(df_subset, "main_hypertension_dataset.csv", row.names=TRUE)
saveRDS(df_subset, "main_hypertension_dataset.rds")