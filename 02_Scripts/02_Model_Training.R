
rm(list = ls())
library(randomForest)

load("03_Clean_Data/GSE51808_Cleaned.RData")

x_data <- train_data[, colnames(train_data) != "Condition"]
y_label <- train_data$Condition

set.seed(123) 
initial_rf <- randomForest(x = x_data, y = y_label, ntree = 500, importance = TRUE)

var_imp <- importance(initial_rf)
top_genes_df <- var_imp[order(var_imp[, "MeanDecreaseGini"], decreasing = TRUE), ]
top_20_genes <- rownames(top_genes_df)[1:20]

# حفظ الجينات في ملف CSV داخل مجلد 04_Results
write.csv(data.frame(Gene_Symbol = top_20_genes, Importance_Score = top_genes_df[1:20, "MeanDecreaseGini"]),
          file = "04_Results/Top_20_Biomarkers.csv", row.names = FALSE)

smart_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if(s == 0) return(rep(0, length(x)))
  return((x - mean(x, na.rm = TRUE)) / s)
}

x_train_top20 <- train_data[, top_20_genes]
x_train_scaled <- as.data.frame(lapply(x_train_top20, smart_scale))
rownames(x_train_scaled) <- rownames(train_data)

set.seed(123)
final_rf_model <- randomForest(x = x_train_scaled, y = y_label, ntree = 500, importance = TRUE)
features <- top_20_genes

save(final_rf_model, smart_scale, features, file = "05_Models/Final_RF_Model.RData")

