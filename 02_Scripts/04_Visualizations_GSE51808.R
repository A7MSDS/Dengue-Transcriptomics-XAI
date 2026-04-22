

rm(list = ls())
library(pROC)
library(pheatmap)
library(randomForest)

load("03_Clean_Data/GSE51808_Cleaned.RData")
load("05_Models/Final_RF_Model.RData")

probs <- predict(final_rf_model, type = "prob")[, "DHF"]
roc_train <- roc(train_data$Condition, probs)
auc_val <- round(auc(roc_train), 3)



png("04_Results/ROC_Internal_GSE51808.png", width = 800, height = 600, res = 150)
plot(roc_train, col = "darkgreen", lwd = 5, main = paste("Internal Performance (AUC =", auc_val, ")"))
dev.off()


heatmap_data <- train_data[, features]
heatmap_data_scaled <- t(scale(heatmap_data)) 

annotation_col <- data.frame(Status = train_data$Condition)
rownames(annotation_col) <- rownames(train_data)
ann_colors = list(Status = c(DF = "#3498db", DHF = "#e74c3c"))


png("04_Results/Heatmap_Top20_GSE51808.png", width = 1000, height = 1200, res = 150)
pheatmap(heatmap_data_scaled, 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         main = "Genetic Signature of Dengue Severity (Top 20 Genes)",
         show_colnames = FALSE, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_method = "complete")
dev.off()

