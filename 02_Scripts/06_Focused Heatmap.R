
rm(list = ls())
if(!require(pheatmap)) install.packages("pheatmap")
if(!require(randomForest)) install.packages("randomForest")
library(pheatmap)
library(randomForest)

load("03_Clean_Data/GSE51808_Cleaned.RData")
load("05_Models/Final_RF_Model.RData")

imp_matrix <- importance(final_rf_model)
imp_df <- data.frame(
  Gene = rownames(imp_matrix),
  Gini_Score = imp_matrix[, "MeanDecreaseGini"]
)
imp_df <- imp_df[order(imp_df$Gini_Score, decreasing = TRUE), ]
top_genes <- as.character(imp_df$Gene[1:20])

heat_data <- t(train_data[, top_genes])

annotation_col <- data.frame(Condition = train_data$Condition)
rownames(annotation_col) <- rownames(train_data)

ann_colors <- list(
  Condition = c(DHF = "firebrick", DF = "navy")
)

png("04_Results/Focused_XAI_Heatmap.png", width = 800, height = 600, res = 150)
pheatmap(heat_data,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         show_colnames = FALSE,
         cluster_rows = FALSE,       
         main = "Directional Impact of Top 20 Predictors (Explainable AI)")
dev.off()

