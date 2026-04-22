

rm(list = ls())
library(ggplot2)
library(randomForest)

load("03_Clean_Data/GSE51808_Cleaned.RData")
load("05_Models/Final_RF_Model.RData")

imp_matrix <- importance(final_rf_model)
imp_df <- data.frame(
  Gene = rownames(imp_matrix),
  Gini_Score = imp_matrix[, "MeanDecreaseGini"]
)

imp_df <- imp_df[order(imp_df$Gini_Score, decreasing = TRUE), ]
top_20_imp <- imp_df[1:20, ]

png("04_Results/RF_Feature_Importance.png", width = 1000, height = 800, res = 150)
ggplot(top_20_imp, aes(x = reorder(Gene, Gini_Score), y = Gini_Score)) +
  geom_bar(stat = "identity", fill = "#c0392b", color = "black") +    
  coord_flip() +
  labs(title = "Biomarker Importance: Why the Model Predicted DHF?",
       subtitle = "Based on Mean Decrease Gini (Internal Decision Power)",
       x = "Biomarker Genes", 
       y = "Importance Score (Mean Decrease Gini)") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 10, color = "black")
  )
dev.off()


print(head(top_20_imp, 5))
 
