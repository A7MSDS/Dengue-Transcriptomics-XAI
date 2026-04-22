

if(!require(DALEX)) install.packages("DALEX")
library(DALEX)
library(ggplot2)
library(randomForest)

load("03_Clean_Data/GSE51808_Cleaned.RData")
load("05_Models/Final_RF_Model.RData")

model_features <- rownames(final_rf_model$importance)
X_df <- as.data.frame(train_data[, model_features])

y_target <- ifelse(train_data$Condition == "DHF", 1, 0)

rf_explainer <- explain(
  model = final_rf_model,
  data = X_df,
  y = y_target,
  predict_function = function(model, newdata) {
    
    as.numeric(predict(model, newdata, type = "prob")[, "DHF"])
  },
  label = "Random Forest (DHF Predictor)",
  colorize = FALSE
)

cat("⏳ جاري اختبار تأثير كل جين على دقة الموديل (قد يستغرق نصف دقيقة)...\n")
set.seed(123)
rf_parts <- model_parts(rf_explainer, B = 50) 

png("04_Results/XAI_DALEX_Importance.png", width = 1000, height = 800, res = 150)
plot(rf_parts, show_boxplots = FALSE) +
  ggtitle("Explainable AI: Global Feature Importance",
          subtitle = "How dropping each gene affects the model's predictive power") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(face = "bold"),
    plot.title = element_text(color = "darkblue")
  )
dev.off()

