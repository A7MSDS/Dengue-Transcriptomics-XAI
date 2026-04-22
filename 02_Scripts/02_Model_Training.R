# ============================================================
# سكريبت (2): تدريب خوارزمية Random Forest واستخراج أهم 20 جين
# ============================================================

# 1. تنظيف واجهة العمل وتحميل المكتبات
rm(list = ls())
library(randomForest)

cat("1. جاري تحميل البيانات النظيفة من مجلد 03_Clean_Data...\n")
load("03_Clean_Data/GSE51808_Cleaned.RData")

# 2. فصل الجينات عن التشخيص للتدريب الآمن (Matrix Interface)
x_data <- train_data[, colnames(train_data) != "Condition"]
y_label <- train_data$Condition

cat("2. جاري تدريب النموذج الأولي لاكتشاف الجينات المؤثرة (قد يستغرق دقيقة)...\n")
set.seed(123) 
initial_rf <- randomForest(x = x_data, y = y_label, ntree = 500, importance = TRUE)

# 3. استخراج أهم 20 جين
var_imp <- importance(initial_rf)
top_genes_df <- var_imp[order(var_imp[, "MeanDecreaseGini"], decreasing = TRUE), ]
top_20_genes <- rownames(top_genes_df)[1:20]

cat("3. تم استخراج الجينات! جاري حفظها في ملف Excel في مجلد النتائج...\n")
# حفظ الجينات في ملف CSV داخل مجلد 04_Results
write.csv(data.frame(Gene_Symbol = top_20_genes, Importance_Score = top_genes_df[1:20, "MeanDecreaseGini"]),
          file = "04_Results/Top_20_Biomarkers.csv", row.names = FALSE)

# 4. المعايرة الذكية (Z-score Scaling) لتجهيز الموديل النهائي للاختبارات القادمة
cat("4. جاري معايرة المقاييس وتدريب الموديل النهائي...\n")
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

cat("5. جاري حفظ الموديل النهائي في مجلد 05_Models...\n")
# حفظ الموديل النهائي، ودالة المعايرة، وأسماء الجينات كـ "حزمة جاهزة" للسكريبت القادم
save(final_rf_model, smart_scale, features, file = "05_Models/Final_RF_Model.RData")

cat("✅ تمت العملية بنجاح! راجع مجلد 04_Results و 05_Models.\n")