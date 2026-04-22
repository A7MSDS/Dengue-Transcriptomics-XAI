# ============================================================
# سكريبت (4): رسم الـ AUC والـ Heatmap لحجر الأساس (GSE51808)
# ============================================================

rm(list = ls())
library(pROC)
library(pheatmap)
library(randomForest)

# 1. تحميل الداتا والموديل
load("03_Clean_Data/GSE51808_Cleaned.RData")
load("05_Models/Final_RF_Model.RData")

# 2. حساب الـ AUC الداخلي (Internal Performance)
# ملاحظة: سنستخدم الـ Out-of-Bag (OOB) ليكون الرقم واقعياً
probs <- predict(final_rf_model, type = "prob")[, "DHF"]
roc_train <- roc(train_data$Condition, probs)
auc_val <- round(auc(roc_train), 3)

cat("📊 الـ AUC الداخلي لحجر الأساس هو:", auc_val, "\n")

# 3. حفظ رسمة الـ ROC
png("04_Results/ROC_Internal_GSE51808.png", width = 800, height = 600, res = 150)
plot(roc_train, col = "darkgreen", lwd = 5, main = paste("Internal Performance (AUC =", auc_val, ")"))
dev.off()

# 4. تجهيز الداتا للـ Heatmap (الـ 20 جين فقط)
# سنقوم بعمل Z-score scaling لتوضيح الفروقات
heatmap_data <- train_data[, features]
heatmap_data_scaled <- t(scale(heatmap_data)) # قلب المصفوفة لتكون الجينات في الصفوف

# تجهيز عمود الألوان (التصنيف)
annotation_col <- data.frame(Status = train_data$Condition)
rownames(annotation_col) <- rownames(train_data)
ann_colors = list(Status = c(DF = "#3498db", DHF = "#e74c3c")) # أزرق للـ DF وأحمر للـ DHF

# 5. رسم وحفظ الـ Heatmap
png("04_Results/Heatmap_Top20_GSE51808.png", width = 1000, height = 1200, res = 150)
pheatmap(heatmap_data_scaled, 
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         main = "Genetic Signature of Dengue Severity (Top 20 Genes)",
         show_colnames = FALSE, # لإخفاء أسماء العينات الكثيرة
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_method = "complete")
dev.off()

cat("✅ تم حفظ الـ AUC والـ Heatmap في مجلد 04_Results بنجاح!\n")