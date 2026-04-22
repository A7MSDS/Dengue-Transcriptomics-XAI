# ============================================================
# سكريبت (1): تجهيز وتنظيف حجر الأساس (GSE51808) - النسخة المستقرة
# ============================================================

# تنظيف الذاكرة وتحميل المكتبات
rm(list = ls())
library(GEOquery)

cat("1. جاري قراءة الملف الخام من مجلد 01_Raw_Data...\n")
gse51808 <- getGEO(filename = "01_Raw_Data/GSE51808_series_matrix.txt.gz")

exprs_data <- exprs(gse51808)
pheno_data <- pData(gse51808)

cat("2. جاري تصنيف المرضى وعزل الحالات المطلوبة (DF vs DHF)...\n")
condition <- rep(NA, nrow(pheno_data))
# تحديد الحالات بناءً على العناوين السريرية
condition[grepl("DHF", pheno_data$title, ignore.case = TRUE)] <- "DHF"
condition[grepl("DF", pheno_data$title, ignore.case = TRUE) & !grepl("DHF", pheno_data$title, ignore.case = TRUE)] <- "DF"

pheno_data$Condition <- condition

# استبعاد أي عينات أخرى (مثل الأصحاء أو فترة النقاهة)
valid_samples <- !is.na(pheno_data$Condition)
pheno_clean <- pheno_data[valid_samples, ]
exprs_clean <- exprs_data[, valid_samples]

cat("3. جاري تنظيف أسماء الجينات وتجميع البيانات (الطريقة الآمنة)...\n")
fdata <- fData(gse51808)
gene_symbols <- sapply(strsplit(as.character(fdata$`Gene Symbol`), " /// "), `[`, 1)

# استبعاد الصفوف التي لا تحتوي على اسم جين صريح
valid_genes <- gene_symbols != "" & !is.na(gene_symbols)
exprs_genes <- exprs_clean[valid_genes, ]
symbols_valid <- gene_symbols[valid_genes]

df <- as.data.frame(exprs_genes)

# تجميع القراءات المتعددة باستخدام نظام القوائم (لتجنب إيرور المعادلات)
df_agg <- aggregate(df, by = list(Symbol = symbols_valid), FUN = mean)

# تحويل أسماء الجينات إلى أسماء صفوف
rownames(df_agg) <- df_agg$Symbol
df_agg$Symbol <- NULL

# قلب المصفوفة لتناسب خوارزمية التدريب (المرضى في الصفوف والجينات في الأعمدة)
train_data <- as.data.frame(t(df_agg))
train_data$Condition <- as.factor(pheno_clean$Condition)

cat("4. جاري حفظ البيانات النظيفة في مجلد 03_Clean_Data...\n")
# حفظ الجدول النهائي كملف RData جاهز للاستخدام في السكريبت القادم
save(train_data, file = "03_Clean_Data/GSE51808_Cleaned.RData")

cat("✅ اكتمل التجهيز بنجاح! يمكنك الآن إغلاق هذا السكريبت.\n")