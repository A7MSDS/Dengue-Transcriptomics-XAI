

rm(list = ls())
library(GEOquery)
library(pROC)
library(randomForest)

gse18090 <- getGEO(filename = "01_Raw_Data/GSE18090_series_matrix.txt.gz")
exprs_test <- exprs(gse18090)
pheno_test <- pData(gse18090)

condition_test <- rep(NA, nrow(pheno_test))

condition_test[grepl("DHF", pheno_test$title, ignore.case = TRUE)] <- "DHF"
condition_test[grepl("DF", pheno_test$title, ignore.case = TRUE) & !grepl("DHF", pheno_test$title, ignore.case = TRUE)] <- "DF"

pheno_test$Condition <- condition_test

valid_samples <- !is.na(pheno_test$Condition)
exprs_test <- exprs_test[, valid_samples]
pheno_test <- pheno_test[valid_samples, ]


print(table(pheno_test$Condition))
cat("=========================================\n\n")


load("05_Models/Final_RF_Model.RData")


fdata_test <- fData(gse18090)
gene_symbols_test <- sapply(strsplit(as.character(fdata_test$`Gene Symbol`), " /// "), `[`, 1)

keep_rows <- gene_symbols_test %in% features
exprs_test_small <- exprs_test[keep_rows, , drop = FALSE]
symbols_small <- gene_symbols_test[keep_rows]

test_df <- as.data.frame(exprs_test_small)
test_df_agg <- aggregate(test_df, by = list(Symbol = symbols_small), FUN = mean)
rownames(test_df_agg) <- test_df_agg$Symbol
test_df_agg$Symbol <- NULL

test_data <- as.data.frame(t(test_df_agg))
test_data$Condition <- as.factor(pheno_test$Condition)


for(g in setdiff(features, colnames(test_data))) {
  test_data[[g]] <- 0
}
x_test <- test_data[, features]


x_test_scaled <- as.data.frame(lapply(x_test, smart_scale))
rownames(x_test_scaled) <- rownames(x_test)

preds <- predict(final_rf_model, newdata = x_test_scaled, type = "prob")[, "DHF"]
roc_val <- roc(test_data$Condition, preds)
auc_score <- round(auc(roc_val), 3)


png("04_Results/ROC_Validation_GSE18090.png", width = 800, height = 600, res = 150)
plot(roc_val, col = "red", lwd = 5, main = "External Validation (GSE18090)")
legend("bottomright", legend = paste("AUC =", auc_score), col = "red", lwd = 5, cex = 1.2)
dev.off() 

