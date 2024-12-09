setwd("D:/R/项目/AD/实验")

library(mice)
library(nortest)
library(VIM)
library(tidyverse)

####处理缺失值####
# 加载数据
ndata <- read.csv("alzheimer_combined.csv", row.names = 1)
# 筛选需要处理的数据
data <- ndata[, c(4:ncol(ndata))]
# # 可视化缺失值
# aggr(data, numbers = FALSE, prop = TRUE, cex.axis = 0.7)
# dev.off()
# 计算缺失比例
missing_ratio <- colMeans(is.na(data))
missing_data <- data.frame(missing_ratio)
write.csv(missing_data, "missing_data.csv")
# 剔除缺失比例 > 40%的变量
data <- data[, missing_ratio <= 0.4]
# 定义函数检查正态性
is_normal <- function(x) {
  ad_test <- ad.test(x[!is.na(x)]) # 使用Anderson-Darling正态性检验
  return(ad_test$p.value > 0.05)
}
# 处理缺失比例 < 5%的变量
for (col in names(data)) {
  if (missing_ratio[col] < 0.05) {
    if (is.numeric(data[[col]])) {
      if (is_normal(data[[col]])) {
        # 正态分布，使用均值填补
        data[[col]][is.na(data[[col]])] <- mean(data[[col]], na.rm = TRUE)
      } else {
        # 偏态分布，使用中位数填补
        data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)
      }
    }
  }
}
# 处理缺失比例 > 5%的变量
data <- mice(data, m = 5, method = 'pmm', maxit = 50, seed = 500)
# 完成多重插补并获取完整数据集
completed_data <- complete(data)
# 拼接完整的数据
combined_data <- cbind(ndata[, c(1:2)], completed_data)
# 整理数据
colnames(combined_data)[11] <- "CCI"
#输出处理后的数据
write.csv(combined_data, "alzheimer_combined_imputed.csv")

####划分数据集####
combined_data <- read.csv("alzheimer_combined_imputed.csv", row.names = 1)
table(combined_data$CCI <= 5)
combined_data$Group <- ifelse(combined_data$CCI <= 5, "CCI-low", "CCI-high")
write.csv(combined_data, "alzheimer_combined_imputed.csv")
# 设置随机种子以保证结果可重复
set.seed(123)
# 获取数据集的行数
n <- nrow(combined_data)
# 计算训练集的行数（70%的比例）
train_index <- sample(1:n, size = round(0.7 * n))
# 划分训练集和验证集
train_data <- combined_data[train_index, ]
validation_data <- combined_data[-train_index, ]
#保存数据集
write.csv(train_data, "train_data.csv")
write.csv(validation_data, "validation_data.csv")

####基线表#####
library(compareGroups)
combined_data <- read.csv("alzheimer_combined_imputed.csv", row.names = 1)
# 整理基线表的变量
vars <- colnames(combined_data)
# 使用sapply生成一个命令字符串，每个列名对应一个 = NA
method_string <- paste(vars, "= NA", sep = "", collapse = ", ")
# 将Group按照"CCI-low"和"CCI-high"的顺序因子化，这样"CCI-low"就会在"CCI-high"之前
combined_data$Group <- factor(combined_data$Group, levels = c("aCCI-low", "aCCI-high"))
# 生成基线表，分组变量为Group
datatab <- descrTable(Group~., 
                      data = combined_data,
                      method = method_string,
                      show.all = TRUE,
                      show.p.overall = TRUE  # 禁用P值显示
)
datatab
# 打印表格，显示全部变量水平
export2word(datatab, "datatab.docx")
export2csv(datatab, "datatab.csv")

####多因素逻辑回归####
library(rms)
library(pROC)

train <- read.csv("train_data.csv", row.names = 1)
validation <- read.csv("validation_data.csv", row.names = 1)
var <- c("age",	"resprate",	"MCHC",	"Base.Excess",	"Anion.Gap", "Glucose",
         "RDW",	"Alkaline.Phosphatase",	"Potassium..Whole.Blood",	"Hematocrit",
         "Urea.Nitrogen", "Phosphate", "Creatinine", "Hemoglobin",	"MCH", "Group"
)
train <- train[, var]
validation <- validation[, var]

# 检查数据类型
str(train)
train$Group <- factor(train$Group, levels = c("aCCI-low", "aCCI-high"))
validation$Group <- factor(validation$Group, levels = c("aCCI-low", "aCCI-high"))

set.seed(1234)
logistic_model <- glm(Group ~ ., data = train,
                      family = binomial())

summary(logistic_model)

train_pre_lg <- predict(logistic_model, newdata = train, type = "response")
validation_pre_lg <- predict(logistic_model, newdata = validation, type = "response")

train_roc_lg <- roc(train$Group, train_pre_lg)
validation_roc_lg <- roc(validation$Group, validation_pre_lg)

train_auc_lg <- auc(train_roc_lg)
validation_auc_lg <- auc(validation_roc_lg)

plot(train_roc_lg, col = "#FD8D62", lwd = 3, main = "ROC Curve for Training and Validation Sets")
lines(validation_roc_lg, col = "#66C3A5", lwd = 3)
legend_text <- c(paste0("Training (AUC = ", round(train_auc_lg, 3), ")"),
                 paste0("Validation (AUC = ", round(validation_auc_lg, 3), ")"))
legend("bottomright", legend = legend_text, col = c("#FD8D62", "#66C3A5"),lwd=2
       # ,cex = 0.8, box.lwd = 0, inset = c(0.01, 0.01)
)
dev.off()

####LASSO####
library(glmnet)
library(foreign)

# 读取文件
train <- read.csv("train_data.csv", row.names = 1)
validation <- read.csv("validation_data.csv", row.names = 1)
var <- c("age",	"resprate",	"MCHC",	"Base.Excess",	"Anion.Gap", "Glucose",
         "RDW",	"Alkaline.Phosphatase",	"Potassium..Whole.Blood",	"Hematocrit",
         "Urea.Nitrogen", "Phosphate", "Creatinine", "Hemoglobin",	"MCH", "Group"
)
train <- train[, var]
validation <- validation[, var]

# 检查数据类型
str(train)
train$Group <- ifelse(train$Group=='aCCI-low',0,1)
validation$Group <- ifelse(validation$Group=='aCCI-low',0,1)


# 设置随机种子，K折交叉验证，每次的数据都是随机的，随机数种子一致，就结果给固定住。
set.seed(1234)
y_train <- as.matrix(train[, ncol(train)])
x_train <- as.matrix(train[,c(1:ncol(train)-1)])
# 构建模型
fit = glmnet(x_train, y_train, family="binomial", alpha=1) #这里alpha=1为LASSO回归，如果等于0就是岭回归

plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit=cv.glmnet(x_train, y_train, family="binomial", nfolds = 10)
plot(cvfit)
dev.off()

# 输出筛选的特征基因
best_lambda <- cvfit$lambda.min
coef <- coef(fit, s = best_lambda)
index <- which(coef != 0)
actCoef <- coef[index]
lassovar  <- row.names(coef)[index]
varCoef <- cbind(Feature=lassovar,Coef=actCoef)
varCoef <- varCoef[-1, ]
varCoef 
# write.csv(, file = "lasso_varcoff.csv")

# 训练集上的预测概率
train_pre_lasso <- predict(fit, newx = x_train, s = best_lambda, type = "response")

# 验证集上的预测概率
x_validation <- as.matrix(validation[, c(1:(ncol(validation)-1))])
y_validation <- as.matrix(validation[, ncol(validation)])
validation_pre_lasso <- predict(fit, newx = x_validation, s = best_lambda, type = "response")

# 计算训练集的ROC曲线
train_roc_lasso <- roc(y_train, as.vector(train_pre_lasso))

# 计算验证集的ROC曲线
validation_roc_lasso <- roc(y_validation, as.vector(validation_pre_lasso))

# 计算AUC值
train_auc_lasso <- auc(validation_roc_lasso)
validation_auc_lasso <- auc(validation_roc_lasso)

# 绘制训练集和验证集的ROC曲线
plot(train_roc_lasso, col = "#FD8D62", lwd = 3, main = "ROC Curve for Training and Validation Sets")
lines(validation_roc_lasso, col = "#66C3A5", lwd = 3)
# 在图例中添加AUC值
legend_text <- c(paste0("Training (AUC = ", round(train_auc_lasso, 3), ")"),
                 paste0("Validation (AUC = ", round(validation_auc_lasso, 3), ")"))
legend("bottomright", legend = legend_text, col = c("#FD8D62", "#66C3A5"), lty=1, lwd=2
       # ,cex = 0.8, box.lwd = 0, inset = c(0.01, 0.01)
)
dev.off()

####随机森林####
# 加载必要的包
library(randomForestSRC)
library(pROC)

# 读取数据
train <- read.csv("train_data.csv", row.names = 1)
validation <- read.csv("validation_data.csv", row.names = 1)
var <- c("age",	"resprate",	"MCHC",	"Base.Excess",	"Anion.Gap", "Glucose",
         "RDW",	"Alkaline.Phosphatase",	"Potassium..Whole.Blood",	"Hematocrit",
         "Urea.Nitrogen", "Phosphate", "Creatinine", "Hemoglobin",	"MCH", "Group"
)
train <- train[, var]
validation <- validation[, var]

# 检查数据类型
str(train)
train$Group <- factor(train$Group, levels = c("aCCI-low", "aCCI-high"))
validation$Group <- factor(validation$Group, levels = c("aCCI-low", "aCCI-high"))

# 参数 cv = TRUE 进行十折交叉验证，cv = FALSE 不进行交叉验证
set.seed(1234)
rf_model <- rfsrc(
  Group ~ ., 
  data = train, 
  ntree = 1000,
  mtry = sqrt(ncol(train) - 1) / 2, # 减少 mtry，减少模型复杂度
  nodesize = 15, # 增大 nodesize，减少过拟合
  nsplit = 2,
  samptype = "swr", # 采用有放回采样
  splitrule = "gini",
  importance = TRUE
)

print(rf_model)
plot(rf_model)
dev.off()

# 查看变量重要性(VIMP+min_depth)
library(ggRandomForests)
gg_dta <- gg_minimal_vimp(rf_model)
plot(gg_dta)
# write.csv(gg_dta, "RF_importance.csv")
dev.off()

# 训练集上的预测概率
train_pre_rf <- predict(rf_model, newdata = train)$predicted
# 验证集上的预测概率
validation_pre_rf <- predict(rf_model, newdata = validation)$predicted

# 计算训练集的ROC曲线
train_roc_rf <- roc(train$Group, train_pre_rf[, "aCCI-high"])

# 计算验证集的ROC曲线
validation_roc_rf <- roc(validation$Group, validation_pre_rf[, "aCCI-high"])

# 计算AUC值
train_auc_rf <- auc(train_roc_rf)
validation_auc_rf <- auc(validation_roc_rf)

# 绘制训练集和验证集的ROC曲线
plot(train_roc_rf, col = "#FD8D62", lwd = 3, main = "ROC Curve for Training and Validation Sets")
lines(validation_roc_rf, col = "#66C3A5", lwd = 3)
# 在图例中添加AUC值
legend_text <- c(paste0("Training (AUC = ", round(train_auc_rf, 3), ")"),
                 paste0("Validation (AUC = ", round(validation_auc_rf, 3), ")"))
legend("bottomright", legend = legend_text, col = c("#FD8D62", "#66C3A5"), lty=1, lwd=2
       # ,cex = 0.8, box.lwd = 0, inset = c(0.01, 0.01)
)
dev.off()

####XGBoost####
library(xgboost)
library(caret)
library(pROC)

# 读取数据
train <- read.csv("train_data.csv", row.names = 1)
validation <- read.csv("validation_data.csv", row.names = 1)
var <- c("age",	"resprate",	"MCHC",	"Base.Excess",	"Anion.Gap", "Glucose",
         "RDW",	"Alkaline.Phosphatase",	"Potassium..Whole.Blood",	"Hematocrit",
         "Urea.Nitrogen", "Phosphate", "Creatinine", "Hemoglobin",	"MCH", "Group"
)
train <- train[, var]
validation <- validation[, var]

# 检查数据类型
str(train)
train$Group <- factor(make.names(train$Group))
validation$Group <- factor(make.names(validation$Group))

# # 转换标签形式
# train$Group <- as.numeric(train$Group)-1
# Validation$Group <- as.numeric(Validation$Group)-1

# 构建XGBoost模型
set.seed(1234)
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE)
xgb_grid <- expand.grid(
  nrounds = 100,
  max_depth = 4,
  eta = 0.05,
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

xgb_model_caret <- train(
  Group ~ ., data = train,
  method = "xgbTree",
  trControl = train_control,
  tuneGrid = xgb_grid,
  metric = "ROC",
  lambda = 10,       # 增加 L2 正则化力度
  alpha = 1,       # 增加 L1 正则化力度
  objective = "binary:logistic" 
)


# 查看模型结果
print(xgb_model_caret)

# 获取特征重要性
xgb_model <- xgb_model_caret$finalModel
importance <- xgb.importance(feature_names = colnames(train[, -ncol(train)]), model = xgb_model)
write.csv(importance, "XGBoost_importance.csv")

# 绘制特征重要性
xgb.plot.importance(importance)

# 训练集上的预测概率
train_pre_xgb <- predict(xgb_model, newdata = as.matrix(train[, -ncol(train)]))
# 验证集上的预测概率
validation_pre_xgb <- predict(xgb_model, newdata = as.matrix(validation[, -ncol(validation)]))

# 计算训练集的ROC曲线
train_roc_xgb <- roc(train$Group, train_pre_xgb)

# 计算验证集的ROC曲线
validation_roc_xgb <- roc(validation$Group, validation_pre_xgb)

# 计算AUC值
train_auc_xgb <- auc(train_roc_xgb)
validation_auc_xgb <- auc(validation_roc_xgb)

# 绘制训练集和验证集的ROC曲线
plot(train_roc_xgb, col = "#FD8D62", lwd = 3, main = "ROC Curve for Training and Validation Sets")
lines(validation_roc_xgb, col = "#66C3A5", lwd = 3)
# 在图例中添加AUC值
legend_text <- c(paste0("Training (AUC = ", round(train_auc_xgb, 3), ")"),
                 paste0("Validation (AUC = ", round(validation_auc_xgb, 3), ")"))
legend("bottomright", legend = legend_text, col = c("#FD8D62", "#66C3A5"), lty=1, lwd=2
       # ,cex = 0.8, box.lwd = 0, inset = c(0.01, 0.01)
)
dev.off()

####SVM####
library(e1071)
library(caret)
library(pROC)

# 读取数据
train <- read.csv("train_data.csv", row.names = 1)
validation <- read.csv("validation_data.csv", row.names = 1)
var <- c("age",	"resprate",	"MCHC",	"Base.Excess",	"Anion.Gap", "Glucose",
         "RDW",	"Alkaline.Phosphatase",	"Potassium..Whole.Blood",	"Hematocrit",
         "Urea.Nitrogen", "Phosphate", "Creatinine", "Hemoglobin",	"MCH", "Group"
)
train <- train[, var]
validation <- validation[, var]

# 调整数据
train[1:3, ]
train <- train[, c(16, 1:15)]
validation <- validation[, c(16, 1:15)]

# 转换标签形式
str(train)
train$Group <- factor(train$Group, levels = c("aCCI-low", "aCCI-high"))
validation$Group <- factor(validation$Group, levels = c("aCCI-low", "aCCI-high"))

# 构建SVM模型
set.seed(1234)
svm_model <- svm(
  Group ~ .,
  data = train,
  kernel = "radial",
  method = "svmRadial",  # 使用 RBF 核函数
  metric = "ROC",  # 使用 AUC 作为评估标准
  probability = TRUE  # 启用概率预测
)

# 获取模型的权重系数
w <- t(svm_model$coefs) %*% svm_model$SV

# 计算每个特征的权重绝对值
feature_importance <- abs(w)
names(feature_importance) <- colnames(train)[-which(colnames(train) == "Group")]

# 计算每个特征的相对贡献（绝对值除以所有权重绝对值的总和）
relative_importance <- feature_importance / sum(feature_importance)
# write.csv(relative_importance, "SVM_importeance.csv")

# 将数据转换为数据框
importance_df <- data.frame(Feature = colnames(relative_importance), 
                            Importance = as.vector(relative_importance[1,]))

# 对特征按重要性排序
importance_df <- importance_df[order(importance_df$Importance, decreasing = FALSE), ]

# 将特征名称转为因子，并按重要性排序
importance_df$Feature <- factor(importance_df$Feature, levels = importance_df$Feature)

# 使用ggplot2绘制条形图
library(ggplot2)
ggplot(importance_df, aes(x = Feature, y = Importance)) +
  geom_bar(stat = "identity", fill = "#66C3A5") +
  coord_flip() +  # 横向显示特征
  labs(title = "Feature Importance from SVM Model",
       x = "Features",
       y = "Relative Importance") +
  theme_minimal()  # 使用简洁的主题
dev.off()

# 设置相对贡献的阈值 (例如，相对贡献小于 7% 的特征将被排除)
threshold <- 0.07
best_features <- colnames(relative_importance)[relative_importance >= threshold]

# 输出保留的特征
print(best_features)

# 训练集上的预测概率
train_pre_svm <- predict(svm_model, newdata = train, probability = TRUE)
train_prob_svm <- attr(train_pre_svm, "probabilities")[, "aCCI-high"]

# 验证集上的预测概率
validation_pre_svm <- predict(svm_model, newdata = validation, probability = TRUE)
validation_prob_svm <- attr(validation_pre_svm, "probabilities")[, "aCCI-high"]

# 计算训练集的ROC曲线
train_roc_svm <- roc(train$Group, train_prob_svm)

# 计算验证集的ROC曲线
validation_roc_svm <- roc(validation$Group, validation_prob_svm)

# 计算AUC值
train_auc_svm <- auc(train_roc_svm)
validation_auc_svm <- auc(validation_roc_svm)

# 绘制训练集和验证集的ROC曲线
plot(train_roc_svm, col = "#FD8D62", lwd = 3, main = "ROC Curve for Training and Validation Sets")
lines(validation_roc_svm, col = "#66C3A5", lwd = 3)
# 在图例中添加AUC值
legend_text <- c(paste0("Training (AUC = ", round(train_auc_svm, 3), ")"),
                 paste0("Validation (AUC = ", round(validation_auc_svm, 3), ")"))
legend("bottomright", legend = legend_text, col = c("#FD8D62", "#66C3A5"),lwd=2
       # ,cex = 0.8, box.lwd = 0, inset = c(0.01, 0.01)
)
dev.off()

####验证模型####
library(rms)
library(pROC)

# 读取数据
data <- read.csv("alzheimer_combined_imputed.csv", row.names = 1)
data$Group <- factor(data$Group, levels = c("aCCI-low", "aCCI-high"))

# 提取特征
varCoef <- read.csv("lasso_varcoff.csv", row.names = 1)
gg_dta <- read.csv("RF_importance.csv", row.names = 1)

features_lg <- c("age","Glucose","Phosphate")
features_lasso <- varCoef$Feature
features_rf <- gg_dta$names[gg_dta$col != "-"]

data_lg <- data[, c("Group", features_lg)]
data_lasso <- data[, c("Group", features_lasso)]
data_rf <- data[, c("Group", features_rf)]

set.seed(1234)
# 构建模型公式
model_formula_lg <- as.formula(paste("Group ~", paste(features_lg, collapse = " + ")))
model_formula_lg
model_formula_lasso <- as.formula(paste("Group ~", paste(features_lasso, collapse = " + ")))
model_formula_rf <- as.formula(paste("Group ~", paste(features_rf, collapse = " + ")))

# 配置 datadist
dd <- datadist(data)
options(datadist = "dd")
# 构建多因素逻辑回归模型
model_glm <- lrm(model_formula_lg, data = data_lg, x = TRUE, y = TRUE)

model_lasso <- lrm(model_formula_lasso, data = data_lasso, x = TRUE, y = TRUE)

model_rf <- lrm(model_formula_rf, data = data_rf, x = TRUE, y = TRUE)

#####ROC#####
# 预测概率
pred_glm <- predict(model_glm, newdata = data_lg, type = "fitted")
pred_lasso <- predict(model_lasso, newdata = data_lasso, type = "fitted")
pred_rf <- predict(model_rf, newdata = data_rf, type = "fitted")

# 计算 ROC 曲线
roc_glm <- roc(data_lg$Group, pred_glm)
roc_lasso <- roc(data_lasso$Group, pred_lasso)
roc_rf <- roc(data_rf$Group, pred_rf)

# 设置颜色
colors <- c("#1f77b4", "#FD8D62", "#66C3A5")
# 绘制 ROC 曲线
plot(roc_glm, col = colors[1], lwd = 2, main = "ROC Curves for Different Models")
lines(roc_lasso, col = colors[2], lwd = 2)
lines(roc_rf, col = colors[3], lwd = 2)

# 添加图例并显示 AUC 值
legend("bottomright", legend = c(paste0("GLM (AUC = ", round(auc(roc_glm), 3), ")"),
                                 paste0("LASSO (AUC = ", round(auc(roc_lasso), 3), ")"),
                                 paste0("RF (AUC = ", round(auc(roc_rf), 3), ")")),
       col = colors, lty = 1, lwd = 2)
dev.off()

#####校准曲线####
# 绘制校准曲线
cal_glm <- calibrate(model_glm, method = "boot", B = 1000)
cal_lasso <- calibrate(model_lasso, method = "boot", B = 1000)
cal_rf <- calibrate(model_rf, method = "boot", B = 1000)

plot(cal_glm, col = colors[1], lwd = 2, main = "Calibration Curves")
lines(cal_lasso, col = colors[2], lwd = 2)
lines(cal_rf, col = colors[3], lwd = 2)
legend("topleft", legend = c("GLM", "LASSO", "RF"),
       col = colors, lty = 1, lwd = 2)
dev.off()

#####DCA曲线#####
library(rmda)

# 将预测概率添加到对应的数据框中
data_lg$pred_glm <- pred_glm
data_lasso$pred_lasso <- pred_lasso
data_rf$pred_rf <- pred_rf

# 将 Group 变量转换为数值类型（0 和 1）
data_lg$Group <- as.numeric(data_lg$Group) - 1  # 假设原始 Group 为因子 "CCI-low" 和 "CCI-high"
data_lasso$Group <- as.numeric(data_lasso$Group) - 1
data_rf$Group <- as.numeric(data_rf$Group) - 1

# 计算 DCA 曲线
dca_glm <- decision_curve(Group ~ pred_glm, data = data_lg, family = binomial, thresholds = seq(0, 1, by = 0.01))
dca_lasso <- decision_curve(Group ~ pred_lasso, data = data_lasso, family = binomial, thresholds = seq(0, 1, by = 0.01))
dca_rf <- decision_curve(Group ~ pred_rf, data = data_rf, family = binomial, thresholds = seq(0, 1, by = 0.01))

# 绘制 DCA 曲线
plot_decision_curve(list(dca_glm, dca_lasso, dca_rf), curve.names = c("GLM", "LASSO", "RF"),
                    col = colors, lty = 1, lwd = 2, 
                    main = "Decision Curve Analysis")
dev.off()

####列线图模型####
final_model <- model_lasso

# 获取模型的系数
coefficients <- coef(final_model)
print(coefficients)

# 提取变量名称和系数
variables <- names(coefficients)
formula_string <- paste(coefficients, "*", variables, collapse = " + ")

# 构建最终的公式字符串
final_formula <- paste("Probability of High aCCI =", formula_string)
print(final_formula)

nomogram_model <- nomogram(final_model, 
                           fun = plogis, 
                           lp = FALSE, 
                           fun.at = seq(0.1, 0.9, by = 0.1),  # 修改这里的刻度间隔
                           funlabel = "Probability of High CCI")

# 绘制列线图
plot(nomogram_model)
dev.off()
