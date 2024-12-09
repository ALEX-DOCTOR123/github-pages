library(xgboost)
library(caret)
library(tidyverse)
library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)

set.seed(123)
data<-read.table(file="diffGeneExp.txt",sep = "\t",header = T,check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

TrainControl <- trainControl( method = "repeatedcv", number = 10, repeats = 4)
model<- train(x=data,y=as.factor(group),  method = "xgbTree", trControl = TrainControl,verbose = FALSE)


plot(varImp(model))
importance <- varImp(model)
head(importance)
important <- as.data.frame(importance$importance) 
a<-important
varimpdf <- data.frame(var = row.names(a),
                       impor = a[,1])


ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
  geom_col(colour = "lightblue",fill = "lightblue")+
  labs(title="Feature gene importance (XGBoost)", x="",y = "importance")+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  theme(axis.text.x = element_text(size = 3))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))



