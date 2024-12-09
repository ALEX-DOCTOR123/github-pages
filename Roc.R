
library(glmnet)  
library(survival)
library(tidyverse)
library(caret)
library(pROC)
library(limma)
rt=read.table("TCGA_TPM.txt", header=T, sep="\t", check.names=F, row.names=1)
inTrain<-createDataPartition(y=1:length(colnames(rt)),p=0.7,list=F)
train<-rt[,inTrain]
test<-rt[,-inTrain]
data=t(train)
x=as.matrix(data)
y=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)", "\\4",rownames(data))
y=sapply(strsplit(y,""), "[", 1)
fit=glmnet(x, y, family = "binomial", alpha=1)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()


coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

rt1=t(data)
lassoexp=rt1[lassoGene,,drop=F]
lassoexp=as.data.frame(lassoexp)
lassoexptest=as.data.frame(test[lassoGene,,drop=F])

clinical=read.table("clinical.txt",header = T,sep = "\t",row.names = 1)
clinical=na.omit(clinical)
clinical$futime=clinical$futime/365
tumordata= lassoexp%>% dplyr::select(str_which(colnames(.), "-01A"))%>%as.data.frame()#仅保留肿瘤数据
tumordatatest=lassoexptest%>%dplyr::select(str_which(colnames(.), "-01A"))%>%as.data.frame()#仅保留肿瘤数据
colnames(tumordata)=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)", "\\1\\-\\2\\-\\3\\",colnames(tumordata))
colnames(tumordatatest)=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)", "\\1\\-\\2\\-\\3\\",colnames(tumordatatest))
sameid=intersect(colnames(tumordata),rownames(clinical))
sameidtest=intersect(colnames(tumordatatest),rownames(clinical))
tumordatatest=tumordatatest[,sameidtest]
tumordata=tumordata[,sameid]
clinicaltest=clinical[sameidtest,]
clinical=clinical[sameid,]
tumordata=t(tumordata)%>%as.data.frame()
tumordatatest=t(tumordatatest)%>%as.data.frame()
rt=data.frame(futime=clinical$futime,fustat=clinical$fustat,tumordata)
rttest=data.frame(futime=clinicaltest$futime,fustat=clinicaltest$fustat,tumordatatest)
cox <- coxph(Surv(futime, fustat) ~ ., data = rt)
cox=step(cox,direction = "both")
riskScore=as.data.frame(predict(cox,type="risk",newdata=rt))
colnames(riskScore)="risk"
risk=data.frame(riskScore,fustat=clinical$fustat)
roc1=roc(risk$fustat, as.numeric(risk$risk),smooth=T)
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC_train.pdf", width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Train")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()
riskScoretest=as.data.frame(predict(cox,type="risk",newdata=rttest))
colnames(riskScoretest)="risk"
risktest=data.frame(riskScoretest,fustat=clinicaltest$fustat)
roc1=roc(risktest$fustat, as.numeric(risktest$risk),smooth=T)
ci1=ci.auc(roc1, method="bootstrap")
ciVec=as.numeric(ci1)
pdf(file="ROC_test.pdf", width=5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Test")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()
