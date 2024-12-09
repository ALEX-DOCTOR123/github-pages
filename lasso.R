
set.seed(123)
library(glmnet)                   #引用包
inputFile="diffGeneExp.txt"       #输入文件


#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#构建模型
x=as.matrix(rt)
y=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#输出筛选的特征基因
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)

rt1=t(rt)
lassoexp=rt1[lassoGene,,drop=F]
lassoexp=as.data.frame(lassoexp)
colnames(lassoexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(lassoexp))
write.table(lassoexp, file="LASSO.geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)
