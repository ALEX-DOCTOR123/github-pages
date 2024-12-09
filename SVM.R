
#引用包
library(e1071)
library(kernlab)
library(caret)

set.seed(123)
inputFile="diffGeneExp.txt"        #输入文件

#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\-(.*)\\-(.*)\\-(.*)\\-(.*)", "\\5", row.names(data))

#SVM-RFE分析
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")

#绘制图形
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#标注交叉验证误差最小的点
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

#输出选择的基因
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)
rt1=t(data)
svmexp=rt1[lassoGene,,drop=F]
svmexp=as.data.frame(svmexp)
colnames(svmexp)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", colnames(svmexp))
write.table(lassoexp, file="SVM.geneExp.txt", sep="\t", quote=F, row.names=T, col.names=T)
