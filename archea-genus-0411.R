library(gplots)

class <- read.csv("archaea-genus0411.csv",sep=";")
row.names(class) <- class$Groups

#只保留数据列
class <- class[,2:ncol(class)]/100
class <- data.matrix(class)

#定义36个colorkey的分割点,分割点的单位是不同的
bk <- c(0,0.0001,0.0007,0.002,0.006,0.01,0.05,0.09,0.1,0.2,0.3,0.5,0.6,0.8,0.9,1) 
#定义渐变颜色:white blue green orange red
mycol <- colorRampPalette(colors=c("gray","orange","red","purple"))
mycol <- mycol(15)
#使用自己改编(heatmap.2)的函数绘制
classmap <- myheatmap3(class,scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",
                       col=mycol,breaks=bk,margin=c(5,12))
