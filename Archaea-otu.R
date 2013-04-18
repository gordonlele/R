library(gplots)

class <- read.csv("Archaea-otu-phylum-0409.csv",sep=";")
row.names(class) <- class$Groups
#删除数据全部为0的行或者缺少数据的行
class <- delete_all_zero_rows(class,3)

#按照phylum分组,即class的第2列,分组包含的行多的在前
class <- order_by_groups(class,2)

#如果是新的数据集,只要修改class$phylum为需要排序的列即可
#获得自定义的颜色值集合,以定义左下角的图标样例的颜色
vec.dphylum <- GetTopColors(nlevels(class$phylum))
#用因子的方法获得图形索引id,用于之后的显示
f.dphylum <- factor(class$phylum)
dphylum.color <- rep(0,length(f.dphylum))
for(i in 1:length(dphylum.color))
{
  dphylum.color[i] <- vec.dphylum[ f.dphylum[i]==unique(f.dphylum) ]
}
#只保留数据列
class <- class[,3:ncol(class)]
class <- set_data_to_percent(class)  #数据处理为每个数占每列和的百分百
class <- data.matrix(class)

#定义36个colorkey的分割点,分割点的单位是不同的
bk <- c(0,0) #0值用白色标示,特殊处理
bk <- append(bk,seq(0.0001,0.001,by=0.0001))
bk <- append(bk,seq(0.002,0.01,by=0.001))
bk <- append(bk,seq(0.02,0.1,by=0.01))
bk <- append(bk,seq(0.2,1,by=0.1))
#定义渐变颜色:white blue green orange red
mycol <- colorRampPalette(colors=c("blue","green","orange","red"))
mycol <- append(c("white"),mycol(37)) 
#使用自己改编(heatmap.2)的函数绘制
classmap <- myheatmap2(class,scale="none",key=TRUE, symkey=FALSE, density.info="none", trace="none",
                       RowSideColors=dphylum.color,Rowv=NA,labRow=NA,
                       col=mycol,breaks=bk,
                       sample_legend=unique(f.dphylum),sample_fill=vec.dphylum)
