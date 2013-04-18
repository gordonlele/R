source('~/Documents/个人数据中心/git/R/my_heatmap.R')

#把一个数据框的每一列的值初始化为这个数除以这一列总和*100,即这个数在列和的百分百
set_data_to_percent <- function (x)
{
  sum.x <- apply(x,2,sum)
  for (i in 1:ncol(x))
  {
    x[,i] <- x[,i]/sum.x[i]*100
  }
  x
}


#删除列全部为0的行,无效数据
delete_all_zero_rows <- function (x,startCol=1,endCol=ncol(x))
{
  
  s.row <- rep(TRUE,nrow(x))
  s.value <- apply(x[,startCol:endCol],1,sum)
  for(i1 in 1:length(s.value))
  {
    if (s.value[i1] == 0 | is.na(s.value[i1]))
      s.row[i1] <- FALSE
  }
  x <- x[which(s.row),]
}
