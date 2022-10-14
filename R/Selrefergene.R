Selrefergene <- function(data1,data2,threshold = 0.25){
  cv.de <- function(test,num.sample){
    data <- matrix(unlist(test[,-1]),byrow = F,ncol = num.sample)
    id <- test[,1]
    sd <- apply(data,1,function(x)
      sd(x))
    mean <- apply(data,1,function(x)
      mean(x))
    cv <- sd/mean
    all <- cbind.data.frame(id,sd,mean,cv)
    all <- all[order(all[,4],decreasing = F),]
    return(all)}

  inter.id <- intersect(data1[,1],data2[,1])
  data1 <- data1[match(inter.id,data1[,1]),]
  data2 <- data2[match(inter.id,data2[,1]),]
  data1.cv <- cv.de(data1,dim(data1[,-1])[2])
  data2.cv <- cv.de(data2,dim(data2[,-1])[2])
  data1.stable <- data1.cv[1:ceiling(dim(data1.cv)[1]*threshold),]
  data2.stable <- data2.cv[1:ceiling(dim(data2.cv)[1]*threshold),]
  refergene <- intersect(data1.stable[,1],data2.stable[,1])
  return(refergene)
}
