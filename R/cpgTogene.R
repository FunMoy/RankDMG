cpgTogene <- function(data,platform = "integrated"){
  if (platform == "450k") {
    probeid_geneid = cpg.450
  }
  if (platform == "850k") {
    probeid_geneid = cpg.850
  }
  if (platform == "integrated") {
    probeid_geneid = cpg.integrated
  }
  colname = c('gene',colnames(data))
  reid=match(probeid_geneid[,1],rownames(data))
  index=which(is.na(reid)==F)
  reid=reid[index]
  probeid_geneid=probeid_geneid[index,]
  data=data[reid,]
  geneid=probeid_geneid[,2]
  exp.matrix=data
  uni.geneid=unique(geneid)
  matrix.mean <- lapply(uni.geneid,function(x){
    index.i=which(geneid==x)
    if(length(index.i)>1){
      row.i=colMeans(exp.matrix[index.i,])
    }
    else{
      row.i=exp.matrix[index.i,]
    }
    row.i=c(as.matrix(x),as.matrix(row.i))
    return(row.i)
  })
  matrix.mean = do.call('rbind',matrix.mean)
  exp = matrix.mean[,-1]
  exp = matrix(as.numeric(exp),byrow = F,ncol = dim(exp)[2])
  data = cbind.data.frame(matrix.mean[,1],exp)
  colnames(data)=colname
  return(data)
}
