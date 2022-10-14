simudata <- function(data,refergene){
  gene = data[,1]
  methy = data[,-1]
  can.gene = gene[-match(refergene,gene,nomatch = 0)]
  can.methy = methy[-match(refergene,gene,nomatch = 0),]
  samplenum=dim(methy)[2]
  for(j in 1:samplenum){
    can.select.up1 = which(can.methy[,j]<0.75)###hyper 0.2-0.25
    select.dmg.up1.index1=sample(can.select.up1,size=300)##select 300 genes

    can.select.up2 = which(can.methy[,j]<0.70)###hyper 0.25-0.30
    leftid1=setdiff(can.select.up2,select.dmg.up1.index1)
    select.dmg.up2.index2=sample(leftid1,size=300)##select 300 genes

    can.select.up3 = which(can.methy[,j]<0.65)###hyper 0.30-0.35
    union12=union(select.dmg.up1.index1,select.dmg.up2.index2)
    leftid2=setdiff(can.select.up3,union12)
    select.dmg.up3.index3=sample(leftid2,size=300)##select 300 genes

    can.select.up4 = which(can.methy[,j]<0.60)###hyper 0.35-0.40
    union123=union(union12,select.dmg.up3.index3)
    leftid3=setdiff(can.select.up4,union123)
    select.dmg.up4.index4=sample(leftid3,size=300)##select 300 genes

    can.select.down1 = which(can.methy[,j]>0.25)###hypo 0.20-0.25
    union1234=union(union123,select.dmg.up4.index4)
    leftid4=setdiff(can.select.down1,union1234)
    select.dmg.down1.index1=sample(leftid4,size=300)##select 300 genes

    can.select.down2 = which(can.methy[,j]>0.30)###hypo 0.25-0.30
    union12345=union(union1234,select.dmg.down1.index1)
    leftid5=setdiff(can.select.down2,union12345)
    select.dmg.down2.index2=sample(leftid5,size=300)##select 300 genes

    can.select.down3 = which(can.methy[,j]>0.35)###hypo 0.30-0.35
    union123456=union(union12345,select.dmg.down2.index2)
    leftid6=setdiff(can.select.down3,union123456)
    select.dmg.down3.index3=sample(leftid6,size=300)##select 300 genes

    can.select.down4 = which(can.methy[,j]>0.40)###hypo 0.35-0.40
    union1234567=union(union123456,select.dmg.down3.index3)
    leftid7=setdiff(can.select.down4,union1234567)
    select.dmg.down4.index4=sample(leftid7,size=300)##select 300 genes

    select.upg1=can.gene[select.dmg.up1.index1]
    select.upg2=can.gene[select.dmg.up2.index2]
    select.upg3=can.gene[select.dmg.up3.index3]
    select.upg4=can.gene[select.dmg.up4.index4]

    select.downg1=can.gene[select.dmg.down1.index1]
    select.downg2=can.gene[select.dmg.down2.index2]
    select.downg3=can.gene[select.dmg.down3.index3]
    select.downg4=can.gene[select.dmg.down4.index4]

    up.index1 = match(select.upg1,gene)
    up.index2 = match(select.upg2,gene)
    up.index3 = match(select.upg3,gene)
    up.index4 = match(select.upg4,gene)

    down.index1=match(select.downg1,gene)
    down.index2=match(select.downg2,gene)
    down.index3=match(select.downg3,gene)
    down.index4=match(select.downg4,gene)
    for(i in 1:300){
      methy[up.index1[i],j]=methy[up.index1[i],j]+runif(1,0.20,0.25)
      methy[up.index2[i],j]=methy[up.index2[i],j]+runif(1,0.25,0.30)
      methy[up.index3[i],j]=methy[up.index3[i],j]+runif(1,0.30,0.35)
      methy[up.index4[i],j]=methy[up.index4[i],j]+runif(1,0.35,0.40)
      methy[down.index1[i],j]=methy[down.index1[i],j]-runif(1,0.20,0.25)
      methy[down.index2[i],j]=methy[down.index2[i],j]-runif(1,0.25,0.30)
      methy[down.index3[i],j]=methy[down.index3[i],j]-runif(1,0.30,0.35)
      methy[down.index4[i],j]=methy[down.index4[i],j]-runif(1,0.35,0.40)
    }

  }###end for j
  methy = cbind.data.frame(gene,methy)
  simu=methy
  return(simu)
}
