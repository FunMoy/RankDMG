RankDMG <- function(normal,case,refergene,freq = 0.99,threshold = 0.05){
  ref.gene = intersect(normal[,1],case[,1])
  index.nor = match(ref.gene,normal[,1])
  index.case = match(ref.gene,case[,1])
  normal = normal[index.nor,]
  case = case[index.case,]
  gene=as.data.frame(c(1:length(ref.gene)))
  index = na.omit(match(refergene,ref.gene))
  ingene = index
  control_exp<-normal[,-1]
  case_exp<-case[,-1]
  print("differential methylated analysis... \n");
  outlier_dir<-NULL;
  outlier_pvalue<-NULL;
  Lgene<-dim(gene)[1];
  for (k in 1:Lgene){
    Nnorm=as.matrix(control_exp[k,])
    colN<-dim(control_exp)[2]
    colC<-dim(case_exp)[2]
    ingene = intersect(ingene,gene[,1])
    lstable = length(ingene)
    stable.matrix = control_exp[match(ingene,gene[,1]),]
    N_tmp=matrix(rep(Nnorm,lstable),ncol=colN,byrow=T)-stable.matrix
    Nloc_up=ingene[which((rowSums(N_tmp>0)/colN)>freq)]
    Nloc_down=ingene[which((rowSums(N_tmp<0)/colN)>freq)]
    Nloc_up=match(Nloc_up,gene[,1])
    Nloc_down=match(Nloc_down,gene[,1])
    reverse=matrix(0,colC,4)
    reverse[,1]=rep(length(Nloc_up),colC)
    reverse[,2]=rep(length(Nloc_down),colC)
    Tcanc=as.matrix(case_exp[k,])
    if (length(Nloc_up)>0){
      N_tmp=matrix(rep(Tcanc,length(Nloc_up)),ncol=colC,byrow=T)-case_exp[Nloc_up,]
      case_p=colSums(N_tmp<0)
      reverse[,3]=case_p
    }
    if (length(Nloc_down)>0){
      N_tmpp=matrix(rep(Tcanc,length(Nloc_down)),ncol=colC,byrow=T)-case_exp[Nloc_down,]
      case_pp=colSums(N_tmpp>0)
      reverse[,4]=case_pp;
    }
    GenePair_sig=NULL
    GenePair=rep(0,colC)
    GenePair[which(reverse[,3]>reverse[,4])]<--1
    GenePair[which(reverse[,3]<reverse[,4])]<-1
    tmp=matrix(c(reverse[,1],reverse[,2],reverse[,1]-reverse[,3]+reverse[,4], reverse[,2]-reverse[,4]+reverse[,3]),ncol=4)
    GenePair_sig<-apply(tmp,1,function(x) fisher.test(matrix(x,ncol=2,byrow=T))$p.value)
    outlier_dir=rbind(outlier_dir,GenePair)
    outlier_pvalue=rbind(outlier_pvalue,GenePair_sig)
  }
  fdr<-apply(outlier_pvalue,2,function(x) p.adjust(x,method="fdr",length(x)))
  Methout<-list(ref.gene,outlier_dir,fdr);
  names(Methout)<-c('gene','dir','fdr')

  id = Methout$gene
  dire = Methout$dir
  fdr = Methout$fdr
  for(i in 1:dim(fdr)[2])
  {
    index<-which(fdr[,i]<=0.05);
    dire[-index,i]<-0;
  }
  dire_all<-abs(dire);
  P<-matrix(0,1,dim(fdr)[2]);
  for(j in 1:dim(fdr)[2])
  {
    P[j]<-sum(dire_all[,j])/length(dire[,1]);
  }
  P_all=rowSums(P)/length(dire[1,]);

  dire_up=dire;
  index1<-matrix(0,1,dim(fdr)[2]);
  index2<-matrix(0,1,dim(fdr)[2]);
  for(i in 1:dim(fdr)[2])
  {
    index1[i]<-length(which(dire_up[,i]==-1));
    index2[i]<-length(which(dire_up[,i]==1));
  }
  P<-matrix(0,1,dim(fdr)[2]);
  for(j in 1:dim(fdr)[2])
  {
    P[j]<-index2[j]/(index1[j]+index2[j])
  }
  P_up=rowSums(P)/dim(fdr)[2];
  P0_up<-P_all*P_up;

  p_value<-matrix(0,dim(dire_up)[1],1);
  for(x in 1:dim(dire_up)[1])
  {
    k<-length(which(dire_up[x,]==1));
    n2<-dim(fdr)[2];
    p_value[x]<-1-pbinom(k-1,n2,P0_up);
  }
  fdr_up<-p.adjust(p_value,method="fdr",length(p_value));
  result_index<-which(fdr_up<=threshold);
  result_id_up<-id[result_index];

  dire_down=dire;
  index1<-matrix(0,1,dim(fdr)[2]);
  index2<-matrix(0,1,dim(fdr)[2]);
  for(i in 1:dim(fdr)[2])
  {
    index1[i]<-length(which(dire_down[,i]==-1));
    index2[i]<-length(which(dire_down[,i]==1));
  }
  P<-matrix(0,1,dim(fdr)[2]);
  for(j in 1:dim(fdr)[2])
  {
    P[j]<-index1[j]/(index1[j]+index2[j])
  }
  P_down=rowSums(P)/dim(fdr)[2];

  P0_down<-P_all*P_down;

  p_value<-matrix(0,dim(dire_down)[1],1);
  for(x in 1:dim(dire_down)[1])
  {
    k<-length(which(dire_down[x,]==-1));
    n2<-dim(fdr)[2];
    p_value[x]<-1-pbinom(k-1,n2,P0_down);
  }
  fdr_down<-p.adjust(p_value,method="fdr",length(p_value));
  result_index<-which(fdr_down<=threshold);
  result_id_down<-id[result_index];
  inter.id = intersect(result_id_down,result_id_up)
  if (length(inter.id)>0) {
    result_id_up = result_id_up[-match(inter.id,result_id_up)]
    result_id_down = result_id_down[-match(inter.id,result_id_down)]
  }
  up = result_id_up
  down = result_id_down
  result = list(hyper=up,hypo=down)
  out = list(individual=Methout,population=result)
  return(out)
}
