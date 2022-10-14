\name{RankDMG}
\alias{RankDMG}
\title{Identification of differential methylation genes(DMGs) in methylation array data.}
\usage{
RankDMG(normal,case,refergene,freq = 0.99,threshold = 0.05)
}
\arguments{
   \item{normal}{The methylation beta value dataframe of normal samples.The first column is the gene symbols and the remaining columns are the gene methylation values of the normal samples.}
   \item{case}{The methylation beta value dataframe of disease samples.The first column is the gene symbols and the remaining columns are the gene methylation values of the disease samples.}
   \item{refergene}{The gene set with the smallest coefficient of variation in certain disease.}
   \item{freq}{The criteria for identifying stable gene pairs. The default setting of freq is 0.99.}
   \item{threshold}{The threshold for identifying DMGs.The default setting of threshold is 0.05.}
}
\value{
A list with the following two items.

individual
A list with individual-level result.It contains all gene symbols,individual-level DMGs in each case sample (value is 1 for hypermethylated,value is -1 for hypomethylated) and adjusted pvalue
using the Benjamini and Hochberg approach.

population
A list with population-level DMGs.
}
\description{
RankDMG can identify of individual-level and population-level DMGs in methylation array data.
}
\examples{
data(example)
data(beta)
methy <- cpgTogene(beta,platform = "850k")
case <- methy[,c(1,match(clininf[grep("cancer",clininf[,2]),1],colnames(methy)))]
refergene <- Selrefergene(data1,data2,threshold = 0.25)
DMGs <- RankDMG(normal = normal,case =case,refergene = refergene,freq = 0.99,threshold = 0.05)
}
