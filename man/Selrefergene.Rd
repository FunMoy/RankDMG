\name{Selrefergene}
\alias{Selrefergene}
\title{Selected a gene set with the smallest coefficient of variation in certain disease.}
\usage{
Selrefergene(data1,data2,threshold = 0.25)
}
\arguments{
   \item{data1}{A methylation beta value dataframe for certain disease (contains control and case samples).The first columns is the gene symbols and the remaining columns are the gene methylation values.}
   \item{data2}{Another methylation beta value dataframe for same disease as data1 (contains control and case samples).The first columns is the gene symbols and the remaining columns are the gene methylation values.}
   \item{threshold}{The threshold to selected a gene set with the smallest coefficient of variation.The default setting of threshold is 0.25.}
}
\value{Character.A gene set.}
\description{
The function of \code{Selrefergene} is used to selected a gene set with the smallest coefficient of variation in certain disease.
}
\examples{
data(example)
refergene <- Selrefergene(data1,data2,threshold = 0.25)
}
