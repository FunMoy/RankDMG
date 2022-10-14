\name{cpgTogene}
\alias{cpgTogene}
\title{Data translator}
\usage{
cpgTogene(data,platform = "integrated")
}
\arguments{
   \item{data}{A beta value matrix with rows labeling the CpGs and columns labeling samples.}
   \item{platform}{The platform to be used("450k", "850k" or "integrated").}
   }
\value{The methylation beta value dataframe.The first columns is the gene symbols and the remaining columns are the gene methylation values.}

\description{A beta value matrix with rows labeling the CpGs and columns labeling samples is transformed into a beta value dataframe with gene symbols and gene methylation beta values.
}
\examples{
data(beta)
methy <- cpgTogene(beta,platform = "850k")
}
