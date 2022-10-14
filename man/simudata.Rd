\name{simudata}
\alias{simudata}
\title{Generated simulated data.}
\usage{
simudata(data,refergene)
}
\arguments{
   \item{data}{The methylation beta value dataframe of control samples.The first column is the gene symbols and the remaining columns are the gene methylation values of the control samples.}
   \item{refergene}{The gene set with the smallest coefficient of variation in certain disease.}
}
\value{A methylation beta value dataframe for certain disease.The first column is the gene symbols and the remaining columns are the gene methylation values of certain disease.}
\description{
The function of simudata is used to generate simulated data for simulation experiment.
}
\examples{
data(example)
data(beta)
methy <- cpgTogene(beta,platform = "850k")
control <- methy[,c(1,match(clininf[grep("normal",clininf[,2]),1],colnames(methy)))]
refergene <- Selrefergene(data1,data2,threshold = 0.25)
simu <- simudata(data = control,refergene = refergene)
}
