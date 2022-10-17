# Introduction
The average methylation level of CpG sites within a promoter represented its gene methylation level. RankDMG is an algorithm to identify individualized differentially methylated genes (DMGs) of disease samples, and could be applied in a data with limited normal samples (small data) or no normal control (one-phenotype data). RankDMG could analyzed not only samples measured by the 450K or EPIC, but also the integration data of the two platforms. 
# Install
To install the RankDMG, install from github using devtools:
```
library(devtools)
install_github("FunMoy/RankDMG")
```
# Calculate the methylation levels of genes
Here, the methylation level of a gene was represented by the average methylation level of CpG sites within its promoter. The program of **cpgTogene** is to calculate the methylation level of each gene. We used the dataset, GSE149282, as an example to illustrate how to used the program of **cpgTogene**.
```
library(RankDMG)
data(beta)
methy <- cpgTogene(beta,platform = "850k")
```
The detailed description of **cpgTogene** was shown in the manaul page of **cpgTogene**. 
The data of **beta** is the outputted matrix from ChAMP. A value of the matrix represented the methylation level of a CpG site in a sample. The row of the matrix represented the CpG probes, and the column of the matrix represented the samples.


# Identify reference genes
The reference genes are those with the smallest coefficient of variation across cancer and normal samples. we used two datasets, with the disease and normal samples, to evaluate the reproducibility of reference genes, and the overlapped reference genes of the two datasets were finally defined as reference genes for a disease analysis. The program of **Selrefergene** was used to identify reference genes. The following is an example:
```
library(RankDMG)
data(example)
refergene <- Selrefergene(data1,data2,threshold = 0.25)
```
The detailed of **Selrefergene** was shown in the manual page of **Selrefergene**. The parameters of data1 and data2 were the two datasets used to identify reference genes.
# How to generate simulated data
The simulated data was used to evaluate the performance of RankDMG. Here, to ensure the inherent characteristics of normal samples, we used the normal samples of an independent dataset to simulate the corresponding disease samples. The program of **simudata** was used to generate simulated data, and the detailed description was shown in the manual page of **simudata**. We used an example to explain the simulation experiment:
```
library(RankDMG)
data(example)
data(beta)
methy <- cpgTogene(beta,platform = "850k")
control <- methy[,c(1,match(clininf[grep("normal",clininf[,2]),1],colnames(methy)))]
refergene <- Selrefergene(data1,data2,threshold = 0.25)
simu <- simudata(data = control,refergene = refergene)
```
The parameter of data was the methylation levels of genes in normal samples of an independent dataset, and refergene was the reference genes.

# Identify individualized and population-level DMGs
The program of **RankDMG** was the main program in our algorithm. **RankDMG** could identify individualized DMGs of disease samples. After identifying individualized DMGs in a one-phenotype data or small data, **RankDMG** not only can infer genes hypermethylated (hypomethylated) in a nonrandom high percentage of samples, namely population-level hypermethylated (hypomethylated) genes, using binomial test, but also can provide the individualized DMGs to explain the inter-individual heterogeneity. This may meet the different analysis needs of users.


# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to:
Qi Fan <FunMoy@163.com>; Haidan Yan <Joyan168@126.com>
