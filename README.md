# Introduction
RankDMG is an algorithm to identification of differential methylation genes(DMGs) in methylation array data. It can identify of individual-level and population-level DMGs in methylation array data.
# Install
To install the RankDMG, install from github using devtools:
```
library(devtools)
install_github("FunMoy/RankDMG")
```
# How to convert a beta value matrix with CpG sites as rownames into a beta value dataframe with gene symbols 
The purpose of our algorithm is to identify DMGs. So our first step is to convert beta value matrix with rows labeling the CpGs and columns labeling samples into beta value dataframe with gene symbols. We use the code **cpgTogene** as translator. There is a example to show how to use it:
```
library(RankDMG)
data(beta)
methy <- cpgTogene(beta,platform = "850k")
```
The detailed of **cpgTogene** in manaul page of **cpgTogene**. The **beta** is a beta value matrix for the GSE149282 dataset.

# How to select refergene
The refergene is a gene set with the smallest coefficient of variation in certain disease, which is important for our 
algorithm to identify of DMGs for certain disease. We could use the code **Selrefergene** to select refergene. There is a example:
```
library(RankDMG)
data(example)
refergene <- Selrefergene(data1,data2,threshold = 0.25)
```
The detailed of **Selrefergene** in manaul page of **Selrefergene**. And the **example** contains several data sets, which uesd to demonstrate the code.These data sets are **clininf**,**data1**,**data2** and **normal**. The detailed of **example** in manaul page of **example**. 

# How to generate simulated data
We used simulated data to test the performance of our algorithm. And there is a example to generate simulated data:
```
library(RankDMG)
data(example)
data(beta)
methy <- cpgTogene(beta,platform = "850k")
control <- methy[,c(1,match(clininf[grep("normal",clininf[,2]),1],colnames(methy)))]
refergene <- Selrefergene(data1,data2,threshold = 0.25)
simu <- simudata(data = control,refergene = refergene)
```
The code **simudata** is uesd to generate simulated data, its specific usage in manaul page of **simudata**.

# How to identify of individual-level and population-level DMGs
The code **RankDMG** is the main program in our algorithm, which is used to identify of individual-level and population-level DMGs.  Its specific usage in manaul page of **RankDMG**.


# Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to:
Qi Fan <FunMoy@163.com>; Haidan Yan <Joyan168@126.com>
