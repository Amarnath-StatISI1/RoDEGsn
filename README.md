## RoDEGsn
Robust modeling of noisy miRNA expression data based on expression values for finding deferentially expressed genes

## Introduction
RoDEGsn- A parametric approach to find the deferentially expressed miRNAs where the distributions of the expression values are moderately skewed. Here, the inherent model used is based on the famous skew normal (SN) distribution, which includes the both symmetric and skewed classes. This function takes the log transformed and normalized expression values along with the tuning parameter $\alpha$ and return the top $d$ miRNAs along with the corresponding adjusted p-values.

## Installation of required packages

Here we have used R 3.6.1. First we have to download packages from Bioconductor. In this connection, we have to install "BiocManager" at first and then proceed as follows to install "limma","Biobase" and "Geoquery". Note that "BiocManager" is available for R version 3.6 and above, so we must have R version 3.6 and above.

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("Biobase")
BiocManager::install("GEOquery")
```
Besides these, we have to install the following package also:
```r
utils::install.packages("sn")
```

## Loading of required packages
Before running the program, we have to load the following packages
```r
library(limma)
library(Biobase)
library(GEOquery)
library(sn)
library(MASS)
library(doParallel)
```

## Loading the required data set

Here we have used Mice tumor data (Accession number- GSE53388) as example. We have run the following code to load the data from GEO website:

```r
#load series and platform data from GEO

gset=getGEO("GSE53388", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx =grep("GPL1261", attr(gset, "names")) 
else idx=1
gset=gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset)=make.names(fvarLabels(gset))

# group names for all samples

gsms <- "111111111111111111000000000000000000"
sml=c()
for (i in 1:nchar(gsms)) 
{ sml[i]=substr(gsms,i,i) 
}
```

## Data preprocessing
Here we have used log transformation and quantile normalization for filtering the data.
```r
ex=exprs(gset)
qx=as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=(qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)]=NaN
exprs(gset)=log2(ex) }

```
## Application of RoDEGsn
Here, we have applied our method to the first 500 genes of the data and summarized the adjusted p-value of the top 10 genes. We have considered tuning parameter $\alpha=0.5$. 

```r
gd=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
ex_gp1=ex[,which(gd==1)]
ex_gp2=ex[,which(gd==0)]
p_val05=rep(0,500)
result=RoDEGsn(ex[1:500,],gd,10,0.5)
```
## Output
```r
 result
                     p_adj
1415737_at    1.749096e-39
1415777_at    8.665615e-38
1415783_at    5.471823e-36
1415806_at    2.901688e-34
1415817_s_at  6.738432e-33
1415862_at    3.750261e-31
1415740_at   4.457271e-28
1415676_a_at 8.229201e-24
1415696_at   2.442977e-18
```
