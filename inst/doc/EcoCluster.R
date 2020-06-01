## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----eval=TRUE----------------------------------------------------------------
rm(list=ls(all=TRUE))
set.seed(1)

#general settings
nloc=300 #number of locations
nspp=50  #number of species
ngroup=5 #number of groups

#set parameters
z=sample(1:ngroup,size=nloc,replace=T)         #cluster assignment of each location
phi=matrix(rbeta(ngroup*nspp,1,1),ngroup,nspp) #species composition of each cluster  

#generate data
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  phi1=phi[z[i],]
  y[i,]=rbinom(nspp,size=1,prob=phi1)
}
colnames(y)=paste0('spp',1:nspp)
rownames(y)=paste0('loc',1:nloc)

head(y[,1:13])

## ----eval=TRUE,results='hide'-------------------------------------------------
library('EcoCluster')

