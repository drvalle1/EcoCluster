## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ---- eval=TRUE----------------------------------------------------------
library(tsbp)
library(mvtnorm)
data(data_mixture)
data(data_sbm)
data(data_sam) 

## ---- eval=TRUE----------------------------------------------------------
# fake data mixtures
data_mixture[1:5, 1:10]
# fake data SAM model
y[1:5, 1:10]
head(xmat)
# fake data SBM
head(data_sbm)

## ---- eval=TRUE----------------------------------------------------------
ngibbs <- 10
PosMix <- mixture.gibbs.main.func(dat=data_mixture,ngroup=50,ngibbs=ngibbs,burnin=ngibbs/2)

## ---- eval=TRUE----------------------------------------------------------
PosMix$phi[1:5,1:10]
PosMix$theta[1:5,1:5]

## ---- eval=TRUE----------------------------------------------------------
ngibbs <- 10
ngroup.loc <- 10
ngroup.spp <- 10
PosSBM <- SBM(dat=data_sbm,ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp,ngibbs=ngibbs,burnin=ngibbs/2)
names(PosSBM)
PosSBM$gamma

## ---- eval=TRUE----------------------------------------------------------
ngibbs <- 10
ngroups <- 50
PosSAM <- gibbs.SAM(y=y,xmat=xmat,ngroups=ngroups,ngibbs=ngibbs,burnin=ngibbs/2)
names(PosSAM)

